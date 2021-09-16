import re
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from BCBio import GFF

from os import mkdir
from os import environ
from os.path import join
from os.path import exists
from sys import argv

from binascii import crc32
from subprocess import PIPE
from subprocess import Popen
from collections import namedtuple

from log import logging
from cli import parse_args
from cache import Cache
from tools import check_all_tools
from tools import Tool
from utilities import IOWrapper
from utilities import ToolExecutionError


log = logging.getLogger('sudan')


FORMAT_FASTA = 'fasta'
MIN_CONTIG_LEN = 200
MAX_TRNA_LEN = 500



def deeplen(iterable):
    return sum(map(len, iterable.values()))


class Context:
    def __init__(self):
        self.CONTIG = {}
        self.TOTAL_BP = -1
        self.CPU = 8
        self.KINGDOM = 'Bacteria'
        self.args = None
        self.all_rna = []

    def add_features(self, features):
        """accepts features in a dict {contig_id: [features]}
        appends given features to contig with associated id
        """
        for cid, fs in features.items():
            if cid in self.CONTIG:
                self.CONTIG[cid].features += fs


context = Context()


def prepare_input(seq_iterator):
    ncontig = 0
    total_bp = 0
    ids = set()

    for contig in seq_iterator:
        if len(contig.seq) < MIN_CONTIG_LEN:
            log.debug(f'Skipping short contig {contig.id}')

        ncontig += 1

        if '|' in contig.id:
            old_id = contig.id
            contig.id = old_id.replace('|', '_')
            log.debug(f'Changing illegal name {old_id} to {contig.id}')

        if contig.id in ids:
            raise ValueError(
                f'Error: duplicate id in sequence file {contig.id}')
        ids.add(contig.id)

        seq = str(contig.seq.upper())
        seq, n = re.subn(r'[*-]', '', seq)
        if n > 0:
            log.debug(f'Removed {n} gaps/pads')

        seq, n = re.subn(r'[^ACTG]', 'N', seq)
        if n > 0:
            log.debug(f'Replaced {n} wacky IUPAC with N')

        # todo: contig id rename

        contig.seq = Seq(seq)
        contig.description = ''
        total_bp += len(contig.seq)

        contig.annotations['molecule_type'] = 'DNA'
        context.CONTIG[contig.id] = contig
        yield contig

    if ncontig <= 0:
        raise ValueError(f'FASTA file {input} contains no suitable sequence entries')

    log.info(f'Total {total_bp} BP read')
    context.TOTAL_BP = total_bp


def cleanup_input(input):
    output = join(context.args.output_dir, 'after_cleanup.fna')

    data = SeqIO.parse(input, format=FORMAT_FASTA)
    prepared = prepare_input(data)

    SeqIO.write(prepared, output, format=FORMAT_FASTA)
    log.info(f'File {output} written')

    return output


def _run_tool(tool_name, cmd, env=None):
    popen_args = {
        'stdout': PIPE,
        'stderr': PIPE,
        'text': 'utf-8',
        'env': env
    }
    with Popen(cmd.split(), **popen_args) as process:
        for line in process.stdout:
            yield line

        retcode = process.wait(3)
        error = process.stderr.read().strip()

        if retcode != 0 or error:
            raise ToolExecutionError(tool_name, error)


def run_tool(tool_name, cmd, env=None):
    """ see _run_tool Generates line by line tool output, 
    throws ToolExecutionError if errcode != 0 or stderr is not empty
    this wrapper is for setup cache
    """
    if not context.args.cache:
        return _run_tool(tool_name, cmd, env)

    suffix = hex(crc32(cmd.encode()))
    cache_name = f'{tool_name}-{suffix}.cached'
    return context.cache.cachelines(cache_name, _run_tool(tool_name, cmd, env))


TrnaRecord = namedtuple('TrnaRecord', ['num', 'name', 'pos', 'nums', 'codon'])


def t_rna(input):
    aragorn = context.Tool.aragorn

    log.info(f'Running: ' + aragorn.name)
    cmd = f"aragorn -l -gc11 -w {input}"

    trna_output = run_tool(aragorn.name, cmd)

    re_coords = re.compile(r'(c)?\[-?(\d+),(\d+)\]')

    contig_id = None
    features = {}

    for line in map(str.strip, trna_output):
        if line.startswith('>'):
            contig_id = line[1:]
            features[contig_id] = []
            continue

        data = line.split()
        if len(data) != 5 or not re.match(r'\d+', data[0]):
            continue

        record = TrnaRecord(*data)

        if '?' in record.name:
            log.debug(f'tRNA {record.pos} is a pseudo/wacky gene - skipping')
            continue

        match = re_coords.match(record.pos)
        if not match:
            log.debug(f'Invalid position format {record}')
        revcom, start, end = match.groups()

        start = int(start)
        end = int(end)
        strand = -1 if revcom else +1

        if start > end:
            log.debug(f'tRNA {record.pos} has start > end - skipping.')
            continue

        start = max(start, 1)
        min(end, len(context.CONTIG[contig_id].seq))

        if abs(end-start) > MAX_TRNA_LEN:
            log.debug(f'tRNA {record.pos} is too big (>{MAX_TRNA_LEN}) - skipping.')
            continue

        tool = 'Aragorgn:' + aragorn.version

        ftype = 'tRNA'
        product = record.name + record.codon
        gene = {}

        if 'tmRNA' in record.name:
            ftype = 'tmRNA'
            product = 'transfer-messenger RNA, SsrA'
            gene = {'gene': 'ssrA'}

        qualifiers = {
            'product': product,
            'inference': f'COORDINATES:profile:{tool}'
        }
        qualifiers.update(gene)

        features[contig_id].append(SeqFeature(
            FeatureLocation(start, end, strand),
            id=f'{record.num} {product}',
            type=ftype,
            strand=strand,
            qualifiers=qualifiers
        ))

    log.info(f'Found {deeplen(features)} tRNAs')
    return features


def r_rna(input):
    barrnap = context.Tool.barrnap

    tool_out = run_tool(
        barrnap.name,
        f'barrnap --kingdom bac --threads 8 --quiet {input}',
        env=dict(environ, LC_ALL='C')
    )

    features = {}
    for record in GFF.parse(IOWrapper(map(str.strip, tool_out))):
        features[record.id] = record.features

    log.info(f'Found {deeplen(features)} rRNAs')

    # Q: do we need rnammer alternative rrna annotation?
    return features


def nc_rna(input):
    infernal = context.Tool.cmscan
    toolname = f"infernal:" + infernal.version

    icpu = context.CPU | 1
    dbsize = context.TOTAL_BP * 2 / 10e6
    cmdb = f'./prokka/db/cm/{context.KINGDOM}'

    cmd = f"cmscan -Z {dbsize} --cut_ga --rfam --nohmmonly --fmt 2 --cpu {icpu}" + \
          f" --tblout /dev/stdout -o /dev/null --noali {cmdb} {input}"

    features = {}
    for line in map(str.strip, run_tool('infernal', cmd)):
        if line.startswith('#'):
            continue

        data = line.split()
        if len(data) < 2 or not data[2].startswith('RF'):
            continue

        contig_id = data[3]
        if not contig_id in context.CONTIG:
            continue

        # Overlaps with a higher scoring match
        if len(data) > 19 and data[19] == '=':
            continue

        if contig_id not in features:
            features[contig_id] = []

        qualifiers = {
            'product': data[1],
            'inference': f'COORDINATES:profile:{toolname}',
            'accession': data[2],
            'note': ' '.join(data[26:]),
            'score': float(data[16])
        }

        start = min(map(int, [data[9], data[10]]))
        end = max(map(int, [data[9], data[10]]))
        strand = +1 if data[11] == '-' else +1

        features[contig_id].append(SeqFeature(
            FeatureLocation(start, end, strand),
            id=f'{data[0]} {data[1]}',
            type='misc_RNA',
            strand=strand,
            qualifiers=qualifiers,
        ))

    log.info(f'Found {deeplen(features)} ncRNAs')
    return features


def crisprs(input):
    minced = context.Tool.minced

    if not minced.have:
        log.warning('Skipping search for CRISPR repeats. Install minced to enable.')
        return

    tool_out = run_tool(minced.name, f'minced -gff {input}')

    features = {}
    for record in GFF.parse(IOWrapper(map(str.strip, tool_out))):
        features[record.id] = []
        for feature in record.features:
            n_repears = feature.qualifiers['score'][0]
            feature.qualifiers['note'] = f'CRISPR with {n_repears} repeat units'
            if 'rpt_family' not in feature.qualifiers:
                features.qualifiers['rpt_family'] = 'CRISPR'
            if 'rpt_type' not in feature.qualifiers:
                features.qualifiers['rpt_type'] = 'direct'

            features[record.id].append(feature)
            context.all_rna.append(feature)

    log.info(f'Found {deeplen(features)} CRISPR repeats')
    return features


def cds(input):
    prodigal = context.Tool.prodigal

def main(args):
    context.args = args
    log.setLevel(max(0, 20 - args.verboose*10))

    check_all_tools()
    context.Tool = Tool

    if args.cache:
        context.cache = Cache(join(args.output_dir, '.cache'))

    out = cleanup_input(args.input_file)

    trna_features = t_rna(out)
    context.add_features(trna_features)

    rrna_features = r_rna(out)
    context.add_features(rrna_features)

    ncrna_features = nc_rna(out)
    context.add_features(ncrna_features)

    crisprs_features = crisprs(out)
    context.add_features(crisprs_features)

    SeqIO.write(context.CONTIG.values(), 'out.gb', 'gb')
    


if __name__ == '__main__':
    main(parse_args())
