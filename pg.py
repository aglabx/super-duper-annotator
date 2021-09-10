from contextlib import contextmanager
import re

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

from hashlib import sha256
from subprocess import PIPE
from subprocess import Popen
from collections import namedtuple

from cache import Cache
from utilities import IOWrapper
from utilities import ToolExecutionError


FORMAT_FASTA = 'fasta'
MIN_CONTIG_LEN = 200
MAX_TRNA_LEN = 500


cache = Cache('temp')


def deeplen(iterable):
    return sum(map(len, iterable.values()))


class Context:
    def __init__(self):
        self.CONTIG = {}
        self.TOTAL_BP = -1
        self.CPU = 8
        self.KINGDOM = 'Bacteria'

    def add_features(self, features):
        """accepts features in a dict {contig_id: [features]}
        appends given features to contig with associated id
        """
        for cid, fs in features.items():
            if cid in self.CONTIG:
                self.CONTIG[cid].features += fs


CONTEXT = Context()


def prepare_input(seq_iterator):
    ncontig = 0
    total_bp = 0
    ids = set()

    for contig in seq_iterator:
        if len(contig.seq) < MIN_CONTIG_LEN:
            print(f'Skipping short contig {contig.id}')

        ncontig += 1

        if '|' in contig.id:
            old_id = contig.id
            contig.id = old_id.replace('|', '_')
            print(f'Changing illegal name {old_id} to {contig.id}')

        if contig.id in ids:
            raise ValueError(f'Error: duplicate id in sequence file {contig.id}')
        ids.add(contig.id)

        seq = str(contig.seq.upper())
        seq, n = re.subn(r'[*-]', '', seq)
        if n > 0:
            print(f'Removed {n} gaps/pads')
        seq, n = re.subn(r'[^ACTG]', 'N', seq)
        if n > 0:
            print(f'Replaced {n} wacky IUPAC with N')

        # todo: contig id rename

        contig.seq = Seq(seq)
        contig.description = ''
        total_bp += len(contig.seq)

        CONTEXT.CONTIG[contig.id] = contig
        yield contig

    if ncontig <= 0:
        raise ValueError(
            f'FASTA file {input} contains no suitable sequence entries')
    print(f'Total {total_bp} BP read')

    CONTEXT.TOTAL_BP = total_bp


def first(input):
    if not exists(input):
        raise FileExistsError(input)

    output = './pg/after_cleanup.fna'

    data = SeqIO.parse(input, format=FORMAT_FASTA)
    prepared = prepare_input(data)

    SeqIO.write(prepared, output, format=FORMAT_FASTA)
    print(f'File {output} written')

    return output


def _run_tool(tool_name, cmd, env=None):
    with Popen(cmd.split(), text='utf-8', stdout=PIPE, stderr=PIPE, env=env) as process:
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
    suffix = sha256(cmd.encode()).hexdigest()
    cache_name = f'{tool_name}{suffix}.cached'
    yield from cache.cachelines(cache_name, _run_tool(tool_name, cmd, env))


TrnaRecord = namedtuple('TrnaRecord', ['num', 'name', 'pos', 'nums', 'codon'])


def t_rna(input):
    if not exists(input):
        raise FileExistsError(input)

    print(f'Running: aragorn')
    cmd = f"aragorn -l -gc11 -w {input}"

    trna_output = run_tool('aragorn', cmd)

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
            print(f'tRNA {record.pos} is a pseudo/wacky gene - skipping')
            continue

        match = re_coords.match(record.pos)
        if not match:
            print(f'Invalid position format {record}')
        revcom, start, end = match.groups()

        start = int(start)
        end = int(end)
        strand = -1 if revcom else +1

        if start > end:
            print(f'tRNA {record.pos} has start > end - skipping.')
            continue

        start = max(start, 1)
        min(end, len(CONTEXT.CONTIG[contig_id].seq))

        if abs(end-start) > MAX_TRNA_LEN:
            print(f'tRNA {record.pos} is too big (>{MAX_TRNA_LEN}) - skipping.')
            continue

        tool = 'Aragorgn'  # todo: + tool.aragorn.version

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

    print(f'Found {deeplen(features)} tRNAs')
    return features


def r_rna(input):
    if not exists(input):
        raise FileExistsError(input)

    tool_out = run_tool(
        'barrnap',
        f'barrnap --kingdom bac --threads 8 --quiet {input}',
        env=dict(environ, LC_ALL='C')
    )

    features = {}
    for record in GFF.parse(IOWrapper(map(str.strip, tool_out))):
        features[record.id] = record.features

    print(f'Found {deeplen(features)} rRNAs')

    # Q: do we need rnammer alternative rrna annotation?
    return features


def nc_rna(input):
    if not exists(input):
        raise FileExistsError(input)

    icpu = CONTEXT.CPU | 1
    dbsize = CONTEXT.TOTAL_BP * 2 / 10e6
    cmdb = f'./prokka/db/cm/{CONTEXT.KINGDOM}'
    tool = f"Infernal"  # todo: add version

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
        if not contig_id in CONTEXT.CONTIG:
            continue

        # Overlaps with a higher scoring match
        if len(data) > 19 and data[19] == '=':
            continue

        if contig_id not in features:
            features[contig_id] = []
            
        qualifiers = {
            'product': data[1],
            'inference': f'COORDINATES:profile:{tool}',
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

    print(f'Found {deeplen(features)} ncRNAs')
    return features



def main():
    out = first(argv[1])
    trna_features = t_rna(out)
    CONTEXT.add_features(trna_features)

    rrna_features = r_rna(out)
    CONTEXT.add_features(rrna_features)

    ncrna_features = nc_rna(out)
    CONTEXT.add_features(ncrna_features)

    print(f"Total {deeplen(trna_features) + deeplen(rrna_features)} tRNA + rRNA features")
    print('we\'re done for now, come back later')

if __name__ == '__main__':
    main()
