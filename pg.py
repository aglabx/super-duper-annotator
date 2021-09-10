import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from BCBio import GFF

from sys import argv
from os import environ
from os.path import exists
from io import StringIO, TextIOWrapper

from subprocess import Popen
from subprocess import PIPE
from collections import namedtuple


FORMAT_FASTA = 'fasta'
MIN_CONTIG_LEN = 200
MAX_TRNA_LEN = 500


class ToolExecutionError(Exception):

    def __init__(self, name, error):
        self.name = name
        self.error = error
        super().__init__(f'{self.name} error: {self.error}')


class Context:
    def __init__(self):
        self.CONTIG = {}

    def add_features(self, features):
        for cid, fs in features.items():
            if cid in self.CONTIG:
                self.CONTIG[cid].features += fs


CONTEXT = Context()


def first(input):
    if not exists(input):
        raise FileExistsError(input)

    output = './pg/after_cleanup.fna'

    data = SeqIO.parse(input, format=FORMAT_FASTA)

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
                raise ValueError(
                    f'Error: duplicate id in sequence file {contig.id}')
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

    prepared = prepare_input(data)

    SeqIO.write(prepared, output, format=FORMAT_FASTA)
    print(f'File {output} written')

    return output


def run_tool(tool_name, cmd, env=None):
    """ Generates line by line tool output, 
    throws ToolExecutionError if errcode != 0 or stderr is not empty
    """
    with Popen(cmd.split(), text='utf-8', stdout=PIPE, stderr=PIPE, env=env) as process:
        for line in process.stdout:
            yield line.strip()

        retcode = process.wait(3)
        error = process.stderr.read().strip()

        if retcode != 0 or error:
            raise ToolExecutionError(tool_name, error)


def tRna(input):
    if not exists(input):
        raise FileExistsError(input)

    output = './pg/after_aragorn.fna'

    cmd = f"aragorn -l -gc11 -w {input}"
    print(f'Running: aragorn')
    trna_output = run_tool('aragorn', cmd)

    TrnaRecord = namedtuple(
        'TrnaRecord', ['num', 'name', 'pos', 'nums', 'codon'])

    re_coords = re.compile(r'(c)?\[-?(\d+),(\d+)\]')

    num_trna = 0
    contig_id = None
    features = {}
    for line in trna_output:
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

        num_trna += 1

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

    print(f'Found {num_trna} tRNAs')
    return features
    #   There is no point to write out new file since
    #   feature info will not be presented in fasta format
    #   we'll store it in memory for now

    # SeqIO.write(CONTEXT.CONTIG.values(), output, format=FORMAT_FASTA)
    # print(f'File {output} written')
    # return output


class IOWrapper:

    def __init__(self, obj):
        self._iter = obj

    def read(self):
        pass

    def readline(self):
        try:
            return next(self._iter)
        except StopIteration:
            return ""

    def close(self):
        pass


def rRna(input):
    if not exists(input):
        raise FileExistsError(input)

    output = './pg/after_barrnap.fna'

    tool_out = run_tool(
        'barrnap',
        f'barrnap --kingdom bac --threads 8 --quiet {input}',
        env=dict(environ, LC_ALL='C')
    )

    num_rrna = 0
    features = {}
    for record in GFF.parse(IOWrapper(tool_out)):
        features[record.id] = record.features
        num_rrna += len(record.features)

    print(f'Found {num_rrna} rRNAs')
    return features


out = first(argv[1])
trna_features = tRna(out)
rrna_features = rRna(out)

CONTEXT.add_features(trna_features)
CONTEXT.add_features(rrna_features)
