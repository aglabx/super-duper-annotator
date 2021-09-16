import re

from os import environ
from os import popen
from os.path import join
from os.path import exists
from os.path import abspath
from logging import getLogger
from distutils.spawn import find_executable

log = getLogger('sudan.tools')
BIDEC = r'(\d+\.\d+)'


PARALLEL = 'parallel'
ARAGORN = 'aragorn'
RNAMMER = 'rnammer'
BARRNAP = 'barrnap'
PRODIGAL = 'prodigal'
SIGNALP = 'signalp'
MINCED = 'minced'
CMSCAN = 'cmscan'
CMPRESS = 'cmpress'
HMMSCAN = 'hmmscan'
HMMPRESS = 'hmmpress'
BLASTP = 'blastp'
MAKEBLASTDB = 'makeblastdb'
TBL2ASN = 'tbl2asn'


class _Tool:

    def __init__(self, name, getver, regexp, minver, needed):
        self.name = name
        self.needed = needed
        self.getver = getver
        self.regexp = regexp
        self.minver = minver

    @classmethod
    def standard_tool(cls, name):
        return cls(name, None, None, None, needed=True)


class Tool:
    tools_checked = False

    def __init__(self, name, version, have):
        self.name = name
        self.version = version
        self.have = have


__tools = [
    # just the standard unix tools we need
    # yes, we need this before we can test versions :-/
    _Tool.standard_tool('grep'),
    _Tool.standard_tool('egrep'),
    _Tool.standard_tool('sed'),
    _Tool.standard_tool('find'),
    _Tool.standard_tool('java'),
    # for the new --proteins option ability to take .gbk or .embl files
    _Tool.standard_tool('prokka-genbank_to_fasta_db'),

    _Tool(
        PARALLEL,
        "parallel --version | grep -E 'parallel 2[0-9]{7}$'",
        re.compile(r'parallel (\d+)'),
        "20130422",
        True
    ),
    _Tool(
        ARAGORN,
        "aragorn -h 2>&1 | grep -i '^ARAGORN v'",
        re.compile(rf'{BIDEC}'),
        "1.2",
        True
    ),
    _Tool(
        RNAMMER,
        "rnammer -V 2>&1 | grep -i 'rnammer [0-9]'",
        re.compile(rf'{BIDEC}'),
        "1.2",
        False
    ),
    _Tool(
        BARRNAP,
        "LC_ALL=C barrnap --version 2>&1",
        re.compile(rf'{BIDEC}'),
        "0.4",
        False
    ),
    _Tool(
        PRODIGAL,
        "prodigal -v 2>&1 | grep -i '^Prodigal V'",
        re.compile(rf'{BIDEC}'),
        "2.6",
        True
    ),
    _Tool(
        SIGNALP,
        # this is so long-winded as -v changed meaning (3.0=version, 4.0=verbose !?)
        "if [ \"`signalp -version 2>&1 | grep -Eo '[0-9]+\.[0-9]+'`\" != \"\" ]; then echo `signalp -version 2>&1 | grep -Eo '[0-9]+\.[0-9]+'`; else signalp -v < /dev/null 2>&1 | egrep ',|# SignalP' | sed 's/^# SignalP-//'; fi",
        re.compile(rf'^{BIDEC}'),
        "3.0",
        False
    ),
    _Tool(
        MINCED,
        "minced --version | sed -n '1p'",
        re.compile(rf'minced\s+\d+\.{BIDEC}'),
        "2.0",
        False  # really 0.3 as we skip first number.
    ),
    _Tool(
        CMSCAN,
        "cmscan -h | grep '^# INFERNAL'",
        re.compile(rf'INFERNAL\s+{BIDEC}'),
        "1.1",
        False
    ),
    _Tool(
        CMPRESS,
        "cmpress -h | grep '^# INFERNAL'",
        re.compile(rf'INFERNAL\s+{BIDEC}'),
        "1.1",
        False
    ),
    _Tool(
        HMMSCAN,
        "hmmscan -h | grep '^# HMMER'",
        re.compile(rf'HMMER\s+{BIDEC}'),
        "3.1",
        True
    ),
    _Tool(
        HMMPRESS,
        "hmmpress -h | grep '^# HMMER'",
        re.compile(rf'HMMER\s+{BIDEC}'),
        "3.1",
        True
    ),
    _Tool(
        BLASTP,
        "blastp -version",
        re.compile(rf'blastp:\s+{BIDEC}'),
        "2.8",
        True
    ),
    _Tool(
        MAKEBLASTDB,
        "makeblastdb -version",
        re.compile(rf'makeblastdb:\s+{BIDEC}'),
        "2.8",
        False
    ),
    _Tool(
        TBL2ASN,
        "tbl2asn - | grep '^tbl2asn'",
        re.compile(rf'tbl2asn\s+{BIDEC}'),
        "24.3",
        True
    )
]


def parse_version(s_version):
    return tuple(map(int, s_version.split('.')))


def check_tool(tool):
    executable = find_executable(tool.name)
    if not executable and tool.needed:
        raise Exception(
            f'Required executable {tool.name} was not found in PATH')
    if not executable:
        log.debug(f'Tool {tool.name} not found. skipped.')
        return False, 0
    if not tool.getver:
        return True, 1

    ver_output = popen(tool.getver).read().strip()

    match = tool.regexp.search(ver_output)
    if not match:
        raise Exception(f'Cound not determine version of {tool.name}. ' +
                        f'Please install version {tool.minver} or higher')

    s_version = match.group(1)
    version = parse_version(s_version)
    min_version = parse_version(tool.minver)

    if version < min_version:
        raise Exception(
            f'Tool {tool.name} version required {min_version} or higher')

    log.debug(f'Found {tool.name}:{s_version}')

    return True, s_version


def check_all_tools():
    """
    Initializes Tool class with static fields which contains info about tool
    """
    log.info('Checking tools')
    for tool in __tools:
        have, version = check_tool(tool)
        tool = Tool(tool.name, version, have)
        setattr(Tool, tool.name, tool)

    Tool.tools_checked = True


def main():
    list(check_all_tools())


if __name__ == '__main__':
    main()
