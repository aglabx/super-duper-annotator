import argparse as ap

from os.path import isfile


def file(arg):
    if isfile(arg):
        return arg
    raise ValueError(arg)


def parse_args():
    arp = ap.ArgumentParser(
        'sudadn',
        usage='use it wisely',
        description='annotates marvelously',
    )

    def flag(*name_or_flags, help=None):
        return arp.add_argument(*name_or_flags, action='store_true', help=help)

    arp.add_argument('input_file', type=file)
    arp.add_argument('--output-dir', '-o', default='sudan_output')

    flag('--cache', help='use cache to speed up multiple runs with same argss')
    arp.add_argument('--verboose', '-v', action='count', default=0)

    return arp.parse_args()


def main():
    a = parse_args()
    print(a)



if __name__ == '__main__':
    main()