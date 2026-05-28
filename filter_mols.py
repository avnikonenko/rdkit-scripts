#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse


def filter_mols(input_fname, output_fname, names_fname, order):

    with open(names_fname)as f:
        names = []
        names_set = set()
        for line in f:
            name = line.strip()
            if name not in names_set:
                names.append(name)
                names_set.add(name)

    with open(input_fname) as f_in, open(output_fname, 'wt') as f_out:
        molstr = []
        mols_by_name = {} if order else None
        for line in f_in:
            molstr.append(line)
            if line.strip() == '$$$$':
                name = molstr[0].strip()
                if name in names_set:
                    if order:
                        mols_by_name.setdefault(name, []).append(molstr)
                    else:
                        f_out.writelines(molstr)
                molstr = []

        if order:
            for name in names:
                for molstr in mols_by_name.get(name, []):
                    f_out.writelines(molstr)


def main():
    parser = argparse.ArgumentParser(description='Filter input molecules by names.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output SDF file.')
    parser.add_argument('-n', '--names', metavar='FILENAME', required=True, type=str,
                        help='text file with molecule names to keep.')
    parser.add_argument('--ordered', action='store_true',
                        help='return molecules in the order specified in the names file. By default input order is '
                             'preserved.')
    args = parser.parse_args()

    filter_mols(args.input, args.output, args.names, args.order)


if __name__ == '__main__':
    main()
