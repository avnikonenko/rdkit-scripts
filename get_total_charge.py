#!/usr/bin/env python3
import argparse
import sys
from rdkit.Chem import rdmolops
from read_input import read_input


def main():
    parser = argparse.ArgumentParser(description='''Returns the formal charge for the molecule using RDKiT''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, default=None,
                        help='input file with structures. If omitted STDIN will be read as SMILES.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='Output text file. If omitted output will be in stdout.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    input_format = 'smi' if args.input is None else None
    with open(args.output, 'wt') as fout:
        i = 0
        for i, (mol, mol_name) in enumerate(read_input(args.input, input_format=input_format), 1):
            charge = rdmolops.GetFormalCharge(mol)
            fout.write('\t'.join([mol_name, str(charge)]) + '\n')
            if args.verbose and i % 1000 == 0:
                sys.stderr.write(f'\r{i} molecules were processed')
        sys.stderr.write(f'\r{i} molecules were processed\n')


if __name__ == '__main__':
    main()
