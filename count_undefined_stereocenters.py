#!/usr/bin/env python3

__author__ = 'pavel'

import argparse
from gen_stereo_rdkit import get_unspec_double_bonds
from rdkit import Chem
from read_input import read_input


def main():

    parser = argparse.ArgumentParser(description='Count stereoisomers with RDKit.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=False, default=None,
                        help='input file in SDF or SMILES format. SMILES input should have no header, '
                             'the first column is SMILES string and the second column with ID is optional. '
                             'If omitted STDIN will be read as SMILES.')
    parser.add_argument('--no_tetrahedral', required=False, action='store_true', default=False,
                        help='do not count stereoisomers for unspecified tetrahedral centers.')
    parser.add_argument('--no_double_bond', required=False, action='store_true', default=False,
                        help='do not count stereoisomers for unspecified double bonds.')

    args = parser.parse_args()

    input_format = 'smi' if args.input is None else None
    for mol, mol_name in read_input(args.input, input_format=input_format):
        undef = 0
        if not args.no_tetrahedral:
            undef += sum(i[1] == '?' for i in Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        if not args.no_double_bond:
            undef += len(get_unspec_double_bonds(mol))
        print(f'{mol_name}\t{undef}')


if __name__ == '__main__':
    main()
