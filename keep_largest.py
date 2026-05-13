#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import sys
from multiprocessing import Pool
from rdkit import Chem
from read_input import read_input


def get_largest(mol, mol_name):
    output_frag = None
    try:
        max_hac = 0
        if mol is not None:
            frags = Chem.GetMolFrags(mol, asMols=True)
            for frag in frags:
                if frag.GetNumHeavyAtoms() > max_hac:
                    output_frag = frag
                    max_hac = frag.GetNumHeavyAtoms()
    except Exception as e:
        sys.stderr.write(f'Error: {e}\n')
        sys.stderr.write(f'Molecule {mol_name} is skipped\n')
    return output_frag, mol_name


def get_largest_mp(items):
    return get_largest(*items)


def main_params(in_fname, out_fname, ncpu, verbose):

    pool = Pool(ncpu)

    input_format = 'smi' if in_fname is None else None

    fo = open(out_fname, 'wt') if out_fname is not None else sys.stdout
    try:
        for i, (mol, mol_name) in enumerate(pool.imap(get_largest_mp, read_input(in_fname, input_format=input_format, sanitize=True)), 1):
            if mol:
                fo.write(f'{Chem.MolToSmiles(mol, isomericSmiles=True)}\t{mol_name}\n')
            if verbose and i % 1000 == 0:
                sys.stderr.write(f'\r{i} records were processed')
        if verbose:
            sys.stderr.write('\n')
    finally:
        if out_fname is not None:
            fo.close()


def main():
    parser = argparse.ArgumentParser(description='Keep the largest fragment by the number of heavy atoms '
                                                 'in each compound record. If components have the same number a random '
                                                 'one will be selected.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, default=None,
                        help='input file in SDF or SMILES format. SMILES input should have no header, '
                             'the first column is SMILES string and the second column with ID is optional. '
                             'If omitted STDIN will be read as SMILES.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='output file in SMILES format. If omitted output will be redirected to STDOUT.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=1, type=int,
                        help='Number of CPU cores to use. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = parser.parse_args()

    if args.input == "/dev/stdin":
        args.input = None

    main_params(args.input, args.output, args.ncpu, args.verbose)


if __name__ == '__main__':
    main()
