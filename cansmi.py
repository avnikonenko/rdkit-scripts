#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 26-07-2019
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2019
# license         : 
#==============================================================================

import sys
import argparse
from read_input import read_input
from multiprocessing import Pool, cpu_count
from rdkit import Chem


def calc(items):
    mol, mol_name = items
    # items is a tuple (mol, mol_name)
    Chem.AssignStereochemistryFrom3D(mol)
    smi = Chem.MolToSmiles(mol, isomericSmiles=True)
    return smi, mol_name


def main():

    parser = argparse.ArgumentParser(description='Conversion of input file to canonical SMILES with RDKit.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=False, default=None,
                        help='input file with compounds. Supported formats: SMILES (*.smi), '
                             'SDF (*.sdf, *.sdf.gz), Python pickled (*.pkl). '
                             'If omitted STDIN will be read as SMILES.')
    parser.add_argument('-o', '--output', metavar='output.smi', required=False, default=None,
                        help='output file with canonical SMILES. If omitted output will be redirected to STDOUT.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=1, type=int,
                        help='number of CPUs to use for computation.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = parser.parse_args()

    input_fname = None if args.input == '/dev/stdin' else args.input
    input_format = 'smi' if input_fname is None else None
    f = open(args.output, 'wt') if args.output is not None else sys.stdout
    try:
        if args.ncpu > 1:
            p = Pool(max(1, min(args.ncpu, cpu_count())))
            iterator = p.imap(calc, read_input(input_fname, input_format=input_format))
        else:
            iterator = (calc(line) for line in read_input(input_fname, input_format=input_format))
        for i, res in enumerate(iterator, 1):
            if res[0]:
                f.write('\t'.join(res) + '\n')
            if args.verbose and i % 10000 == 0:
                sys.stderr.write(f'\r{i} molecules passed')
    finally:
        if args.output is not None:
            f.close()


if __name__ == '__main__':
    main()
