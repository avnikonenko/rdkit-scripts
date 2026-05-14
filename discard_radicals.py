#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import sys
from multiprocessing import Pool, cpu_count

from rdkit import Chem

from read_input import read_input


def process_mol(items):
    mol, mol_name = items
    radical_electrons = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    smi = Chem.MolToSmiles(mol, isomericSmiles=True)
    return f'{smi}\t{mol_name}\n', radical_electrons > 0


def filter_radicals(input_fname, output_fname, radicals_fname, ncpu, verbose):
    input_format = 'smi' if input_fname is None else None
    ncpu = max(1, min(ncpu, cpu_count()))

    fout = open(output_fname, 'wt') if output_fname is not None else sys.stdout
    frad = open(radicals_fname, 'wt') if radicals_fname is not None else None

    processed = 0
    kept = 0
    discarded = 0

    if ncpu > 1:
        pool = Pool(ncpu)
        iterator = pool.imap(process_mol, read_input(input_fname, input_format=input_format), chunksize=100)
    else:
        pool = None
        iterator = (process_mol(items) for items in read_input(input_fname, input_format=input_format))

    try:
        for processed, (line, is_radical) in enumerate(iterator, 1):
            if is_radical:
                discarded += 1
                if frad is not None:
                    frad.write(line)
            else:
                kept += 1
                fout.write(line)

            if verbose and processed % 1000 == 0:
                sys.stderr.write(f'\rProcessed: {processed}; kept: {kept}; discarded: {discarded}')
                sys.stderr.flush()
    finally:
        if pool is not None:
            pool.close()
            pool.join()
        if output_fname is not None:
            fout.close()
        if frad is not None:
            frad.close()

    if verbose:
        sys.stderr.write(f'\rProcessed: {processed}; kept: {kept}; discarded: {discarded}\n')
        sys.stderr.flush()


def main():
    parser = argparse.ArgumentParser(description='Discard structures containing radical electrons.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, default=None,
                        help='input SDF, SMI, SDF.GZ or PKL file. If omitted STDIN will be read as SMILES.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='output file with non-radical structures in SMILES format. '
                             'If omitted output will be redirected to STDOUT.')
    parser.add_argument('-d', '--discarded', metavar='FILENAME', required=False, default=None,
                        help='optional output file with discarded radical structures in SMILES format.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=1, type=int,
                        help='number of CPUs to use for computation.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    input_fname = None if args.input == '/dev/stdin' else args.input
    filter_radicals(input_fname, args.output, args.discarded, args.ncpu, args.verbose)


if __name__ == '__main__':
    main()
