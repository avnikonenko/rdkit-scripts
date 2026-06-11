#!/usr/bin/env python3

import argparse
import sys
from read_input import read_input
from rdkit import Chem


def remove_isotopes_from_mol(mol):
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    mol = Chem.RemoveHs(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def remove_isotopes(input_fname, output_fname, verbose):
    input_format = 'smi' if input_fname is None else None
    f = open(output_fname, 'wt') if output_fname is not None else sys.stdout
    i = 0
    try:
        for i, (mol, mol_name) in enumerate(read_input(input_fname, input_format=input_format), 1):
            if mol:
                mol = remove_isotopes_from_mol(mol)
                f.write(f'{Chem.MolToSmiles(mol, isomericSmiles=True)}\t{mol_name}\n')
            if verbose and i % 1000 == 0:
                sys.stderr.write(f'\rProcessed {i} molecules')
    finally:
        if output_fname is not None:
            f.close()
    if verbose:
        sys.stderr.write(f'\rProcessed {i} molecules\n')


def main():
    parser = argparse.ArgumentParser(description='Remove isotope labels and removable explicit hydrogens from input molecules.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, default=None, type=str,
                        help='input SDF or SMILES file. If omitted STDIN will be read as SMILES.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None, type=str,
                        help='output SMILES file. If omitted output will be redirected to STDOUT.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    remove_isotopes(args.input, args.output, args.verbose)


if __name__ == '__main__':
    main()
