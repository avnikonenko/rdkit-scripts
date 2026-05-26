#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Scaffolds import rdScaffoldNetwork
from read_input import read_input


def print_counts(scaffold_col, counts):
    print(f"{scaffold_col}\tnum_compounds")
    for scaffold, num_compounds in counts:
        print(f"{scaffold}\t{num_compounds}")


def murcko_smiles(mol):
    if mol is None:
        return None
    scaf = MurckoScaffold.GetScaffoldForMol(mol)
    if scaf is None or scaf.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(scaf, isomericSmiles=False)


def build_network(scaffold_smiles):
    mols = [Chem.MolFromSmiles(s) for s in scaffold_smiles if s]
    mols = [m for m in mols if m is not None]

    params = rdScaffoldNetwork.ScaffoldNetworkParams()
    params.includeGenericScaffolds = False
    params.includeScaffoldsWithAttachments = False
    params.includeScaffoldsWithoutAttachments = True
    params.includeNames = False

    net = rdScaffoldNetwork.CreateScaffoldNetwork(mols, params)
    return net


def network_parents(net):
    """
    Returns adjacency from child scaffold to parent scaffolds.

    RDKit scaffold network edges encode relationships between scaffold nodes.
    In most use cases, smaller/parent scaffolds are upstream of larger/child
    scaffolds, but we defensively orient edges by atom count.
    """
    nodes = list(net.nodes)
    parents = defaultdict(set)

    for edge in net.edges:
        a = edge.beginIdx
        b = edge.endIdx
        sa = nodes[a]
        sb = nodes[b]

        ma = Chem.MolFromSmiles(sa)
        mb = Chem.MolFromSmiles(sb)

        if ma is None or mb is None:
            continue

        na = ma.GetNumAtoms()
        nb = mb.GetNumAtoms()

        if na <= nb:
            parent, child = sa, sb
        else:
            parent, child = sb, sa

        parents[child].add(parent)

    return parents


def ancestors_by_level(scaffold, parents, max_level=5):
    """
    level 0: exact scaffold
    level 1: immediate parent(s)
    level 2: parents of parents
    etc.

    A scaffold can have multiple parents in a scaffold network, so this returns
    a tuple of ancestor scaffold strings per level.
    """
    levels = {0: {scaffold}}
    current = {scaffold}

    for level in range(1, max_level + 1):
        nxt = set()
        for s in current:
            nxt.update(parents.get(s, []))
        if not nxt:
            break
        levels[level] = nxt
        current = nxt

    return levels


def main():
    parser = argparse.ArgumentParser(description='Cluster input molecules by Murcko scaffold network levels.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, default=None, type=str,
                        help='input file in SDF, SDF.GZ, SMI, or PKL formats. If omitted STDIN '
                             'will be read as SMILES.')
    parser.add_argument('-o', '--output', metavar='PREFIX', required=True, type=str,
                        help='output file prefix. Assignments will be saved to PREFIX_assignments.csv and '
                             'cluster tables to PREFIX_largest_clusters_level_N.csv.')
    parser.add_argument('-m', '--max_level', '--max-level', dest='max_level', metavar='INTEGER', type=int, default=5,
                        help='maximum number of scaffold ancestor levels to report. Default: 5.')
    parser.add_argument('-t', '--top', metavar='INTEGER', type=int, default=20,
                        help='number of largest clusters to print for each scaffold level. Default: 20.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    if args.max_level < 0:
        raise SystemExit("Maximum scaffold level should be zero or greater.")
    if args.top < 1:
        raise SystemExit("Number of clusters to print should be greater than zero.")

    molecules = []
    input_format = 'smi' if args.input is None else None
    for i, (mol, name) in enumerate(read_input(args.input, input_format=input_format), 1):
        molecules.append({
            'Name': name,
            'SMILES': Chem.MolToSmiles(mol, isomericSmiles=True),
            'murcko_scaffold': murcko_smiles(mol),
        })
        if args.verbose and i % 1000 == 0:
            sys.stderr.write(f'\r{i} molecules passed')
            sys.stderr.flush()

    if args.verbose and molecules:
        sys.stderr.write(f'\r{len(molecules)} molecules passed\n')

    if not molecules:
        raise SystemExit("No molecules found.")

    valid = [row for row in molecules if row["murcko_scaffold"] is not None]

    if not valid:
        raise SystemExit("No valid Murcko scaffolds found.")

    unique_scaffolds = sorted({row["murcko_scaffold"] for row in valid})

    net = build_network(unique_scaffolds)
    parents = network_parents(net)

    # Map each exact scaffold to its ancestors at each level.
    assignments = []

    for row in valid:
        name = row["Name"]
        smi = row["SMILES"]
        exact = row["murcko_scaffold"]

        levels = ancestors_by_level(exact, parents, max_level=args.max_level)

        out = {
            "Name": name,
            "SMILES": smi,
            "murcko_scaffold": exact,
        }

        for level in range(args.max_level + 1):
            ancestors = sorted(levels.get(level, []))
            # Multiple ancestors are possible in a scaffold network.
            # Join them so the assignment is explicit.
            out[f"scaffold_level_{level}"] = ".".join(ancestors) if ancestors else ""

        assignments.append(out)

    assignment_fields = ["Name", "SMILES", "murcko_scaffold"]
    assignment_fields.extend(f"scaffold_level_{level}" for level in range(args.max_level + 1))
    with open(f"{args.output}_assignments.csv", "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=assignment_fields)
        writer.writeheader()
        writer.writerows(assignments)

    # Count largest clusters at each level.
    for level in range(args.max_level + 1):
        col = f"scaffold_level_{level}"

        scaffold_names = defaultdict(set)
        for assignment in assignments:
            for scaffold in assignment[col].split("."):
                if scaffold:
                    scaffold_names[scaffold].add(assignment["Name"])
        counts = sorted(
            ((scaffold, len(names)) for scaffold, names in scaffold_names.items()),
            key=lambda item: (-item[1], item[0])
        )

        with open(f"{args.output}_largest_clusters_level_{level}.csv", "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow([col, "num_compounds"])
            writer.writerows(counts)

        print(f"\nLargest clusters at scaffold level {level}")
        print_counts(col, counts[:args.top])

    print(f"\nWrote: {args.output}_assignments.csv")
    print(f"Wrote largest-cluster tables for levels 0..{args.max_level}")


if __name__ == "__main__":
    main()
