#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import sys
import warnings
from collections import deque

warnings.filterwarnings(
    'ignore',
    message='Deprecated in version 2.8.0.*',
    category=DeprecationWarning,
)

import prolif as plf
from rdkit import Chem


def read_contacts(contact_args, contacts_file):
    contacts = set()
    for contact in contact_args or []:
        contact = contact.strip().lower()
        if contact:
            contacts.add(contact)

    if contacts_file is not None:
        with open(contacts_file) as f:
            for line in f:
                line = line.split('#', 1)[0].strip().lower()
                if line:
                    contacts.add(line)

    return contacts


def residues_from_contacts(contacts):
    residues = set()
    for contact in contacts:
        if '.' not in contact:
            continue
        residue = contact.rsplit('.', 1)[0]
        if residue:
            residues.add(residue.upper())
    return sorted(residues)


def load_protein(protein_fname, parser):
    parent_mol = Chem.MolFromPDBFile(
        protein_fname,
        removeHs=False,
        sanitize=True,
    )
    if parent_mol is None:
        raise ValueError('RDKit failed to read protein PDB: %s' % protein_fname)

    if parser == 'rdkit':
        return parent_mol, plf.Molecule(parent_mol)

    if parser == 'mda':
        import MDAnalysis as mda

        universe = mda.Universe(protein_fname)
        return parent_mol, plf.Molecule.from_mda(universe)

    raise ValueError('Unknown protein parser: %s' % parser)


def get_mol_id(mol, entry_index, field_name):
    if field_name:
        if mol.HasProp(field_name):
            return mol.GetProp(field_name)
        return ''
    if mol.HasProp('_Name'):
        name = mol.GetProp('_Name').strip()
        if name:
            return name
    return str(entry_index)


def contact_key(protein_residue, interaction_name):
    return ('%s.%s' % (protein_residue, interaction_name)).lower()


def ligand_indices_from_metadata(metadata):
    indices = metadata.get('parent_indices', {}).get('ligand')
    if indices is None:
        indices = metadata.get('indices', {}).get('ligand', ())
    return indices


def protein_indices_from_metadata(metadata):
    indices = metadata.get('parent_indices', {}).get('protein')
    if indices is None:
        indices = metadata.get('indices', {}).get('protein', ())
    return tuple(indices)


def protein_hbond_parent_matches(protein_mol):
    hbacceptor = plf.interactions.HBAcceptor()
    return {
        # ProLIF HBAcceptor means ligand acceptor, protein donor.
        'HBAcceptor': set(protein_mol.GetSubstructMatches(hbacceptor.prot_pattern)),
        # ProLIF HBDonor means ligand donor, protein acceptor.
        'HBDonor': set(protein_mol.GetSubstructMatches(hbacceptor.lig_pattern)),
    }


def valid_parent_protein_hbond(interaction_name, metadata, parent_matches):
    matches = parent_matches.get(interaction_name)
    if matches is None:
        return True
    return protein_indices_from_metadata(metadata) in matches


def heavy_atom_index(mol, atom_idx):
    if atom_idx < 0 or atom_idx >= mol.GetNumAtoms():
        return None
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetAtomicNum() != 1:
        return atom_idx
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() != 1:
            return neighbor.GetIdx()
    return None


def collect_matching_interactions(mol, ifp, requested_contacts, parent_matches):
    found_contacts = set()
    interacting_atoms = set()

    for residue_pair, interactions in ifp.items():
        protein_residue = str(residue_pair[1])
        for interaction_name, occurrences in interactions.items():
            key = contact_key(protein_residue, interaction_name)
            if key not in requested_contacts:
                continue

            found_occurrence = False
            for metadata in occurrences:
                if not valid_parent_protein_hbond(
                    interaction_name,
                    metadata,
                    parent_matches,
                ):
                    continue
                found_occurrence = True
                for atom_idx in ligand_indices_from_metadata(metadata):
                    heavy_idx = heavy_atom_index(mol, atom_idx)
                    if heavy_idx is not None:
                        interacting_atoms.add(heavy_idx)
            if found_occurrence:
                found_contacts.add(key)

    return found_contacts, interacting_atoms


def ring_systems(mol):
    systems = []
    for ring in mol.GetRingInfo().AtomRings():
        ring = set(ring)
        merged = []
        for i, system in enumerate(systems):
            if ring & system:
                ring.update(system)
                merged.append(i)
        for i in reversed(merged):
            systems.pop(i)
        systems.append(ring)
    return systems


def connected_components(mol, atoms):
    atoms = set(atoms)
    components = []

    while atoms:
        start = atoms.pop()
        component = {start}
        queue = deque([start])

        while queue:
            atom_idx = queue.popleft()
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in atoms:
                    atoms.remove(neighbor_idx)
                    component.add(neighbor_idx)
                    queue.append(neighbor_idx)

        components.append(component)

    return components


def shortest_path_between_sets(mol, starts, targets, allowed_atoms):
    starts = set(starts) & allowed_atoms
    targets = set(targets) & allowed_atoms
    if not starts or not targets:
        return None

    parents = {}
    queue = deque()
    for atom_idx in starts:
        parents[atom_idx] = None
        queue.append(atom_idx)

    while queue:
        atom_idx = queue.popleft()
        if atom_idx in targets:
            path = []
            while atom_idx is not None:
                path.append(atom_idx)
                atom_idx = parents[atom_idx]
            return list(reversed(path))

        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in allowed_atoms or neighbor_idx in parents:
                continue
            parents[neighbor_idx] = atom_idx
            queue.append(neighbor_idx)

    return None


def add_connector_atoms(mol, scaffold_atoms, allowed_atoms):
    scaffold_atoms = set(scaffold_atoms)

    while True:
        components = connected_components(mol, scaffold_atoms)
        if len(components) <= 1:
            break

        best_path = None
        for i, component_a in enumerate(components[:-1]):
            for component_b in components[i + 1:]:
                path = shortest_path_between_sets(
                    mol, component_a, component_b, allowed_atoms
                )
                if path is None:
                    continue
                if best_path is None or len(path) < len(best_path):
                    best_path = path

        if best_path is None:
            break

        before = len(scaffold_atoms)
        scaffold_atoms.update(best_path)
        if len(scaffold_atoms) == before:
            break

    return scaffold_atoms


def expand_touched_ring_systems(mol, atoms):
    atoms = set(atoms)
    for system in ring_systems(mol):
        if system & atoms:
            atoms.update(system)
    return atoms


def expand_double_bond_neighbors(mol, atoms):
    atoms = set(atoms)
    expanded = set(atoms)

    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in atoms:
            expanded.add(end)
        if end in atoms:
            expanded.add(begin)

    return expanded


def extract_contact_scaffold_atoms(mol, interacting_atoms):
    interacting_atoms = set(interacting_atoms)
    if not interacting_atoms:
        return set()

    heavy_atoms = {
        atom.GetIdx() for atom in mol.GetAtoms()
        if atom.GetAtomicNum() != 1
    }

    scaffold_atoms = set(interacting_atoms) & heavy_atoms
    interacting_atoms = sorted(scaffold_atoms)

    for i, start in enumerate(interacting_atoms[:-1]):
        for end in interacting_atoms[i + 1:]:
            path = shortest_path_between_sets(
                mol,
                {start},
                {end},
                heavy_atoms,
            )
            if path is not None:
                scaffold_atoms.update(path)

    scaffold_atoms = expand_touched_ring_systems(mol, scaffold_atoms)
    scaffold_atoms = expand_double_bond_neighbors(mol, scaffold_atoms)
    return expand_touched_ring_systems(mol, scaffold_atoms)



def scaffold_smiles(mol, scaffold_atoms, isomeric_smiles):
    if not scaffold_atoms:
        return ''
    return Chem.MolFragmentToSmiles(
        mol,
        atomsToUse=sorted(scaffold_atoms),
        isomericSmiles=isomeric_smiles,
        canonical=True,
    )


def write_results(ligand_fname, protein_mol, protein, requested_contacts, residue_filter,
                  output_fname, field_name, sep, isomeric_smiles, verbose):
    supplier = Chem.SDMolSupplier(ligand_fname, removeHs=False, sanitize=True)
    fingerprint = plf.Fingerprint(interactions='all', count=True)
    parent_matches = protein_hbond_parent_matches(protein_mol)
    out = open(output_fname, 'wt') if output_fname is not None else sys.stdout

    try:
        out.write(sep.join([
            'entry_index',
            'molecule_id',
            'scaffold_smiles',
            'found_contacts_count',
            'found_contacts',
        ]) + '\n')

        for entry_index, mol in enumerate(supplier, start=1):
            if mol is None:
                sys.stderr.write(
                    'Entry %i cannot be read from SDF and was skipped\n'
                    % entry_index
                )
                continue

            mol_id = get_mol_id(mol, entry_index, field_name)

            try:
                ligand = plf.Molecule(mol)
                ifp = fingerprint.generate(
                    ligand,
                    protein,
                    residues=residue_filter or None,
                    metadata=True,
                )
            except Exception as exc:
                sys.stderr.write(
                    'Entry %i (%s) failed in ProLIF and was skipped: %s\n'
                    % (entry_index, mol_id, exc)
                )
                continue

            found_contacts, interacting_atoms = collect_matching_interactions(
                mol, ifp, requested_contacts, parent_matches
            )
            if not found_contacts:
                continue

            scaffold_atoms = extract_contact_scaffold_atoms(
                mol, interacting_atoms
            )
            smi = scaffold_smiles(mol, scaffold_atoms, isomeric_smiles)
            out.write(sep.join([
                str(entry_index),
                mol_id,
                smi,
                str(len(found_contacts)),
                ';'.join(sorted(found_contacts)),
            ]) + '\n')

            if verbose and entry_index % 1000 == 0:
                sys.stderr.write('\r%i SDF entries processed' % entry_index)
                sys.stderr.flush()
    finally:
        if output_fname is not None:
            out.close()


def main():
    parser = argparse.ArgumentParser(
        description='Extract contact-anchored ligand scaffolds from docked SDF poses.'
    )
    parser.add_argument('-i', '--input', metavar='poses.sdf', required=True,
                        help='input SDF file with ligand poses.')
    parser.add_argument('-p', '--protein', metavar='protein.pdb', required=True,
                        help='protein PDB file with all atoms.')
    parser.add_argument('-o', '--output', metavar='output.tsv', default=None,
                        help='output TSV file. If omitted, STDOUT is used.')
    parser.add_argument('-c', '--contacts', metavar='CONTACT', nargs='*',
                        default=[],
                        help='requested contacts, e.g. leu83.a.hbdonor.')
    parser.add_argument('--contacts-file', metavar='contacts.txt', default=None,
                        help='optional file with one requested contact per line.')
    parser.add_argument('positional_contacts', metavar='CONTACT', nargs='*',
                        help='requested contacts, e.g. leu83.a.hbdonor.')
    parser.add_argument('-f', '--field-name', metavar='FIELD_NAME',
                        default=None,
                        help='SDF field containing molecule id. Default: mol title.')
    parser.add_argument('--protein-parser', choices=['rdkit', 'mda'],
                        default='rdkit',
                        help='protein PDB parser. Default: rdkit.')
    parser.add_argument('-s', '--sep', metavar='CHAR', default='\t',
                        help='output separator. Default: tab.')
    parser.add_argument('--isomeric-smiles', action='store_true',
                        default=False,
                        help='include stereochemistry in scaffold SMILES.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = parser.parse_args()
    contacts = read_contacts(
        args.contacts + args.positional_contacts,
        args.contacts_file,
    )
    if not contacts:
        parser.error('at least one contact must be given via --contacts or --contacts-file')

    protein_mol, protein = load_protein(args.protein, args.protein_parser)
    residue_filter = residues_from_contacts(contacts)

    write_results(
        args.input,
        protein_mol,
        protein,
        contacts,
        residue_filter,
        args.output,
        args.field_name,
        args.sep,
        args.isomeric_smiles,
        args.verbose,
    )


if __name__ == '__main__':
    main()
