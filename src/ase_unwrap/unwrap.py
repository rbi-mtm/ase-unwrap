"""Unwrap structure across periodic boundaries"""

import sys

import numpy as np

from scipy import sparse

import ase.io

from ase.neighborlist import natural_cutoffs
from ase.neighborlist import NeighborList

def reorder_mol(neighbor_list, mol_indices, center_atom=None, reordered=None):
    """Reorder a list of atom indices starting from the ``center_atom``

    The reordering is done recursively. The ``center_atom`` is the
        first one added to the list. Then, one of its neighbors is added.
        Then, one of the neighbors of the neighbor is added...

    fix_mol() then iterates over the reordered list and translates the atoms
        for the molecule to be unwrapped.

    Parameters
    ----------
    neigbor_list : ase.neighborlist.Neighborlist
        ASE neighborlist
    mol_indices : list of int
        List of atomic indices defining the molecule
    center_atom : None | int
        Atom index from which to start recursively finding closest neighbors;
            if None, mol_indices[0] is chosen

    Returns
    -------
    reordered : list of int
        Reordered list of atomic indices

    """

    if not reordered:
        reordered = []

    if not center_atom:
        center_atom = mol_indices[0]

    if center_atom not in reordered:
        reordered.append(center_atom)

    nb_indices, _ = neighbor_list.get_neighbors(center_atom)

    for nb_index in nb_indices:

        if nb_index not in reordered:

            reordered = reorder_mol(neighbor_list,
                                    mol_indices,
                                    center_atom=nb_index,
                                    reordered=reordered)

        else:
            pass

    return reordered

def get_center_atom(atoms, mol_indices):
    """Get the atom closest to the center of the box for nicer unwrapping

    Parameters
    ----------
    atoms : ase.atoms.Atoms
        An ASE atoms object
    mol_indices : list of int
        List of atom indices

    Returns
    -------
    center_atom : int
        Index of the atom closest to the center of the box

    """

    scaled_positions = atoms[mol_indices].get_scaled_positions()

    distances = np.linalg.norm(scaled_positions - 0.5, axis=1)

    center_atom = mol_indices[np.argmin(distances)]

    return center_atom

def fix_mol(atoms, neighbor_list, mol_indices):
    """Unwraps a connected collection of atoms (molecule)

    Parameters
    ----------
    atoms : ase.atoms.Atoms
        An ASE atoms object
    neigbor_list : ase.neighborlist.Neighborlist
        ASE neighborlist
    mol_indices : list of int
        List of atomic indices defining the molecule

    """

    center_atom = get_center_atom(atoms=atoms, mol_indices=mol_indices)

    reordered = reorder_mol(neighbor_list=neighbor_list,
                            mol_indices=mol_indices,
                            center_atom=center_atom)

    cell = atoms.cell

    for atom_index in reordered:

        nb_indices, offsets = neighbor_list.get_neighbors(atom_index)

        for n, nb_index in enumerate(nb_indices):

            offset_n = offsets[n]

            for j, offset in enumerate(offset_n):

                atoms[nb_index].position += offset * cell[j]
                neighbor_list.update(atoms)

def unwrap(atoms, cutoff_correction=1.0):
    """Unwrap an ASE-readable structure

    Parameters
    ----------
    atoms : ase.atoms.Atoms
        An ASE atoms object
    cutoff_correction : float
        Parameter with which the ASE natural cutoffs will be multiplied

    """

    cutoffs = natural_cutoffs(atoms, mult=cutoff_correction)

    neighbor_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
    neighbor_list.update(atoms)

    matrix = neighbor_list.get_connectivity_matrix()
    n_components, component_list = sparse.csgraph.connected_components(matrix)

    mol_idxs = []

    for molecule_index in range(n_components):
        mol_idxs.append([i for i in range(len(component_list))
                   if component_list[i] == molecule_index])

    for molecule_indices in mol_idxs:
        fix_mol(atoms, neighbor_list, molecule_indices)

def unwrap_structure_file(input_file, cutoff_correction=1.0, output_file=None):
    """Reads a structure, unwraps it and saves it to a file

    Parameters
    ----------
    input_file : str
        Filepath to a structure to be unwrapped
    cutoff_correction : float
        Parameter with which the ASE natural cutoffs will be multiplied
    output : str
        Name of the output file; if None, will be automatically generated

    """

    atoms = ase.io.read(input_file)

    unwrap(atoms=atoms, cutoff_correction=cutoff_correction)

    if not output_file:

        prefix, suffix = input_file.split('.')
        output_file = prefix + '_unwrapped.' + suffix

    ase.io.write(filename=output_file, images=atoms)

def _main_cli():

    try:
        cutoff_corr = float(sys.argv[2])

    except IndexError:
        cutoff_corr = 1.0 #pylint: disable=invalid-name

    try:
        output = sys.argv[3]

    except IndexError:
        output = None #pylint: disable=invalid-name

    try:
        inputf = sys.argv[1]

    except IndexError:
        print('Input file must be provided.')

    else:
        unwrap_structure_file(input_file=inputf,
               cutoff_correction=cutoff_corr,
               output_file=output)

if __name__ == '__main__':
    _main_cli()
