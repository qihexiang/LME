from uuid import UUID
from libs.Atom import Atom
from libs.UUIDPair import UUIDPair
from pydash import py_

def atoms_bonds_loader(atoms_dict, bonds_dict):
    atoms = {
        UUID(atom_id): Atom(
            atoms_dict[atom_id]["element"], atoms_dict[atom_id]["position"]
        )
        if atoms_dict[atom_id] is not None
        else None
        for atom_id in atoms_dict.keys()
    }

    bonds = {
        UUIDPair(tuple(py_.map(bond_id.split(" "), UUID))): bonds_dict[bond_id]
        for bond_id in bonds_dict.keys()
    }

    return atoms, bonds