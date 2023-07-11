from lib import EPS, AtomPair
import numpy as np
from pydash import py_
from copy import deepcopy


class EditableLayer:
    def __init__(self, eps=EPS) -> None:
        self.eps = eps
        self.atoms = []
        self.bonds = {}
        self.subscribers = []

    def add_subscribers(self, subscriber) -> None:
        self.subscribers.append(subscriber)

    def remove_subscriber(self, subscriber) -> None:
        self.subscribers.remove(subscriber)

    def update_state(self, update_fn):
        atoms, bonds = update_fn(self.atoms, self.bonds)
        self.atoms = atoms
        self.bonds = bonds
        for subscriber in self.subscribers:
            subscriber(atoms, bonds)

    def add_atom(self, new_atom):
        duplicated = py_.find(
            self.atoms, lambda atom: atom.get_id() == new_atom.get_id()
        )
        if duplicated is not None:
            raise ValueError(f"Duplicated atom id {new_atom.get_id()}.")

        too_closed = py_.find(
            self.atoms,
            lambda atom: np.linalg.norm(atom.position - new_atom.position) <= self.eps,
        )
        if too_closed is not None:
            raise ValueError(
                f"Atoms too closed: exsited {too_closed}, new {new_atom}, eps: {self.eps}"
            )

        def update(atoms, bonds):
            atoms = atoms + [new_atom]
            return atoms, bonds

        self.update_state(update)
        return new_atom.get_id()

    def remove_atom(self, atom_id):
        def update(atoms, bonds):
            atoms = py_.filter(atoms, lambda atom: atom.get_id() != atom_id)
            return atoms, bonds

        self.update_state(update)

    def replace_atom(self, atom_id, element):
        target = py_.find_index(self.atoms, lambda atom: atom.get_id() == atom_id)
        if target == -1:
            raise LookupError(f"No such atom {atom_id} found.")

        def update(atoms, bonds):
            atoms = (
                atoms[0:target] + [atoms[target].replace(element)] + atoms[target + 1 :]
            )
            return atoms, bonds

        self.update_state(update)

    def move_atom(self, atom_id, vector):
        target = py_.find_index(self.atoms, lambda atom: atom.get_id() == atom_id)
        if target is None:
            raise LookupError(f"No such atom {atom_id} found.")

        def update(atoms, bonds):
            atoms = atoms[0:target] + [atoms[target].move(vector)] + atoms[target + 1 :]
            return atoms, bonds

        self.update_state(update)

    def set_bonds(self, atom_ids, bond_type):
        target_pair = AtomPair(atom_ids)

        def update(atoms, bonds):
            bonds[target_pair] = bond_type
            return atoms, deepcopy(bonds)

        self.update_state(update)

    def get_bonds(self, from_atom):
        atom_pairs = self.bonds.keys()
        atom_pairs = py_.filter(
            atom_pairs,
            lambda atom_pair: atom_pair.has_atom(from_atom)
            and self.bonds[atom_pair] is not None,
        )
        return atom_pairs

    def get_bond(self, from_atom, to_atom):
        target_pair = AtomPair((from_atom, to_atom))
        return self.bonds[target_pair]
