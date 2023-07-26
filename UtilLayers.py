from pydash import py_
import numpy as np
from copy import deepcopy
from libs.Atom import AtomWithId
from libs.UUIDPair import UUIDPair

from libs.constants import EPS


class DedupLayer:
    def __init__(self, eps=EPS) -> None:
        self.eps = eps

    def __call__(self, atoms, bonds):
        atoms = AtomWithId.from_atoms(atoms)

        def find_dups(atom):
            dups = (
                py_.chain(atoms)
                .filter(
                    lambda target: target.element == atom.element
                    and np.linalg.norm(target.position - atom.position) < self.eps
                )
                .map(lambda target: UUIDPair((target.get_id(), atom.get_id())))
                .value()
            )
            dups = set(dups)
            dups.remove(UUIDPair((atom.get_id(), atom.get_id())))
            return list(dups)

        dups = py_.map(atoms, find_dups)
        dups = set(py_.flatten(dups))

        def split_groups(acc, ne):
            a, b = ne.a, ne.b
            idx = py_.find_index(acc, lambda s: a in s or b in s)
            if idx == -1:
                return acc + [{a, b}]
            target = acc[idx] | {a, b}
            return acc[0:idx] + [target] + acc[idx + 1 :]

        dup_groups = py_.reduce(dups, split_groups, [])
        atoms = deepcopy(atoms)
        bonds = deepcopy(bonds)
        for group in dup_groups:
            group = list(group)
            kept, removed = group[0], group[1:]
            connected = {}
            for atom in removed:
                targets = py_.filter(bonds.keys(), lambda ap: ap.has_uuid(atom))
                if len(targets) == 0:
                    pass
                else:
                    bond_type = bonds[targets[0]]
                    anothers = py_.map(targets, lambda ap: ap.another_atom(atom))
                    connected = connected | {another: bond_type for another in anothers}
                    for target in targets:
                        bonds.pop(target)
            bonds = bonds | {
                UUIDPair((kept, another)): connected[another]
                for another in connected.keys()
            }
            atoms = py_.filter(atoms, lambda atom: atom.get_id() not in removed)
        atoms = AtomWithId.to_atoms_dict(atoms)
        return atoms, bonds

    @property
    def export(self):
        return {"type": "dedup", "eps": self.eps}


class AutoBondLayer:
    def __init__(self, radius_table) -> None:
        elements = radius_table.keys()
        self.elements_table = py_.flatten(
            py_.map(elements, lambda e1: py_.map(elements, lambda e2: (e1, e2)))
        )
        radiuses = py_.map(elements, lambda element: radius_table[element])
        self.bond_length_matrix = py_.flatten(
            py_.map(radiuses, lambda r1: py_.map(radiuses, lambda r2: r1 + r2))
        )
        self.max_bond_length = np.array(self.bond_length_matrix, dtype="float64").max()

    def max_bond_length(self, target):
        idx = py_.find_index(self.elements_table, target)
        if idx == -1:
            return None
        return self.bond_length_matrix[idx]

    # def __call__(self, atoms, bonds) -> Any:
