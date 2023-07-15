import numpy as np
from lib import UUIDPair, StateContainer, rotate_matrix
from uuid import uuid4 as uuid
from pydash import py_
import json
from uuid import UUID


def molecule_text(target) -> str:
    output = ""
    output += "atoms:\n"
    atom_ids = py_.filter(
        target.atom_ids, lambda atom_id: target.atoms[atom_id] is not None
    )
    atoms = target.atoms
    atoms = py_.map(atom_ids, lambda atom_id: atoms[atom_id])
    for i, atom in enumerate(atoms):
        output += f"{i+1} {atom}\n"

    bond_ids = target.bond_ids
    bonds = target.bonds
    output += "bonds:\n"
    for bond_id in bond_ids:
        a_idx, b_idx = py_.find_index(
            atom_ids, lambda atom_id: atom_id == bond_id.a
        ), py_.find_index(atom_ids, lambda atom_id: atom_id == bond_id.b)
        if -1 not in [a_idx, b_idx]:
            output += f"{a_idx + 1} {b_idx + 1} {bonds[bond_id]}\n"

    return output


def atoms_bonds_loader(atoms_dict, bonds_dict):
    atoms = {
        UUID(atom_id): Atom(atoms_dict[atom_id]["element"], atoms_dict[atom_id]["position"])
        if atoms_dict[atom_id] is not None
        else None
        for atom_id in atoms_dict.keys()
    }

    bonds = {
        UUIDPair(tuple(py_.map(bond_id.split(" "), UUID))): bonds_dict[bond_id]
        for bond_id in bonds_dict.keys()
    }

    return atoms, bonds


class StaticLayer:
    def __init__(self, atoms=dict(), bonds=dict(), contains=None, load=None) -> None:
        if load is not None:
            if load["type"] != "static":
                raise ValueError("Not a StaticLayer dict")
            self.__contains = load["contains"]
            self.__atoms, self.__bonds = atoms_bonds_loader(
                load["atoms"], load["bonds"]
            )
            return None
        if contains is not None:
            (editable, transfomers) = contains
            self.__contains = {
                "editable": editable.export,
                "transfomers": py_.map(
                    transfomers, lambda transformer: transformer.export
                ),
            }
            atoms, bonds = editable.atoms, editable.bonds
            for transformer in transfomers:
                atoms, bonds = transformer(atoms, bonds)
        else:
            self.__contains = None
        atom_ids = py_.filter(atoms.keys(), lambda atom_id: atoms[atom_id] is not None)
        self.__atoms = {atom_id: atoms[atom_id] for atom_id in atom_ids}
        bond_ids = py_.filter(
            bonds.keys(),
            lambda bond_id: bonds[bond_id] is not None
            and bond_id.a in atom_ids
            and bond_id.b in atom_ids,
        )
        self.__bonds = {bond_id: bonds[bond_id] for bond_id in bond_ids}

    @property
    def atoms(self):
        return self.__atoms | {}

    @property
    def bonds(self):
        return self.__bonds | {}

    @property
    def atom_ids(self):
        return self.atoms.keys()

    @property
    def bond_ids(self):
        return self.bonds.keys()

    def __repr__(self) -> str:
        return molecule_text(self)

    @property
    def export(self):
        atoms = self.atoms
        bonds = self.bonds
        atoms = {
            str(atom_id): atoms[atom_id].export if atoms[atom_id] is not None else None
            for atom_id in atoms.keys()
        }
        bonds = {bond_id.export: bonds[bond_id] for bond_id in bonds.keys()}
        return {
            "type": "static",
            "atoms": atoms,
            "bonds": bonds,
            "contains": self.__contains,
        }


class EditableLayer(StateContainer):
    def __init__(self, base=StaticLayer(), load=None) -> None:
        if load is not None:
            if load["type"] != "editable":
                raise ValueError("Not a EditableLayer dict")
            atoms, bonds = atoms_bonds_loader(load["atoms"], load["bonds"])
            super().__init__((atoms, bonds, set()))
            self.base = StaticLayer(load=load["base"])
        else:
            super().__init__((dict(), dict(), set()))
            self.base = base

    @property
    def atoms(self):
        return self.base.atoms | self.state[0]

    @property
    def bonds(self):
        return self.base.bonds | self.state[1]

    @property
    def selected(self):
        return self.state[2]

    def __patch_to_atoms(self, patch):
        def updator(state):
            a, b, s = state
            return a | patch, b, s

        self.update(updator)

    def __patch_to_bonds(self, patch):
        def updator(state):
            a, b, s = state
            return a, b | patch, s

        self.update(updator)

    def select(self, atom_ids):
        if (
            py_.chain(atom_ids)
            .map(lambda atom_id: self.atoms.get(atom_id))
            .every(lambda atom: atom is not None)
            .value()
        ):

            def updator(state):
                a, b, s = state
                return a, b, (s | set(atom_ids))

            self.update(updator)
            return 0
        raise KeyError("Some of atoms to select are not existed")

    def deselect(self, atom_ids):
        def updator(state):
            a, b, s = state
            return a, b, set(py_.filter(s, lambda selected: selected not in atom_ids))

        self.update(updator)
        return 0

    def deselect_all(self):
        def updator(state):
            a, b, _ = state
            return a, b, set()

        self.update(updator)
        return 0

    def select_all(self):
        def updator(state):
            a, b, s = state
            return a, b, set(a.keys())

        self.update(updator)
        return 0

    @property
    def atom_ids(self):
        return set(self.atoms.keys())

    @property
    def bond_ids(self):
        return set(self.bonds.keys())

    def add_atoms(self, atoms):
        # to_add = {uuid(): atom for atom in atoms}
        to_add = py_.map(atoms, lambda atom: (uuid(), atom))
        patch = {atom_id: atom for atom_id, atom in to_add}
        self.__patch_to_atoms(patch)
        return py_.map(to_add, lambda id_atom: id_atom[0])

    def __remove_atoms(self, atom_ids):
        if self.atom_ids.issuperset(atom_ids):
            self.__patch_to_atoms({atom_id: None for atom_id in atom_ids})
            return 0
        raise KeyError("Non-existed id found.")

    def set_bond(self, atom_a_id, atom_b_id, bond):
        atom_a, atom_b = self.atoms.get(atom_a_id), self.atoms.get(atom_b_id)
        if None not in [atom_a, atom_b]:
            self.__patch_to_bonds({UUIDPair((atom_a_id, atom_b_id)): bond})
            return 0
        raise KeyError("At least one of the atoms not existed.")

    def remove_selected(self):
        self.__remove_atoms(self.selected)
        bonds_remove = py_.filter(
            self.bond_ids,
            lambda atom_pair: atom_pair.a in self.selected
            or atom_pair.b in self.selected,
        )
        self.__patch_to_bonds({bond_id: None for bond_id in bonds_remove})
        self.deselect_all()
        return 0

    def set_element_selected(self, element):
        updated = {
            atom_id: self.atoms[atom_id].replace(element) for atom_id in self.selected
        }
        self.__patch_to_atoms(updated)
        return 0

    def translation_selected(self, vector):
        translated = {
            atom_id: self.atoms[atom_id].move(vector) for atom_id in self.selected
        }
        self.__patch_to_atoms(translated)
        return 0

    def rotation_selected(self, axis, center, angle):
        center = np.array(center, dtype="float64")
        matrix = rotate_matrix(axis, angle)

        def calculate_target_position(origin):
            return np.matmul((origin - center), matrix) + center

        atom_ids = list(self.selected)
        rotated_atoms = (
            py_.chain(atom_ids)
            .map(lambda atom_id: self.atoms[atom_id])
            .map(lambda atom: atom.move_to(calculate_target_position(atom.position)))
            .value()
        )
        rotated = {atom_id: atom for (atom_id, atom) in zip(atom_ids, rotated_atoms)}
        self.__patch_to_atoms(rotated)
        return 0

    def __repr__(self):
        return molecule_text(self)

    @property
    def export(self):
        atoms, bonds, _ = self.state
        atoms = {
            str(atom_id): atoms[atom_id].export if atoms[atom_id] is not None else None
            for atom_id in atoms.keys()
        }
        bonds = {bond_id.export: bonds[bond_id] for bond_id in bonds.keys()}
        return {
            "type": "editable",
            "atoms": atoms,
            "bonds": bonds,
            "base": self.base.export,
        }


if __name__ == "__main__":
    from lib import Atom
    from symmetry_layers import RotationLayer
    import json

    editor = EditableLayer()
    C3_axis = RotationLayer([0, 0, 1], 3, "C")
    C2_axis = RotationLayer([0, 1, 0], 2, "C")
    [C, H] = editor.add_atoms([Atom("C", [0, 0, 1]), Atom("H", [1, 0, 1.2])])
    editor.set_bond(C, H, 1.0)
    locked = StaticLayer(contains=(editor, [C3_axis, C2_axis]))
    editor = EditableLayer(base=locked)
    [C1, C2] = py_.filter(
        editor.atom_ids, lambda atom_id: editor.atoms[atom_id].element == "C"
    )
    editor.set_bond(C1, C2, 1.0)
    print(editor)
    with open("./export.json", "w") as f:
        f.write(json.dumps(editor.export))
    
    with open("./export.json", "r") as f:
        data = f.read()
        data = json.loads(data)
        print(EditableLayer(load=data))
