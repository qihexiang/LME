from copy import deepcopy
import uuid
import numpy as np
from pydash import py_
import json
from libs.StateContainer import StateContainer
from libs.UUIDPair import UUIDPair
from libs.atoms_bonds_loader import atoms_bonds_loader
from libs.constants import EPS
from libs.matrix import rotate_matrix
from libs.molecule_text import atoms_bonds_from_mol2, mol2_to_atom, mol2_to_bond, molecule_text
from layers.StaticLayer import StaticLayer
import re

class DiffLayer:
    def __init__(self, atoms, bonds) -> None:
        self.__atoms = deepcopy(atoms)
        self.__bonds = deepcopy(bonds)

    def __call__(self, atoms, bonds):
        return atoms | self.__atoms, bonds | self.__bonds

    @property
    def export(self):
        return {
            "type": "diff",
            "atoms": self.__atoms,
            "bonds": self.__bonds
        }

class EditableLayer(StateContainer):
    @staticmethod
    def from_mol2(text):
        editable = EditableLayer()
        atoms_lines, bonds_lines = atoms_bonds_from_mol2(text)
        atoms = py_.map(atoms_lines, mol2_to_atom)
        atoms_ids = py_.map(atoms, lambda atom: atom[0])
        atom_uuid_list = editable.add_atoms(py_.map(atoms, lambda atom: atom[1]))
        atom_uuids = {
            atom_id: atom_uuid for (atom_id, atom_uuid) in zip(atoms_ids, atom_uuid_list)
        }
        bonds = py_.map(bonds_lines, mol2_to_bond)
        for (a, b, order) in bonds:
            uuid_a, uuid_b = atom_uuids[a], atom_uuids[b]
            editable.set_bond(uuid_a, uuid_b, order)

        return editable
        
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
        overlayed = self.base.atoms | self.state[0]
        existed = {
            atom_id: overlayed[atom_id]
            for atom_id in py_.filter(
                overlayed.keys(), lambda atom_id: overlayed[atom_id] is not None
            )
        }
        return existed

    @property
    def bonds(self):
        existed_atom_ids = self.atom_ids
        overlayed = self.base.bonds | self.state[1]
        existed = {
            bond_id: overlayed[bond_id]
            for bond_id in py_.filter(
                overlayed.keys(),
                lambda bond_id: bond_id.a in existed_atom_ids
                and bond_id.b in existed_atom_ids
                and overlayed[bond_id] is not None,
            )
        }
        return existed

    @property
    def selected(self):
        return self.state[2]

    @property
    def atom_ids(self):
        return set(self.atoms.keys())

    @property
    def bond_ids(self):
        return set(self.bonds.keys())

    def find_with_classname(self, class_name):
        return py_.filter(self.atom_ids, lambda atom_id: class_name in self.atoms[atom_id].class_name.split(" "))

    def detach_diff_layer(self):
        return DiffLayer(self.state[0], self.state[1])

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
            return a, b, set(self.atom_ids)

        self.update(updator)
        return 0

    def import_atoms_bonds(self, atoms, bonds):
        self.deselect_all()
        self.__patch_to_atoms(atoms)
        self.__patch_to_bonds(bonds)
        self.select(atoms.keys())

    def add_atoms(self, atoms):
        # to_add = {uuid(): atom for atom in atoms}
        to_add = py_.map(atoms, lambda atom: (uuid.uuid4(), atom))
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
    
    def add_substitute(self, substitute, center_id, entry_id):
        center = self.atoms[center_id]
        entry = self.atoms[entry_id]
        bond = self.bonds[UUIDPair((center_id, entry_id))]
        direction = center.position - entry.position
        direction = direction / np.linalg.norm(direction)
        atoms, bonds, center_idx, entry_idx = substitute.output(direction)
        self.import_atoms_bonds(atoms, bonds)
        new_center = self.atoms[center_idx]
        to_move = center.position - new_center.position
        self.translation_selected(to_move)
        self.deselect_all()
        self.select([center_id, entry_idx])
        self.remove_selected()
        self.set_bond(center_idx, entry_id, bond)

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


substitute_re = re.compile("Substitute\(entry=(?P<entry>.*),center=(?P<center>.*)\)")
class Substitute(StaticLayer):
    @staticmethod
    def from_mol2(text: str):
        layer = EditableLayer.from_mol2(text)
        entry_name, center_name = substitute_re.search(text).groups()
        entry_idx = layer.find_with_classname(entry_name)[0]
        center_idx = layer.find_with_classname(center_name)[0]
        return Substitute.build_from(layer, entry_idx, center_idx)


    @staticmethod
    def build_from(layer, entry_idx, center_idx):
        return Substitute(layer.atoms, layer.bonds, entry_idx, center_idx)

    def __init__(self, atoms, bonds, entry_idx, center_idx, eps = EPS) -> None:
        super().__init__(atoms, bonds)
        self._vector = self.atoms[center_idx].position - self.atoms[entry_idx].position
        self._vector = self._vector / np.linalg.norm(self._vector)
        self._entry_idx = entry_idx
        self._center_idx = center_idx
        self.eps = eps

    def output(self, direction):
        updated_atom_ids = {
            origin_id: uuid.uuid4() for origin_id in self.atom_ids
        }

        updated_atoms_table = {
            updated_atom_ids[origin_id]: self.atoms[origin_id] for origin_id in self.atom_ids
        }

        updated_bond_ids = {
            uuid_pair: UUIDPair(tuple(py_.map(uuid_pair, lambda origin_uuid: updated_atom_ids[origin_uuid]))) for uuid_pair in self.bond_ids
        }

        updated_bonds_table = {
            updated_bond_ids[origin_uuid_pair]: self.bonds[origin_uuid_pair] for origin_uuid_pair in self.bond_ids
        }

        updated_center_idx = updated_atom_ids[self._center_idx]
        updated_entry_idx = updated_atom_ids[self._entry_idx]

        direction = np.array(direction) / np.linalg.norm(direction)
        axis = np.cross(direction, self._vector)
        if(np.linalg.norm(axis) == 0.):
            [x, y, z] = direction
            if x == 0 and y == 0:
                axis = np.array([1, 1, 0], dtype="float64")
            else:
                axis = np.array([y, -x, 0], dtype="float64")
        else:
            axis = axis/np.linalg.norm(axis)
        angle = np.arccos(np.dot(self._vector, direction)) / (2*np.pi) * 360
        editable_layer = EditableLayer()
        editable_layer.import_atoms_bonds(updated_atoms_table, updated_bonds_table)
        editable_layer.select_all()
        editable_layer.rotation_selected(axis, editable_layer.atoms[updated_entry_idx].position, angle)
        return editable_layer.atoms, editable_layer.bonds, updated_center_idx, updated_entry_idx

if __name__ == "__main__":
    from libs.Atom import Atom
    from SymmetryLayers import RotationLayer
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
    exported = editor.__repr__()
    print(editor)
    with open("./export.json", "w") as f:
        f.write(json.dumps(editor.export))

    with open("./export.json", "r") as f:
        data = f.read()
        data = json.loads(data)
        layer = EditableLayer(load=data)
        imported = layer.__repr__()
        assert exported == imported
        editable, transformers = layer.base.extract()
        print(editable, transformers)
        # print(layer.base.contains)

