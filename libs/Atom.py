import numpy as np
from pydash import py_
import uuid

class Atom:
    def __init__(self, element, position, class_name = "") -> None:
        self.__element = element
        self.__position = np.array(position, dtype="float64")
        self.class_name = class_name

    @property
    def element(self):
        return self.__element

    @property
    def position(self):
        return self.__position

    def replace(self, element):
        return Atom(element, self.position, self.class_name)

    def move_to(self, target) -> None:
        return Atom(self.element, target, self.class_name)

    def move(self, vector) -> None:
        target = self.position + vector
        return self.move_to(target)

    def copy(self):
        return Atom(self.element, self.position, self.class_name)

    def __repr__(self) -> str:
        x, y, z = self.position
        return f"{self.element} {x} {y} {z} {self.class_name}"

    @property
    def export(self):
        return {"element": self.element, "position": list(self.position), "class_name": self.class_name}


class AtomWithId(Atom):
    @staticmethod
    def from_atoms(atoms_dict):
        return py_.filter(
            [
                AtomWithId(atom_id, atoms_dict[atom_id])
                if atoms_dict[atom_id] is not None
                else None
                for atom_id in atoms_dict.keys()
            ],
            lambda atom_with_id: atom_with_id is not None,
        )

    @staticmethod
    def to_atoms_dict(atoms_list):
        return {
            atom_with_id.get_id(): atom_with_id.__to_atom()
            for atom_with_id in atoms_list
        }

    def __init__(self, atom_id, atom):
        super().__init__(atom.element, atom.position)
        self.__id = atom_id

    def get_id(self):
        return self.__id

    def replace(self, element):
        return AtomWithId(self.get_id(), super().replace(element))

    def move_to(self, target) -> None:
        return AtomWithId(self.get_id(), super().move_to(target))

    def copy(self):
        return AtomWithId(uuid.uuid4(), super().copy())

