import uuid
import numpy as np
from pydash import py_

EPS = 1e-8

class Atom:
    def __init__(self, element, position, atom_id=None) -> None:
        if atom_id == None:
            self.__atom_id = uuid.uuid4()
        else:
            self.__atom_id = atom_id

        self.element = element
        self.position = np.array(position, dtype="float64")

    def get_id(self):
        return self.__atom_id

    def replace(self, element):
        return Atom(element, self.position, self.get_id())

    def move_to(self, target) -> None:
        return Atom(self.element, target, self.get_id())
    
    def move(self, vector) -> None:
        target = self.position + vector
        return self.move_to(target)
    
    def copy(self):
        return Atom(self.element, self.position, uuid.uuid4())
    
    def __repr__(self) -> str:
        return f"{self.element} {self.position}"

class AtomPair:
    def __init__(self, ids) -> None:
        self.a, self.b = ids

    def has_atom(self, atom_id) -> bool:
        return self.another_atom(atom_id) is not None
    
    def another_atom(self, atom_id):
        if self.a == atom_id:
            return self.b
        
        if self.b == atom_id:
            return self.a
        
        return None

    def __eq__(self, another) -> bool:
        return (another.a == self.a and another.b == self.b) or (another.a == self.b and another.b == self.a)
    
    def __hash__(self) -> int:
        sorted_id = py_.sort([self.a, self.b])
        return hash(tuple(sorted_id))
    
    def __repr__(self) -> str:
        return f"{self.a} {self.b}"