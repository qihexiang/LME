from copy import deepcopy
import uuid
import numpy as np
from pydash import py_
import os
from scipy.spatial.transform import Rotation as R

PRODUCTION = os.getenv("LME_PRODUCTION") in ["Y", "y", "YES", "yes", 1]

EPS = 1e-8


class StateContainer:
    def __init__(self, init_state=None) -> None:
        self.__state__ = init_state
        self.subscribers = set()

    @property
    def state(self):
        # In prodcution mode, return the __state__ directly
        if PRODUCTION:
            return self.__state__
        # In development mode, return an copied state to avoid change the state from outside.
        return deepcopy(self.__state__)

    def add_subscriber(self, subscriber):
        self.subscribers.add(subscriber)

    def remove_subscriber(self, subscriber):
        self.subscribers.remove(subscriber)

    def __broadcast__(self):
        for subscriber in self.subscribers:
            subscriber(self.state)

    def update(self, updator):
        self.__state__ = updator(self.state)
        self.__broadcast__()


class Atom:
    def __init__(self, element, position) -> None:
        self.__element = element
        self.__position = np.array(position, dtype="float64")

    @property
    def element(self):
        return self.__element
    
    @property
    def position(self):
        return self.__position

    def replace(self, element):
        return Atom(element, self.position)

    def move_to(self, target) -> None:
        return Atom(self.element, target)

    def move(self, vector) -> None:
        target = self.position + vector
        return self.move_to(target)

    def copy(self):
        return Atom(self.element, self.position)

    def __repr__(self) -> str:
        x, y, z = self.position
        return f"{self.element} {x} {y} {z}"


class UUIDPair:
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
        return (another.a == self.a and another.b == self.b) or (
            another.a == self.b and another.b == self.a
        )

    def __hash__(self) -> int:
        sorted_id = py_.sort([self.a, self.b])
        return hash(tuple(sorted_id))

    def __repr__(self) -> str:
        return f"{self.a} {self.b}"


def mirror_matrix(law_vector):
    [x, y, z] = law_vector / np.linalg.norm(law_vector)
    return np.array(
        [
            [1 - 2 * (x**2), -2 * x * y, -2 * x * z],
            [-2 * x * y, 1 - 2 * (y**2), -2 * y * z],
            [-2 * x * z, -2 * y * z, 1 - 2 * (z**2)],
        ],
        dtype="float64",
    )


def rotate_matrix(axis, times):
    axis = np.array(axis, dtype="float64")
    axis = axis / np.linalg.norm(axis)
    return R.from_rotvec(axis * (360 / times), degrees=True).as_matrix()


if __name__ == "__main__":
    a1 = Atom("C", [1, 1, 0])
    a2 = a1.move_to([1, 2, 0])
    s = {a1: 1, a2: 2}
    print(s)
