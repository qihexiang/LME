import numpy as np
from copy import deepcopy
from pydash import py_
from lib import AtomPair, EPS
from scipy.spatial.transform import Rotation as R

def compose_layers(layers):
    def composed(a, b):
        for layer in layers:
            a, b = layer(a, b)
        return a, b
    return composed

class SymmetryLayer:
    """
    对称层基类。
    """
    def __init__(self) -> None:
        pass

    def should_ignore_on_copy(self, atom):
        """
        该函数用于在调用`copy_atoms`方法时判断是否需要对原子进行复制(即生成新的ID),对每个子类,应该手动实现该方法。
        """
        raise NotImplementedError("Must implemented before using default copy_atoms")

    def copy_atoms(self, origin_atoms):
        """
        调用该方法从原始原子数组中生成一份拷贝，如果该原子被`should_ignore_on_copy`返回为`True`,则使用Atom.copy(),否则使用deepcopy()
        """
        return py_.map(
            origin_atoms,
            lambda atom: atom.copy()
            if not self.should_ignore_on_copy(atom)
            else deepcopy(atom),
        )

    def generate_bonds(self, origin_atoms, new_atoms, bonds):
        """
        根据原始原子数组和键连接的情况，生成新的原子数组中原子的键连接情况
        """
        origin_ids = py_.map(origin_atoms, lambda atom: atom.get_id())
        new_ids = py_.map(new_atoms, lambda atom: atom.get_id())
        id_refs = {}
        for origin_id, new_ids in zip(origin_ids, new_ids):
            id_refs[origin_id] = new_ids

        new_bonds = {}
        for atom_pair in bonds.keys():
            a, b = atom_pair.a, atom_pair.b
            a, b = id_refs[a], id_refs[b]
            new_pair = AtomPair((a, b))
            new_bonds[new_pair] = bonds[atom_pair]

        return new_bonds


class InverseLayer(SymmetryLayer):
    """
    中心对称层
    """
    def __init__(self, center=[0.0, 0.0, 0.0], eps=EPS) -> None:
        super().__init__()
        self.eps = eps
        self.center = np.array(center, dtype="float64")

    def on_center(self, position):
        return np.linalg.norm(position - self.center) < self.eps
    
    def should_ignore_on_copy(self, atom):
        return self.on_center(atom.position)
    
    def inverse(self,positions):
        centered = positions - self.center
        inversed = -1 * centered
        return inversed + self.center

    def __call__(self, atoms, bonds):
        inversed_positions = self.inverse(py_.map(atoms, lambda atom: atom.position))
        inversed = (
            py_.chain(zip(self.copy_atoms(atoms), inversed_positions))
            .map(lambda ap: ap[0].move_to(ap[1]))
            .value()
        )
        atoms = py_.uniq_by(atoms + inversed, lambda atom: atom.get_id())
        bonds = bonds | self.generate_bonds(atoms, inversed, bonds)
        return atoms, bonds


class MirrorLayer(SymmetryLayer):
    def __init__(self, law_vector, center=[0.0, 0.0, 0.0], eps=EPS) -> None:
        super().__init__()
        law_vector = np.array(law_vector, dtype="float64")
        self.law_vector = law_vector / np.linalg.norm(law_vector)
        self.center = np.array(center, dtype="float64")
        self.eps = eps

    def on_mirror(self, position):
        OP = position - self.center
        return abs(np.dot(OP, self.law_vector)) < self.eps

    def mirror(self, position):
        OP = position - self.center
        mirror_P = np.dot(OP, self.law_vector) * self.law_vector
        target = position - 2 * mirror_P
        return target
    
    def should_ignore_on_copy(self, atom):
        return self.on_mirror(atom.position)

    def __call__(self, atoms, bonds):
        mirrored = (
            py_.chain(self.copy_atoms(atoms))
            .map(lambda atom: atom.move_to(self.mirror(atom.position)))
            .value()
        )
        atoms = py_.uniq_by(atoms + mirrored, lambda atom: atom.get_id())
        bonds = bonds | self.generate_bonds(atoms, mirrored, bonds)
        return atoms, bonds


class RotationLayer(SymmetryLayer):
    def __init__(self, axis, times, center=[0.0, 0.0, 0.0], eps=EPS) -> None:
        super().__init__()
        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)
        self.center = np.array(center)
        self.times = int(times)
        self.eps = eps
        self.axis = axis
        self.matrices = [
            R.from_rotvec(axis * (360 / times * i), degrees=True).as_matrix()
            for i in range(0, times)
        ]

    def rotate(self, positions):
        positions = positions - self.center
        return [np.matmul(positions, matrix) + self.center for matrix in self.matrices]

    def on_axis(self, position):
        OP = position - self.center
        # self.axis is a unit vector
        delta = OP - np.dot(OP, self.axis) * self.axis
        delta = np.linalg.norm(delta)
        return delta < self.eps

    def should_ignore_on_copy(self, atom):
        return self.on_axis(atom.position)

    def __call__(self, atoms, bonds):
        input_positions = py_.map(atoms, lambda atom: atom.position)
        # Take only rotated groups
        rotated_position_groups = self.rotate(input_positions)[1:]
        rotated_atoms = [
            py_.map(zip(self.copy_atoms(atoms), group), lambda ap: ap[0].move_to(ap[1]))
            for group in rotated_position_groups
        ]
        rotated_bonds = [
            self.generate_bonds(atoms, new_atoms, bonds) for new_atoms in rotated_atoms
        ]
        atoms = (
            py_.chain([atoms] + rotated_atoms)
            .flatten()
            .uniq_by(lambda atom: atom.get_id())
            .value()
        )
        bonds = py_.reduce([bonds] + rotated_bonds, lambda acc, next: acc | next)
        return atoms, bonds


if __name__ == "__main__":
    from lib import Atom
    from editable_layer import EditableLayer

    C3 = RotationLayer([0,0,1], 3)

    layer = EditableLayer()
    layer.add_subscribers(lambda a,b: print(C3(a,b)))
    N = layer.add_atom(Atom("N", [0,0,0]))
    H = layer.add_atom(Atom("H", [1,0,-1]))
    layer.set_bonds((N,H), 1.0)