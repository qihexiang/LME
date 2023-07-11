import numpy as np
from copy import deepcopy
from pydash import py_
from lib import AtomPair, EPS
from scipy.spatial.transform import Rotation as R


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
    对称元素i
    """
    def __init__(self, center=[0.0, 0.0, 0.0], eps=EPS) -> None:
        super().__init__()
        self.eps = eps
        self.center = np.array(center, dtype="float64")

    def on_center(self, position):
        return np.linalg.norm(position - self.center) < self.eps

    def should_ignore_on_copy(self, atom):
        return self.on_center(atom.position)

    def inverse(self, positions):
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


class MirrorLayer(SymmetryLayer):
    def __init__(self, law_vector, center=[0.0, 0.0, 0.0], eps=EPS) -> None:
        super().__init__()
        law_vector = np.array(law_vector, dtype="float64")
        self.law_vector = law_vector / np.linalg.norm(law_vector)
        self.matrix = mirror_matrix(law_vector)
        self.center = np.array(center, dtype="float64")
        self.eps = eps

    def on_mirror(self, position):
        OP = position - self.center
        return abs(np.dot(OP, self.law_vector)) < self.eps

    def mirror(self, positions):
        positions = positions - self.center
        targets = np.matmul(positions, self.matrix)
        return targets + self.center

    def should_ignore_on_copy(self, atom):
        return self.on_mirror(atom.position)

    def __call__(self, atoms, bonds):
        mirrored_positions = self.mirror(py_.map(atoms, lambda atom: atom.position))
        mirrored = (
            py_.chain(zip(self.copy_atoms(atoms), mirrored_positions))
            .map(lambda ap: ap[0].move_to(ap[1]))
            .value()
        )
        atoms = py_.uniq_by(atoms + mirrored, lambda atom: atom.get_id())
        bonds = bonds | self.generate_bonds(atoms, mirrored, bonds)
        return atoms, bonds


class RotationLayer(SymmetryLayer):
    """
    Cn, In, Sn对称元素
    """

    def __init__(self, axis, times, mode="C", center=[0.0, 0.0, 0.0], eps=EPS) -> None:
        """
        初始化一个旋转操作层
        ---
        parameters:

        - axis: 旋转轴坐标
        - times: 旋转轴次数n
        - mode: 可选值为"C", "S"或"I", 分别对应旋转轴、映轴和反轴, 设置为以外的任何值时，作为旋转轴处理，默认为"C"
        - center: 轴中点，默认为[0,0,0]
        - eps: 容差范围, 用于判断是否与轴重合, 默认值为1E-8
        """
        super().__init__()
        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)
        self.center = np.array(center)
        self.times = int(times)
        self.eps = eps
        self.axis = axis
        before = np.eye(3)
        if mode == "I":
            before = before * -1.0
        if mode == "S":
            before = mirror_matrix(axis)
        self.matrices = [
            np.matmul(
                before ** (i if i > 0 else 1),
                R.from_rotvec(axis * (360 / times * i), degrees=True).as_matrix(),
            )
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

    I4 = RotationLayer([0, 0, 1], 4, "S")

    layer = EditableLayer()
    layer.add_subscribers(lambda a, b: print(I4(a, b)))
    C = layer.add_atom(Atom("C", [0, 0, 0]))
    H = layer.add_atom(Atom("H", [1, 0, 1]))
    layer.set_bonds((C, H), 1.0)
