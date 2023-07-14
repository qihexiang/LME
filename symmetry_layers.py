import numpy as np
from copy import deepcopy
from pydash import py_
from editable_layer import StaticLayer, molecule_output
from lib import UUIDPair, EPS, Atom, mirror_matrix, rotate_matrix
from scipy.spatial.transform import Rotation as R
from util_layers import DedupLayer
from uuid import uuid4 as uuid


def compose_layers(layers):
    def composed(a, b):
        for layer in layers:
            a, b = layer(a, b)
        return a, b

    return composed


class AtomWithId(Atom):
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
        return AtomWithId(uuid(), super().copy())


def atoms_format_transformer(fn):
    """
    该装饰器会将参数atoms的结构{UUID: Atom}转换为[AtomWithId]后传给被装饰函数进行计算

    计算结束后, 该装饰器会将返回的atoms从[AtomWithId]转换回{UUID: Atom}
    """

    def wrapped(self, atoms, bonds):
        atom_ids = py_.filter(atoms.keys(), lambda atom_id: atoms[atom_id] is not None)
        atoms = py_.map(atom_ids, lambda uid: atoms[uid])
        atom_with_ids = [
            AtomWithId(atom_id, atom) for atom_id, atom in zip(atom_ids, atoms)
        ]
        atoms, bonds = fn(self, atom_with_ids, bonds)
        atoms = {atom.get_id(): Atom(atom.element, atom.position) for atom in atoms}
        return atoms, bonds

    return wrapped


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
            new_pair = UUIDPair((a, b))
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

    @atoms_format_transformer
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

    @atoms_format_transformer
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

    当旋转操作层是反轴或映轴时, 输出结果会通过`DedupLayer`进行过滤
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
        self.matrix = np.matmul(before, rotate_matrix(axis, times))
        self.dedup = DedupLayer(self.eps) if mode in ["I", "S"] else None

    def rotate(self, positions):
        positions = positions - self.center
        return np.matmul(positions, self.matrix)

    def on_axis(self, position):
        OP = position - self.center
        # self.axis is a unit vector
        delta = OP - np.dot(OP, self.axis) * self.axis
        delta = np.linalg.norm(delta)
        return delta < self.eps

    def on_center(self, position):
        return np.linalg.norm(position - self.center) < self.eps

    def should_ignore_on_copy(self, atom):
        if self.dedup is None:
            return self.on_axis(atom.position)
        return self.on_center(atom.position)

    @atoms_format_transformer
    def __call__(self, atoms, bonds):
        rotated_atoms = deepcopy(atoms)
        rotated_bonds = deepcopy(bonds)
        current_positions = py_.map(atoms, lambda atom: atom.position)
        for _ in range(
            1, self.times * (1 if self.dedup is None or self.times % 2 == 0 else 2)
        ):
            current_positions = self.rotate(current_positions)
            new_atoms = self.copy_atoms(atoms)
            atom_positions = zip(new_atoms, current_positions)
            rotated_atoms += py_.map(atom_positions, lambda ap: ap[0].move_to(ap[1]))
            rotated_bonds = rotated_bonds | self.generate_bonds(atoms, new_atoms, bonds)
        rotated_atoms = py_.uniq_by(rotated_atoms, lambda atom: atom.get_id())
        if self.dedup is None:
            return rotated_atoms, rotated_bonds
        return self.dedup(rotated_atoms, rotated_bonds)


if __name__ == "__main__":
    from lib import Atom
    from editable_layer import EditableLayer

    symmetry = RotationLayer([0, 0, 1], 6, "S")

    output_layer = StaticLayer()

    def output(state):
        global output_layer
        a, b, _ = state
        a, b = symmetry(a, b)
        output_layer = StaticLayer(a, b)
        # print(molecule_output(output_layer))

    layer = EditableLayer()
    layer.add_subscriber(output)
    [C, H] = layer.add_atoms([Atom("C", [0, 0, 1]), Atom("H", [1, 0, 1.2])])
    layer.set_bond(C, H, 1.0)
    layer.select([H])
    layer.set_element_selected("Cl")

    layer = EditableLayer(output_layer)
    [C1, C2] = py_.filter(
        layer.atom_ids, lambda atom_id: layer.atoms[atom_id].element == "C"
    )
    layer.set_bond(C1, C2, 1.0)
    print(molecule_output(layer))
