from copy import deepcopy
from libs.atoms_bonds_loader import atoms_bonds_loader
from libs.constants import PRODUCTION
from libs.molecule_text import molecule_text
from layers.SymmetryLayers import SymmetryLayer
from layers.UtilLayers import DedupLayer
from pydash import py_


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
                "transformers": py_.map(
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
    def contains(self):
        if PRODUCTION:
            return self.__contains
        else:
            return deepcopy(self.__contains)

    @property
    def atoms(self):
        return {
            atom_id: self.__atoms[atom_id]
            for atom_id in py_.filter(
                self.__atoms.keys(), lambda atom_id: self.__atoms[atom_id] is not None
            )
        }

    @property
    def bonds(self):
        existed_atoms = self.atom_ids
        return {
            bond_id: self.__bonds[bond_id]
            for bond_id in py_.filter(
                self.__bonds.keys(),
                lambda bond_id: bond_id.a in existed_atoms
                and bond_id.b in existed_atoms
                and self.__bonds[bond_id] is not None,
            )
        }

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
    
    @staticmethod
    def __extract_transformer(data):
        if(data["type"].startswith("symmetry.")):
            return SymmetryLayer.from_dict(data)
        if(data["type"] == "dedup"):
            return DedupLayer(data["eps"])
        raise ValueError("Invalid transformer inforamtion")

    def extract(self):
        from EditableLayer import EditableLayer
        editable = EditableLayer(load=self.contains["editable"])
        transformers = py_.map(self.contains["transformers"], StaticLayer.__extract_transformer)
        return editable, transformers

    def to_editable_layer(self):
        from layers.EditableLayer import EditableLayer
        editable = EditableLayer()
        editable.import_atoms_bonds(self.atoms, self.bonds)
        editable.deselect_all()
        return editable
