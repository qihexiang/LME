from typing import Any
from pydash import py_
from fs.osfs import OSFS
from layers.EditableLayer import Substitute, EditableLayer
from libs.molecule_text import atoms_bonds_to_mol2
from posixpath import join, isabs, relpath, dirname
from openbabel import openbabel
from openbabel import pybel
import re

atom_section_re = re.compile("@<TRIPOS>ATOM\n((.*\n)*)\n@")

class AddSubsititute:
    def __init__(self, options, metas) -> None:
        self.replace = options["replace"]
        self.substitutes_lib = py_.map(metas["substitutes"], lambda libpath: join(metas["rootDirectory"], libpath)) + [join(dirname(__file__.replace("\\", "/")), "..", "Substitutes")]
        self.substitutes_lib = [OSFS(directory) for directory in self.substitutes_lib]
        self.substitutes = options["substitutes"]
        if self.substitutes == "all":
            self.substitutes = self.__all_substitutes()


    def __all_substitutes(self):
        return (
            py_
                .chain(self.substitutes_lib)
                .map(lambda directory: directory.listdir("."))
                .flatten()
                .filter(lambda filename: filename.endswith(".mol2"))
                .uniq()
                .map(lambda filename: filename[0:-5])
                .value()
        )
    
    def load_substitute(self, substitute_name):
        target_filename = f"{substitute_name}.mol2"
        for directory in self.substitutes_lib:
            try:
                mol2 = directory.readtext(target_filename)
                return Substitute.from_mol2(mol2)        
            except Exception:
                pass
        raise FileNotFoundError(f"No {target_filename} found in given directories")
    
    def __call__(self, target, names) -> Any:
        def generate_for_subsitite(sub_entry):
            editable = EditableLayer(target)
            indexes = py_.map(self.replace, lambda names: py_.map(names, lambda name: editable.find_with_classname(name)[0]))
            if type(sub_entry) == str:
                tag_name = sub_entry
                sub_name = sub_entry
            else:
                tag_name = sub_entry["name"]
                sub_name = sub_entry["substitute"]
            substitute = self.load_substitute(sub_name)
            for [entry_idx,center_idx] in indexes:
                editable.add_substitute(substitute, center_idx, entry_idx)
            return editable.to_static_layer(), tag_name
        return [generate_for_subsitite(sub_entry) for sub_entry in self.substitutes]

class AtomModify:
    def __init__(self, options, metas) -> None:
        self.tasks = options

    def __call__(self, target, names) -> Any:
        editable = EditableLayer(target)
        for task in self.tasks:
            editable.deselect_all()
            editable.select(editable.find_with_classname(task["select"]))
            if "element" in task:
                editable.set_element_selected(task["element"])
            if "translation" in task:
                editable.translation_selected(task["translation"])
        return editable.to_static_layer(), ""

class BondModify:
    def __init__(self, options, metas) -> None:
        self.tasks = options

    def __call__(self, target, names) -> Any:
        editable = EditableLayer(target)
        connects = [[editable.find_with_classname(n1)[0], editable.find_with_classname(n2)[0], bond_type] for [n1, n2, bond_type] in self.tasks]
        
        for [a, b, bond] in connects:
            editable.set_bond(a, b, bond)
        
        return editable.to_static_layer(), ""

class ImportStructure:
    def __init__(self, options, metas) -> Any:
        rootDirectory = metas["rootDirectory"]
        targetPath = join(*options["filepath"])
        targetPath = relpath(targetPath, rootDirectory) if isabs(rootDirectory) else targetPath
        editable = EditableLayer.from_mol2(OSFS(rootDirectory).readtext(targetPath))
        self.atoms = editable.atoms
        self.bonds = editable.bonds
        self.rotation = options.get("rotation")
        self.translation = options.get("translation")

    def __call__(self, target, names) -> Any:
        editable = EditableLayer(target)
        editable.import_atoms_bonds(self.atoms, self.bonds)
        if self.rotation is not None and self.rotation is not False:
            editable.rotation_selected(self.rotation["axis"], self.rotation["center"], self.rotation["angle"])

        if self.translation is not None and self.translation is not False:
            editable.translation_selected(self.translation)

        return editable.to_static_layer(), ""

def debug_output(a):
    print(a)
    return a

class Output:
    def __init__(self, options, metas) -> None:
        self.rootDirectory = OSFS(metas["rootDirectory"])
        self.filenamePattern = join(*options["pattern"])
        self.clean = options.get("clean")
        self.freeze = self.clean.get("freeze")
        self.freeze = self.freeze if self.freeze is not None else []
    
    def __clean__(self, mol2, tags):
        atoms_section = atom_section_re.findall(mol2)
        if self.clean is None:
            pass
        rules = self.clean["rules"]
        for [tag, value] in rules:
            if(tags.get(tag) == value):
                print(tags.get(tag))
                obmol = pybel.readstring("mol2", mol2).OBMol
                constraints = openbabel.OBFFConstraints()
                to_freezes = (py_.chain(atoms_section[0][0].split("\n"))
                             .map(lambda line: py_.filter(line.split(" "), lambda token: token != ""))
                             .filter(lambda line: len(line) != 0)
                            #  .map(debug_output)
                             .filter(lambda line: line[1] in self.freeze)
                             .map(lambda line: int(line[0]))
                             .value()
                            )
                for to_freeze in to_freezes:
                    print(to_freeze)
                    constraints.AddAtomConstraint(to_freeze)

                ff = openbabel.OBForceField.FindForceField("UFF")
                ff.Setup(obmol, constraints)
                print("Init:", ff.Energy())
                ff.ConjugateGradients(self.clean["steps"])
                print("Optimed", ff.Energy())
                ff.GetConformers(obmol)
                output = pybel.Molecule(obmol)
                output = output.write("mol2")
                return output
        pass

    def __write__(self, target, content):
        target_dir = dirname(target)
        if not self.rootDirectory.exists(target_dir):
            self.rootDirectory.makedirs(target_dir)
        return self.rootDirectory.writetext(target, content)

    
    def __call__(self, item, names) -> Any:
        mol2 = atoms_bonds_to_mol2(item.atoms, item.bonds)
        cleaned = self.__clean__(mol2, names)
        filename = self.filenamePattern
        for stage_name in names:
            filename = filename.replace(f"{{{stage_name}}}", names[stage_name])
        self.__write__(filename, mol2)
        if cleaned is not None:
            self.__write__(f"{filename}.cleaned.mol2", cleaned)
        return item, "output"

default_runners = {
    "add_substitute": (AddSubsititute, True),
    "modify_atom": (AtomModify, False),
    "modify_bond": (BondModify, False),
    "import": (ImportStructure, False),
    "output": (Output, False)
}