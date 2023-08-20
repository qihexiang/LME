from typing import Any
from pydash import py_
from fs.osfs import OSFS
from os.path import join, dirname
from layers.EditableLayer import Substitute, EditableLayer
from libs.molecule_text import atoms_bonds_to_mol2, atoms_bonds_from_mol2
from posixpath import join, isabs, relpath

class AddSubsititute:
    def __init__(self, options, metas) -> None:
        self.replace = options["replace"]
        self.substitutes_lib = py_.map(metas["substitutes"], lambda libpath: join(metas["rootDirectory"], libpath)) + [join(dirname(__file__), "..", "Substitutes")]
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

class BondModify:
    def __init__(self, options, metas) -> None:
        self.tasks = options

    def __call__(self, target, names) -> Any:
        editable = EditableLayer(target)
        connects = py_.map(self.tasks, lambda names: py_.map(names, lambda target: [editable.find_with_classname(target[0])[0], editable.find_with_classname(target[1])[0], target[2]]))
        
        for [a, b, bond] in connects:
            editable.set_bond(a, b, bond)
        
        return editable.to_static_layer(), ""

class ImportStructure:
    def __call__(self, options, metas) -> Any:
        rootDirectory = metas["rootDirectory"]
        targetPath = join(*options["filepath"])
        targetPath = relpath(targetPath, rootDirectory) if isabs(rootDirectory) else targetPath
        self.atoms, self.bonds = atoms_bonds_from_mol2(OSFS(rootDirectory).readtext(targetPath))
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

class Output:
    def __init__(self, options, metas) -> None:
        self.rootDirectory = OSFS(metas["rootDirectory"])
        self.filenamePattern = options["pattern"]

    def __write__(self, target, content):
        target_dir = dirname(target)
        if not self.rootDirectory.exists(target_dir):
            self.rootDirectory.makedirs(target_dir)
        return self.rootDirectory.writetext(target, content)

    
    def __call__(self, item, names) -> Any:
        mol2 = atoms_bonds_to_mol2(item.atoms, item.bonds)
        filename = self.filenamePattern
        for stage_name in names:
            filename = filename.replace(f"{{{stage_name}}}", names[stage_name])
        self.__write__(filename, mol2)
        return item, "output"

default_runners = {
    "add_substitute": (AddSubsititute, True), 
    "modify_bond": (BondModify, False),
    "import": (ImportStructure, False),
    "output": (Output, False)
}