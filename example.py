from libs.molecule_text import atoms_bonds_to_mol2
from EditableLayer import EditableLayer, Polymer
from fs.osfs import OSFS
from pydash import py_

text = OSFS("example_data").readtext("2_NH3_AV.mol2")
mole = EditableLayer.from_mol2(text)
N_id = py_.find(mole.atom_ids, lambda atom_id: mole.atoms[atom_id].element == "N")
H_id = py_.find(mole.atom_ids, lambda atom_id: mole.atoms[atom_id].element == "H")
polymer = Polymer.build_from(mole, H_id, N_id)
mole.add_polymer(polymer, H_id, N_id)
mol2 = atoms_bonds_to_mol2(mole.atoms, mole.bonds)

OSFS("example_data").writetext("editable.mol2", mol2)