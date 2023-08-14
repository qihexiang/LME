from libs.molecule_text import atoms_bonds_to_mol2
from layers.EditableLayer import EditableLayer, Substitute
from fs.osfs import OSFS
from pydash import py_

text = OSFS("example_data").readtext("NH3.mol2")
mole = EditableLayer.from_mol2(text)
N_id = py_.find(mole.atom_ids, lambda atom_id: mole.atoms[atom_id].element == "N")
H_id = py_.find(mole.atom_ids, lambda atom_id: mole.atoms[atom_id].element == "H")
substitute = Substitute.build_from(mole, H_id, N_id)
mole.add_substitute(substitute, H_id, N_id)
mol2 = atoms_bonds_to_mol2(mole.atoms, mole.bonds)

OSFS("output").writetext("editable.mol2", mol2)