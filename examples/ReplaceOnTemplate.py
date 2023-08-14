from layers.EditableLayer import EditableLayer, Substitute
from fs.osfs import OSFS
from pydash import py_

from layers.StaticLayer import StaticLayer
from libs.molecule_text import atoms_bonds_to_mol2

directory = OSFS("Substitutes")

substitutes = [(entry, Substitute.from_mol2(directory.readtext(entry))) for entry in directory.listdir("./")]

template = OSFS("example_data").readtext("template.mol2")
template_layer = StaticLayer(contains=(EditableLayer.from_mol2(template), []))

to_replace = [
    ("P11", ["H33", "H36"]),
    ("P12", ["H34", "H35"])
]

diff_layer_groups = {}

for (entry_name, center_names) in to_replace:
    diff_layers = []
    for (_, substitute) in substitutes:
        layer = EditableLayer(template_layer)
        entry_index = layer.find_with_classname(entry_name)[0]
        for center_name in center_names:      
            center_index = layer.find_with_classname(center_name)[0]
            layer.add_substitute(substitute, center_index, entry_index)
        diff_layers.append(layer.detach_diff_layer())
    diff_layer_groups[entry_name] = diff_layers

alldiff = py_.flatten(py_.map(diff_layer_groups["P11"], lambda diff1: py_.map(diff_layer_groups["P12"], lambda diff2: [diff1, diff2])))

output = OSFS("output")

for i, diffs in enumerate(alldiff):
    layer = StaticLayer(contains=(template_layer, diffs))
    mol2 = atoms_bonds_to_mol2(layer.atoms, layer.bonds)
    output.writetext(f"{i}.mol2", mol2)