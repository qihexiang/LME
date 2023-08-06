from fs.osfs import OSFS
from EditableLayer import EditableLayer

directory = OSFS("example_data")
text = directory.readtext("2_NH3_AV.mol2")

editable = EditableLayer.from_mol2(text)

print(editable)