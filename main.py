from editable_layer import EditableLayer
from symmetry_layers import InverseLayer, MirrorLayer
from lib import Atom

if __name__ == "__main__":
    layer = EditableLayer()
    inv_layer = InverseLayer()
    mirror_layer = MirrorLayer([1.0, 0.0, 0.0])

    def layers(a, b):
        a, b = inv_layer(a, b)
        a, b = mirror_layer(a, b)
        return a, b

    layer.add_subscribers(lambda a, b: print(layers(a, b)))
    # C = layer.add_atom(Atom("C", [0,0,0]))
    # H = layer.add_atom(Atom("H", [1,0,0]))
    # layer.set_bonds((C,H), 1.0)
    # H = layer.add_atom(Atom("H", [1,1,0]))
    # layer.set_bonds((C,H), 1.0)
    # H = layer.add_atom(Atom("H", [0,0,1]))
    # layer.set_bonds((C,H), 1.0)
    # layer.replace_atom(C, "N")
    C = layer.add_atom(Atom("C", [1, 1, 0]))
