from pydash import py_

def molecule_text(target) -> str:
    output = ""
    output += "atoms:\n"
    atom_ids = py_.filter(
        target.atom_ids, lambda atom_id: target.atoms[atom_id] is not None
    )
    atoms = target.atoms
    atoms = py_.map(atom_ids, lambda atom_id: atoms[atom_id])
    for i, atom in enumerate(atoms):
        output += f"{i+1} {atom}\n"

    bond_ids = target.bond_ids
    bonds = target.bonds
    output += "bonds:\n"
    for bond_id in bond_ids:
        a_idx, b_idx = py_.find_index(
            atom_ids, lambda atom_id: atom_id == bond_id.a
        ), py_.find_index(atom_ids, lambda atom_id: atom_id == bond_id.b)
        if -1 not in [a_idx, b_idx]:
            output += f"{a_idx + 1} {b_idx + 1} {bonds[bond_id]}\n"

    return output
