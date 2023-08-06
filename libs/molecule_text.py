
from pydash import py_

from libs.Atom import Atom

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

def mol2_get_section(lines, section_name):
    section_start = py_.find_index(lines, lambda line: line == f"@<TRIPOS>{section_name}")
    section_offset = py_.find_index(lines[section_start + 1:], lambda line: line.startswith("@<TRIPOS>")) 
    if section_offset == -1:
        return lines[section_start + 1:]
    return lines[section_start + 1:section_start + 1 + section_offset]

def mol2_to_atom(line):
    [atom_id, atom_name, x, y, z, element] = py_.filter(line.split(" "), lambda item: item != "")[0:6]
    atom = Atom(element, [float(x), float(y), float(z)], atom_name)
    return int(atom_id), atom

def mol2_to_bond(line):
    [a, b, order] = py_.filter(line.split(" "), lambda item: item != "")[1:4]
    return int(a), int(b), float(order)

def atoms_bonds_from_mol2(text):
    lineEnd = "\r\n" if "\r" in text else "\n"
    lines = text.split(lineEnd)
    
    atoms_lines = py_.filter(mol2_get_section(lines, "ATOM"), lambda line: line != "")
    bonds_lines = py_.filter(mol2_get_section(lines, "BOND"), lambda line: line != "")

    return atoms_lines, bonds_lines    
