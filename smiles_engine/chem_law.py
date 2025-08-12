"""chem_law.py — Infers chemical details from a raw molecular graph."""

# Simple valence dictionary: typical max valence electrons for common elements
VALENCE = {
    "H": 1,
    "C": 4,
    "N": 3,
    "O": 2,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
}

def infer_implicit_hydrogens(mol):
    """
    Calculate implicit hydrogens for each atom in the molecule.
    Implicit hydrogens are hydrogens not explicitly shown but
    implied by valence rules and bonding.

    Args:
        mol: CompoundGraph object containing atoms and bonds

    Sets:
        atom.implicit_h attribute for each atom, number of implicit hydrogens.
    """
    for atom in mol.atoms:
        # Maximum valence for this atom type (default 0 if unknown)
        max_valence = VALENCE.get(atom.symbol, 0)
        # Sum total bond order of bonds attached to this atom
        current_bond_order = sum(
            bond.order for bond in mol.bonds
            if bond.a1 == atom.idx or bond.a2 == atom.idx
        )
        # Implicit hydrogens = max valence minus bonded electrons
        atom.implicit_h = max_valence - current_bond_order

def get_default_bond_length(atom1, atom2):
    """
    Return a default bond length between two atoms.
    This is a placeholder using a fixed length (1.5 units),
    but could be improved by using atom types or bond order.

    Args:
        atom1, atom2: Atom objects

    Returns:
        float: bond length in arbitrary units
    """
    return 1.5

def get_bond_angle_sequence(num_bonds):
    """
    Given the number of bonds around an atom, return a list of angles
    (in degrees) to position child atoms around it for 2D layout.

    Args:
        num_bonds (int): number of bonds (neighbors) to place

    Returns:
        List[float]: angles in degrees for each bond placement
    """
    if num_bonds == 0:
        return []
    elif num_bonds == 1:
        return [0]
    elif num_bonds == 2:
        # Slightly spread out for two bonds
        return [-30, 30]
    elif num_bonds == 3:
        # Triangular-ish distribution
        return [-60, 0, 60]
    elif num_bonds == 4:
        # Approximate tetrahedral planar projection
        return [-90, -30, 30, 90]
    else:
        # For more bonds, spread evenly in 360 degrees
        step = 360 / num_bonds
        return [i * step for i in range(num_bonds)]

# Color lookup for atom types in rendering
ATOM_COLORS = {
    "H": "lightgray",
    "C": "black",
    "N": "blue",
    "O": "red",
    "F": "green",
    "Cl": "green",
    "Br": "brown",
    "I": "purple",
}

def get_atom_color(symbol):
    """
    Return a color string for the given atom symbol.
    Defaults to gray for unknown elements.

    Args:
        symbol (str): atom element symbol, e.g. "C", "O"

    Returns:
        str: color name or hex string
    """
    return ATOM_COLORS.get(symbol, "gray")

def get_bond_style(bond):
    """
    Return a style dictionary for bond rendering based on bond properties.

    Args:
        bond: Bond object with 'order' and 'aromatic' attributes

    Returns:
        dict: style keys 'width' (line thickness) and 'dash' (line style)
    """
    if bond.aromatic:
        # Aromatic bonds are dotted and thin
        return {"width": 2, "dash": "dot"}
    elif bond.order == 1:
        # Single bond, normal thickness solid line
        return {"width": 2, "dash": "solid"}
    elif bond.order == 2:
        # Double bond thicker solid line
        return {"width": 4, "dash": "solid"}
    elif bond.order == 3:
        # Triple bond even thicker solid line
        return {"width": 6, "dash": "solid"}
    else:
        # Fallback style
        return {"width": 2, "dash": "solid"}
