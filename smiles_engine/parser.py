"""
A minimal SMILES (Simplified Molecular Input Line Entry System) string parser 
that converts chemical structures into an internal graph representation without 
using external chemistry libraries.
"""

import re

class Atom:
    def __init__(self, symbol, idx, isotope=None, charge=0, explicit_h=0):
        """
        Represents an atom in the molecule.

        Args:
            symbol (str): Chemical symbol, e.g. 'C', 'O', 'N'.
            idx (int): Unique index of the atom in the molecule.
            isotope (int, optional): Isotope number, e.g. 13 for Carbon-13.
            charge (int): Formal charge on the atom.
            explicit_h (int): Number of explicit hydrogens attached.
        """
        self.symbol = symbol
        self.idx = idx
        self.isotope = isotope
        self.charge = charge
        self.explicit_h = explicit_h

class Bond:
    def __init__(self, atom1, atom2, order=1, aromatic=False):
        """
        Represents a bond between two atoms.

        Args:
            atom1 (int): Index of the first atom.
            atom2 (int): Index of the second atom.
            order (int): Bond order (1=single, 2=double, 3=triple).
            aromatic (bool): Whether the bond is aromatic.
        """
        self.a1 = atom1
        self.a2 = atom2
        self.order = order
        self.aromatic = aromatic

class CompoundGraph:
    def __init__(self):
        """
        Represents the full compound as a graph of atoms and bonds.
        """
        self.atoms = []       # List of Atom objects
        self.bonds = []       # List of Bond objects
        self.components = []  # List of connected components (lists of atom indices)
        self.rings = []       # List of sets of atom indices forming rings

def parse_smiles_fragment(smiles_fragment, start_idx=0):
    """
    Parse a fragment of a SMILES string (without '.' separators).

    Args:
        smiles_fragment (str): SMILES string fragment.
        start_idx (int): Index offset for atoms in this fragment.

    Returns:
        atoms (list): List of Atom objects parsed.
        bonds (list): List of Bond objects parsed.
        last_idx (int): Next available atom index after this fragment.
        ring_atoms_map (dict): Maps ring numbers to sets of atoms forming rings.
    """
    atoms = []
    bonds = []

    prev_atom = None              # Tracks the previous atom for bonding
    branch_stack = []            # Stack to handle branching '(' and ')'
    ring_map = {}                # Maps ring number to atom index where ring opened
    ring_atoms_map = {}          # Maps ring number to set of atom indices in that ring

    bond_order_map = {"-": 1, "=": 2, "#": 3, ":": 1}  # Bond symbols to order/aromatic

    pending_bond_order = 1       # Bond order for next bond (default single)
    pending_aromatic = False     # Whether next bond is aromatic

    i = 0
    length = len(smiles_fragment)
    while i < length:
        ch = smiles_fragment[i]

        if ch == "(":
            # Start of branch: push current atom to stack
            branch_stack.append(prev_atom)
            i += 1
            continue

        if ch == ")":
            # End of branch: pop atom from stack to continue from
            prev_atom = branch_stack.pop()
            i += 1
            continue

        if ch in bond_order_map:
            # Bond symbol encountered; set pending bond order/aromatic flag
            pending_bond_order = bond_order_map[ch]
            pending_aromatic = (ch == ":")
            i += 1
            continue

        # Handle ring closure digits and multi-digit rings with '%'
        if ch.isdigit() or ch == '%':
            if ch == '%':
                # Multi-digit ring number after '%', e.g. '%12'
                ring_num_str = smiles_fragment[i+1:i+3]
                ring_num = int(ring_num_str)
                i += 3
            else:
                # Single digit ring number
                ring_num = int(ch)
                i += 1

            if ring_num in ring_map:
                # Closing a ring: create bond between prev_atom and stored atom
                other_idx = ring_map.pop(ring_num)
                bonds.append(Bond(prev_atom.idx, other_idx, pending_bond_order, aromatic=pending_aromatic))

                # Track atoms involved in this ring number
                if ring_num not in ring_atoms_map:
                    ring_atoms_map[ring_num] = set()
                ring_atoms_map[ring_num].add(prev_atom.idx)
                ring_atoms_map[ring_num].add(other_idx)
            else:
                # Opening a ring: store current atom index for this ring number
                ring_map[ring_num] = prev_atom.idx

            # Reset bond state after ring closure/opening
            pending_bond_order = 1
            pending_aromatic = False
            continue

        if ch == "[":
            # Parsing an atom specification inside brackets, e.g. '[13CH3+]' or '[O-]'
            j = smiles_fragment.index("]", i)
            bracket_content = smiles_fragment[i + 1:j]

            # Regex to parse isotope, element symbol, hydrogens, charge
            match = re.match(r"(?:(\d+))?([A-Z][a-z]?)(H\d?)?([+-]\d?)?", bracket_content)
            isotope = int(match.group(1)) if match.group(1) else None
            symbol = match.group(2)
            h_part = match.group(3)
            charge_part = match.group(4)

            explicit_h = 0
            if h_part:
                explicit_h = 1 if len(h_part) == 1 else int(h_part[1:])

            charge = 0
            if charge_part:
                # Handle simple +/- or numeric charges
                charge = 1 if len(charge_part) == 1 and charge_part[0] == "+" else \
                         -1 if len(charge_part) == 1 and charge_part[0] == "-" else \
                         int(charge_part)

            # Create atom and add to list
            atom = Atom(symbol, start_idx + len(atoms), isotope=isotope, charge=charge, explicit_h=explicit_h)
            atoms.append(atom)

            # Create bond from previous atom if exists
            if prev_atom is not None:
                bonds.append(Bond(prev_atom.idx, atom.idx, pending_bond_order, aromatic=pending_aromatic))

            prev_atom = atom
            pending_bond_order = 1
            pending_aromatic = False
            i = j + 1
            continue

        if ch.islower():
            # Aromatic atoms are lowercase in SMILES, convert to uppercase symbol and mark aromatic bonds
            symbol = ch.upper()
            atom = Atom(symbol, start_idx + len(atoms))
            atoms.append(atom)
            if prev_atom is not None:
                bonds.append(Bond(prev_atom.idx, atom.idx, pending_bond_order, aromatic=True))
            prev_atom = atom
            pending_bond_order = 1
            pending_aromatic = False
            i += 1
            continue

        if ch.isupper():
            # Normal atoms, possibly two-letter element symbol if next char is lowercase
            if i + 1 < length and smiles_fragment[i + 1].islower():
                symbol = ch + smiles_fragment[i + 1]
                i += 1
            else:
                symbol = ch

            atom = Atom(symbol, start_idx + len(atoms))
            atoms.append(atom)

            # Create bond from previous atom, non-aromatic by default
            if prev_atom is not None:
                bonds.append(Bond(prev_atom.idx, atom.idx, pending_bond_order, aromatic=False))

            prev_atom = atom
            pending_bond_order = 1
            pending_aromatic = False
            i += 1
            continue

        # If character unrecognized, just skip
        i += 1

    last_idx = start_idx + len(atoms)
    return atoms, bonds, last_idx, ring_atoms_map

def parse_smiles(smiles: str) -> CompoundGraph:
    """
    Parse a full SMILES string possibly containing multiple disconnected fragments separated by '.'.

    Args:
        smiles (str): Full SMILES string.

    Returns:
        CompoundGraph: Parsed compound graph with atoms, bonds, components, and rings.
    """
    compound = CompoundGraph()

    fragments = smiles.split(".")
    current_index = 0

    combined_ring_map = {}

    for frag in fragments:
        atoms, bonds, next_index, ring_atoms_map = parse_smiles_fragment(frag, current_index)

        # Append parsed atoms and bonds
        compound.atoms.extend(atoms)
        compound.bonds.extend(bonds)

        # Track connected component indices for this fragment
        component_indices = list(range(current_index, next_index))
        compound.components.append(component_indices)

        # Merge rings from this fragment into combined ring map
        for ring_num, atom_set in ring_atoms_map.items():
            if ring_num not in combined_ring_map:
                combined_ring_map[ring_num] = set()
            combined_ring_map[ring_num].update(atom_set)

        current_index = next_index

    # Convert combined ring map to list of ring atom sets
    compound.rings = list(combined_ring_map.values())

    return compound
