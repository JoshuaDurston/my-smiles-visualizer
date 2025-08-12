"""
Generates 2D coordinates for atoms in a molecule for visualization purposes,
based on the molecular connectivity described in a CompoundGraph.

Enhancements:
- Places aromatic and other rings as regular polygons with equalized bond lengths.
- Uses stronger spring forces for ring bonds to preserve their regularity during relaxation.
"""

import math
from collections import deque
from .chem_law import get_default_bond_length


def get_bond_angle_sequence(num_bonds, hybridization="sp3"):
    """
    Return angles (degrees) evenly spaced around a central atom based on its bonding geometry.
    Uses common idealized angles for sp, sp2, sp3 hybridizations; otherwise divides 360° evenly.
    
    Args:
        num_bonds (int): Number of bonds attached to the atom.
        hybridization (str): Type of hybridization ("sp", "sp2", "sp3", etc.)

    Returns:
        List[float]: Angles in degrees for placing bonds around the atom.
    """
    if num_bonds == 0:
        return []
    elif num_bonds == 1:
        return [0]  # Single bond points arbitrarily at 0°
    elif num_bonds == 2:
        # Linear (sp) prefers 180°, otherwise ~120°
        return [0, 180 if hybridization == "sp" else 120]
    elif num_bonds == 3:
        # Trigonal planar (sp2/sp3) or right-angle (default fallback)
        return [0, 120, 240] if hybridization in ("sp2", "sp3") else [0, 90, 180]
    elif num_bonds == 4:
        # Tetrahedral angles approx (109.5° apart)
        return [0, 109.5, 219, 328.5]
    else:
        # For >4 bonds, spread evenly around 360°
        step = 360 / num_bonds
        return [i * step for i in range(num_bonds)]


def place_ring_atoms(compound_graph, ring_atom_indices, x_offset=0, y_offset=0):
    """
    Place atoms in a ring evenly spaced on a circle, creating a regular polygon.

    Args:
        compound_graph: Molecule graph containing atoms and bonds.
        ring_atom_indices (List[int]): Indices of atoms forming the ring.
        x_offset (float): Horizontal offset to position this ring.
        y_offset (float): Vertical offset to position this ring.

    Returns:
        dict: Mapping of atom indices to (x, y) coordinates.
    """
    positions = {}
    n = len(ring_atom_indices)
    if n == 0:
        return positions  # Empty ring, nothing to place

    # Calculate average ideal bond length for all bonds in the ring
    lengths = []
    for i in range(n):
        a1 = compound_graph.atoms[ring_atom_indices[i]]
        a2 = compound_graph.atoms[ring_atom_indices[(i + 1) % n]]  # next atom in ring, wrap around
        lengths.append(get_default_bond_length(a1, a2))
    bond_length = sum(lengths) / len(lengths)

    # Calculate radius of circumscribed circle for regular polygon
    radius = bond_length / (2 * math.sin(math.pi / n))

    # Position ring center with offset + radius (to keep positive coords)
    center_x = x_offset + radius
    center_y = y_offset + radius

    # Place each ring atom at equal angles around circle
    for i, atom_idx in enumerate(ring_atom_indices):
        angle = 2 * math.pi * i / n  # Angle in radians
        x = center_x + radius * math.cos(angle)
        y = center_y + radius * math.sin(angle)
        positions[atom_idx] = (x, y)

    return positions


def generate_2d_layout_for_component(compound_graph, component_atom_indices, x_offset=0, y_offset=0):
    """
    Generate 2D layout for a connected component (fragment) of the molecule.
    First places any ring atoms using regular polygon placement,
    then places remaining atoms via BFS traversal positioning.

    Args:
        compound_graph: Molecule graph with atoms, bonds, rings.
        component_atom_indices (List[int]): Atoms in this connected component.
        x_offset, y_offset (float): Offsets for positioning this component.

    Returns:
        dict: Atom index -> (x, y) coordinates for this component.
    """
    positions = {}
    visited = set()

    # Build neighbors map restricted to this component for quick lookups
    neighbors_map = {idx: [] for idx in component_atom_indices}
    for bond in compound_graph.bonds:
        if bond.a1 in component_atom_indices and bond.a2 in component_atom_indices:
            neighbors_map[bond.a1].append(bond.a2)
            neighbors_map[bond.a2].append(bond.a1)

    # Place all ring atoms first using place_ring_atoms, then mark as visited
    for ring_atoms in compound_graph.rings:
        if all(atom in component_atom_indices for atom in ring_atoms):
            ring_pos = place_ring_atoms(compound_graph, list(ring_atoms), x_offset, y_offset)
            positions.update(ring_pos)
            visited.update(ring_atoms)

    queue = deque()

    # Start BFS from ring atoms: enqueue their unvisited neighbors with initial bond angles
    for root_atom_idx in visited:
        neighbors = [n for n in neighbors_map[root_atom_idx] if n not in visited]
        bond_angles = get_bond_angle_sequence(len(neighbors))
        for i, nb in enumerate(neighbors):
            queue.append((nb, root_atom_idx, bond_angles[i]))

    # If no rings in component, pick first atom and start layout there
    if not visited:
        root_atom_idx = component_atom_indices[0]
        positions[root_atom_idx] = (x_offset, y_offset)
        visited.add(root_atom_idx)

        neighbors = [n for n in neighbors_map[root_atom_idx] if n not in visited]
        bond_angles = get_bond_angle_sequence(len(neighbors))
        for i, nb in enumerate(neighbors):
            queue.append((nb, root_atom_idx, bond_angles[i]))

    # BFS loop: place atoms relative to their parents using bond lengths and angles
    while queue:
        atom_idx, parent_idx, parent_angle = queue.popleft()
        if atom_idx in visited:
            continue  # Already positioned

        # Calculate ideal bond length between parent and current atom
        bond_length = get_default_bond_length(
            compound_graph.atoms[parent_idx],
            compound_graph.atoms[atom_idx]
        )

        px, py = positions[parent_idx]
        rad = math.radians(parent_angle)

        # Position this atom based on parent's position + bond vector
        ax = px + bond_length * math.cos(rad)
        ay = py + bond_length * math.sin(rad)

        positions[atom_idx] = (ax, ay)
        visited.add(atom_idx)

        # Enqueue children neighbors with angles offset to avoid overlap
        child_neighbors = [n for n in neighbors_map[atom_idx] if n not in visited]
        bond_angles = get_bond_angle_sequence(len(child_neighbors))

        # Start angles offset by 210° to spread out bonds nicely
        start_angle = (parent_angle + 180 + 30) % 360

        for i, nb in enumerate(child_neighbors):
            angle = (start_angle + bond_angles[i]) % 360
            queue.append((nb, atom_idx, angle))

    return positions


def relax_positions(positions, compound_graph, iterations=50, spring_k=0.1, repulsion_k=0.05):
    """
    Apply simple 2D force-directed relaxation to smooth atom positions:
    - Spring forces try to keep bonds at ideal lengths.
    - Repulsive forces avoid atom overlap.
    - Ring bonds get stronger springs and equal target lengths for regularity.

    Args:
        positions (dict): Initial atom positions {atom_idx: (x, y)}.
        compound_graph: Molecule graph to access bonds and atoms.
        iterations (int): Number of relaxation steps.
        spring_k (float): Spring force constant.
        repulsion_k (float): Repulsion force constant.

    Returns:
        dict: Relaxed atom positions {atom_idx: (x, y)}.
    """
    # Convert tuples to mutable lists for position updates
    pos = {k: list(v) for k, v in positions.items()}

    # Identify ring bonds and calculate average ideal bond length per ring
    ring_bonds = set()
    ring_target_lengths = {}
    for ring in compound_graph.rings:
        if len(ring) < 3:
            continue  # Not a proper ring
        # Create pairs of adjacent atoms in ring (with wrap-around)
        ring_pairs = [(ring[i], ring[(i + 1) % len(ring)]) for i in range(len(ring))]
        lengths = [
            get_default_bond_length(compound_graph.atoms[a1], compound_graph.atoms[a2])
            for a1, a2 in ring_pairs
        ]
        avg_len = sum(lengths) / len(lengths)
        for a1, a2 in ring_pairs:
            key = tuple(sorted((a1, a2)))  # Normalize bond key
            ring_bonds.add(key)
            ring_target_lengths[key] = avg_len  # Use same target length for all ring bonds

    for _ in range(iterations):
        forces = {i: [0.0, 0.0] for i in pos}

        # Calculate spring forces for all bonds
        for bond in compound_graph.bonds:
            if bond.a1 in pos and bond.a2 in pos:
                key = tuple(sorted((bond.a1, bond.a2)))
                if key in ring_bonds:
                    ideal_len = ring_target_lengths[key]
                    k_val = spring_k * 5.0  # Stronger spring for ring bonds to keep ring shape
                else:
                    ideal_len = get_default_bond_length(
                        compound_graph.atoms[bond.a1],
                        compound_graph.atoms[bond.a2]
                    )
                    k_val = spring_k

                dx = pos[bond.a2][0] - pos[bond.a1][0]
                dy = pos[bond.a2][1] - pos[bond.a1][1]
                dist = math.sqrt(dx*dx + dy*dy) or 0.01  # Avoid division by zero
                diff = dist - ideal_len
                force = k_val * diff
                fx = force * dx / dist
                fy = force * dy / dist

                # Apply equal and opposite forces on bonded atoms
                forces[bond.a1][0] += fx
                forces[bond.a1][1] += fy
                forces[bond.a2][0] -= fx
                forces[bond.a2][1] -= fy

        # Calculate repulsive forces to avoid atom overlap
        atom_list = list(pos.keys())
        for i in range(len(atom_list)):
            for j in range(i + 1, len(atom_list)):
                a, b = atom_list[i], atom_list[j]
                dx = pos[b][0] - pos[a][0]
                dy = pos[b][1] - pos[a][1]
                dist_sq = dx*dx + dy*dy
                if dist_sq < 0.01:
                    dist_sq = 0.01  # Prevent extreme forces on very close atoms
                force = repulsion_k / dist_sq
                dist = math.sqrt(dist_sq)
                fx = force * dx / dist
                fy = force * dy / dist

                # Repulsive force pushes atoms apart
                forces[a][0] -= fx
                forces[a][1] -= fy
                forces[b][0] += fx
                forces[b][1] += fy

        # Update positions by applying forces
        for i in pos:
            pos[i][0] += forces[i][0]
            pos[i][1] += forces[i][1]

    # Convert back to tuples for immutability and return
    return {k: tuple(v) for k, v in pos.items()}


def generate_2d_layout(compound_graph):
    """
    Generate 2D coordinates for all atoms in the molecule.
    If multiple disconnected components exist, space them apart horizontally.

    Args:
        compound_graph: Molecule graph containing components.

    Returns:
        dict: Atom index -> (x, y) coordinates for the whole molecule.
    """
    all_positions = {}
    x_cursor = 0.0  # Tracks horizontal offset for next component
    spacing = 4.0   # Fixed horizontal gap between components

    for component in compound_graph.components:
        # Generate layout for this connected fragment, offset horizontally
        comp_positions = generate_2d_layout_for_component(compound_graph, component, x_offset=x_cursor, y_offset=0)

        # Relax positions slightly to improve bond lengths and avoid overlaps
        comp_positions = relax_positions(comp_positions, compound_graph)

        # Merge positions into global dictionary
        all_positions.update(comp_positions)

        # Calculate width of this component to offset the next one
        xs = [pos[0] for pos in comp_positions.values()]
        width = (max(xs) - min(xs)) if xs else 0
        x_cursor += width + spacing  # Move cursor for next component

    return all_positions
