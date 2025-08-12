"""
Molecule rendering module that snaps molecular structures onto a triangular lattice grid.

- Parses SMILES strings safely and generates 2D coordinates.
- Aligns atoms onto lattice points preserving molecular geometry.
- Builds aromatic rings and highlights some functional groups (NH2, OH, SH).
- Produces a Plotly figure visualizing the molecule on a triangular lattice.

Depends on RDKit for chemistry, NumPy for math, and custom lattice & drawing utilities.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rd_render.triangular_lattice import draw_triangular_lattice, build_lattice_neighbors
from collections import deque
from rd_render.utils import safe_mol_from_smiles
import numpy as np
import streamlit as st
from . import draw_skeleton

def render_molecule(smiles: str):
    # Parse SMILES safely
    mol = safe_mol_from_smiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    # Add explicit hydrogens to detect NH2, OH, SH groups
    mol_with_hs = Chem.AddHs(mol)

    # Compute 2D coordinates for molecule (no explicit H)
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()

    # Extract raw atom 2D positions as numpy array
    raw_positions = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                              for i in range(mol.GetNumAtoms())])

    # Center and scale molecule coordinates to fit lattice nicely
    centroid = raw_positions.mean(axis=0)
    scaled_positions = (raw_positions - centroid) * (1 / 1.4)

    # Lattice offset (adjust origin)
    dx = np.sqrt(3) / 2
    offset = np.array([(30 * dx) / 2, 15])
    adjusted_positions = scaled_positions + offset

    # Draw lattice and get lattice points + neighbors
    fig, lattice_nodes = draw_triangular_lattice()
    lattice_neighbors = build_lattice_neighbors(lattice_nodes)

    # Map atoms to lattice nodes (snapping)
    snapped_positions = {}
    used_nodes = set()

    # Find aromatic rings for special positioning
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    aromatic_rings = [r for r in rings if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]

    def unit_vector(v):
        norm = np.linalg.norm(v)
        return v / norm if norm > 1e-8 else v

    def hexagon_nodes(center, radius=1.0):
        # Return coordinates of hexagon corners around center
        angles_deg = np.arange(0, 360, 60)
        return [center + np.array([np.cos(np.deg2rad(a)) * dx * 2 * radius,
                                  np.sin(np.deg2rad(a)) * radius])
                for a in angles_deg]

    # Place aromatic ring atoms on hex lattice nodes near ring center
    if aromatic_rings:
        ring = aromatic_rings[0]
        ring_pos = np.mean(adjusted_positions[list(ring)], axis=0)
        center_node_idx = int(np.argmin(np.sum((lattice_nodes - ring_pos)**2, axis=1)))
        center_node_pos = lattice_nodes[center_node_idx]

        neighbors = lattice_neighbors.get(center_node_idx, [])
        if len(neighbors) >= 6:
            # Sort neighbors by angle around center to get hex order
            def angle_sort(idx):
                vec = lattice_nodes[idx] - center_node_pos
                return np.arctan2(vec[1], vec[0])
            ring_hex_nodes = sorted(neighbors, key=angle_sort)[:6]
        else:
            # Fallback: find closest unused nodes around center forming hex
            ring_hex_nodes = []
            for coord in hexagon_nodes(center_node_pos):
                dist = np.sum((lattice_nodes - coord)**2, axis=1)
                for idx in np.argsort(dist):
                    if idx not in used_nodes:
                        ring_hex_nodes.append(idx)
                        break

        # Assign ring atoms to snapped lattice nodes
        for atom_idx, node_idx in zip(ring, ring_hex_nodes):
            snapped_positions[atom_idx] = node_idx
            used_nodes.add(node_idx)
    else:
        # No aromatic ring: pick central atom and snap it to closest lattice node
        center_atom_idx = int(np.argmin(np.sum((adjusted_positions - offset)**2, axis=1)))
        start_node_idx = int(np.argmin(np.sum((lattice_nodes - adjusted_positions[center_atom_idx])**2, axis=1)))
        snapped_positions[center_atom_idx] = start_node_idx
        used_nodes.add(start_node_idx)

    # BFS to snap remaining atoms to lattice nodes near connected atoms
    queue = deque(snapped_positions.keys())
    while queue:
        current = queue.popleft()
        current_node_idx = snapped_positions[current]
        current_pos = lattice_nodes[current_node_idx]
        current_rdkit_pos = adjusted_positions[current]

        for bond in mol.GetAtomWithIdx(current).GetBonds():
            neighbor = bond.GetOtherAtomIdx(current)
            if neighbor in snapped_positions:
                continue

            desired_dir = unit_vector(adjusted_positions[neighbor] - current_rdkit_pos)
            candidate_nodes = lattice_neighbors.get(current_node_idx, [])

            # Pick best candidate node with direction closest to desired_dir
            best_node, best_cos = None, -2
            for node_idx in candidate_nodes:
                if node_idx in used_nodes:
                    continue
                candidate_dir = unit_vector(lattice_nodes[node_idx] - current_pos)
                cos_sim = np.dot(desired_dir, candidate_dir)
                if cos_sim > best_cos:
                    best_cos = cos_sim
                    best_node = node_idx

            if best_node is not None:
                snapped_positions[neighbor] = best_node
                used_nodes.add(best_node)
                queue.append(neighbor)
            else:
                st.warning(f"Could not place atom {neighbor} near atom {current}")

    # Final atom positions on lattice
    positions = {atom_idx: lattice_nodes[node_idx] for atom_idx, node_idx in snapped_positions.items()}

    # Create special labels for functional groups (NH2, SH, OH)
    labels = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        idx = atom.GetIdx()
        atom_hs = mol_with_hs.GetAtomWithIdx(idx)

        if sym == "N":
            hydrogens = [nbr for nbr in atom_hs.GetNeighbors() if nbr.GetSymbol() == "H"]
            if len(hydrogens) == 2:
                labels[idx] = ("NH₂", positions.get(idx, np.array([0,0])) + np.array([0.3, 0.3]))
        elif sym == "S":
            hydrogens = [nbr for nbr in atom_hs.GetNeighbors() if nbr.GetSymbol() == "H"]
            if len(hydrogens) >= 1:
                labels[idx] = ("SH", positions.get(idx, np.array([0,0])) + np.array([0.3, 0.3]))
        elif sym == "O":
            hydrogens = [nbr for nbr in atom_hs.GetNeighbors() if nbr.GetSymbol() == "H"]
            if len(hydrogens) == 1:
                labels[idx] = ("OH", positions.get(idx, np.array([0,0])) + np.array([0.3, 0.3]))

    # Draw all parts: bonds, atoms, labels, aromatic rings
    draw_skeleton.draw_bonds(fig, positions, mol)
    draw_skeleton.draw_atoms(fig, positions, mol, labels=labels)
    if labels:
        draw_skeleton.draw_labels(fig, labels)
    draw_skeleton.draw_aromatic_rings(fig, positions, aromatic_rings)

    return fig
