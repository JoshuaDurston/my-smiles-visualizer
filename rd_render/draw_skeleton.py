"""
Plotting utilities for rendering molecules using Plotly.

- Draw atoms with color-coded markers and element labels (excluding carbons by default).
- Render bonds with correct single, double, and triple bond styling.
- Render aromatic rings with gray circles.
- Add some functional group text labels (NH2, OH, SH).

Uses RDKit molecule objects and NumPy for coordinate calculations.
"""

from rdkit import Chem
import plotly.graph_objects as go
import numpy as np

# Atom colors for plotting
ATOM_COLORS = {
    'O': 'red',
    'S': 'yellow',
    'N': 'blue',
    'H': 'white',
    'P': 'orange',
}

ATOM_SIZE = 12
BOND_COLOR = 'black'
BOND_WIDTH = 2

def draw_atoms(fig, positions, mol, labels=None):
    labels = labels or {}
    for atom_idx, pos in positions.items():
        # Always draw base atoms, even if special label exists

        atom = mol.GetAtomWithIdx(atom_idx)
        symbol = atom.GetSymbol()

        # Show "OH" label explicitly on oxygen with 1 H as fallback
        if symbol == 'O':
            hydrogens = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'H']
            if len(hydrogens) == 1:
                symbol = 'OH'

        # Skip carbons by default (no label)
        if symbol == 'C':
            continue

        color = ATOM_COLORS.get(symbol[0], 'blue')

        fig.add_trace(go.Scatter(
            x=[pos[0]], y=[pos[1]],
            mode='markers+text',
            text=[symbol],
            textposition='top center',
            marker=dict(size=ATOM_SIZE, color=color),
            textfont=dict(size=16, color='black'),
            hoverinfo='text',
            hovertext=[f"{symbol} (Atom {atom_idx})"],
            showlegend=False
        ))

def offset_lines_centered(x0, y0, x1, y1, n=2, spacing=0.05, inner_ratio=1.0, outer_ratio=0.7):
    # Calculate offsets for multiple bond lines (double, triple)
    dx, dy = x1 - x0, y1 - y0
    norm = np.hypot(dx, dy)
    if norm == 0:
        return []
    # Perpendicular offset vector
    ox, oy = -dy / norm, dx / norm
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    segments = []
    for i in range(-(n//2), n//2 + 1):
        if n % 2 == 0 and i == 0:
            continue  # skip center for even n
        ratio = inner_ratio if i == 0 else outer_ratio
        half_len = norm * ratio / 2
        dxn, dyn = dx / norm * half_len, dy / norm * half_len
        fx, fy = ox * i * spacing, oy * i * spacing
        segments.append(([cx - dxn + fx, cx + dxn + fx], [cy - dyn + fy, cy + dyn + fy]))
    return segments

def draw_bonds(fig, positions, mol):
    # Draw bonds with single, double, triple lines accordingly
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a1 not in positions or a2 not in positions:
            continue
        p1, p2 = positions[a1], positions[a2]
        btype = bond.GetBondType()

        if btype == Chem.rdchem.BondType.SINGLE:
            segments = [([p1[0], p2[0]], [p1[1], p2[1]])]
        elif btype == Chem.rdchem.BondType.DOUBLE:
            segments = offset_lines_centered(p1[0], p1[1], p2[0], p2[1], n=2)
        elif btype == Chem.rdchem.BondType.TRIPLE:
            segments = offset_lines_centered(p1[0], p1[1], p2[0], p2[1], n=3)
        else:
            segments = [([p1[0], p2[0]], [p1[1], p2[1]])]

        for xs, ys in segments:
            fig.add_trace(go.Scatter(
                x=xs, y=ys,
                mode='lines',
                line=dict(color=BOND_COLOR, width=BOND_WIDTH),
                hoverinfo='skip',
                showlegend=False
            ))

def draw_labels(fig, labels):
    # Draw text labels for functional groups (NH2, OH, SH)
    for _, (text, pos) in labels.items():
        fig.add_annotation(
            x=pos[0],
            y=pos[1],
            text=text,
            showarrow=False,
            font=dict(size=12, color="black"),
            align="center",
            bgcolor="white",
            borderpad=2,
            opacity=0.85,
        )

def draw_aromatic_rings(fig, positions, aromatic_rings):
    # Draw gray circles around aromatic rings
    for ring in aromatic_rings:
        coords = np.array([positions[i] for i in ring if i in positions])
        if len(coords) >= 6:
            center = coords.mean(axis=0)
            radius = 0.75
            fig.add_shape(
                type='circle',
                x0=center[0] - radius, x1=center[0] + radius,
                y0=center[1] - radius, y1=center[1] + radius,
                line=dict(color='gray', width=2, dash='solid')
            )
