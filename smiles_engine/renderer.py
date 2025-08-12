"""
renderer.py — Converts positioned atoms & bonds into a visual diagram.

Uses:
- Atom positions from layout.py
- Bond orders, aromatic flags from parser.py
- Atom colours & styles from chem_law.py

Output:
- A Plotly Figure object
"""

import plotly.graph_objects as go
import numpy as np
from .chem_law import get_atom_color, get_bond_style

def offset_lines_centered(x0, y0, x1, y1, n=2, spacing=0.1, inner_ratio=1.0, outer_ratio=0.7):
    """
    Calculate offset line segments centered on the bond line for multi-line bonds.

    Args:
        x0, y0 (float): Coordinates of first atom.
        x1, y1 (float): Coordinates of second atom.
        n (int): Number of parallel lines (e.g. 2 for double bond).
        spacing (float): Distance between lines.
        inner_ratio (float): Relative length ratio for center line.
        outer_ratio (float): Relative length ratio for outer lines.

    Returns:
        List of tuples: Each tuple contains two lists ([x_start, x_end], [y_start, y_end])
        representing one line segment's coordinates.
    """
    dx, dy = x1 - x0, y1 - y0
    norm = np.hypot(dx, dy)
    if norm == 0:
        return []
    # Perpendicular unit vector for offset
    ox, oy = -dy / norm, dx / norm
    # Center point of the bond
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    segments = []

    # Iterate over offsets to generate parallel lines
    for i in range(-(n//2), n//2 + 1):
        if n % 2 == 0 and i == 0:
            continue  # Skip center line for even number of lines (e.g. double bond)
        ratio = inner_ratio if i == 0 else outer_ratio
        half_len = norm * ratio / 2
        dxn, dyn = dx / norm * half_len, dy / norm * half_len
        fx, fy = ox * i * spacing, oy * i * spacing
        segments.append(([cx - dxn + fx, cx + dxn + fx], [cy - dyn + fy, cy + dyn + fy]))

    return segments

def render_compound(compound_graph, positions):
    """
    Create a Plotly figure of the compound.

    Args:
        compound_graph: CompoundGraph from parser.py
        positions: dict of {atom_idx: (x, y)} from layout.py
    """

    bond_traces = []
    for bond in compound_graph.bonds:
        x0, y0 = positions[bond.a1]
        x1, y1 = positions[bond.a2]
        style = get_bond_style(bond)

        # Number of lines to draw based on bond order (single=1, double=2, triple=3)
        n_lines = bond.order if bond.order in (1, 2, 3) else 1

        # Generate parallel line segments for multi-line bonds
        segments = offset_lines_centered(x0, y0, x1, y1, n=n_lines, spacing=0.12)

        # Add each line segment as a separate trace
        for xs, ys in segments:
            bond_traces.append(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="lines",
                    line=dict(color="black", width=1, dash=style["dash"]),
                    hoverinfo="skip",
                    showlegend=False
                )
            )

    atom_traces = []
    for atom in compound_graph.atoms:
        x, y = positions[atom.idx]
        atom_traces.append(
            go.Scatter(
                x=[x],
                y=[y],
                mode="markers+text",
                marker=dict(
                    color=get_atom_color(atom.symbol),
                    size=20
                ),
                text=atom.symbol,
                textposition="top center",
                hovertext=f"{atom.symbol} (idx {atom.idx})",
                hoverinfo="text",
                showlegend=False
            )
        )

    fig = go.Figure(data=bond_traces + atom_traces)
    fig.update_layout(
        xaxis=dict(visible=False, scaleanchor="y", scaleratio=1),
        yaxis=dict(visible=False),
        showlegend=False,
        plot_bgcolor="white",
        dragmode='pan',
        height=800
    )
    return fig
