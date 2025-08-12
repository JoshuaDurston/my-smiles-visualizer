"""
Provides utilities for creating and working with a 2D triangular (hexagonal-style) lattice.

Functions:
    draw_triangular_lattice(cols=30, rows=30, spacing=1.0)
        Generates a triangular lattice of points and returns both a Plotly figure
        for visualization and the lattice node coordinates.

    build_lattice_neighbors(nodes, spacing=1.0, tolerance=0.05)
        Given lattice node coordinates, builds a neighbor adjacency list based
        on Euclidean distance with a configurable tolerance.
"""

import plotly.graph_objects as go
import numpy as np

def draw_triangular_lattice(cols=30, rows=30, spacing=1.0):
    # Horizontal and vertical spacing:
    dy = spacing                          # Vertical spacing between nodes.
    dx = spacing * np.sqrt(3) / 2         # Horizontal offset between columns (60° triangle basis).

    x_coords = []  # X positions of all nodes
    y_coords = []  # Y positions of all nodes

    # Generate lattice node positions
    for col in range(cols):
        for row in range(rows):
            x = col * dx
            # Stagger every other column to form the triangle pattern
            y = row * dy + (dy / 2 if col % 2 else 0)
            x_coords.append(x)
            y_coords.append(y)

    # Combine x and y into a (N, 2) numpy array of node positions
    nodes = np.column_stack((x_coords, y_coords))

    # Compute center of the lattice (for zooming/centering view)
    x_center = (cols * dx) / 2
    y_center = (rows * dy) / 2

    zoom_size = 6            # Controls how much of the grid is shown
    circle_radius = 0.05     # Radius of the grey background dots (in plot units)

    # Create grey circle shapes for each node for visual effect
    node_shapes = [
        dict(
            type="circle",
            xref="x", yref="y",  # Position relative to x/y axes
            x0=x - circle_radius,
            x1=x + circle_radius,
            y0=y - circle_radius,
            y1=y + circle_radius,
            line=dict(color='rgba(128, 128, 128, 0.1)'),     # Circle outline
            fillcolor='rgba(128, 128, 128, 0.2)',            # Circle fill
            layer="below"                # Place under other plot layers
        )
        for x, y in zip(x_coords, y_coords)
    ]

    # Initialize the Plotly figure
    fig = go.Figure()

    # Add an invisible scatter plot so you can hover nodes
    fig.add_trace(go.Scatter(
        x=x_coords,
        y=y_coords,
        mode='markers',
        marker=dict(size=20, color='rgba(0,0,0,0)'),  # Invisible marker
        hoverinfo='text',
        text=[f'Node {i} ({round(x,3)}, {round(y,3)})' 
              for i, (x, y) in enumerate(zip(x_coords, y_coords))],
        showlegend=False
    ))

    # Configure layout (zoom level, axis visibility, etc.)
    fig.update_layout(
        shapes=node_shapes,
        xaxis=dict(
            range=[x_center - zoom_size / 2, x_center + zoom_size / 2],
            visible=False
        ),
        yaxis=dict(
            range=[y_center + zoom_size / 2, y_center - zoom_size / 2],
            visible=False,
            scaleanchor='x'  # Equal aspect ratio for x and y
        ),
        width=700,
        height=700,
        margin=dict(l=20, r=20, t=20, b=20),
        plot_bgcolor='white',
        dragmode='pan',  # Allow user to drag the plot
    )

    return fig, nodes


def build_lattice_neighbors(nodes, spacing=1.0, tolerance=0.05):
    neighbors = {}
    threshold = spacing * (1 + tolerance)  # Maximum allowed neighbor distance

    for i, node in enumerate(nodes):
        # Compute distance from node i to all other nodes
        dists = np.linalg.norm(nodes - node, axis=1)

        # Filter: neighbors are nodes within threshold distance, excluding itself
        neighbor_idxs = np.where((dists > 0) & (dists <= threshold))[0]

        neighbors[i] = neighbor_idxs.tolist()

    return neighbors
