"""
smiles_engine package — tools for parsing SMILES strings and rendering schematic 2D diagrams.

Public API:
- render_from_smiles(smiles: str) -> Plotly Figure
- parse_only(smiles: str) -> CompoundGraph
"""

from .parser import parse_smiles
from .layout import generate_2d_layout
from .renderer import render_compound

def render_from_smiles(smiles_str):
    """
    Given a SMILES string, return a ready-to-display Plotly Figure.
    """
    compound_graph = parse_smiles(smiles_str)
    positions = generate_2d_layout(compound_graph)
    fig = render_compound(compound_graph, positions)
    return fig

def parse_only(smiles_str):
    """
    Given a SMILES string, return a CompoundGraph without layout or rendering.
    Useful for analysis/testing.
    """
    return parse_smiles(smiles_str)
