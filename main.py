# MyChemStudy - SMILES Visualizer
# Copyright (c) 2025 Joshua Durston
# Licensed under the MIT License (see LICENSE file for details)
"""
MyChemStudy Streamlit App: Molecular Visualization and Property Lookup

This app allows users to input chemical compounds via names, formulas, or SMILES strings
and visualizes their structures using two different methods side-by-side:

- Left panel: RDKit-based rendering and property display, including molecular mass,
  lipophilicity (LogP), polarity (TPSA), hydrogen bond donors, rotatable bonds,
  molecular formula, and synonyms fetched from PubChem.

- Right panel: Custom SMILES parser and 2D layout engine implemented from scratch,
  demonstrating chemical structure inference without external chemistry libraries.

Input handling includes validation of SMILES strings and PubChem lookup for non-SMILES inputs.
The UI maintains synchronized input and visualization states using Streamlit session state.

Examples supported in input include:
- Compound names (e.g. "Acetylsalicylic acid")
- SMILES strings (e.g. "c1ccccc1N")
- Molecular formulas (e.g. "CHCl3")
- Color additives (e.g. "FD&C Blue 1")
- Elemental allotropes (e.g. "S8")

Developed as a chemical sketch and analysis tool for study purposes.
"""
import streamlit as st
from rd_render import is_valid_smiles, process_smiles
from rd_render import get_smiles_from_input
from smiles_engine import render_from_smiles

st.set_page_config(layout="wide")
st.title("MyChemStudy: SMILES Compound Visualizer")

# Initialize global smiles variable in session state (default: ethanol)
if "smiles" not in st.session_state:
    st.session_state.smiles = "CCO"  # Default starting SMILES
if "input_str" not in st.session_state:
    st.session_state.input_str = st.session_state.smiles

def update_smiles():
    user_input = st.session_state.input_str.strip()
    if user_input == st.session_state.smiles:
        return  # no change

    if is_valid_smiles(user_input):
        st.session_state.smiles = user_input
    else:
        try:
            with st.spinner("Searching..."):
                fetched_smiles = get_smiles_from_input(user_input)
                if fetched_smiles:
                    st.session_state.smiles = fetched_smiles
                else:
                    st.error("Could not interpret input as a valid SMILES, name, or formula.")
        except Exception as e:
            st.error(f"An error occurred while querying PubChem: {e}")

left_col, right_col = st.columns([1, 1])

with left_col:
    st.markdown(
        """
        ### RDKit-based Molecular Visualization & Properties
        This section uses RDKit to generate a 2D structure of the molecule.  
        The molecule is placed on a triangular lattice by analyzing neighbors and snapping atoms to grid points.  
        It also fetches molecular properties and synonyms from PubChem to give a fuller chemical overview.
        """
    )

    user_input = st.text_input(
        "Enter a valid compound name, formula, or SMILES string:",
        key="input_str",
        on_change=update_smiles,
    )
    st.markdown(
        "Complex formulae may take a long time to load or simply timeout. Examples: "
        "'Acetylsalicylic acid', 'c1ccccc1N', 'CHCl3', 'FD&C Blue 1', 'S8'"
    )

    img = mol_mass = logp = tpsa = hbd = rot_bonds = mol_formula = names_display = lattice_fig = None

    results = process_smiles(st.session_state.smiles)
    if results:
        (
            img,
            mol_mass,
            logp,
            tpsa,
            hbd,
            rot_bonds,
            mol_formula,
            names_display,
            lattice_fig,
        ) = results
    else:
        st.error("Failed to process the SMILES string.")

    if img:
        st.image(img, caption="RDKit Compound Structure", width=300)

    if mol_mass is not None:
        st.metric("Molar Mass", f"{mol_mass:.2f} g/mol")

    if mol_formula:
        st.write(f"**Formula:** {mol_formula}")

    prop_cols = st.columns(2)
    with prop_cols[0]:
        if logp is not None:
            st.write(f"**LogP (Lipophilicity):** {logp:.2f}")
        if hbd is not None:
            st.write(f"**H-Bond Donors:** {hbd}")

    with prop_cols[1]:
        if tpsa is not None:
            st.write(f"**TPSA (Polarity):** {tpsa:.2f}")
        if rot_bonds is not None:
            st.write(f"**Rotatable Bonds:** {rot_bonds}")

    if names_display:
        st.subheader(names_display[0])
        if len(names_display) > 1:
            with st.expander("Also known as..."):
                for alt_name in names_display[1:]:
                    st.write(f"- {alt_name}")

    if lattice_fig:
        st.plotly_chart(lattice_fig, use_container_width=True, height=300)

with right_col:
    st.header("SMILES Engine Graph")

    st.markdown(
        """
        ### Custom SMILES Parsing & Layout Engine  
        This section does not rely on RDKit or any external chemistry libraries. Instead, it parses SMILES strings using a custom-built Python parser that interprets atoms, bonds, branches, and ring closures through explicit rules and heuristics.

        It reconstructs the molecule’s atomic connectivity and bond orders into an internal graph representation, then generates 2D coordinates for visualization by placing rings as regular polygons and arranging other atoms via force-directed layout algorithms.

        This experimental approach is not perfect or precise by any means, but demonstrates how chemical structures can be inferred and visualized from SMILES alone, without using third-party cheminformatics tools.
        """
    )

    fig = render_from_smiles(st.session_state.smiles)
    st.plotly_chart(fig, use_container_width=True)
