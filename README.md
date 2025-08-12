# MyChemStudy

MyChemStudy is a work-in-progress molecular visualization and property analysis tool built with Python and Streamlit. It parses SMILES strings to generate 2D structural layouts, visualize molecules, and fetch chemical properties using PubChem.

This project is designed primarily as a learning tool to help me practice chemistry concepts while improving my Python coding skills, and to demonstrate integration of RDKit and PubChemPy. Please note it is **not intended for commercial use**.

## Installation
- Python 3.11.9 is recommended for best compatibility.

1. Clone the repo  
2. Install dependencies: `pip install -r requirements.txt rdkit-pypi`  
3. Run: `streamlit run main.py`

## Features

- Parses SMILES strings into molecular graphs with custom chemical rules.
- Generates 2D layouts for molecule visualization.
- Retrieves molecular properties and synonyms from PubChem.
- Interactive Streamlit app for easy exploration.

## Current Limitations

- Handling of disconnected fragments in SMILES is still basic.
- Complex molecules may not always render perfectly.
- Support for ionic bonds and charges is partial.
- UI and chemical logic remain under active development.

## Future Plans

- Improve fragment connection logic.    
- Optimize layout algorithms for clarity.
- Enhanced atom visualisation for specific compounds.
- Some properties generated purely from SMILES.

Feedback is welcome!