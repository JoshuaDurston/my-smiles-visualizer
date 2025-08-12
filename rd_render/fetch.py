"""
Handles retrieval of SMILES strings from user input via PubChemPy:
- Accepts user input as chemical formula or compound name
- Queries PubChem for matching compounds
- Retrieves SMILES strings from compound properties
- Validates SMILES strings before returning
- Uses Streamlit for debug/status messages
"""

from pubchempy import get_compounds, get_properties, BadRequestError
from rd_render.compound import is_valid_smiles
from rd_render.formula import is_valid_formula
import streamlit as st

def get_smiles_from_input(user_input):
    results = []

    if is_valid_formula(user_input):
        try:
            results = get_compounds(user_input, 'formula')
        except BadRequestError:
            results = []

    if not results:
        try:
            results = get_compounds(user_input, 'name')
        except BadRequestError:
            results = []

    if results:
        for compound in results:
            cid = compound.cid
            st.write(f"Checking compound CID {cid}")
            props = get_properties(['CanonicalSMILES', 'IsomericSMILES', 'ConnectivitySMILES'], str(cid), 'cid')

            if props:
                prop = props[0]
                smiles = prop.get('CanonicalSMILES') or prop.get('IsomericSMILES') or prop.get('ConnectivitySMILES')
                if smiles and is_valid_smiles(smiles):
                    st.write(f"Using valid SMILES: {smiles}")
                    return smiles

    return None
