"""
Handles SMILES string processing and visualization for molecular data.
- Parses SMILES strings into RDKit molecule objects safely (suppresses RDKit errors).
- Returns molecule images, properties, and lattice visualizations.
- Retrieves compound names, synonyms, and molecular formula from PubChem.
- Calculates molecular weight and key descriptors.
"""

from rdkit.Chem import Draw, Descriptors, Crippen, rdMolDescriptors
from pubchempy import get_compounds
from rd_render.compound_render import render_molecule
from rd_render.utils import safe_mol_from_smiles

def process_smiles(smiles):
    mol = safe_mol_from_smiles(smiles)
    if not mol:
        return None, None, None, None, None, None, None, None

    # RDKit bitmap image (PIL Image object)
    img = Draw.MolToImage(mol)

    # RDKit descriptors
    mol_mass = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # PubChem data: names and molecular formula
    molecular_formula = None
    names_display = []
    try:
        compounds = get_compounds(smiles, namespace='smiles')
        if compounds:
            compound = compounds[0]
            molecular_formula = compound.molecular_formula

            names = []
            if compound.iupac_name:
                names.append(compound.iupac_name)
            for syn in compound.synonyms:
                if syn and syn not in names:
                    names.append(syn)
                if len(names) >= 5:
                    break

            names_display = [name.title() for name in names]
        else:
            names_display = ["Unknown Molecule"]
    except Exception:
        names_display = ["Lookup failed"]

    # Lattice rendering figure (plotly)
    fig = None
    try:
        fig = render_molecule(smiles)
    except Exception:
        fig = None

    return img, mol_mass, logp, tpsa, hbd, rotatable_bonds, molecular_formula, names_display, fig

def is_valid_smiles(smiles_string):
    mol = safe_mol_from_smiles(smiles_string)
    return mol is not None
