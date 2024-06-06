import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit import Chem
from rdkit.Chem import Draw
import requests
from io import BytesIO

def read_sdf_file(file):
    sdf_data = file.read().decode('utf-8')
    sdf_supplier = Chem.SDMolSupplier()
    sdf_supplier.SetData(sdf_data)
    mol = next(iter(sdf_supplier), None)
    if mol:
        smiles = Chem.MolToSmiles(mol)
        return smiles
    return None

def fetch_molecule_from_chembl(chembl_id):
    url = f'https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json'
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        smiles = data.get('molecule_structures', {}).get('canonical_smiles', None)
        return smiles
    return None

def ligand_viewer_page():
    st.sidebar.title("Ligand Viewer")
    DEFAULT_MOL = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    SMILE = st.sidebar.text_input("Enter SMILES string", DEFAULT_MOL)
    chembl_id = st.sidebar.text_input("Enter ChEMBL ID")
    mol_file = st.sidebar.file_uploader("Upload mol or sdf file", type=['mol', 'sdf'])

    smiles = SMILE

    if mol_file:
        if mol_file.name.endswith('.mol'):
            mol = Chem.MolFromMolFile(mol_file)
            if mol:
                smiles = Chem.MolToSmiles(mol)
            else:
                st.sidebar.error("Invalid mol file.")
        elif mol_file.name.endswith('.sdf'):
            sdf_smiles = read_sdf_file(mol_file)
            if sdf_smiles:
                smiles = sdf_smiles
            else:
                st.sidebar.error("Invalid sdf file.")

    elif chembl_id:
        fetched_smiles = fetch_molecule_from_chembl(chembl_id)
        if fetched_smiles:
            smiles = fetched_smiles
        else:
            st.sidebar.error("Invalid ChEMBL ID or unable to fetch molecule.")
    else:
        smiles = SMILE

    st.markdown(
        f"""
        <style>
        .title-container {{
            position: fixed;
            top: 0;
            width: 100%;
            background-color: var(--background-color);
            z-index: 1000;
            padding: 40px;
        }}
        .content-container {{
            margin-top: 10px; 
        }}
        </style>
        <div class="title-container">
            <h2> Molecule SMILES: {smiles}</h2>
        </div>
        <div class="content-container">
        """,
        unsafe_allow_html=True
    )

    # Set the desired width and height for the visualization
    page_width = 1440 
    height = 600

    col1, spacer, col2, col3 = st.columns([1, 0.1, 250, 1])

    with col2:
        smile_code = st_ketcher(smiles, height=height)

    st.markdown(
        """
        <style>
        body, html {
            overflow: hidden;
        }
        .main .block-container {
            max-width: 1440px;
            padding-left: 1rem;
            padding-right: 1rem;
        }
        .css-18e3th9 {
            padding-top: 0rem;
            padding-bottom: 0rem;
            padding-left: 1rem;
            padding-right: 1rem;
        }
        .css-1d391kg {
            padding-top: 0rem;
            padding-bottom: 0rem;
            padding-left: 1rem;
            padding-right: 1rem;
        }
        .css-1lcbmhc {
            display: none;
        }
        .css-1d391kg {
            pointer-events: none;
        }
        </style>
        </div>
        """,
        unsafe_allow_html=True
    )