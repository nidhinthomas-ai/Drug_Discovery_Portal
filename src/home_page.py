import streamlit as st

def home_page():
    st.title("Welcome to Drug Discovery Portal!")

    # Inject CSS to justify text on both ends
    st.markdown("""
    <style>
    .justified-text {
        text-align: justify;
    }
    </style>
    """, unsafe_allow_html=True)

    st.markdown("""
    ## Introduction

    <div class="justified-text">
    This application is a portal to multiple tools relevant to Drug Dscovery. Tools available in this application can help you learn more about your molecule of interest, visualize those proteins and ligands. More tools will be added to this application in the future. Stay tuned!
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    ## Molecule Chat

    <div class="justified-text">
    The Molecule Chat page is an interactive chatbot powered by GPT-4o. This tool is designed to provide detailed information about molecule based on either the molecule's name, its PDB ID or CheMBL ID. Simply input the molecule name, PDB ID, or ChEMBL ID, and Molecule Chat will deliver relevant details about the Molecule's structure, function, and other key characteristics. Additionally, the chatbot retains the history of your conversation, allowing for a more fluid and organic query experience. Whether you're a student looking to learn more or a researcher needing specific details, Molecule Chat is here to assist you.
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    ## Protein Viewer

    <div class="justified-text">
    The Mol Viewer page is a robust molecular visualization tool based on Stmol. It allows users to visualize and interact with 3D models of proteins. You can input a PDB ID to load the protein structure and then explore various aspects of the molecule. The tool provides options to change the chain ID, adjust colors, and highlight specific residues, giving you a customized view of the protein. This visualization capability is invaluable for researchers and educators needing to illustrate or analyze protein structures.
    </div>
    """, unsafe_allow_html=True)

    # st.markdown("""
    # ## Binding Domain

    # <div class="justified-text">
    # The Binding Domain page is a predictive tool designed to identify the binding domains of proteins from their amino acid sequences. Leveraging a fine-tuned version of the Protein Language Model (ESM2), this application can predict which residues in a protein are likely to be involved in binding interactions. Simply input the protein sequence, and the tool will analyze it to identify potential binding sites. This feature is particularly useful for drug discovery and research, as understanding binding domains is crucial for developing therapeutic interventions.
    # </div>
    # """, unsafe_allow_html=True)

    st.markdown("""
    ## Ligand Viewer

    <div class="justified-text">
    The Ligand Viewer tool allows you to visualize 2D structures of small molecules. You can enter a SMILES string, a ChEMBL ID, or upload a mol or sdf file, and the tool will display the corresponding molecule. This tool leverages the streamlit-ketcher package for rendering 2D structures, providing an intuitive way to visualize and analyze molecular structures for drug discovery and research.
    </div>
    """, unsafe_allow_html=True)