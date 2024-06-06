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
    This application is a portal to multiple tools relevant to Drug Discovery. Tools available in this application can help you learn more about your protein of interest, visualize those proteins and ligands. More tools will be added to this application in the future. Stay tuned!
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    ## Protein Chat

    <div class="justified-text">
    The Protein Chat page is an interactive chatbot powered by GPT-4o. This tool is designed to provide detailed information about the protein based on either the protein's name or its PDB ID. Simply, input the protein name or PDB ID, and Protein Chat will deliver relevant details about the Protein's structure, function, and other key characteristics. Additionally, the chatbot retains the history of your conversation, allowing for a more fluid and organic query experience. Whether you're a student looking to learn more or a researcher needing specific details, Protein Chat is here to assist you.
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    ## Protein Viewer

    <div class="justified-text">
    The Mol Viewer page is a robust molecular visualization tool based on Stmol. It allows users to visualize and interact with 3D models of proteins. You can input a PDB ID to load the protein structure and then explore various aspects of the molecule. The tool provides options to change the chain ID, adjust colors, and highlight specific residues, giving you a customized view of the protein. This visualization capability is invaluable for researchers and educators needing to illustrate or analyze protein structures.
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    ## Ligand Viewer

    <div class="justified-text">
    The Ligand Viewer tool allows you to visualize 2D structures of small molecules. You can enter a SMILES string, a ChEMBL ID, or upload a mol or sdf file, and the tool will display the corresponding molecule. This tool leverages the streamlit-ketcher package for rendering 2D structures, providing an intuitive way to visualize and analyze molecular structures for drug discovery and research.
    </div>
    """, unsafe_allow_html=True)