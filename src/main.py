import streamlit as st
from home_page import home_page
from prot_viewer_page import prot_viewer_page
from protein_chat_page import protein_chat_page
from ligand_viewer_page import ligand_viewer_page

st.sidebar.title("Drug Discovery Portal")
page = st.sidebar.selectbox("Select Page", ["Home", "Protein Chat", "Protein Viewer", "Ligand Viewer"])

# Handle page navigation
if page == "Home":
    home_page()
elif page == "Protein Viewer":
    prot_viewer_page()
elif page == "Protein Chat":
    protein_chat_page()
elif page == "Ligand Viewer":
    ligand_viewer_page()