import streamlit as st
from home_page import home_page
from prot_viewer_page import prot_viewer_page
from binding_domain_page import binding_domain_page
from protein_chat_page import protein_chat_page
from ligand_viewer_page import ligand_viewer_page

# Sidebar for page navigation
st.sidebar.title("Drug Discovery Portal")
# page = st.sidebar.selectbox("Select Page", ["Home", "Protein Chat", "Protein Viewer", "Ligand Viewer", "Binding Domain"])
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
# elif page == "Binding Domain":
#     binding_domain_page()

# # Additional handling for URL hash navigation
# if 'nav_page' not in st.session_state:
#     st.session_state['nav_page'] = 'Home Page'

# query_params = st.query_params
# if 'page' in query_params:
#     st.session_state['nav_page'] = query_params['page'][0]
