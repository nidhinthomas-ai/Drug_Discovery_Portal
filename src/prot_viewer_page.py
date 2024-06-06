import streamlit as st
from tempfile import NamedTemporaryFile
from streamlit_molstar import st_molstar, st_molstar_rcsb

def prot_viewer_page():
    st.sidebar.title("Protein Viewer")

    pdb_code = st.sidebar.text_input(
        label="Enter PDB Code",
        value="6VN7",
    )

    # Slider for surface transparency
    # surf_transp = st.sidebar.slider("Surface Transparency", min_value=0.0, max_value=1.0, value=0.0)

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
            <h2> PDB <a href="https://www.rcsb.org/structure/{pdb_code}" target="_blank">{pdb_code.upper()}</a></h2>
        </div>
        <div class="content-container">
        """,
        unsafe_allow_html=True
    )

    pdb_file = st.sidebar.file_uploader("Upload PDB file", type=['pdb'])

    displaycol = st.columns([1, 0.1, 250, 1])[2]

    with displaycol:
        if pdb_file:
            with NamedTemporaryFile(delete=False, suffix=".pdb") as temp_file:
                temp_file.write(pdb_file.getvalue())
                temp_file_path = temp_file.name
            st_molstar(temp_file_path, height=600, settings={"backgroundColor": "white", "surfaceTransparency": surf_transp})
        elif pdb_code:
            st_molstar_rcsb(pdb_code, height=600)
        else:
            st.write("Please enter a PDB code or upload a PDB file.")

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
        /* Hide the close button for the sidebar */
        .css-1lcbmhc {
            display: none;
        }
        /* Ensure the sidebar remains open */
        .css-1d391kg {
            pointer-events: none;
        }
        </style>
        </div> <!-- Close the content-container div -->
        """,
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    prot_viewer_page()