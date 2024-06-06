import streamlit as st
import requests
from langchain_community.llms import OpenAI as LangChainOpenAI

def fetch_protein_info_from_rcsb(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        return None

def fetch_protein_info_from_uniprot_by_name(name):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={name}&format=json"
    response = requests.get(url)
    if response.status_code == 200:
        results = response.json()
        if results['results']:
            return results['results'][0]
        else:
            return None
    else:
        return None

def fetch_pdb_ids_from_uniprot(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}/databaseReferences?database=pdb"
    response = requests.get(url)
    if response.status_code == 200:
        results = response.json()
        pdb_ids = [entry['id'] for entry in results.get('databaseReferences', []) if entry['database'] == 'PDB']
        return pdb_ids
    else:
        return None

def format_protein_info(protein_info):
    if not protein_info:
        return "No information found.", None
    
    title = protein_info.get('struct', {}).get('title', 'N/A')
    polymer = protein_info.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_count_protein', 'N/A')
    release_date = protein_info.get('rcsb_accession_info', {}).get('initial_release_date', 'N/A')
    
    details = f"**Title:** {title}\n**Polymer Entity Count (Protein):** {polymer}\n**Release Date:** {release_date}"
    return details, title

def format_protein_info_uniprot(protein_info):
    if not protein_info:
        return "No information found.", None
    
    primary_accession = protein_info.get('primaryAccession', 'N/A')
    protein_name = protein_info.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'N/A')
    organism = protein_info.get('organism', {}).get('scientificName', 'N/A')
    
    details = f"**Primary Accession:** {primary_accession}\n**Protein Name:** {protein_name}\n**Organism:** {organism}"
    return details, primary_accession

def protein_chat_page():
    st.title("Welcome to Protein Chatbot!")

    openai_api_key = st.secrets["OPENAI_API_KEY"]
    llm = LangChainOpenAI(api_key=openai_api_key, temperature=0.1, max_tokens=2048)
    
    if "openai_model" not in st.session_state:
        st.session_state["openai_model"] = llm

    if "temperature" not in st.session_state:
        st.session_state["temperature"] = 0.1  # Default temperature

    if "messages" not in st.session_state:
        st.session_state.messages = [
            {"role": "system", "content": """
            I am here to help you learn about your protein of interest.
            You can provide me with the name of a protein, or its UniProt ID.
            I will search the UniProt and RCSB databases for information about the protein.
            Please provide the name or ID, and I will give you detailed information about the protein, including its structure, function, and relevant references for further learning.
            """}
        ]

    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    if prompt := st.chat_input("Enter the name of the protein or UniProt ID"):
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)

        if prompt.startswith("PDB"):
            pdb_id = prompt.split()[-1]
            protein_info = fetch_protein_info_from_rcsb(pdb_id)
            response_content, title = format_protein_info(protein_info)
        else:
            protein_info = fetch_protein_info_from_uniprot_by_name(prompt)
            if protein_info:
                response_content, primary_accession = format_protein_info_uniprot(protein_info)
                pdb_ids = fetch_pdb_ids_from_uniprot(primary_accession)
                if pdb_ids:
                    pdb_details = "\n".join(pdb_ids)
                    response_content += f"\n**PDB IDs:** {pdb_details}"
                title = primary_accession
            else:
                response_content = "No information found for the given protein name or UniProt ID."
                title = None

        if title:
            detailed_prompt = f"Provide detailed information about the protein with the following details: {response_content}. Include references."
            detailed_response = llm(detailed_prompt)
            response_content += f"\n\n{detailed_response}"

        st.session_state.messages.append({"role": "assistant", "content": response_content})
        with st.chat_message("assistant"):
            st.markdown(response_content)

    temperature = st.sidebar.slider("Temperature", min_value=0.0, max_value=1.0, value=st.session_state["temperature"], step=0.1)
    if temperature != st.session_state["temperature"]:
        st.session_state["temperature"] = temperature
        llm.temperature = temperature  # Update the temperature in LangChain LLM

if __name__ == "__main__":
    protein_chat_page()