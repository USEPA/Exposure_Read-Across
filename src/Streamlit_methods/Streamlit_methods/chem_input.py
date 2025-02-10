import streamlit as st
from streamlit_ketcher import st_ketcher
from ctxpy import Chemical
from rdkit import Chem


class first_section:
    
    def ketcher_smiles():

        with st.container(border=True):
            st.markdown("## Existing Chemical Input")
            st.markdown("### Structure")
            st.markdown("Draw the structure in the provided space below and then click the "
                        "`Apply` button to retrieve information on the chemical.")
            smiles_code = False
            smiles_code = st_ketcher() #This is where the smiles code is returned 
            
            if smiles_code:
                compute = True
                if len(smiles_code) == 0:
                    compute = False
                    st.warning("Missing SMILES string, did you forget to draw a chemical structure or not click `Apply` after drawing?")

            return smiles_code
        
    def initial_details(smiles_code):
        chem = Chemical(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')
        has_dtxsid = False 

        with st.container(border=True):

            inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles_code))
            hits = chem.search(by='equals',word=inchikey) 

            if len(hits) > 1 and type(hits)==list:
                #Alert if multiple entries are found
                st.warning("Multiple entries returned from the dashboard. Information displayed below is from the return with the smallest rank.")
                #Select the dictionary with the highest rank
                hits = min(hits, key=lambda d: d['rank'])
                #code below is expecting hits to be in the form of a list of dictionaries
                hits = [hits]
            #All subsequent information is extracted by DTXSID, so the program canot proceed if a DTXSID cannot be found
            elif len(hits)<1:
                st.warning("No information can be found on a chemical with this structure. Check that the structure you drew is correct.")
                raise Exception 

            try:
                #Lets the user know if the structure being used is a generic structure for a class of compounds 
                if hits[0]['isMarkush']:
                    st.warning("Many isomeric structures available for this compound")
                structure_dtxsid = hits[0]['dtxsid']
                has_dtxsid=True
                preferred_name = hits[0]['preferredName']

                st.markdown(f"#### Name of Chemical: {preferred_name}")
                st.markdown(f"#### Substance SMILES Identifier: {smiles_code}")
                st.markdown(f"#### Substance DSSTox ID: {structure_dtxsid}")
            except:
                #Catches structures for which a DTXSID cannot be assigned
                has_dtxsid = False
                st.warning("No information can be found on a chemical with this structure. Check that the structure you drew is correct.")  
        
        return structure_dtxsid, has_dtxsid



