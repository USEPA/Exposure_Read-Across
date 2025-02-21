import streamlit as st
from streamlit_ketcher import st_ketcher
from ctxpy import Chemical
from rdkit import Chem


class first_section:
    
    def substance_input():

        with st.container(border=True):
            st.markdown("## Chemical Input")
            st.markdown("### Input option 1: Structure")
            st.markdown("Draw the structure in the provided space below and then click the "
                        "`Apply` button to retrieve information on the chemical.")

            smiles_code = False
            smiles_code = st_ketcher() #This is where the smiles code is returned 
            
            st.markdown("### Input option 2: Text ID input")
            st.markdown("Input the CAS-RN or DTXSID or exact name of your chemical of interest and then press your 'enter' keyboard key")
            
            #Will need to see if the interpreter reads past the ketcher call above 
            id_input=False
            id_input = st.text_input('CAS-RN, DTXSID, or exact name')


            return smiles_code, id_input
        
    def initial_details(smiles_code, id_entered):
        chem = Chemical(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')
        has_dtxsid = False 
    
        with st.container(border=True):
            if smiles_code:
                inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles_code))
                hits = chem.search(by='equals',word=inchikey) 
            
            if id_entered:
                hits = chem.search(by='equals',word=id_entered) 

            if smiles_code or id_entered:

                if len(hits) > 1 and type(hits)==list:
                    #Alert if multiple entries are found
                    st.warning("Multiple entries returned from the dashboard. Information displayed below is from the return with the smallest rank.")
                    #Select the dictionary with the highest rank
                    hits = min(hits, key=lambda d: d['rank'])
                    #code below is expecting hits to be in the form of a list of dictionaries
                    hits = [hits]
                #All subsequent information is extracted by DTXSID, so the program canot proceed if a DTXSID cannot be found
                elif len(hits)<1 and smiles_code:
                    st.warning("No information can be found on a chemical with this structure. Check that the structure you drew is correct.")
                    raise Exception 
                elif len(hits)<1 and id_entered:
                    st.warning("No information can be found on a chemical matching your entry. Check that the information you entered is correct.")
                    raise Exception 

                try:
                    #Lets the user know if the structure being used is a generic structure for a class of compounds 
                    if hits[0]['isMarkush']:
                        st.warning("Many isomeric structures available for this compound")
                    
                    preferred_name = hits[0]['preferredName']
                    hits_smiles = hits[0]['smiles']
                    chem_dtxsid = hits[0]['dtxsid']
                    has_dtxsid=True
                    

                    st.markdown(f"#### Name of Chemical: {preferred_name}")
                    st.markdown(f"#### Substance SMILES Identifier: {hits_smiles}")
                    st.markdown(f"#### Substance DSSTox ID: {chem_dtxsid}")
                except:
                    #Catches structures for which a DTXSID cannot be assigned
                    has_dtxsid = False
                    st.warning("No information can be found on a chemical with this structure. Check that the structure you drew is correct.")  

        return chem_dtxsid, has_dtxsid, hits_smiles



