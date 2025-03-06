import streamlit as st
from streamlit_ketcher import st_ketcher
from ctxpy import Chemical
from rdkit import Chem


class FirstSection:
    
    
    @st.dialog('Choose your input method', width='large')
    def id_dialog(self, id_choice):

        if id_choice == "Structure":
            st.session_state.smiles_code = False
            st.markdown("Draw the structure in the provided space below and then click the "
                        "`Apply` button to retrieve information on the chemical.")
            st.session_state.smiles_code = st_ketcher() #This is where the smiles code is returned 
            st.session_state.id_input = False
            if st.session_state.smiles_code:
                 st.rerun()
        
        if id_choice == "DTXSID or CAS-RN":
            st.session_state.id_input = False
            st.markdown("Input the DTXSID or CAS-RN of your chemical of interest and then press the 'enter' keyboard key")
            st.session_state.id_input = st.text_input('DTXSID or CAS-RN', value=None)
            st.session_state.smiles_code = False
            if  st.session_state.id_input:
                st.rerun()
        


    def substance_input(self):
        st.markdown('### Choose your desired method for chemical identity input '
                    'from the drop-down menu. To retrieve information on an additional chemical, '
                    'click the `Reset` button, then choose the desired input type')
        
        searcher=False

        searcher = st.selectbox("Chemical Search Method",("Structure","DTXSID or CAS-RN"), index=None)
        
        #Button must be below selection box for correct logical flow 
        if st.button("Reset"):
            del st.session_state.smiles_code
            del st.session_state.id_input
            searcher = False

        #Condition executes when script has not been run or the "reset" button has been pressed
        if ("smiles_code" not in st.session_state) and ("id_input" not in st.session_state) and searcher:
            self.id_dialog(searcher)
            
        #Only one of the two session_state variables is assigned per run,
        # requiring the set of conditional returns below 
        if "smiles_code" in st.session_state:
            if st.session_state.smiles_code:
                return st.session_state.smiles_code, False
        if "id_input" in st.session_state:   
            if st.session_state.id_input:
                return False, st.session_state.id_input    
    
        return False, False
        
        #return st.session_state.smiles_code, st.session_state.id_input
    

    def initial_details(self, smiles_code, id_entered):
        chem = Chemical()
        has_dtxsid = False 
        #if smiles_code and id_entered:
            #id_entered

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
                    st.warning("No information can be found on a chemical with this structure or ID. Check that the structure or ID you entered is correct.")  

        #Need to return SMILES for later analog calculation
        return chem_dtxsid, has_dtxsid, hits_smiles



