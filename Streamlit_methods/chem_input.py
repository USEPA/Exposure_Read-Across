import streamlit as st
from streamlit_ketcher import st_ketcher
from ctxpy import Chemical
from rdkit import Chem


class FirstSection:

    @st.dialog('Chemical Input', width='large')
    def id_dialog(self, id_choice):

        if id_choice == "Structure":
            st.session_state.smiles_code = False
            st.markdown("Draw the structure in the provided space "
                        "below and then click "
                        "the `Apply` button to retrieve information "
                        "on the chemical.")
            
            # This is where the smiles code is returned
            st.session_state.smiles_code = st_ketcher() 
            st.session_state.id_input = False
            if st.session_state.smiles_code:
                st.rerun()
        
        if id_choice == "DTXSID or CAS-RN":
            st.session_state.id_input = False
            st.markdown("Input a valid DTXSID or CAS-RN and "
                        "then press the `enter` "
                        "keyboard key.")
            st.session_state.id_input = st.text_input('DTXSID or CAS-RN', value=None)
            st.session_state.smiles_code = False
            if st.session_state.id_input:
                st.rerun()
        
    def substance_input(self):
        
        st.markdown('Choose chemical search method from the '
                    'drop-down menu below. To retrieve information '
                    'on an different '
                    'chemical (or if the input box is not popping up '
                    'after the input '
                    'type is selected), click the `Reset` '
                    'button, then choose the '
                    'desired search method.')
        
        searcher = False

        searcher = st.selectbox("Chemical Search Method",
                                ("Structure", "DTXSID or CAS-RN"),
                                index=None)
        
        # Button must be below selection box for correct logical flow
        if st.button("Reset"):
            del st.session_state.smiles_code
            del st.session_state.id_input
            searcher = False

        # Condition executes when script has not been run
        # or the "reset" button has been pressed
        if (("smiles_code" not in st.session_state) and
            ("id_input" not in st.session_state) and
            (searcher)):

            self.id_dialog(searcher)
            
        # Only one of the two session_state variables is assigned per run,
        # requiring the set of conditional returns below
        if "smiles_code" in st.session_state:
            if st.session_state.smiles_code:
                return st.session_state.smiles_code, False
        if "id_input" in st.session_state:
            if st.session_state.id_input:
                return False, st.session_state.id_input
    
        return False, False
            
    def initial_details(self, smiles_code, id_entered):
        chem = Chemical()
        has_dtxsid = False

        with st.container(border=True):
            if smiles_code:
                inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles_code))
                hits = chem.search(by='equals', word=inchikey)
            
            if id_entered:
                hits = chem.search(by='equals', word=id_entered)

            if smiles_code or id_entered:

                if len(hits) > 1 and isinstance(hits, list):
                    # Alert if multiple entries are found
                    st.warning("Multiple chemical matches are returned by CTX."
                               " Information displayed below is from "
                               "the best possible "
                               "match (as determined by CTX).")
                    # Select the dictionary with the highest rank
                    hits = min(hits, key=lambda d: d['rank'])
                    # code below is expecting hits to be in the form of 
                    # a list of dictionaries
                    hits = [hits]
                # All subsequent information is extracted by DTXSID, 
                # so the program cannot proceed if a DTXSID cannot be found
                elif len(hits) < 1:
                    st.warning("No information can be found on a chemical "
                                "with this structure or ID. Check that the "
                                "structure or ID you entered is correct.")  
                    raise Exception

                try:
                    # Lets the user know if the structure being used is 
                    # a generic structure for a class of compounds 
                    if hits[0]['isMarkush']:
                        st.warning("This chemical has a Markush SMILES string."
                                   " This implies that there may be one or "
                                   "more isomeric matches to this structure."
                                   " If you used the Structure search method "
                                   "to identify this compound, "
                                   "please ensure that the returned chemical "
                                   "information is for the correct chemical; "
                                   "otherwise, use the CompTox Chemicals "
                                   "Dashboard to "
                                   "search for your chemical and find "
                                   "the correct "
                                   "DTXSID to use for searching here.")
                    
                    preferred_name = hits[0]['preferredName']
                    hits_smiles = hits[0]['smiles']
                    chem_dtxsid = hits[0]['dtxsid']
                    has_dtxsid = True

                    st.markdown(f'**Name of Chemical**: {preferred_name}')
                    st.markdown(f'**Substance SMILES Identifier**: {hits_smiles}')
                    st.markdown(f'**Substance DSSTox ID**: {chem_dtxsid}')
                except:
                    # Catches structures for which a DTXSID cannot be assigned
                    has_dtxsid = False
                    st.warning("No information can be found on a chemical with"
                               " this "
                               "structure or ID. Check that the structure or "
                               "ID you "
                               "entered is correct.")  

        # Need to return SMILES for later analog calculation
        return chem_dtxsid, has_dtxsid, hits_smiles
