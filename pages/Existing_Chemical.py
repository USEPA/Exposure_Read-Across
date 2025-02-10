#Although streamlit pages have the option of starting the page titles with a number 
# to set their rank, doing this prevents importing functions into unit tests

import streamlit as st
import pathlib
import Streamlit_methods 

script_location = pathlib.Path(__file__).parent.parent.resolve()

with st.container(border=True):
        st.markdown('### The information presented here is useful for new chemical evaluation in two scenarios:' )
        st.markdown('###     When the new substance of interest is in the DSSTox database, but '
                    'not in the TSCA active inventory for the proposed use it is being evaluated for. ')
        st.markdown('###     Or, when the substance of interest is not in the DSSTox database '
                    'and information is desired on the structural analogues of the substance.')
        
#The streamlit_ketcher method does not run if the entry box is empty, 
#so writing a way to catch when it is empty is not necessary

smiles_code = Streamlit_methods.first_section.ketcher_smiles()

if smiles_code:
        
        structure_dtxsid, has_dtxsid = Streamlit_methods.first_section.initial_details(smiles_code)            

        if has_dtxsid:
            Streamlit_methods.cpdat_displays(structure_dtxsid)
            
            Streamlit_methods.usis_info(structure_dtxsid, script_location)
            
            Streamlit_methods.predicted_info(structure_dtxsid, script_location)
                                        

           





                    
