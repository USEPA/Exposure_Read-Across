#Although streamlit pages have the option of starting the page titles with a number 
# to set their rank, doing this prevents importing functions into unit tests

import streamlit as st
import pathlib
import Streamlit_methods 

script_location = pathlib.Path(__file__).parent.parent.resolve()


smiles_code = Streamlit_methods.first_section.ketcher_smiles()

if smiles_code:
        
        structure_dtxsid, has_dtxsid = Streamlit_methods.first_section.initial_details(smiles_code)            

        if has_dtxsid:
            Streamlit_methods.cpdat_displays(structure_dtxsid)
            
            Streamlit_methods.usis_info(structure_dtxsid, script_location)
            
            Streamlit_methods.predicted_info(structure_dtxsid, script_location)

if smiles_code:   
    analog_instance = Streamlit_methods.analog_operations()
    analog_instance.analog_retrieve(script_location, smiles_code)