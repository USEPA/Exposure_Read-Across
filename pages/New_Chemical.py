#Although streamlit pages have the option of starting the page titles with a number 
# to set their rank, doing this prevents importing functions into unit tests

import streamlit as st
import pathlib
import Streamlit_methods 

script_location = pathlib.Path(__file__).parent.parent.resolve()

smiles_code, user_input = Streamlit_methods.first_section.substance_input()


if smiles_code or user_input:
        
    target_dtxsid, has_dtxsid, ctxpy_smiles = Streamlit_methods.first_section.initial_details(smiles_code, user_input)            

    if has_dtxsid:
        Streamlit_methods.cpdat_displays(target_dtxsid)
        usis_full = Streamlit_methods.usis_info(target_dtxsid, script_location)
        Streamlit_methods.predicted_info(target_dtxsid, script_location)

    analog_instance = Streamlit_methods.analog_class()
    #Am now generating analogs using the SMILES code returned by ctx-py, 
    # rather than the SMILES returned by Ketcher  
    analog_table = analog_instance.analog_retrieve(script_location, ctxpy_smiles)                            

    if not usis_full.empty:
        #summary_figure = analog_instance.usis_summary_fig(analog_table, usis_full)
        new_fig = analog_instance.osha_cpdat_cdr_seem(analog_table, usis_full, script_location)

#else:
#    st.error('No data present for calculating exposure summary')








