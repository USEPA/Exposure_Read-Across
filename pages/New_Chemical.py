# Although streamlit pages have the option of starting the page titles with a number 
#  to set their rank, doing this prevents importing functions into unit tests

import Streamlit_methods

# This function runs the code that is shared between the two pages 
entered, osha_full, function_location, ctx_smiles = Streamlit_methods.shared_content()

if entered:

    analog_instance = Streamlit_methods.AnalogClass()
    # Am now generating analogs using the SMILES code returned by ctx-py, 
    #  rather than the SMILES returned by Ketcher  
    analog_table = analog_instance.analog_retrieve(function_location, ctx_smiles)                            

    # OSHA, CPDat, and CDR data display
    analog_instance.osha_cpdat_cdr(analog_table, osha_full, function_location)
    # Exposure pathway probabilities
    analog_instance.exp_pthwy(analog_table, function_location)
    # Model predictions display
    analog_instance.model_prediction(analog_table, function_location)





 






