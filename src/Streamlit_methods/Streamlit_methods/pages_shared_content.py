import pathlib
import Streamlit_methods 

def shared_content():
    #Gives the location for the the function location 
    script_location = pathlib.Path(__file__).parent.parent.parent.parent.resolve()
    smiles_code, user_input = False, False
    smiles_code, user_input = Streamlit_methods.first_section.substance_input()

    go_ahead = False
    if smiles_code or user_input:
            
        go_ahead = True

        target_dtxsid, has_dtxsid, ctxpy_smiles = Streamlit_methods.first_section.initial_details(smiles_code, user_input)            

        if has_dtxsid:
            Streamlit_methods.cpdat_displays(target_dtxsid)
            osha_set = Streamlit_methods.osha_info(target_dtxsid, script_location)
            Streamlit_methods.predicted_info(target_dtxsid, script_location)

        return go_ahead, osha_set, script_location, ctxpy_smiles
    else:

        #Avoids the initial error message before a substance is input 
        return False, False, False, False
    
    




