import pathlib
import Streamlit_methods


def shared_content():
    # Gives the location for the the function location
    script_location = pathlib.Path(__file__).parent.parent.resolve()
    smiles_code, id_entered = False, False

    first_section_instance = Streamlit_methods.FirstSection()

    # Takes the user input of substance and returns 
    # either the SMILES code or the ID entered by the user
    smiles_code, id_entered = first_section_instance.substance_input()

    go_ahead = False
    if smiles_code or id_entered:

        go_ahead = True
        # Displays the initial results returned 
        # from lookup of entered information
        target_dtxsid, has_dtxsid, ctxpy_smiles = (first_section_instance
                                                   .initial_details(smiles_code,
                                                                    id_entered))

        if has_dtxsid:
            reported_data = Streamlit_methods.ReportedInfo(target_dtxsid,
                                                           script_location)
            # Displays functional-use data
            reported_data.functional_use()
            # Displays list presence data
            reported_data.list_presence()
            # Displays product-inclusion data
            reported_data.product_inclusion()

            # osha_set is the OSHA info before any cleaning  
            osha_set = reported_data.osha_info()
            modeled_data = Streamlit_methods.ModeledData(target_dtxsid,
                                                         script_location,
                                                         osha_set)
            # Displays Quantitative Structure-Use Relationship data
            modeled_data.qsur_display()
            # Output from Jeff's model
            modeled_data.minucci_model()
            # Displays seem3 data
            # along with selected minnuci model data
            modeled_data.seem3_and_minucci()

        return go_ahead, osha_set, script_location, ctxpy_smiles
    else:

        # Avoids the initial error message before a substance is input
        return False, False, False, False
     
    




