import streamlit as st
import pandas as pd
from ctxpy import Chemical
import altair as alt

@st.cache_data
def import_usis(script_location):
    
    """
    Reads in the USIS dataset from Jacob Kvasnika

    Returns:
        Dataframe: The Dataframe of USIS data
    """
    #Looks like the interpreter is able to find the file 
    usis = pd.read_feather(script_location/"data"/"usis_2023_with_dtxsid.feather")

    return usis

#Converts measurements in mg/m^3 into units of ppm,
# assuming pressure of 1 atm and temperature of 25 C=77F
def ppm_to_mg (value, unit, mw)->float:

    """
    Converts multiple starting units to ppm

    Args:
        value (float or int): The numerical value of the air concentration of the chemical of interest.
        unit (str): The unit of the measurement.
        air_volume (float or int): The volume of air sampled
        dtxsid (str): The DTXSID of the target substance.

    Returns:
        float: The air concentration.
    """
    value = float(value) 
    
    #Converts ppm measurements into mg/m^3 
    if unit=='P':
        #Conversion formula from https://www.cdc.gov/niosh/docs/2004-101/calc.html                       
        mg_per_m3_value  = (value*mw)/(24.45)
        return mg_per_m3_value  
    elif unit=='M':
        return value
    else:
        return '0'


def usis_info(structure_dtxsid, script_location):

    chem = Chemical(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')
    usis = import_usis(script_location)


    with st.container(border=True):
                usis_of_interest = usis[usis['dtxsid']==structure_dtxsid]

                if usis_of_interest.empty:
                    st.markdown('### No exposure data available')
                else:
                    #Include a note that all zero-values have been dropped
                    usis_of_interest.drop('index', axis=1, inplace=True)
                    
                    csv_usis_of_interest=usis_of_interest.to_csv(index=False)
                    st.markdown('### USIS data on the entered substance')
                    st.download_button(
                                    label="Download measured exposure data as CSV",
                                    data=csv_usis_of_interest,
                                    file_name="usis_data_on_chemical_of_interest.csv"
                                    )
                    #This dataframe will be displayed as a box and whisker plot

                    #Dropping entries with a value of zero, as these do not contain useful information
                    usis_of_interest=usis_of_interest[usis_of_interest['exposure_level'] != '0']
                    usis_of_interest=usis_of_interest[usis_of_interest['exposure_level'] != '0.0']
                    #Dropping "Fibers per cubic centimeter" and empty rows.
                    #Is not necessary to have a separate clause for null entries, 
                    # as those are classified as not being in the list.
                    usis_of_interest = usis_of_interest[usis_of_interest['measure_unit_id'].isin(['M', 'P'])]

                    deets = chem.details(by='dtxsid', word=structure_dtxsid)   
                    molar_mass = deets['monoisotopicMass']    
                    usis_of_interest['exposure_level'] = usis_of_interest.apply(lambda x: ppm_to_mg (x.exposure_level, x.measure_unit_id, molar_mass), axis=1) 
                    
                    #Change the labeled unit type when the conversion is done
                    usis_of_interest['measure_unit_id'] = 'M'
                    usis_of_interest['measure_unit_name'] = 'Milligrams per cubic meter'
                   
                    #The units that are convertable to mg/m^3:
                        #P - mg/m3, #X- micrograms, Y - milligrams,
                    #Not convertable:
                    # F - fibers/cc, % - percentage of bulk material    

                    #Dropping any rows where the ppm_to_mg function returns a value of zero, 
                    # as these cannot be placed on a log scale 
                    usis_of_interest = usis_of_interest[usis_of_interest['exposure_level']!='0'] 

                    naics_sector_or_sub = st.selectbox("Which level of NAICS results would you like to see?", ("NAICS Sector", "NAICS Subsector"))
                    if not usis_of_interest.empty:
                        st.markdown('#### Occupational inhalation exposure data on the entered substance from the combined datasets CEHD, OIS and IMIS, converted to units of mg/m^3')
                        if naics_sector_or_sub == "NAICS Subsector":
                            st.markdown('##### Data is separated by NAICS subsector. Any entries with units that cannot be converted to mg/m^3 or which have a non-inhalation sample type are not displayed.')
                            box_and_whisker_mg_m3 = alt.Chart(usis_of_interest).mark_boxplot().encode(
                                x= alt.X('exposure_level:Q', 
                                        scale=alt.Scale(type="log", 
                                                        domain=[(usis_of_interest['exposure_level'].min())/10,
                                                                (usis_of_interest['exposure_level'].max())*10 ])).title('Air concentration (ppm)'),
                                y= alt.Y('naics_2022_subsector_title:N', 
                                         sort='-x',
                                         axis=alt.Axis(title='NAICS Subsector',
                                         titleX=-370))).configure_axis(labelLimit=1000)                                                                 
                            
                        if naics_sector_or_sub == "NAICS Sector":
                            usis_of_interest_nona = usis_of_interest[usis_of_interest['naics_2022_sector_title'].notna()]
                            st.markdown('##### Data is separated by NAICS subsector. Any entries with units that cannot be converted to mg/m^3 or which have a non-inhalation sample type are not displayed.')
                            box_and_whisker_mg_m3 = alt.Chart(usis_of_interest_nona).mark_boxplot().encode(
                                x= alt.X('exposure_level:Q', 
                                        scale=alt.Scale(type="log", 
                                                        domain=[(usis_of_interest_nona['exposure_level'].min())/10,
                                                                (usis_of_interest_nona['exposure_level'].max())*10 ])).title('Air concentration (ppm)'),
                                y= alt.Y('naics_2022_subsector_title:N', sort='-x',
                                                                            axis=alt.Axis(title='NAICS Subsector',
                                                                            titleX=-370))).configure_axis(labelLimit=1000)                                                                 
                            


                        st.altair_chart(box_and_whisker_mg_m3, use_container_width = True)
                        # st.markdown('##### Figure note: Chart type is box-and-whisker, but for some substances only outliers can be seen')
                    else:   
                        st.markdown('### No occupational inhalation data available')
    
    return usis

    








