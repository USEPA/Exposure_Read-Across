import pandas as pd
import numpy as np
#from pandas.io.formats.style import Styler
import streamlit as st
from streamlit_ketcher import st_ketcher
import pathlib
#from typing import Optional, Tuple
from chembl_webresource_client.new_client import new_client as ch
from rdkit import Chem
from ctxpy import Chemical, Exposure
import altair as alt
from sklearn.neighbors import KNeighborsRegressor
from rdkit.Chem import rdFingerprintGenerator

script_location = pathlib.Path(__file__).parent.parent.resolve()
   
@st.cache_data
def import_usis():
    
    """
    Reads in the USIS dataset from Jacob Kvasnika

    Returns:
        Dataframe: The Dataframe of USIS data
    """
    #Must use an absolute path so that the unit-testing script is able to read it 
    usis = pd.read_feather(script_location/"data"/"usis_2023_with_dtxsid.feather")

    return usis


usis = import_usis()


@st.cache_data
def minucci_loader():#Returns -> DataFrame, but declaring this in the function annotation throws an error
    
    """
    Loads in the Minucci model dataset.

    Returns:
        Dataframe: The Dataframe of the Minucci dataset.
    """
    #This version of consolidated_minucci has been compressed with the "snappy" algorithm, which is fast to retrieve.
    #Consolidated_minucci has been split into two parts so that they will be small enough to upload to github

    minucci_path1 = (script_location/"data"/"consolidated_minucci_snappy_1.parq")
    minucci_path2 = (script_location/"data"/"consolidated_minucci_snappy_2.parq")
    minucci_1 = pd.read_parquet(minucci_path1)
    minucci_2 = pd.read_parquet(minucci_path2)

    df = pd.concat([minucci_1, minucci_2])


    return df


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


expo = Exposure(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')
chem = Chemical(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')


minucci=minucci_loader()



with st.container(border=True):
        st.markdown('### The information presented here is useful for new chemical evaluation in two scenarios:' )
        st.markdown('###     When the new substance of interest is in the DSSTox database, but '
                    'not in the TSCA active inventory for the proposed use it is being evaluated for. ')
        st.markdown('###     Or, when the substance of interest is not in the DSSTox database '
                    'and information is desired on the structural analogues of the substance.')
        
with st.form("ExistingChem"):
    st.markdown("## Existing Chemical Input")
    st.markdown("### Structure")
    st.markdown("Draw the structure in the provided space below and then click "
                "`Apply` to retrieve information on the chemical.")
    smiles_ready = False
    smiles_code = st_ketcher() #This is where the smiles code is returned 
    
    exst_chem_submitted = st.form_submit_button("Submit")
    
    if exst_chem_submitted:
        compute = True
        if len(smiles_code) == 0:
            compute = False
            st.warning("Missing SMILES string, did you forget to draw a chemical structure or not click `Apply` after drawing?")
        else:
            smiles_ready=True 

if exst_chem_submitted:
    # "not compute" indicates that nothing was entered into the structure box
    if not compute:
        st.error("Missing information, see input form for more information.")
    else:
        with st.container(border=True):

            inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles_code))
            hits = chem.search(by='equals',word=inchikey) 

            if len(hits) > 1 and type(hits)==list:
                #Alert if multiple entries are found
                st.warning("Multiple hits returned from the dashboard. Information displayed below is from the return with the smallest rank.")
                #Select the dictionary with the highest rank
                hits = min(hits, key=lambda d: d['rank'])
                #code below is expecting hits to be in the form of a list of dictionaries
                hits = [hits]
            #All subsequent information is extracted by DTXSID, so the program canot proceed if a DTXSID cannot be found
            elif len(hits)<1:
                st.warning("No information can be found on a chemical with this structure. Check that the structure you drew is correct.")
                raise Exception 
            
            

            try:
                #Lets the user know if the structure being used is a generic structure for a class of compounds 
                if hits[0]['isMarkush']:
                    st.warning("Many isomeric structures available for this compound")
                structure_dtxsid = hits[0]['dtxsid']
                has_dtxsid=True
                preferred_name = hits[0]['preferredName']

                st.markdown(f"#### Name of Chemical: {preferred_name}")
                st.markdown(f"#### Substance SMILES Identifier: {smiles_code}")
                st.markdown(f"#### Substance DSSTox ID: {structure_dtxsid}")
            except:
                #Catches structures for which a DTXSID cannot be assigned
                st.warning("No information can be found on a chemical with this structure. Check that the structure you drew is correct.")               

        if has_dtxsid:
        
            with st.container(border=True):
                #Pulls reported functional use, predicted functional use, and presence in consumer/industrial formulations or articles 
                #Dataset displays are sequential due to being too wide to display as side-by-side columns 
                
                st.markdown('## Usage Information ')
                st.markdown('##### Mouse-over any entries whose labels are cut off to see the full labels')
                st.markdown('##### Data pulled from CPDAT by DTXSID ')
                #Cannot find a way to control the size of header text 
                st.header('', divider=True)

                #Pulling reported functional use
                function_info = expo.search_cpdat(vocab_name="fc", dtxsid=structure_dtxsid)
                function_info = pd.DataFrame(function_info)
                function_info_csv = function_info.to_csv(index=False)
                st.download_button(
                                    label="Download all substance functional use data as CSV",
                                    data=function_info_csv,
                                    file_name="functional_use_data.csv"
                                    )
                    
                
                if function_info.empty:
                    st.markdown('### No Reported Function information available')
                else:
                    #Dropping id, dtxsid, docid columns 
                    function_info.drop(columns=['id','dtxsid','docid',], inplace=True)
                    
                    
                    #Dropping rows with no function information 
                    function_info=function_info[function_info['reportedfunction'].notna()]
                    

                    function_info.rename(columns={'datatype':'Data Type', 
                                                'doctitle':'Document Title', 
                                                'docdate':'Document Date', 
                                                'reportedfunction':'Reported Function', 
                                                'functioncategory':'Function Category'}, inplace=True)
                    
                    chart_function_info = function_info.drop_duplicates(subset=['Document Title', 
                                                                                'Reported Function',
                                                                                'Function Category'
                                                                                ],
                                                                        ignore_index=True  
                                                                        )# inplace=True  
                    grouped_cfi=chart_function_info.groupby(by=['Reported Function',
                                                    'Function Category']
                                                    ).size(
                                                    ).reset_index(name='Docs per function')
                    #Experimented with making this plot vertical, but decided that it looked worse,
                    #due to the bars being much shorter and running into the legend. 
                    # I don't know if there is a way to move different components of the plot around  
                    docs_per_function_cat = alt.Chart(grouped_cfi).mark_bar().encode(
                        x="Docs per function",
                        y='Function Category',
                        color=alt.Color('Reported Function', scale=alt.Scale(scheme='category20') )
                    ).configure_axis(labelLimit=1000)
                    
                    st.altair_chart(docs_per_function_cat, use_container_width = True)
                    st.markdown("##### Figure note: Collected documents containing entered substance grouped by reported function within a function category. " 
                                "Only data with a reported function has been included. "
                                "This data facilitates the understanding of the uses of the entered substance by displaying the substance's reported functions."
                                )  
        
                #Now pulling predicted functional use
                
                #Now pulling list presence information 
                st.markdown('### List Presence Information')
                st.markdown('##### The lists on which the entered substance is present. ' 
                            'These provide additional information about the contexts in which the substance is used.')
                list_presence = expo.search_cpdat(vocab_name='lpk', dtxsid=structure_dtxsid)
                list_presence = pd.DataFrame(list_presence)
                if list_presence.empty:
                    st.markdown('### No List Presence information available')
                else:
                    list_presence.drop(columns=['id','dtxsid','docid',], inplace=True)
                    list_presence.rename(columns={'doctitle':'Document Title', 
                                                'docdate':'Document Date',
                                                'reportedfunction':'Reported Function',
                                                'functioncategory':'Function Category',
                                                'keywordset':'Key Word Set',
                                                'docsubtitle':'Document Subtitle',
                                                'organization':'Organization',
                                                }, inplace=True)
                    st.dataframe(list_presence) 
                
                
                st.markdown('### Products Reporting This Substance')
                pucs = expo.search_cpdat(vocab_name='puc', dtxsid=structure_dtxsid)
                pucs_available = pd.DataFrame(pucs)
                pucs_available_csv=pucs_available.to_csv(index=False)
                st.download_button(
                                    label="Download product use data as CSV",
                                    data=function_info_csv,
                                    file_name="functional_use_data.csv"
                                    )
                
                if list_presence.empty:
                    st.markdown('### No PUC Information Available')
                else:
                    try:
                        pucs_available.drop(columns=['id','dtxsid','docid',], inplace=True)
                        pucs_available.rename(columns={'doctitle':'Document Title', 
                                                    'docdate':'Document Date', 
                                                    'productname':'Product Name', 
                                                    'gencat':'General Category',
                                                    'prodfam':'Product Family',
                                                    'prodtype':'Product Type',
                                                    'classificationmethod':'Classification Method',
                                                    'rawmincomp':'Raw Minimum \n Composition',
                                                    'rawmaxcomp':'Raw Maximum Composition',
                                                    'rawcentralcomp':'Raw Central Composition',
                                                    'unittype':'Unit Type',
                                                    'lowerweightfraction':'Lower Weight Fraction',
                                                    'upperweightfraction':'Upper Weight Fraction',
                                                    'centralweightfraction':'Central Weight Fraction',
                                                    'weightfractiontype':'Weight Fraction Type',
                                                    'component':'Component',
                                                    }, inplace=True)
                        no_dups_pucs_available=pucs_available.drop_duplicates(subset=['Product Name', 
                                                                                        'General Category',
                                                                                        'Product Family',
                                                                                        ], ignore_index=True)
                        
                        grouped_pucs = no_dups_pucs_available.groupby(by=['General Category', 'Product Family']
                                                                        ).size(
                                                                        ).reset_index(name='Number of product names')
                        grouped_pucs['PUC general category + product family combination'] = grouped_pucs['General Category'] + " + " + grouped_pucs['Product Family']
                        
                        #Dropping empty cells 
                        grouped_pucs=grouped_pucs[grouped_pucs['Product Family']!='']

                        prods_per_catfam = alt.Chart(grouped_pucs).mark_bar().encode(
                        x=alt.X('Number of product names:Q'),
                        y=alt.Y('PUC general category + product family combination:N',  sort='-x',  
                                                                                        axis=alt.Axis(title='PUC general category + product family combination',
                                                                                        titleX=-370),
                        )).configure_axis(labelLimit=1000)# When this method is used, the labels intrude up into the chart and replace it. 
                        #Probably the height of the chart will need to be defined as well

                        
                        st.altair_chart(prods_per_catfam, use_container_width = True)
                        st.markdown('##### Figure note: Number of products associated with each combination of PUC general category and product family')
                    except:
                        st.write("No product use information is available")

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

                    if not usis_of_interest.empty:
                        st.markdown('#### Occupational inhalation exposure data on the entered substance, converted to units of mg/m^3')
                        st.markdown('##### Data is separated by NAICS subsector. Any entries with units that cannot be converted to mg/m^3 or which have a non-inhalation sample type are not displayed.')
                        box_and_whisker_mg_m3 = alt.Chart(usis_of_interest).mark_boxplot().encode(
                            x= alt.X('exposure_level:Q', 
                                    scale=alt.Scale(type="log", 
                                                    domain=[(usis_of_interest['exposure_level'].min())/10,
                                                            (usis_of_interest['exposure_level'].max())*10 ])).title('Air concentration (ppm)'),
                            y= alt.Y('naics_2022_subsector_title:N', sort='-x',
                                                                        axis=alt.Axis(title='NAICS Subsector',
                                                                        titleX=-370))).configure_axis(labelLimit=1000)                                                                 
                        
                        st.altair_chart(box_and_whisker_mg_m3, use_container_width = True)
                        st.markdown('##### Figure note: Chart type is box-and-whisker, but for some substances only outliers can be seen')
                    else:   
                        st.markdown('### No occupational inhalation data available')

            with st.container(border=True):
                st.header(' Predicted Information', divider=True )
                st.markdown('#### To allow more-informed decision-making, '
                            'model predictions can fill information gaps for substances that are in DSSTox ' 
                            'but do not have sufficient exposure information available in the fields above.')
                st.markdown("### Predicted Functional Use")
                st.markdown('#### Data pulled from EPA QSUR models. Paper available [here](https://pubs.rsc.org/en/content/articlelanding/2017/gc/c6gc02744j#!divCitation):')
                predicted_info = expo.search_qsurs(dtxsid=structure_dtxsid)
                predicted_info = pd.DataFrame(predicted_info)
                if predicted_info.empty:
                    st.markdown('### No Reported Function information available')
                else:
                    predicted_info.rename(columns={'harmonizedFunctionalUse':'Harmonized Functional Use'},
                                        inplace=True)
                    predicted_info.sort_values(by=['probability'], ascending=False, inplace=True)
                    st.dataframe(predicted_info)
                
                    #Output from Jeff's model 
                st.markdown('### Predicted Occupational Exposure')
                st.markdown('#### Paper describing model available [here](https://pubs.acs.org/doi/10.1021/acs.est.2c08234):')
                minucci_portion = minucci[minucci['dtxsid']==structure_dtxsid]                  
                # and then display them in the SEEM3 model chart

                if minucci_portion.empty:
                    st.markdown('### No modeled data available')
                else:
                    minucci_portion.sort_values(by='log_mgm3_pred_50th', ascending=False, inplace=True)
                    st.dataframe(minucci_portion)   

                    
                    # Grabs the top 5 highest predicted concentrations
                    # and then averages them
                    #If the list of predictions is shorter than 5, takes all of them
                    if minucci_portion.shape[0]>=5:
                        top_av = (minucci_portion['log_mgm3_pred_50th'][0:5].sum())/5
                    else:
                        top_av = (minucci_portion['log_mgm3_pred_50th'].sum())/minucci_portion.shape[0]
                    
                    #Now converting the concentration predictions to exposure estimates mg/kg/day
                    #Using the formula "ADD : average daily dose (mg/kg-day)" 
                    # from the user_defined_inhalation model. 
                    # It is defined by the following equations 
                    # (the following information is from the package CLOET):

                    #  ADD = (I * ED * EY) / (BW * AT * days_per_year)
                    #Where: 
                    ED = 250 #integer, days exposed per year; 250 is all the weekdays in a year minus two weeks of vacation
                    Cm = top_av #None or float,  mass concentration of chemical in air (mg/m^3)
                    b = 0.636 # volumetric inhalation rate; 0 <= b <=7.9 
                    #(the default value suggested by CLOET is 1.25 m^3/hr, but the ATSDR document "Guidance for Inhalation Exposures" suggests 0.0106 m^3/min (0.636 m^3/hr))
                    #Document link: https://www.atsdr.cdc.gov/pha-guidance/resources/ATSDR-EDG-Inhalation-508.pdf
                    h = 8.5 #: daily exposure duration (hrs/day)
                    EY = 40 #integer, years of occupational exposure; 0 <= EY (default: 40 years)
                    BW = 70 #float, body weight; (0 <= BW) (default: 70 kg)
                    AT = 40 # float, averaging time (EY <= AT <= ATc); (default: 40 years)
                    days_per_year = 365
                    I = Cm * b * h
                    ADD = (I * ED * EY) / (BW * AT * days_per_year)

                    Jeff_data_to_add ={
                        'predictor' : ['ChemSTEER User-Defined Model'],
                        'Predicted_Value' : [ADD],
                        'Quantity Type': ['Average of five highest median concentrations']
                    } 

                    Jeff_df = pd.DataFrame(Jeff_data_to_add)
                                        

            #Displaying figures of SEEM3 component model predictions     
            with st.container(border=True):
                    st.markdown('## SEEM3 Component Model Predictions')
                    st.markdown('##### SEEM3 is the weighted consensus of 13 exposure models (known as its "component models"). Paper describing model available [here](https://doi.org/10.1021/acs.est.8b04056): ')
                    st.markdown('##### Table 2 of the above SEEM3 paper describes the component models displayed in the two figures below, including links to their sources')
                    
                    dat_sources = expo.search_exposures(by="seem", dtxsid=structure_dtxsid)
                    dat_sources = pd.DataFrame(dat_sources)
                    
                    #1)Drop all SEEM2 rows
                    dat_sources = dat_sources[dat_sources['predictor'].str[0:5] != 'SEEM2']

                    #2)Add the label "5th percentile value" to all rows
                    dat_sources['Quantity Type'] = "5th percentile value"
                    dat_sources.rename(columns={'l95':'Predicted_Value'}, inplace=True)
                    
                    #3)Concatenate the datasource label and the median value after renaming 
                    lower_df = dat_sources[['predictor', 'median', 'units']].copy()
                    lower_df['Quantity Type'] = 'Median Value'
                    lower_df.rename(columns={'median':'Predicted_Value'}, inplace=True)
                    dat_sources_plus_median = pd.concat([dat_sources, lower_df], ignore_index=True)
                    
                    #4)Concatenate the datasource label and the upper value 
                    upper_df = dat_sources[['predictor','u95', 'units']].copy()
                    upper_df['Quantity Type'] = '95th Percentile Value'
                    upper_df.rename(columns={'u95':'Predicted_Value'}, inplace=True)
                    dat_sources_plus_medianandupper = pd.concat([dat_sources_plus_median, upper_df], ignore_index=True)
                    
                    #5)drop rows with no information
                    dat_sources_plus_medianandupper=dat_sources_plus_medianandupper[dat_sources_plus_medianandupper['Predicted_Value']!=0]
                    
                    #6)Copy the dataset. For one copy, drop all rows with a unit of "intake fraction" 
                    # and for the other copy, drop all rows with a unit of "mg/kg/day" 
                    #Not displaying any model predictions in other units 
                    mg_data = dat_sources_plus_medianandupper[dat_sources_plus_medianandupper['units']=='mg/kg/day'].copy()
                    #Need to append the average from Jeff's model here to mg_data 
                    intake_fraction_data = dat_sources_plus_medianandupper[dat_sources_plus_medianandupper['units']=='intake fraction'].copy()
                    

                    if not minucci_portion.empty:
                        mg_data = pd.concat([mg_data, Jeff_df], ignore_index=True)

                    mg_chart = alt.Chart(mg_data).mark_circle(filled=True, size=100).encode(
                        x=alt.X('predictor').title('Component Model'),
                        y=alt.Y('Predicted_Value', 
                                scale=alt.Scale(type='log', 
                                                domain=[(mg_data['Predicted_Value'].min())/10, 
                                                        (mg_data['Predicted_Value'].max())*10])).title('Predicted_Value (mg/kg/day)'),
                        color='Quantity Type',
                    ).configure(autosize='fit').configure_axis(labelLimit=1000)

                    intake_fraction_chart = alt.Chart(intake_fraction_data
                                                    ).mark_circle(filled=True, size=100).encode(
                        x=alt.X('predictor').title('Component Model'),
                        y=alt.Y('Predicted_Value', 
                                scale=alt.Scale(type='log', 
                                                domain=[(mg_data['Predicted_Value'].min())/10, 
                                                        (mg_data['Predicted_Value'].max())*10])).title('Predicted_Value (Intake Fraction)'),
                        color='Quantity Type',
                        
                        ).configure_axis(labelLimit=1000)
                    
                    st.header('Predictions in units of mg/kg/day')
                    st.altair_chart(mg_chart, use_container_width=True)
                    st.markdown('#####  Figure note: "ChemSTEER user-defined model" is the average of the top five highest exposure estimates '
                                'derived from the concentrations predicted by the "Predicted Occupational Exposure" model')






                    
