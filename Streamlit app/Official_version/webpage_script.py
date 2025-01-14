import pandas as pd
import numpy as np
#from pandas.io.formats.style import Styler
import streamlit as st
from streamlit_ketcher import st_ketcher
from pathlib import Path
#from typing import Optional, Tuple
from chembl_webresource_client.new_client import new_client as ch
from rdkit import Chem
from ctxpy import Chemical, Exposure
import altair as alt
from sklearn.neighbors import KNeighborsRegressor
from rdkit.Chem import rdFingerprintGenerator
#import sys



counter=0
def generate_ecfp(smiles: str, morg_fing):# -> numpy.ndarray
    """"
    Generates the Morgan fingerprint of the target chemical.

    Args: 
        smiles (str): The smiles string of the target chemical

    Returns:
        numpy.ndarray: The Morgan fingerprint of the target substance
    
    """

    global counter
    print(counter)
    counter+=1

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    return morg_fing.GetFingerprintAsNumPy(molecule)


def analog_finder(MFPs_address, chem_info_address, chem_smiles, num_features):
    #What inputs are needed that I will have to find a match for in the Shah-2016 situation?
        #Morgan fingerprints of all substances in the match-source dataset
        #Dataset of additional information on chemicals 

    #Chose feather as the datatype because it was the fastest datatype to read and write 
    if MFPs_address.suffix==".feather":
        FP1 = pd.read_feather(MFPs_address)
    if MFPs_address.suffix==".parq":
        FP1 = pd.read_parquet(MFPs_address)
    consolidated_dsstox = pd.read_parquet(chem_info_address)
    consolidated_dsstox.rename(columns={'DTXSID':'ID'}, inplace=True)
    mfgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=num_features)
    interest_fp = pd.DataFrame(np.vstack(np.array(generate_ecfp(chem_smiles, mfgen))))
    interest_fp = interest_fp.T

    KNR = KNeighborsRegressor(n_neighbors=10, algorithm='auto', leaf_size=30, p=2, metric='jaccard', 
                        metric_params=None, n_jobs=-1,)

    Y = np.ones(FP1.shape[0])
    KNR.fit(FP1, Y)

    #I2I assigns a numerical index value to the DTXSID index of FP1
    I2I = dict(zip(range(FP1.shape[0]),FP1.index))
    
    #This is the line that actually assigns the nearest neighbors and outputs the 
    #dataframe with ID of substance of interest and list of NN's 
    #whole_Sim is the similarities of the top ten most similar
    #whole_Ind is the indices in FP1 of the nearest neighbors  
    whole_Sim, whole_Ind = KNR.kneighbors(interest_fp)
    
    #While the target substance cannot be referenced by DTXSID, the selected structural analogs can.
    
    #Since whole_Sim and whole_Ind are for only a single chemical, 
    # the relevant target chemical does not need to be pulled out by index
    Sim = whole_Sim.T.flatten() 
    
    #Need to generate Ind by pulling the DTXSID's of the I2I values at the index values of Ind
    Ind = [I2I[index_val] for index_val in whole_Ind.T.flatten()]

    k = 10
    NNi= pd.DataFrame(dict(ID=Ind,sim=Sim)).iloc[:(k+1)]
    NNi = NNi.merge(consolidated_dsstox,on='ID')
    NNi["sim"] = 1-NNi['sim']
    

    return NNi
    
@st.cache_data
def import_cehd():#Returns -> DataFrame, but declaring this in the function annotation throws an error
    #Column datatypes assigned so that read-in will be less likely to cause memory issues
    """
    Reads in the CEHD dataset

    Returns:
        Dataframe: The Dataframe of CEHD data

    """
    #Must use an absolute path so that the unit-testing script is able to read it 
    osha=pd.read_parquet(r"C:\Users\svanambe\OneDrive - Environmental Protection Agency (EPA)\Profile\Work\Projects\Exposure web app\Streamlit app\Official_version\FullCEHD_includingDTXSID.parq") #Dataset generated in the file "CEHD_withDTXSID_generator.ipynb"
                   
    return osha
cehd=import_cehd()


@st.cache_data
def minucci_loader():#Returns -> DataFrame, but declaring this in the function annotation throws an error
    
    """
    Loads in the Minucci model dataset.

    Returns:
        Dataframe: The Dataframe of the Minucci dataset.
    """
    #This version of consolidated_minucci has been compressed with the "snappy" algorithm, which is fast to retrieve
    #consolidated_minucci has been split into two parts so that they will be small enough to upload to github

    #This comment is a test change

    minucci_path1 = Path(r"C:\Users\svanambe\OneDrive - Environmental Protection Agency (EPA)\Profile\Work\Projects\Exposure web app\Streamlit app\Official_version\predictions_with_extrapolation\consolidated_minucci_snappy_1.parq")
    minucci_path2 = Path(r"C:\Users\svanambe\OneDrive - Environmental Protection Agency (EPA)\Profile\Work\Projects\Exposure web app\Streamlit app\Official_version\predictions_with_extrapolation\consolidated_minucci_snappy_2.parq")
    minucci_1 = pd.read_parquet(minucci_path1)
    minucci_2 = pd.read_parquet(minucci_path2)

    df = pd.concat([minucci_1, minucci_2])


    return df


#Converts measurements in mg/m^3 into units of ppm,
# assuming pressure of 1 atm and temperature of 20 C=68F
def mg_to_ppm(value, unit, air_volume, dtxsid)->float:
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

    try:
        air_volume = float(air_volume)
    except:
        return 0

    value = float(value) 
    
    if air_volume == 0:
        return  

    if unit == 'X':
        #Converting air volume in liters into m^3  
        air_volume = float(air_volume) / 1000
        #Converting micrograms to milligrams and then dividing by air volume 
        # to get mg/m^3 
        value = value / (1000 * air_volume)
    elif unit=='Y':
        #Converting air volume in liters into m^3  
        air_volume = air_volume / 1000
        #Dividing by air volume to get mg/m^3 
        value = value / air_volume

    if unit!='P':
        deets=chem.details(by='dtxsid', word=dtxsid)   
        molar_mass=deets['monoisotopicMass']     
        #Will not include derivation of the below formula here, 
        # but I have derived it by hand and compared to online sources to insure accuracy                             
        ppm_value=(24.06*value)/molar_mass
        return ppm_value  
    else:
        return value


expo = Exposure(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')
chem = Chemical(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')

# Need to switch back to using the terminology of "New Chemical" and "Existing Chemical" 

#Am titling the "Substance Structural Analogs Data" tab as such because it will eventually include summaries of the data connected to the analogs
tab1, tab2 = st.tabs(["Substance Structural Analogs Data", "Substance Data for Substances in DSSTox"])#(["New Chemical (not in DSSTox)", "Existing Chemical (in DSSTox)"])

minucci=minucci_loader()

with tab1:
    structure_dtxsid = False
    
    with st.form("NewChem"):
        st.markdown("## New Chemical Input")
        st.markdown("### Structure")
        st.markdown("Draw the structure in the provided space below and then click "
                    "`Apply` to search for structural analogs for the chemical.")
        smiles_ready = False
        smiles_code = st_ketcher() #This is where the smiles code is returned 

        submitted = st.form_submit_button("Submit")
        
        if submitted:
            compute = True
            if len(smiles_code) == 0:
                compute = False
                st.warning("Missing SMILES string, did you forget to draw a chemical structure or not click `Apply` after drawing?")
            else:
                smiles_ready=True 

    if submitted:
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
                elif len(hits)<1 :
                    st.warning("No information can be found on a chemical with this structure. Check that the structure you drew is correct.")
                    raise Exception 
                
                try:
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
    
    if smiles_ready:    
        st.header("This is the section for the structural analogs")
        
        #First thing is to assign the Morgan fingerprints to the smiles code 
        FP_input_path = "./Morgan_fingerprints_of_DSSTox.feather"
        cd_input_path = "./brotli_Consolidated_DSSTox_QSAR_smiles_only.parq"
        final_table=analog_finder(FP_input_path, cd_input_path, smiles_code, 2048)
        final_table.drop(columns=['index'], inplace=True )
        final_table.rename(columns={'ID':'DTXSID', 'sim':'Similarity'}, inplace=True)
        with st.container(border=True):
            st.dataframe(final_table)

with tab2:
    with st.container(border=True):
        st.markdown('### This page is useful in two scenarios:' )
        st.markdown('###     When the new substance of interest is in the DSSTox database, but '
                    'not in the TSCA active inventory. ')
        st.markdown('###     Or, when the substance of interest is not in the DSSTox database '
                    'and information is desired on the structural analogues of the substance.')
        
    with st.form("ExistingChem"):
        st.markdown("## Existing Chemical Input")
        st.markdown("### Structure")
        st.markdown("Draw the structure in the provided space below and then click "
                    "`Apply` to retrieve information on the chemical.")
        smiles_ready = False
        smiles_code = st_ketcher() #This is where the smiles code is returned 
        
        #Needs to be something other than "submitted", so that it is different than the tab1 conditional
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
                #elif hits["title"]=="Bad Request":
                 #   st.warning("No information can be found on a chemical with this structure. Check that the structure you drew is correct.")
                  #  raise Exception 
                
                

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
                    
                    st.markdown('## Measured Usage Information ')
                    #st.markdown('### All entries in CPDAT on the entered substance')
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
                        
                        docs_per_function_cat = alt.Chart(grouped_cfi).mark_bar().encode(
                            x='Function Category',
                            y="Docs per function",
                            color=alt.Color('Reported Function', scale=alt.Scale(scheme='category20') )
                        ).configure_axis(labelLimit=1000)

                        #color=alt.Color('Reported Function:N', scale=alt.Scale(scheme='category20') )
                        
                        
                        st.altair_chart(docs_per_function_cat)
                        st.markdown("##### Figure note: Collected documents containing entered substance grouped by reported function within a function category. " 
                                    "Only data with a reported function has been included. "
                                    "This data facilitates the understanding of the uses of the entered substance by displaying the substance's reported functions."
                                    )

                        #Display Chemical name and DTXSID once 
                        #st.dataframe(function_info, width=1000000)
                    #with st.container(border=True):  
            
                    #Now pulling predicted functional use
                    
                    #Now pulling list presence information 
                    st.markdown('### List Presence Information')
                    st.markdown('##### The lists on which the entered substance is present. ' 
                                'These provide additional information about the contexts in which the substance is used.')
                    #st.markdown('###### Data pulled from CPDAT by DTXSID')
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
                    #st.markdown('###### Data pulled from CPDAT by DTXSID')
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
                            x='PUC general category + product family combination:N',
                            y='Number of product names:Q',
                            )#.configure_axis(labelLimit=1000) When this method is used, the labels intrude up into the chart and replace it. 
                            #Probably the height of the chart will need to be defined as well

                            #.properties( width=1000,height=1000)
                                #color='Reported Function:N'
                            
                            #st.dataframe(grouped_pucs)
                            st.altair_chart(prods_per_catfam)
                            st.markdown('##### Figure note: Number of products associated with each combination of PUC general category and product family')
                           # st.dataframe(pucs_available)
                        except:
                            st.write("No product use information is available")

                with st.container(border=True):
                    #st.markdown('## Measured Exposure Information')
                    cehd_of_interest=cehd[cehd['DTXSID']==structure_dtxsid]

                    if cehd_of_interest.empty:
                        st.markdown('### No exposure data available')
                        empty_cehd=True
                    else:
                        empty_cehd=False
                        #Include a note that all zero-values have been dropped
                        cehd_of_interest.drop('Unnamed: 0', axis=1, inplace=True)
                        df_cehd_of_interest=cehd_of_interest.to_csv(index=False)
                        st.markdown('### OSHA CEHD data on the entered substance')
                        st.download_button(
                                        label="Download measured exposure data as CSV",
                                        data=df_cehd_of_interest,
                                        file_name="CEHD_data_on_chemical_of_interest.csv"
                                        )
                        #This dataframe will be displayed as a box and whisker plot
                            
                        #First three digits of NAICS code is subsector 
                            #Three sectors are represented by a range of 2-digit codes:
                                #Manufacturing (31-33), Retail Trade (44-45) and Transportation and Warehousing (48-49)
                                #Instead of grouping those above categories, just display all subsectors and display the above as a message
                       
                        # Creating new dataset with a column for the subsectors 
                        NAICS_subs = cehd_of_interest.copy()
                        NAICS_subs['naics_subsector'] = NAICS_subs['naics_code'].apply(lambda NAICS: str(NAICS)[0:3] if len(str(NAICS))>=3 else NAICS)
                        
                        #Dropping entries with a subsector of 0.0 or nan
                        NAICS_subs=NAICS_subs[NAICS_subs['naics_subsector'] != '0']
                        NAICS_subs=NAICS_subs[NAICS_subs['naics_subsector'] != '0.0']
                        NAICS_subs=NAICS_subs[NAICS_subs['naics_subsector'] != 'nan']

                        #Dropping entries that are labeled as blanks, as these don't count as samples 
                        NAICS_subs=NAICS_subs[NAICS_subs['qualifier'] != 'BLK']
                        
                        #Choosing to drop rows whose sample_type is not present in the official list of sample types
                        NAICS_subs = NAICS_subs[NAICS_subs['sample_type'].isin(['P', 'A', 'B', 'W'])] 

                        #Dropping entries with a value of zero, as these do not contain useful information
                        NAICS_subs=NAICS_subs[NAICS_subs['sample_result'] != '0']
                        NAICS_subs=NAICS_subs[NAICS_subs['sample_result'] != '0.0']
                        
                        #Dropping rows with units that aren't in the interpretable list or that are missing
                        #Is not necessary to have a sepparate clause for null entries, 
                        # as those are classified as not being in the list.
                        #Leaving out measurements with units of bulk % and fibers/cc
                        NAICS_subs_to_ppm = NAICS_subs[NAICS_subs['unit_of_measurement'].isin(['M', 'P', 'X', 'Y'])]
                
                        #And convert all units to ppm in a function
                        NAICS_subs_to_ppm['sample_result'] = NAICS_subs_to_ppm.apply(lambda x: mg_to_ppm(x.sample_result, x.unit_of_measurement, x.air_volume_sampled, structure_dtxsid), axis=1)
                     
                        #Change the labeled unit type when the conversion is done
                        NAICS_subs_to_ppm['unit_of_measurement'] = 'P'

                        #The units that are convertable to ppm:
                            #M - mg/m3, #X- micrograms, Y - milligrams,
                        #Not convertable:
                        # F - fibers/cc, % - percentage of bulk material    
                                      
                        if not NAICS_subs_to_ppm.empty:
                            st.markdown('#### Occupational inhalation exposure data on the entered substance, converted to units of PPM')
                            st.markdown('##### Data is separated by NAICS subsector. All entries with units that cannot be converted to PPM or which have a non-inhalation sample type are not displayed.')
                            box_and_whisker_ppm = alt.Chart(NAICS_subs_to_ppm).mark_boxplot().encode(
                            x= alt.X('naics_subsector:N').title('NAICS subsector'),
                            y= alt.Y('sample_result:Q').title('Air concentration (ppm)')
                                #tooltip = len('naics_subsector'))
                            )
                            st.altair_chart(box_and_whisker_ppm)
                            st.markdown('##### Figure note: Chart type is box-and-whisker, but for some substances only outliers can be seen')
                        else:   
                            st.markdown('### No data with units of ppm available')
                    
                        #st.dataframe(NAICS_subs_to_ppm)

                with st.container(border=True):
                    st.header(' Predicted Information', divider=True )
                    st.markdown('#### To allow more-informed decision-making, '
                                'model predictions can fill information gaps for substances that are in DSSTox ' 
                                'but do not have sufficient exposure information available in the fields above.')
                    st.markdown("### Predicted Functional Use")
                    #Could try to add footnote
                    st.markdown('#### Data pulled from EPA QSUR models. Paper available here:')
                    st.link_button('QSUR Model Paper', "https://pubs.rsc.org/en/content/articlelanding/2017/gc/c6gc02744j#!divCitation")
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
                    st.markdown('#### Paper describing model available here:')
                    st.link_button('Occupational Exposure Model Paper', "https://pubs.acs.org/doi/10.1021/acs.est.2c08234")

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
                        st.markdown('##### SEEM3 is the weighted consensus of 13 exposure models (known as its "component models") ')
                        st.link_button('SEEM3 Consensus Model Paper', 'https://doi.org/10.1021/acs.est.8b04056')
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
                            # y='Predicted_Value',
                            color='Quantity Type',
                        ).configure(autosize='fit').configure_axis(labelLimit=1000) #.interactive(), tooltip=['Predicted_Value']

                        intake_fraction_chart = alt.Chart(intake_fraction_data
                                                        ).mark_circle(filled=True, size=100).encode(
                            x=alt.X('predictor').title('Component Model'),
                            y=alt.Y('Predicted_Value', 
                                    scale=alt.Scale(type='log', 
                                                    domain=[(mg_data['Predicted_Value'].min())/10, 
                                                            (mg_data['Predicted_Value'].max())*10])).title('Predicted_Value (Intake Fraction)'),
                            color='Quantity Type',
                            
                            ).configure_axis(labelLimit=1000)#.interactive(), tooltip=['Predicted_Value']

                      
                        st.header('Predictions in units of mg/kg/day')
                        st.altair_chart(mg_chart, use_container_width=True)
                        st.markdown('#####  Figure note: "ChemSTEER user-defined model" is the average of the top five highest exposure estimates '
                                    'derived from the concentrations predicted by the "Predicted Occupational Exposure" model')
                        st.header('Predictions in units of Intake Fraction')
                        st.altair_chart(intake_fraction_chart, use_container_width=True)























"""
 dtype={
                        'establishment_name':str,
                        'city':str,
                        'state':str,
                        'zip_code':float,
                        'sic_code':str,
                        'naics_code':str,
                        'sampling_number':str,
                        'office_id':float,
                        'date_sampled':str,
                        'date_reported':str,
                        'eight_hour_twa_calc':str,
                        'lab_number':str,
                        'field_number':str,
                        'sample_type':str,
                        'blank_used':str,
                        'time_sampled':float,
                        'air_volume_sampled':str,
                        'sample_weight':str,
                        'imis_substance_code':str,
                        'substance':str,
                        'sample_result':str,
                        'unit_of_measurement':str,
                        'qualifier':str,
                        'file_name':str,
                        'instrument_type':str,
                        'media':str,
                        'species':str,
                        'country':str,
                        'sample_year': str,
                        'sample_month':str,
                        'QA_flags':str,
                        'FOUND_BY':str,
                        'DTXSID':str,
                        'PREFERRED_NAME':str
                        }"""