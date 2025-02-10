import streamlit as st
import pandas as pd
from ctxpy import Exposure
import altair as alt

@st.cache_data
def minucci_loader(script_location):#Returns -> DataFrame, but declaring this in the function annotation throws an error
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



def predicted_info(structure_dtxsid, script_location):
    
    minucci=minucci_loader(script_location)
    expo = Exposure(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')

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
