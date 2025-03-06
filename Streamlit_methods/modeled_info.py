import streamlit as st
import pandas as pd
from ctxpy import Exposure
import altair as alt


@st.cache_data
def minucci_loader(script_location):#Returns -> DataFrame, but declaring this in the function annotation throws an error
    """
    Loads in the Minucci model dataset.
    Args:
        script_location (str): the location of the script that called the module
    
    Returns:
        Dataframe: The Dataframe of the Minucci dataset.
    """
    #This version of consolidated_minucci has been compressed with the "snappy" algorithm, which is fast to retrieve.
    #Consolidated_minucci has been split into two parts so that they will be small enough to upload to github

    minucci_path1 = (script_location/"data"/"consolidated_minucci_snappy_1.parq")
    minucci_path2 = (script_location/"data"/"consolidated_minucci_snappy_2.parq")
    minucci_1 = pd.read_parquet(minucci_path1)
    minucci_2 = pd.read_parquet(minucci_path2)

    df = pd.concat([minucci_1, minucci_2], ignore_index=True)

    return df


def ChemSTEER_udm(air_conc):
    """
    Converts the concentration predictions to exposure estimates mg/kg/day
    Using the formula for average daily dose (mg/kg-day) from the user_defined_inhalation model. 
    Information is from the package CLOET.

    Args:
        air_conc(float): concentration in air to be converted

    Returns:
        float: predicted exposure dose from the input air concentration
    
    """
    air_conc = float(air_conc)
    #Where: 
    ED = 250 #integer, days exposed per year; 250 is all the weekdays in a year minus two weeks of vacation
    Cm = (10**(air_conc)) #None or float,  mass concentration of chemical in air (mg/m^3)
    b = 1.25 # volumetric inhalation rate; 0 <= b <=7.9. The default value suggested by CLOET is 1.25 m^3/hr
    h = 8.5 #: daily exposure duration (hrs/day)
    EY = 40 #integer, years of occupational exposure; 0 <= EY (default: 40 years)
    BW = 70 #float, body weight; (0 <= BW) (default: 70 kg)
    AT = 40 # float, averaging time (EY <= AT <= ATc); (default: 40 years)
    days_per_year = 365
    I = Cm * b * h
    ADD = (I * ED * EY) / (BW * AT * days_per_year)

    return ADD


def predicted_info(structure_dtxsid, script_location, osha_data):
    
    minucci=minucci_loader(script_location)
    expo = Exposure()

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
        minucci_portion = minucci[minucci['dtxsid']==structure_dtxsid].copy()                
        # and then display them in the SEEM3 model chart

        if minucci_portion.empty:
            st.markdown('### No modeled data available')
        else:
            minucci_portion.sort_values(by='log_mgm3_pred_50th', ascending=False, inplace=True)
            st.dataframe(minucci_portion)   

            
            # Grabs the top 5 highest predicted concentrations
            #If the list of predictions is shorter than 5, takes all of them
            num_preds = minucci_portion.shape[0]
            if num_preds>=5:
                top_jeff_cats = minucci_portion.iloc[0:5]
            else:
                top_jeff_cats = minucci_portion.iloc[0:num_preds]
            
            top_jeff_cats['median'] = top_jeff_cats['log_mgm3_pred_50th'].apply(ChemSTEERUdm)
            top_jeff_cats['upper_confidence_interval'] = top_jeff_cats['log_mgm3_pred_97.5th'].apply(ChemSTEERUdm)
            top_jeff_cats['lower_confidence_interval'] = top_jeff_cats['log_mgm3_pred_2.5th'].apply(ChemSTEERUdm)

            top_jeff_cats.rename(columns={'subsector_name':'predictor'}, inplace=True)
            top_jeff_cats['Quantity Type'] = 'Predicted Occupational Inhalation Exposures'
                
            #st.write('The top jeff data')
           # st.dataframe(top_jeff_cats)

            Jeff_df = pd.DataFrame(top_jeff_cats)

        #Displaying figures of SEEM3 component model predictions     
    with st.container(border=True):
                st.markdown('## SEEM3 Component Model Predictions')
                st.markdown('##### SEEM3 is the weighted consensus of 13 exposure models (known as its "component models"). Paper describing model available [here](https://doi.org/10.1021/acs.est.8b04056): ')
                st.markdown('##### Table 2 of the above SEEM3 paper describes the component models displayed in the two figures below, including links to their sources')
                
                dat_sources= expo.search_exposures(by="seem", dtxsid=structure_dtxsid)
                dat_sources= pd.DataFrame(dat_sources)
                
                dat_sources=dat_sources[dat_sources['predictor']=='SEEM3 Consensus']
                dat_sources.rename(columns={'l95': "lower_confidence_interval", 'u95':'upper_confidence_interval'}, inplace=True)
                dat_sources['Quantity Type'] = 'Predicted Consumer Exposures'

                if not minucci_portion.empty or not dat_sources.empty:
                    if minucci_portion.empty:
                        combined_dat = dat_sources
                    elif dat_sources.empty:
                        combined_dat = Jeff_df
                    else:
                        # Append the predictions from Jeff's model here to dat_sources 
                        combined_dat = pd.concat([dat_sources, Jeff_df], ignore_index=True)
                    # st.write('after merging')
                    # st.dataframe(combined_dat)
                    combined_bars = alt.Chart(combined_dat ).mark_errorbar().encode(
                        alt.X("upper_confidence_interval", scale=alt.Scale(type='log', 
                                                domain=[(combined_dat ['lower_confidence_interval'].min())/10, 
                                                        (combined_dat ['upper_confidence_interval'].max())*10])).title(''),
                        alt.X2("lower_confidence_interval"),
                        alt.Y('predictor:N')
                    )

                    combined_chart = alt.Chart(combined_dat ).mark_circle(filled=True, size=100).encode(
                        x=alt.X('median:Q', 
                                scale=alt.Scale(type='log', 
                                                domain=[(combined_dat ['median'].min())/10, 
                                                        (combined_dat ['median'].max())*10])).title('Predicted_Value (mg/kg/day)'),
                        y=alt.Y('predictor:N').title(''),
                        color='Quantity Type'
                    )#

                    # Choosing to not display intake-fraction data 
                    #intake_fraction_chart = alt.Chart(intake_fraction_data
                    #                               ).mark_circle(filled=True, size=100).encode(
                    #   x=alt.X('predictor').title('Component Model'),
                    #  y=alt.Y('Predicted_Value', 
                    #         scale=alt.Scale(type='log', 
                        #                        domain=[(dat_sources['Predicted_Value'].min())/10, 
                        #                               (dat_sources['Predicted_Value'].max())*10])).title('Predicted_Value (Intake Fraction)'),
                        #color='Quantity Type',
                        
                        #).configure_axis(labelLimit=1000)
                    final_chart = alt.layer(combined_bars, combined_chart).configure(autosize='fit')
                    st.header('Predictions in units of mg/kg/day')
                    st.altair_chart(final_chart, use_container_width=True) #, use_container_width=True
                    st.markdown('#####  Figure note: "ChemSTEER user-defined model" is the average of the top five highest exposure estimates '
                                'derived from the concentrations predicted by the "Predicted Occupational Exposure" model')
                else:
                    st.error('No data is available for figure')




