import streamlit as st
import pandas as pd
from ctxpy import Exposure
import altair as alt


class ModeledData:
    def __init__(self, structure_dtxsid, script_location, osha_set):
        self.structure_dtxsid = structure_dtxsid
        self.script_location = script_location
        self.osha_set = osha_set
        self.modeled_container = st.container(border=True)

    @st.cache_data
    def minucci_loader(_self):  # Returns -> DataFrame, but declaring this in the function annotation throws an error
        """
        Loads in the Minucci model dataset.
        Args:
            script_location (str): the location of the script that called the module
        
        Returns:
            Dataframe: The Dataframe of the Minucci dataset.
        """
        # This version of consolidated_minucci has been compressed with the "snappy" algorithm, which is fast to retrieve.
        # Consolidated_minucci has been split into two parts so that they will be small enough to upload to github

        minucci_path1 = (_self.script_location/"data"/"consolidated_minucci_snappy_1.parq")
        minucci_path2 = (_self.script_location/"data"/"consolidated_minucci_snappy_2.parq")
        minucci_1 = pd.read_parquet(minucci_path1)
        minucci_2 = pd.read_parquet(minucci_path2)

        df = pd.concat([minucci_1, minucci_2], ignore_index=True)

        return df

    def ChemSTEER_udm(self, air_conc):
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

    def qsur_display(self):
        # Displays Quantitative Structure-Use Relationship data
        expo = Exposure()
        self.modeled_container.header(' Predicted Information', divider=True)
        self.modeled_container.markdown('#### To allow more-informed decision-making, '
                    'model predictions can fill information gaps for substances that are in DSSTox ' 
                    'but do not have sufficient exposure information available in the fields above.')
        self.modeled_container.markdown("### Predicted Functional Use")
        self.modeled_container.markdown('#### Data pulled from EPA QSUR models. Paper available [here](https://pubs.rsc.org/en/content/articlelanding/2017/gc/c6gc02744j#!divCitation):')
        predicted_info = expo.search_qsurs(dtxsid=self.structure_dtxsid)
        predicted_info = pd.DataFrame(predicted_info)
        if predicted_info.empty:
            self.modeled_container.markdown('### No Reported Function information available')
        else:
            predicted_info.rename(columns={'harmonizedFunctionalUse':'Harmonized Functional Use'},
                                inplace=True)
            predicted_info.sort_values(by=['probability'], ascending=False, inplace=True)
            self.modeled_container.dataframe(predicted_info)
        
    def minucci_model(self):
        # Output from Jeff's model
        self.modeled_container.markdown('### Predicted Occupational Exposure')
        self.modeled_container.markdown('#### Paper describing model available [here](https://pubs.acs.org/doi/10.1021/acs.est.2c08234):')
        minucci = self.minucci_loader()
        # Need to be able to share minucci_portion with seem3_and_minucci method
        ModeledData.minucci_portion = minucci[minucci['dtxsid'] == self.structure_dtxsid].copy()

        if ModeledData.minucci_portion.empty:
            self.modeled_container.markdown('### No modeled data available')
        else:
            ModeledData.minucci_portion.sort_values(by='log_mgm3_pred_50th', ascending=False, inplace=True)
            self.modeled_container.dataframe(ModeledData.minucci_portion)   

    def seem3_and_minucci(self):
        expo = Exposure()
        # Grabs the top 5 highest predicted concentrations
        # If the list of predictions is shorter than 5, takes all of them
        num_preds = ModeledData.minucci_portion.shape[0]
        if num_preds >= 5:
            top_jeff_cats = ModeledData.minucci_portion.iloc[0:5]
        else:
            top_jeff_cats = ModeledData.minucci_portion.iloc[0:num_preds]
        
        top_jeff_cats['median'] = top_jeff_cats['log_mgm3_pred_50th'].apply(self.ChemSTEER_udm)
        top_jeff_cats['upper_confidence_interval'] = top_jeff_cats['log_mgm3_pred_97.5th'].apply(self.ChemSTEER_udm)
        top_jeff_cats['lower_confidence_interval'] = top_jeff_cats['log_mgm3_pred_2.5th'].apply(self.ChemSTEER_udm)

        top_jeff_cats.rename(columns={'subsector_name':'predictor'}, inplace=True)
        top_jeff_cats['Quantity Type'] = 'Predicted Occupational Inhalation Exposures'

        Jeff_df = pd.DataFrame(top_jeff_cats)
        
        self.modeled_container.markdown('## SEEM3 Component Model Predictions')
        self.modeled_container.markdown('##### SEEM3 is the weighted consensus of 13 exposure models (known as its "component models"). Paper describing model available [here](https://doi.org/10.1021/acs.est.8b04056): ')
        self.modeled_container.markdown('##### Table 2 of the above SEEM3 paper describes the component models displayed in the two figures below, including links to their sources')
        
        dat_sources= expo.search_exposures(by="seem", dtxsid=self.structure_dtxsid)
        dat_sources= pd.DataFrame(dat_sources)
        
        dat_sources=dat_sources[dat_sources['predictor']=='SEEM3 Consensus']
        dat_sources.rename(columns={'l95': "lower_confidence_interval", 'u95':'upper_confidence_interval'}, inplace=True)
        dat_sources['Quantity Type'] = 'Predicted Consumer Exposures'

        if not ModeledData.minucci_portion.empty or not dat_sources.empty:
            if ModeledData.minucci_portion.empty:
                combined_dat = dat_sources
            elif dat_sources.empty:
                combined_dat = Jeff_df
            else:
                # Append the predictions from Jeff's model here to dat_sources 
                combined_dat = pd.concat([dat_sources, Jeff_df], ignore_index=True)
        
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
            )

            final_chart = alt.layer(combined_bars, combined_chart).configure(autosize='fit')
            self.modeled_container.header('Predictions in units of mg/kg/day')
            self.modeled_container.altair_chart(final_chart, use_container_width=True)
            self.modeled_container.markdown('#####  Figure note: "ChemSTEER user-defined model" is the average of the top five highest exposure estimates '
                        'derived from the concentrations predicted by the "Predicted Occupational Exposure" model')
        else:
            self.modeled_container.error('No data is available for figure')




