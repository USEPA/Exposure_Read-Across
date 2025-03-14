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
    def minucci_loader(_self):
        """
        Loads in the Minucci model dataset.
        Args:
            script_location (str): the location of the script that called the module
        
        Returns:
            Dataframe: The Dataframe of the Minucci dataset.
        """
        # This version of consolidated_minucci has been compressed with the "snappy" 
        # algorithm, which is fast to retrieve.
        # Consolidated_minucci has been split into two parts so that they will be 
        # small enough to upload to github

        minucci_path1 = (_self.script_location / "data" /
                         "consolidated_minucci_snappy_1.parq")
        minucci_path2 = (_self.script_location / "data" /
                         "consolidated_minucci_snappy_2.parq")
        minucci_1 = pd.read_parquet(minucci_path1)
        minucci_2 = pd.read_parquet(minucci_path2)

        df = pd.concat([minucci_1, minucci_2], ignore_index=True)

        return df

    def ChemSTEER_udm(self, air_conc):
        """
        Converts the concentration predictions to exposure estimates mg/kg/day
        Using the formula for average daily dose (mg/kg-day) from the 
        user_defined_inhalation model. 
        Information is from the package CLOET.

        Args:
            air_conc(float): concentration in air to be converted

        Returns:
            float: predicted exposure dose from the input air concentration
        
        """
        air_conc = float(air_conc)
        
        #integer, days exposed per year; 250 is all the weekdays in a year minus two
        # weeks of vacation
        ED = 250 
        
        # None or float,  mass concentration of chemical in air (mg/m^3)
        Cm = (10**(air_conc)) 
        
        # volumetric inhalation rate; 0 <= b <=7.9. The default value suggested by 
        # CLOET is 1.25 m^3/hr
        b = 1.25
        
        # daily exposure duration (hrs/day)
        h = 8.5
        
        #integer, years of occupational exposure; 0 <= EY (default: 40 years)
        EY = 40
        
        #float, body weight; (0 <= BW) (default: 70 kg)
        BW = 70
        
        # float, averaging time (EY <= AT <= ATc); (default: 40 years)
        AT = 40
        days_per_year = 365
        I = Cm * b * h
        ADD = (I * ED * EY) / (BW * AT * days_per_year)

        return ADD

    def qsur_display(self):
        # Displays Quantitative Structure-Use Relationship data
        expo = Exposure()
        self.modeled_container.header(' Predicted Information', divider=True)
        
        (self
         .modeled_container
         .markdown("Reported Information is always a gold standard for data. "
                   "However, as is often the case with exposure science, many data are "
                   "lacking for the thousands of chemicals humans interact with. For "
                   "these cases, EPA and others have developed *gap-filling* models to "
                   "predict data when no reported data are available. Use of predicted "
                   "information should acknowledge that these data are not perfect."))
        self.modeled_container.markdown("### Predicted Functional Use")
        qsur_doi = "https://doi.org/10.1039/c6gc02744j"
        (self
         .modeled_container
         .markdown('EPA Quantitative Structure-Use Relationship (QSUR) model '
                   'predictions. These are a set of machine-learning based models '
                   'that predict uses for chemicals based on their structures. See '
                   f'[Phillips et al. 2017]({qsur_doi}) for more information'))

        predicted_info = expo.search_qsurs(dtxsid=self.structure_dtxsid)
        predicted_info = pd.DataFrame(predicted_info)
        
        if predicted_info.empty:
            self.modeled_container.error('No Reported Function information available')
        else:
            (predicted_info
             .rename(columns={'harmonizedFunctionalUse':'Harmonized Functional Use'},
                     inplace=True))
            (predicted_info
             .sort_values(by=['probability'],
                          ascending=False,
                          inplace=True))
            self.modeled_container.dataframe(predicted_info)
        
    def minucci_model(self):
        # Output from Jeff's model
        minnuci_doi = "https://doi.org/10.1021/acs.est.2c08234"
        (self
         .modeled_container
         .markdown('### Predicted Workplace Inhalation Exposure Model'))
        (self
         .modeled_container
         .markdown('This model predicts 1) if a chemical is likely to be measured in a '
                   'particualr industrial sector and 2) if so, what is the likely air '
                   'concentration of the chemical in that industrial sector. More '
                   f'details can be found in [Minucci et al. 2021]({minnuci_doi}).'))
        minucci = self.minucci_loader()
        # Need to be able to share minucci_portion with seem3_and_minucci method
        idx = minucci['dtxsid'] == self.structure_dtxsid
        ModeledData.minucci_portion = minucci[idx].copy()

        if ModeledData.minucci_portion.empty:
            self.modeled_container.error('No modeled data available')
        else:
            (ModeledData
             .minucci_portion
             .sort_values(by='log_mgm3_pred_50th',
                          ascending=False,
                          inplace=True))
            self.modeled_container.dataframe(ModeledData.minucci_portion)

    def seem3_and_minucci(self):
        expo = Exposure()
        # Grabs the top 5 highest predicted concentrations
        # If the list of predictions is shorter than 5, takes all of them
        num_preds = ModeledData.minucci_portion.shape[0]
        if num_preds >= 5:
            top_cats = ModeledData.minucci_portion.iloc[0:5]
        else:
            top_cats = ModeledData.minucci_portion.iloc[0:num_preds]
        
        top_cats['median'] = (top_cats['log_mgm3_pred_50th']
                              .apply(self.ChemSTEER_udm))
        top_cats['upper_confidence_interval'] = (top_cats['log_mgm3_pred_97.5th']
                                                 .apply(self.ChemSTEER_udm))
        top_cats['lower_confidence_interval'] = (top_cats['log_mgm3_pred_2.5th']
                                                 .apply(self.ChemSTEER_udm))

        top_cats.rename(columns={'subsector_name':'predictor'}, inplace=True)
        top_cats['Quantity Type'] = 'Predicted Occupational Inhalation Exposures'

        Jeff_df = pd.DataFrame(top_cats)
        
        ring_doi = "https://doi.org/10.1021/acs.est.8b04056"
        self.modeled_container.markdown('## Comparison of Exposure Models')
        (self
         .modeled_container
         .markdown('SEEM3 uses a consensus modeling approach, combined with Bayesian '
                   'methods to first predict possible exposure pathways for a '
                   'chemical and then predict the exposure of the chemical for a '
                   'general population. Predictions from many different models '
                   '(both developed by EPA and by outside entities) are used to form '
                   'a consensus exposure prediction. For more details see '
                   f'[Ring et al. 2018]({ring_doi}).'))
        
        # (self
        #  .modeled_container
        #  .markdown("Comparison of exposures predicted with exposure doses calculated"
        #           "using the workplace air concentration model from Minucci et al. "
        #           "2021 and consensus general population exposure calculated from "
        #           "Ring et al. 2018. Workplace exposure doses are calculated by"
        #           "using the predicted concentration as input to ChemSTEER's ",
        #           "User-defined Inhalation model with all the default parameters "
        #           "accepted. Only results for the NAICS sectors with the largest"
        #           "predicted concentrations are show for simplicity."))
        
        dat_sources= expo.search_exposures(by="seem", dtxsid=self.structure_dtxsid)
        dat_sources= pd.DataFrame(dat_sources)
        
        dat_sources=dat_sources[dat_sources['predictor']=='SEEM3 Consensus']
        dat_sources.rename(columns={'l95': "lower_confidence_interval",
                                    'u95':'upper_confidence_interval'},
                           inplace=True)
        dat_sources['Quantity Type'] = 'Predicted Consumer Exposures'

        if not ModeledData.minucci_portion.empty or not dat_sources.empty:
            if ModeledData.minucci_portion.empty:
                combined_dat = dat_sources
            elif dat_sources.empty:
                combined_dat = Jeff_df
            else:
                # Append the predictions from Jeff's model here to dat_sources 
                combined_dat = pd.concat([dat_sources, Jeff_df], ignore_index=True)
            bar_range = [(combined_dat['lower_confidence_interval'].min())/10,
                          (combined_dat['upper_confidence_interval'].max())*10]
            combined_bars = (alt
                             .Chart(combined_dat )
                             .mark_errorbar()
                             .encode((alt.X("upper_confidence_interval",
                                            scale=alt.Scale(type='log',
                                                            domain=bar_range))
                                      .title('')),
                                     alt.X2("lower_confidence_interval"),
                                     alt.Y('predictor:N')))

            chart_range = [(combined_dat ['median'].min())/10, 
                           (combined_dat ['median'].max())*10]
            combined_chart = (alt
                              .Chart(combined_dat )
                              .mark_circle(filled=True, size=100)
                              .encode(x=(alt.X('median:Q',
                                               scale=(alt.Scale(type='log',
                                                                domain=chart_range)))
                                               .title('Predicted_Value (mg/kg/day)')),
                                      y=alt.Y('predictor:N').title(''),
                                      color='Quantity Type'))

            final_chart = (alt
                           .layer(combined_bars,
                                  combined_chart)
                           .configure(autosize='fit'))
           # self.modeled_container.header('Predictions in units of mg/kg/day')
            self.modeled_container.altair_chart(final_chart,
                                                use_container_width=True)
            
            # self.modeled_container.markdown('#####  Figure note: "ChemSTEER user-defined model" is the average of the top five highest exposure estimates '
            #             'derived from the concentrations predicted by the "Predicted Occupational Exposure" model')
        else:
            self.modeled_container.error('No data is available for figure')
        (self
         .modeled_container
         .caption("Comparison of exposures predicted with exposure doses calculated "
                  "using the workplace air concentration model from Minucci et al. "
                  "2021 and consensus general population exposure calculated from "
                  "Ring et al. 2018. Workplace exposure doses are calculated by "
                  "using the predicted concentration as input to ChemSTEER's "
                  "User-defined Inhalation model with all the default parameters "
                  "accepted. Only results for the NAICS sectors with the largest "
                  "predicted concentrations are show for simplicity."))




