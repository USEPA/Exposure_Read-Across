import streamlit as st
from rdkit import Chem
from sklearn.neighbors import KNeighborsRegressor
import pandas as pd
from rdkit.Chem import rdFingerprintGenerator
import numpy as np
from ctxpy import Chemical, Exposure
import altair as alt

class analog_operations:

    def ppm_to_mg (self, value, unit, mw)->float:

        """
        Converts multiple starting units to ppm, assuming pressure of 1 atm and temperature of 25 C=77F

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


    def generate_ecfp(self, smiles: str, morg_fing):# -> numpy.ndarray
        """
        Generates the Morgan fingerprint of the target chemical.

        Args: 
            smiles (str): The smiles string of the target chemical

        Returns:
            numpy.ndarray: The Morgan fingerprint of the target substance
        
        """
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return None
        return morg_fing.GetFingerprintAsNumPy(molecule)


    def analog_finder(self, MFPs_address, chem_info_address, chem_smiles=None, unit_test=False, shah_id=None): # -> DataFrame
        """
        Generates the list of analogs for the target chemical.

        Args:
            MFPs_address (str): File address of the set of morgan fingerprints from which analogs are drawn

            chem_info_address (str): File address of the set of additional information on the chemicals from which analogs are drawn

            chem_smiles (str): the smiles string for the chemical of interest 

            unit_test (bool): A flag that denotes whether the function is being accessed by the unit-testing suite 

            shah_id (str): Index value of target chemical in Shah SI dataset 

        Returns:
            DataFrame: the set of analogs of the target chemical 

        """
        #Chose feather as the datatype because it was the fastest datatype to read and write 
        if str(MFPs_address).endswith(".feather"):
            FP1 = pd.read_feather(MFPs_address)
        #The files for unit testing are .parq type, requiring the conditional below
        if str(MFPs_address).endswith(".parq"):
            FP1 = pd.read_parquet(MFPs_address)
        consolidated_dsstox = pd.read_parquet(chem_info_address)
        consolidated_dsstox.rename(columns={'DTXSID':'ID'}, inplace=True)

        #The below two conditionals are necessary to avoid an error with pytest 
        if unit_test==True:
            mfgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
        if unit_test==False:
            mfgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

        if unit_test == False:
            interest_fp = pd.DataFrame(np.vstack(np.array(self.generate_ecfp(chem_smiles, mfgen))))
            interest_fp = interest_fp.T
            print(interest_fp)
        if unit_test == True:
            #Transposition necessary because dataframe is single column, when needs to be a row 
            interest_fp=pd.DataFrame(FP1.loc[shah_id]).T
            print(interest_fp)

        KNR = KNeighborsRegressor(n_neighbors=10, algorithm='auto', leaf_size=30, p=2, metric='jaccard', 
                            metric_params=None, n_jobs=-1,)

        Y = np.ones(FP1.shape[0])
        KNR.fit(FP1, Y)

        #I2I assigns a numerical index value to the DTXSID index of FP1
        I2I = dict(zip(range(FP1.shape[0]),FP1.index))
        
        #This is the line that actually assigns the nearest neighbors and outputs the 
        # dataframe with ID of substance of interest and list of NN's 
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

    def analog_retrieve(self, script_location, smiles_code):  
        
        #Assign the Morgan fingerprints to the smiles code 
        FP_input_path = (script_location/"data"/"Morgan_fingerprints_of_DSSTox.feather")
        cd_input_path = (script_location/"data"/"brotli_Consolidated_DSSTox_QSAR_smiles_only.parq")
        final_table=self.analog_finder(FP_input_path, cd_input_path, chem_smiles=smiles_code)
        final_table.drop(columns=['index'], inplace=True )
        final_table.rename(columns={'ID':'DTXSID', 'sim':'Similarity'}, inplace=True)       
        with st.container(border=True):
            st.header("Structural Analogs")
            st.markdown('##### Analogs are pulled from the set of substances in DSSTox that possess a QSAR-ready SMILES')
            analogs_csv = final_table.to_csv(index=False)
            st.download_button(
                                label="Download analog table",
                                data=analogs_csv,
                                file_name="analogs.csv"
                               ) 
            st.dataframe(final_table)
        return final_table
    

    def usis_summary_fig(self, returned_table, usis_data):
        chem = Chemical(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')
        input_set=pd.DataFrame()
        for dtxsid, chem_name in zip (returned_table['DTXSID'], returned_table['PREFERRED_NAME']):
            usis = usis_data[usis_data['dtxsid']==dtxsid]
            # Cleaning USIS for use in the analog-summarizing program
            usis.drop('index', axis=1, inplace=True)
            usis=usis[usis['exposure_level'] != '0']
            usis=usis[usis['exposure_level'] != '0.0']
            usis = usis[usis['measure_unit_id'].isin(['M', 'P'])]
            deets = chem.details(by='dtxsid', word=dtxsid)   
            molar_mass = deets['monoisotopicMass']    
            usis['exposure_level'] = usis.apply(lambda x: self.ppm_to_mg (x.exposure_level, x.measure_unit_id, molar_mass), axis=1) 
            usis['measure_unit_id'] = 'M'
            usis['measure_unit_name'] = 'Milligrams per cubic meter'
            usis = usis[usis['exposure_level']!='0'] 

            grouped_bysubs = usis.groupby('naics_2022_subsector_title')['exposure_level'].median().reset_index()

            grouped_bysubs['substance_name'] = chem_name
            grouped_bysubs ['Presence'] = '1'
            if not grouped_bysubs.empty:
                input_set = pd.concat([input_set, grouped_bysubs])
            else:
                dummy_dict = {'substance_name':[chem_name]}
                dummy_df =pd.DataFrame(dummy_dict)
                input_set = pd.concat([input_set, dummy_df])

        analog_htmp = alt.Chart(input_set).mark_rect().encode(
             x = alt.X('substance_name:N', 
              axis=alt.Axis(title=''),
              scale=alt.Scale(domain=list(input_set['substance_name'].unique()))), 
            y = alt.Y('naics_2022_subsector_title:N',
                       axis=alt.Axis(title = 'NAICS Subsector',
                                     titleX=-300)),
            color = alt.Color('Presence:Q', legend=None),
            tooltip = 'exposure_level:Q').configure_axis(labelLimit=1000)

        st.altair_chart(analog_htmp)
   
   #Annotations and axis labels provide information 
   # on meaning of variables and structure
    def osha_cpdat_cdr_seem(self, returned_table, usis_data, run_from_location):
       
        expo = Exposure(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')
        input_set=pd.DataFrame()

        with st.container(border=True):
            st.header("Summary of Structural Analog Data ")
            st.markdown("#### Availability of data on analogs from CPDat and USIS.")
            
            #Iterates through the DataFrame of analogs, 
            # adding data from each source on to a common dataframe 
            for dtxsid, chem_name in zip (returned_table['DTXSID'], returned_table['PREFERRED_NAME']):
                #USIS data
                usis = usis_data[usis_data['dtxsid']==dtxsid]
                num_usis_measurements = usis.shape[0]
                usis_addition = pd.DataFrame({'Substance_Name':[chem_name],
                                    'Num_records':[num_usis_measurements],
                                    'record_type':['Record in USIS']
                                    })
                input_set = pd.concat([input_set, usis_addition])
                #List presence data 
                list_presence = expo.search_cpdat(vocab_name='lpk', dtxsid=dtxsid)
                list_presence_df = pd.DataFrame(list_presence)
                num_lpk_measurements = list_presence_df.shape[0]
                lpk_addition = pd.DataFrame({'Substance_Name':[chem_name],
                                    'Num_records':[num_lpk_measurements],
                                    'record_type':['CPDat: List Presence Keywords']
                                    })
                input_set = pd.concat([input_set, lpk_addition])
                #PUC data
                puc_presence = expo.search_cpdat(vocab_name='puc', dtxsid=dtxsid)
                puc_presence_df = pd.DataFrame(puc_presence)
                num_puc_measurements = puc_presence_df.shape[0]
                puc_addition = pd.DataFrame({'Substance_Name':[chem_name],
                                    'Num_records':[num_puc_measurements],
                                    'record_type':['CPDat: PUCS']
                                    })
                input_set = pd.concat([input_set, puc_addition])
                #Function category data 
                fc_presence = expo.search_cpdat(vocab_name='fc', dtxsid=dtxsid)
                fc_presence_df = pd.DataFrame(fc_presence)
                num_fc_measurements = fc_presence_df.shape[0]
                fc_addition = pd.DataFrame({'Substance_Name':[chem_name],
                                    'Num_records':[num_fc_measurements],
                                    'record_type':['CPDat: Functional Use']
                                    })
                input_set = pd.concat([input_set, fc_addition])

            cpdat_htmp = alt.Chart(input_set).mark_rect().encode(
                x = alt.X('Substance_Name:N', 
                axis=alt.Axis(title=''),
                scale=alt.Scale(domain=list(input_set['Substance_Name'].unique()))), 
                y = alt.Y('record_type:N', title=''),
                color = alt.Color('Num_records:Q').scale(scheme='turbo').title("Number of Records")
                ).properties(width=660,height=330).configure_axis(labelLimit=1000)
            st.altair_chart(cpdat_htmp)


            ip_input = pd.DataFrame()
            #Exposure pathway probabilities
            for dtxsid, chem_name in zip (returned_table['DTXSID'], returned_table['PREFERRED_NAME']):

                input_pathways = expo.search_exposures(by='pathways', dtxsid=dtxsid)

                if input_pathways.empty:

                    ip_addition = pd.DataFrame({'Exposure_Type': ['Dietary Exposure Predicted:', 'Residential Exposure Predicted:', 'Pesticide Exposure Predicted:', 'Industrial Exposure Predicted:' ],
                                                'Exposure_Probabilities':[None, None, None, None]
                                                })
                    
                else:
                    diet_prob = input_pathways['probabilityDietary'].iloc[0]
                    res_prob = input_pathways['probabilityResidential'].iloc[0]
                    pest_prob = input_pathways['probabilityPesticde'].iloc[0]
                    industry_prob = input_pathways['probabilityIndustrial'].iloc[0]
                    
                    ip_addition = pd.DataFrame({'Exposure_Type': ['Dietary Exposure Predicted:', 'Residential Exposure Predicted:', 'Pesticide Exposure Predicted:', 'Industrial Exposure Predicted:' ],
                                                'Exposure_Probabilities':[diet_prob, res_prob, pest_prob, industry_prob]
                                                })
                    
                ip_addition['substance_name'] = chem_name
                
                ip_input = pd.concat([ip_input, ip_addition])

            ip_htmp = alt.Chart(ip_input).mark_rect().encode(
                x = alt.X('substance_name:N', 
                    axis=alt.Axis(title=''),
                    scale=alt.Scale(domain=list(ip_input['substance_name'].unique()))),
                y = alt.Y('Exposure_Type:N', title=''),
                    color = alt.Color('Exposure_Probabilities:Q').scale(scheme='turbo')
                    .title("Exposure Predicted")).properties(width=660,height=330).configure_axis(labelLimit=1000)
            
            st.markdown('#### Estimated probability of exposure through common pathways')
            st.altair_chart(ip_htmp)


            nhanes_path = (run_from_location/"data"/"NHANES_inferences.csv")
            nhanes_inferences = pd.read_csv(nhanes_path)


            
            third_fig_input = pd.DataFrame()
            #Model prediction data
            for dtxsid, chem_name in zip (returned_table['DTXSID'], returned_table['PREFERRED_NAME']):
                
                nhanes_presence = nhanes_inferences[nhanes_inferences['dsstox_substance_id'] == dtxsid]
                
                #Cannot make whole dataframe at once with all inputs, 
                # because each input must be filtered independently depending on whether it contains information  

                if not nhanes_presence.empty:
                        
                    nhanes_info = nhanes_presence['mgpkgpday'].iloc[0]
                    nhanes_df = pd.DataFrame({'database':["NHANES inferred exposures"], 
                                            'substance_name':[chem_name],
                                                "exposure_dose":[nhanes_info]})
                    
                else:
                    nhanes_df = pd.DataFrame({'database':["NHANES inferred exposures"],
                                            'substance_name':[chem_name],
                                            "exposure_dose":[None]})
                    
                third_fig_input= pd.concat([third_fig_input, nhanes_df])


                seem = pd.DataFrame(expo.search_exposures(by="seem", dtxsid=dtxsid))

                if not seem.empty:
                    consensus = seem[seem['predictor']=='SEEM3 Consensus']
                    if not consensus.empty:
                        pred_exp = consensus['median'].iloc[0]
                    else:
                        pred_exp = None

                    sheds_direct = seem[seem['predictor']=='SHEDS.Direct']
                    if not sheds_direct.empty:
                        pred_shed_dir = sheds_direct['median'].iloc[0]
                    else:
                        pred_shed_dir = None

                    sheds_indirect = seem[seem['predictor']=='SHEDS.Indirect']
                    if not sheds_indirect.empty:
                        pred_shed_indir = sheds_indirect['median'].iloc[0]
                    else:
                        pred_shed_indir = None

                    raidar_farfield = seem[seem['predictor']=='RAIDAR']
                    if not raidar_farfield.empty:
                        raidar_far_pred = raidar_farfield['median'].iloc[0]
                    else:
                        raidar_far_pred = None


                    raidar_ice = seem[seem['predictor']=='RAIDAR.ICE']
                    if not raidar_ice.empty:
                        raidar_ice_pred = raidar_ice['median'].iloc[0]
                    else:
                        raidar_ice_pred = None

                    seem_df = pd.DataFrame({'database':["SEEM3 consensus", "SHEDS-HT direct", "SHEDS-HT direct", "RAIDAR Far-Field", "RAIDAR Indoor and Consumer" ],                                        
                                            'substance_name':[chem_name, chem_name, chem_name, chem_name, chem_name],
                                            "exposure_dose":[pred_exp, pred_shed_dir, pred_shed_indir, raidar_far_pred, raidar_ice_pred]})
                else:

                    seem_df = pd.DataFrame({'database':["SEEM3 consensus", "SHEDS-HT direct", "SHEDS-HT direct", "RAIDAR Far-Field", "RAIDAR Indoor and Consumer" ],
                                            'substance_name':[chem_name, chem_name, chem_name, chem_name, chem_name],
                                            "exposure_dose":[None, None, None, None, None,]})
                    
                third_fig_input= pd.concat([third_fig_input, seem_df])

            modeled_htmp = alt.Chart(third_fig_input).mark_rect().encode(
            x = alt.X('substance_name:N', 
                axis=alt.Axis(title=''),
                scale=alt.Scale(domain=list(third_fig_input['substance_name'].unique()))),
            y = alt.Y('database:N', title=''),
                color = alt.Color('exposure_dose:Q').scale(scheme='turbo').title("Exposure Dose Predicted (mg/kg/day)")
                ).properties(width=660,height=330).configure_axis(labelLimit=1000)

            st.markdown('####  Predicted exposure doses in units of mg/kg/day')
            st.altair_chart(modeled_htmp)
            









            

                
                




            













