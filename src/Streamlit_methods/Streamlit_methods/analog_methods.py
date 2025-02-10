import streamlit as st
from rdkit import Chem
from sklearn.neighbors import KNeighborsRegressor
import pandas as pd
from rdkit.Chem import rdFingerprintGenerator
import numpy as np
from ctxpy import Chemical

class analog_operations:
    #Converts measurements in mg/m^3 into units of ppm,
    # assuming pressure of 1 atm and temperature of 25 C=77F
    def ppm_to_mg (self, value, unit, mw)->float:

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
        st.header("Structural Analogs")
        #Commenting out while I work on the above 
        
        #First thing is to assign the Morgan fingerprints to the smiles code 
        FP_input_path = (script_location/"data"/"Morgan_fingerprints_of_DSSTox.feather")
        cd_input_path = (script_location/"data"/"brotli_Consolidated_DSSTox_QSAR_smiles_only.parq")
        final_table=self.analog_finder(FP_input_path, cd_input_path, chem_smiles=smiles_code)
        final_table.drop(columns=['index'], inplace=True )
        final_table.rename(columns={'ID':'DTXSID', 'sim':'Similarity'}, inplace=True)       
        with st.container(border=True):
            st.markdown('##### Analogs are pulled from the set of substances in DSSTox that possess a QSAR-ready SMILES')
            st.dataframe(final_table)
        return final_table
    

    def sum_fig(self, returned_table, usis_data):
        chem = Chemical(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')
        for dtxsid, chem_name in zip (returned_table['DTXSID'], returned_table['PREFERRED_NAME']):
            expo_data = usis_data[usis_data['dtxsid']==dtxsid]
            # Cleaning USIS for use in the analog-summarizing program
            usis.drop('index', axis=1, inplace=True)
            usis=usis[usis['exposure_level'] != '0']
            usis=usis[usis['exposure_level'] != '0.0']
            usis = usis[usis['measure_unit_id'].isin(['M', 'P'])]
            deets = chem.details(by='dtxsid', word=dtxsid)   
            molar_mass = deets['monoisotopicMass']    
            usis['exposure_level'] = usis.apply(lambda x: self.ppm_to_mg (x.exposure_level, x.measure_unit_id, molar_mass), axis=1) 

            grouped_bysubs = expo_data.groupby('naics_2022_subsector_title')['exposure_level'].median()
            st.write("HERE'S SOME INFO")
            st.dataframe()










