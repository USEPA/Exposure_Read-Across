import pytest
import sys
import os
from pathlib import Path
import pandas as pd

#print(os.getcwd(  ))
#Comment
#current_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(r'C:\Users\svanambe\OneDrive - Environmental Protection Agency (EPA)\Profile\Work\Projects\Exposure web app\Streamlit app\Official_version')
from webpage_script import analog_finder

chem_info_path = Path (r"C:\Users\svanambe\OneDrive - Environmental Protection Agency (EPA)\Profile\Work\Projects\Exposure web app\Streamlit app\Unit tests\chem_info.parq") 
morgan_fp_path = Path (r"C:\Users\svanambe\OneDrive - Environmental Protection Agency (EPA)\Profile\Work\Projects\Exposure web app\Streamlit app\Unit tests\shah-morgan-fp.parq")

one_chloro_four_nitrobenzene_path = Path (r"C:\Users\svanambe\OneDrive - Environmental Protection Agency (EPA)\Profile\Work\Projects\Exposure web app\Streamlit app\Unit tests\one-Chloro-4-nitrobenzene-NN.parq") 
nitrofen_path = Path (r"C:\Users\svanambe\OneDrive - Environmental Protection Agency (EPA)\Profile\Work\Projects\Exposure web app\Streamlit app\Unit tests\nitrofen-NN.parq")
four_nitroaniline_path = Path (r"C:\Users\svanambe\OneDrive - Environmental Protection Agency (EPA)\Profile\Work\Projects\Exposure web app\Streamlit app\Unit tests\4-Nitroaniline-NN.parq")

one_chloro_four_nitrobenzene_df = pd.read_parquet(one_chloro_four_nitrobenzene_path)
nitrofen_df = pd.read_parquet(nitrofen_path)
four_nitroaniline_df = pd.read_parquet(four_nitroaniline_path)

#chem_info = pd.read_parquet(chem_info_path)
#morgan_fp = pd.read_parquet(morgan_fp_path)

one_chloro_four_nitrobenzene_SMILES =  "C1=CC(=CC=C1[N+](=O)[O-])Cl"
nitrofen_SMILES = "[O-][N+](=O)C1=CC=C(OC2=C(Cl)C=C(Cl)C=C2)C=C1"
four_nitroaniline_SMILES = "NC1=CC=C(C=C1)[N+]([O-])=O"



def test_analogues():
    #Function returns table of analogs and their info 
    #assert analog_finder(morgan_fp_path, chem_info_path, one_chloro_four_nitrobenzene_SMILES, 1024)["dtxsid"].equals(one_chloro_four_nitrobenzene_df["dtxsid"][0:10])
    #print(analog_finder(morgan_fp_path, chem_info_path, one_chloro_four_nitrobenzene_SMILES, 1024))
    print(analog_finder(morgan_fp_path, chem_info_path, nitrofen_SMILES, 1024)[['sim', 'chemical_name']])
    print(nitrofen_df [['sim', 'chemical_name']])

test_analogues()


