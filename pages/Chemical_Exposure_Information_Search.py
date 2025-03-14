
# Although streamlit pages have the option of
# starting the page titles with a number
# to set their rank, doing this prevents importing functions into unit tests
import streamlit as st
import Streamlit_methods

st.title("Chemical Exposure Information Search")
st.divider()
# The shared_content module contains all code in common between the two pages of the 
# web app
entered, osha_full, function_location, ctx_smiles = Streamlit_methods.shared_content()
