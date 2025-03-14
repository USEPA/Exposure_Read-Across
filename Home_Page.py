import streamlit as st

st.set_page_config(layout="wide")
st.title("Exposure Data Viewer")

st.header("Overview", divider=True)
st.markdown("This application retrieves, visualizes, and summarizes "
            "exposure data from publicly-available sources using databases, models, "
            "and tools developed by EPA's Center for Computational Toxicology and "
            "Exposure.")
ctx_url = "https://api-ccte.epa.gov/docs/index.html"
ctxpy_url = "https://pypi.org/project/ctx-python/"
prereqs = (
    "This application accesses EPA data via EPA's "
    "[Computational Toxicology and Exposure (CTX) application program interface (API)]"
    "({ctx_url}). This is done via the [`ctx-python` client]({ctxpy_url}). In order to "
    "successfully use these APIs, an individual-use API key must be accquired by "
    "emailing [ccte_api@epa.gov](mailto:ccte_api@epa.gov) and requesting an API key. "
    "`ctx-python` documentation provides a tool that to save an individual's API key "
    "to a .env file for access by `ctx-python` and this application."
    )

st.header("Prerequisites", divider=True)
st.markdown(prereqs)
st.header("Navigation", divider=True)
st.markdown("The Exposure Data Viewer has two different pages that provide different "
            "functionalities. These pages are the **Chemical Exposure Information "
            "Search** page and the **Chemical Exposure Read-Across** page.")
st.subheader('Chemical Exposure Information Search')
st.markdown('This page allows users to search for reported and predicted information '
            'available for a searched chemical')

st.subheader('Chemical Exposure Read-Across')
st.markdown('This page allows users to search for chemical analogs for a specified '
            'chemical, once the analogs are identified, then relevant exposure '
            'information for those analogs are provided.')

st.header("Usage", divider=True)
st.markdown('On either page, begin by the search method from the drop-down menu.')
st.markdown('For users who prefer to draw the structural representation of a chemical '
            'choose the "Structure" option. For users who prefer to enter a chemical '
            'identifier, choose the "DTXSID or CAS-RN" option.')
st.markdown("Once the structure is drawn or identifier is provided, follow the prompts "
            "to submit the chemical for a search. Note that the search window will "
            "close automatically when the search is complete.")

st.markdown('Results from different data sources and models will appear after the '
            'search window closes. Descriptions, visualizations and/or data '
            'resulting from different data sources and models will appear and '
            'interpretive approaches for the information will be displayed.')


disclaimer = """
Installation, use and modification to this application are be made at the user's own 
risk. Neither the U.S. EPA nor the program author(s) can assume responsibility for 
program modification, content, output, interpretation, or usage.

This application has been extensively tested and verified. However, as for all complex 
software, these programs may not be completely free of errors and may not be applicable 
for all cases. In no event will the U.S. EPA be liable for direct, indirect, special, 
incidental, or consequential damages arising out of the use of the programs and/or 
associated documentation.

This application has been reviewed in accordance with U.S. Environmental Protection 
Agency, Office of Research and Development, and approved for dissemination. Mention of 
trade names or commercial products does not constitute endorsement or recommendation 
for use.
"""
st.header("Disclaimer", divider=True)
st.markdown(disclaimer)


