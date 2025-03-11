import streamlit as st
st.title("Exposure Read-Across Application")

st.header("Overview", divider=True)
st.markdown('This application retrieves, visualizes, and summarizes '
            'exposure data from publicly-available EPA and non-EPA sources.')

prereqs = """
This application accesses EPA data via EPA's 
[Computational Toxicology and Exposure (CTX) application program interface (API)](https://api-ccte.epa.gov/docs/index.html). This is 
done via the [`ctx-python` client](https://pypi.org/project/ctx-python/). In order to 
successfully use these APIs, an individual-use API key must be accquired by emailing 
[ccte_api@epa.gov](mailto:ccte_api@epa.gov) and 
requesting an API key. `ctx-python` documentation provides a tool that to save an 
individual's API key to a .env file for access by `ctx-python` and this application.
"""
st.header("Prerequisites", divider=True)
st.markdown(prereqs)
st.header("Navigation", divider=True)
st.markdown('Choose the "Analog Information Not Included" '
            'page in the right-hand bar to see reported'
            ' and predicted information on your chemical of interest.')

st.markdown('Choose the "Analog Information Included" page '
            'in the left-hand bar to see '
            'the information provided on the "Analog Information Not Included"'
            ' page, as well'
            ' as a summary of the information available on analogs '
            'of your chemical of interest.')
st.header("Usage Instructions", divider=True)
st.markdown('On either page, begin by choosing from the drop-down menu'
            ' whether you would prefer to'
            ' identify your chemical of interest by drawing its structure or'
            ' by entering a text ID (CAS-RN or DTXSID).'
            ' This choice will cause a input field to pop up.'
            ' Once you have submitted your chemical of interest,'
            ' the application will retrieve all available information on it.')

st.markdown('Annotations above each data-source '
            'suggest potential uses and interpretive approaches'
            ' for the information displayed.')


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


