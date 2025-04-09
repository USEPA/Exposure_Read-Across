# Exposure Data Viewer
### Created by: Samuel Van Amberg

A web app for searching and viewing exposure relevant data for chemicals.

## Installation Steps

1. Obtain a personal API Key for the CTX APIs. You can request a key by emailing [ccte_api@epa.gov](mailto:ccte_api@epa.gov).
2. Clone the repository. This can be done by navigating on the command line to your desired location for the installation of this repo and then typing `git clone https://github.com/USEPA/Exposure_Read-Across.git` 
3. Navigate to the `Exposure_Read_Across` directory that you created by cloning above and type `pip install .`
4. After installation of dependencies (namely `ctx-python`) type `ctx_init --x-api-key YourPersonalAPIKey` into the command line. `YourPersonalAPIKey` is the key you obtain from installation step 1.

## Running the app
 
Navigate to the cloned `Exposure_Read-Across` folder in your preferred command-line interface and then run the line `streamlit run Home_Page.py`. This will open up the app in your browser with all three pages available. If the app does not open automatically, click the link provided by the command-line interface to launch the application in a new tab/window in your default browser.

## Usage
The application currently has two use-cases:
1. To look-up exposure-relevant information for a specified chemical (done via the `Chemical_Exposure_Information_Search.py` page)
2. To look-up and summarize exposure-relevant information for chemical analogues of a specified chemical (done via the `Chemical_Exposure_Read_Across.py` page)

## Credit

The file structure for this repo is based on the [andymcdgeo/cookiecutter_streamlit_app](https://github.com/andymcdgeo/cookiecutter-streamlit) Cookiecutter project template.

## Disclaimer
This software/application was developed by the U.S. Environmental Protection Agency (USEPA). No warranty expressed or implied is made regarding the accuracy or utility of the system, nor shall the act of distribution constitute any such warranty. The USEPA has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the USEPA. The USEPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by the USEPA or the United States Government.
