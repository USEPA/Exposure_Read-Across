# Exposure_Read_Across
## Created by: Samuel Van Amberg

A web app for presenting data on analogs of new chemicals to OPPT New Chemical evaluators

## Installation

Clone the repository with 

```
git clone https://github.com/USEPA/Exposure_Read-Across.git

```

Assuming you have pip already installed, install the required packages with `pip install -r requirements.txt`

## Running the app

Navigate to the cloned `Exposure_Read-Across` folder in your preferred command-line interface and then run the line `streamlit run webpage_welcome.py`. This will open up the app in your browser with all three pages available. 

## Navigating the app

Choose the "Existing Chemical" page in the right-hand bar to see measured and predicted information on your chemical of interest. Choose the "New Chemical" page in the right-hand bar to see the information provided on the "Existing Chemical" page, as well as a summary of the information available on analogs of your chemical of interest. On either page, draw the structure of your chemical of interest and click "apply" to retrieve information on the chemical.
The pages contain annotations by each data-source that suggest potential uses for the information displayed.


## Credits

The file structure for this repo derives from the [andymcdgeo/cookiecutter_streamlit_app](https://github.com/andymcdgeo/cookiecutter-streamlit) Cookiecutter project template.
