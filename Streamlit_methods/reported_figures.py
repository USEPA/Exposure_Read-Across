import streamlit as st
from ctxpy import Exposure
import pandas as pd
import altair as alt
from ctxpy import Chemical


class ReportedInfo:
    def __init__(self, retrieved_dtxsid, script_location):
        self.structure_dtxsid = retrieved_dtxsid
        self.reported_container = st.container(border=True)
        self.script_location = script_location
   
    def functional_use(self):
        expo = Exposure()
        
        # Pulls reported functional use, predicted functional use,
        # and presence in consumer/industrial formulations or articles
        # Dataset displays are sequential due to being
        # too wide to display as side-by-side columns

        self.reported_container.markdown('## Reported Information ')
        self.reported_container.markdown('##### Mouse-over any entries whose labels are cut off to see the full labels')
        self.reported_container.markdown('##### Data pulled from CPDAT by DTXSID ')
        # Cannot find a way to control the size of header text;
        # none of the options desrcribed in the documentation work
        self.reported_container.divider()

        # Pulling reported functional use
        function_info = expo.search_cpdat(vocab_name="fc", 
                                            dtxsid=self.structure_dtxsid)
        function_info = pd.DataFrame(function_info)
        function_info_csv = function_info.to_csv(index=False)
        self.reported_container.download_button(
                            label="Download all substance functional use data as CSV",
                            data=function_info_csv,
                            file_name="functional_use_data.csv"
                            )
            
        if function_info.empty:
            self.reported_container.markdown('### No Reported Function information available')
        else:
            #Dropping id, dtxsid, docid columns 
            function_info.drop(columns=['id','dtxsid','docid',], inplace=True)
            
            
            #Dropping rows with no function information 
            function_info=function_info[function_info['reportedfunction'].notna()]
            

            function_info.rename(columns={'datatype':'Data Type', 
                                        'doctitle':'Document Title', 
                                        'docdate':'Document Date', 
                                        'reportedfunction':'Reported Function', 
                                        'functioncategory':'Function Category'}, inplace=True)
            
            chart_function_info = function_info.drop_duplicates(subset=['Document Title', 
                                                                        'Reported Function',
                                                                        'Function Category'
                                                                        ],
                                                                ignore_index=True  
                                                                )# inplace=True  
            
            grouped_cfi=chart_function_info.groupby(by=['Reported Function',
                                            'Function Category']
                                            ).size(
                                            ).reset_index(name='Docs per function')
            #Experimented with making this plot vertical, but decided that it looked worse,
            #due to the bars being much shorter and running into the legend. 
            # I don't know if there is a way to move different components of the plot around  
            docs_per_function_cat = alt.Chart(grouped_cfi).mark_bar().encode(
                x="Docs per function",
                y='Function Category',
                color=alt.Color('Reported Function', scale=alt.Scale(scheme='category20') )
            ).configure_axis(labelLimit=1000)
            self.reported_container.markdown('### Functional Use Information')
            self.reported_container.altair_chart(docs_per_function_cat, use_container_width = True)
            self.reported_container.markdown("##### Figure note: Collected documents containing entered substance grouped by reported function within a function category. " 
                        "Only data with a reported function has been included. "
                        "This data facilitates the understanding of the uses of the entered substance by displaying the substance's reported functions."
                                                )

    # Now pulling list presence information                 
    def list_presence(self):
        expo = Exposure()
    
        self.reported_container.markdown('### List Presence Information')
        self.reported_container.markdown('##### The lists on which the entered substance is present. ' 
                    'These provide additional information about the contexts in which the substance is used.')
        list_presence = expo.search_cpdat(vocab_name='lpk', dtxsid=self.structure_dtxsid)
        list_presence = pd.DataFrame(list_presence)
        if list_presence.empty:
            self.reported_container.markdown('### No List Presence information available')
        else:
            list_presence.drop(columns=['id', 'dtxsid', 'docid',], inplace=True)
            list_presence.rename(columns={'doctitle': 'Document Title', 
                                            'docdate': 'Document Date',
                                            'reportedfunction': 'Reported Function',
                                            'functioncategory': 'Function Category',
                                            'keywordset': 'Key Word Set',
                                            'docsubtitle': 'Document Subtitle',
                                            'organization': 'Organization',
                                        }, inplace=True)
            self.reported_container.dataframe(list_presence) 
        
    # Pulling product-inclusion data 
    def product_inclusion(self):
            expo = Exposure()
        
            self.reported_container.markdown('### Products Reporting This Substance')
            pucs = expo.search_cpdat(vocab_name='puc', dtxsid=self.structure_dtxsid)
            pucs_available = pd.DataFrame(pucs)
            if pucs_available.empty:
                self.reported_container.markdown('### No PUC Information Available')
            else:
                try:
                    pucs_available.rename(columns={'doctitle':'Document Title', 
                                                'docdate':'Document Date', 
                                                'productname':'Product Name', 
                                                'gencat':'General Category',
                                                'prodfam':'Product Family',
                                                'prodtype':'Product Type',
                                                'classificationmethod':'Classification Method',
                                                'rawmincomp':'Raw Minimum \n Composition',
                                                'rawmaxcomp':'Raw Maximum Composition',
                                                'rawcentralcomp':'Raw Central Composition',
                                                'unittype':'Unit Type',
                                                'lowerweightfraction':'Lower Weight Fraction',
                                                'upperweightfraction':'Upper Weight Fraction',
                                                'centralweightfraction':'Central Weight Fraction',
                                                'weightfractiontype':'Weight Fraction Type',
                                                'component':'Component',
                                                }, inplace=True)
                    
                    pucs_available_csv=pucs_available.to_csv(index=False)
                    self.reported_container.download_button(
                                label="Download product use data as CSV",
                                data=pucs_available_csv,
                                file_name="functional_use_data.csv"
                                )
                    pucs_available.drop(columns=['id','dtxsid','docid',], inplace=True)
                    puc_show_opts = self.reported_container.selectbox("Which level of PUC results would you like to see?", ("General Category", 'PUC general category + PUC product family'))

                    if puc_show_opts=='PUC general category + PUC product family':
                        #Duplicates are dropped because the metric of interest is the number of 
                        #distinct products, not the number of records 
                        no_dups_pucs_available=pucs_available.drop_duplicates(subset=['Product Name', 
                                                                                        'General Category',
                                                                                        'Product Family',
                                                                                        ], ignore_index=True)
                        grouped_pucs = no_dups_pucs_available.groupby(by=['General Category', 'Product Family']
                                                                        ).size(
                                                                        ).reset_index(name='Number of product names')
                        #String concatenation to get the combined names 
                        grouped_pucs['PUC general category + PUC product family'] = grouped_pucs['General Category'] + " + " + grouped_pucs['Product Family']
                        #Dropping empty cells
                        grouped_pucs=grouped_pucs[grouped_pucs['Product Family']!='']
                        grouped_pucs=grouped_pucs[grouped_pucs['General Category']!='']
                        
                    
                    if puc_show_opts=="General Category":
                        #Duplicates are dropped because the metric of interest is the number of 
                        #distinct products, not the number of records 
                        no_dups_pucs_available=pucs_available.drop_duplicates(subset=['Product Name', 
                                                                                        'General Category',
                                                                                        ], ignore_index=True)
                        grouped_pucs = no_dups_pucs_available.groupby(by=['General Category']
                                                                        ).size(
                                                                        ).reset_index(name='Number of product names')
                        
                        #Dropping empty cells 
                        grouped_pucs=grouped_pucs[grouped_pucs['General Category']!='']
                    


                    show_all_check = self.reported_container.checkbox("Show all PUC categories")
                    self.reported_container.write("PUC categories represented by only a small number of products" 
                             "have been hidden. Click the checkbox above to show them") 
                    if not show_all_check:
                        grouped_pucs.sort_values(by='Number of product names', 
                                                 ascending=False, 
                                                 ignore_index=True, 
                                                 inplace=True)
                        topval = grouped_pucs['Number of product names'][0]
                        grouped_pucs = grouped_pucs[grouped_pucs['Number of product names'] > (topval*.025)]                            

                    prods_per_catfam = alt.Chart(grouped_pucs).mark_bar().encode(
                        x=alt.X('Number of product names:Q'),
                        y=alt.Y(puc_show_opts + ':N',  
                                sort='-x',  
                                axis=alt.Axis(title=puc_show_opts,
                                titleX=-350),
                        )).configure_axis(labelLimit=1000) 
                    self.reported_container.altair_chart(prods_per_catfam, use_container_width = True)
                
                # Bare "except" used due to wanting to inform user 
                # that error was encountered without
                # program halting
                except:
                    self.reported_container.error("Error encountered when retrieving product use information")

    # Underscores used in method argument names to
    # tell streamlit not to hash the arguments
    # (as recommended by error message)
    @st.cache_data
    def import_osha(_self):
        
        """
        Reads in the osha dataset from Jacob Kvasnika

        Returns:
            Dataframe: The Dataframe of osha data
        """
        osha = pd.read_feather(_self.script_location/"data"/"osha_2023_with_dtxsid.feather")

        return osha

    # Converts measurements in ppm into units of mg/m^3,
    #  assuming pressure of 1 atm and temperature of 25 C=77F
    def ppm_to_mg(self, value, unit, mw) -> float:

        """
        Converts multiple starting units to mg/m^3

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
        if unit == 'P':
            #Conversion formula from https://www.cdc.gov/niosh/docs/2004-101/calc.html                       
            mg_per_m3_value  = (value*mw)/(24.45)
            return mg_per_m3_value  
        elif unit=='M':
            return value
        else:
            return '0'

    def osha_info(self):

        chem = Chemical()
        osha = self.import_osha()

        self.reported_container.container
        osha_of_interest = osha[osha['dtxsid'] == self.structure_dtxsid].copy()

        if osha_of_interest.empty:
            self.reported_container.markdown('### No exposure data available')
        else:
            #Include a note that all zero-values have been dropped
            osha_of_interest.drop('index', axis=1, inplace=True)
            
            csv_osha_of_interest=osha_of_interest.to_csv(index=False)
            self.reported_container.markdown('### OSHA data on the entered substance')
            self.reported_container.download_button(
                            label="Download measured exposure data as CSV",
                            data=csv_osha_of_interest,
                            file_name="osha_data_on_chemical_of_interest.csv"
                            )
            #This dataframe will be displayed as a box and whisker plot

            #Dropping entries with a value of zero, as these do not contain useful information
            osha_of_interest=osha_of_interest[osha_of_interest['exposure_level'] != '0']
            osha_of_interest=osha_of_interest[osha_of_interest['exposure_level'] != '0.0']
            #Dropping "Fibers per cubic centimeter" and empty rows.
            #Is not necessary to have a separate clause for null entries, 
            # as those are classified as not being in the list.
            osha_of_interest = osha_of_interest[osha_of_interest['measure_unit_id'].isin(['M', 'P'])]

            deets = chem.details(by='dtxsid', word=self.structure_dtxsid)   
            molar_mass = deets['averageMass']    
            osha_of_interest['exposure_level'] = osha_of_interest.apply(lambda x: self.ppm_to_mg(x.exposure_level, x.measure_unit_id, molar_mass), axis=1) 
            
            #Change the labeled unit type when the conversion is done
            osha_of_interest['measure_unit_id'] = 'M'
            osha_of_interest['measure_unit_name'] = 'Milligrams per cubic meter'
        
            #The units that are convertable to mg/m^3:
                #P - mg/m3, #X- micrograms, Y - milligrams,
            #Not convertable:
            # F - fibers/cc, % - percentage of bulk material    

            #Dropping any rows where the ppm_to_mg function returns a value of zero, 
            # as these cannot be placed on a log scale 
            osha_of_interest = osha_of_interest[osha_of_interest['exposure_level']!='0'] 

            naics_sector_or_sub = self.reported_container.selectbox("Which level of NAICS results would you like to see?", ("NAICS Sector", "NAICS Subsector"))
            if not osha_of_interest.empty:
                self.reported_container.markdown('#### Occupational inhalation exposure data on the entered substance from the combined datasets CEHD, OIS and IMIS, converted to units of mg/m^3')
                if naics_sector_or_sub == "NAICS Subsector":
                    self.reported_container.markdown('##### Data is separated by NAICS subsector. Any entries with units that cannot be converted to mg/m^3 or which have a non-inhalation sample type are not displayed.')
                    box_and_whisker_mg_m3 = alt.Chart(osha_of_interest).mark_boxplot().encode(
                        x= alt.X('exposure_level:Q', 
                                scale=alt.Scale(type="log", 
                                                domain=[(osha_of_interest['exposure_level'].min())/10,
                                                        (osha_of_interest['exposure_level'].max())*10 ])).title('Air concentration (mg/m^3)'),
                        y= alt.Y('naics_2022_subsector_title:N', 
                                sort='-x',
                                axis=alt.Axis(title='NAICS Subsector',
                                titleX=-370))).configure_axis(labelLimit=1000)                                                                 
                    
                if naics_sector_or_sub == "NAICS Sector":
                    osha_of_interest_nona = osha_of_interest[osha_of_interest['naics_2022_sector_title'].notna()]
                    self.reported_container.markdown('##### Data is separated by NAICS sector. Any entries with units that cannot be converted to mg/m^3 or which have a non-inhalation sample type are not displayed.')
                    box_and_whisker_mg_m3 = alt.Chart(osha_of_interest_nona).mark_boxplot().encode(
                        x= alt.X('exposure_level:Q', 
                                scale=alt.Scale(type="log", 
                                                domain=[(osha_of_interest_nona['exposure_level'].min())/10,
                                                        (osha_of_interest_nona['exposure_level'].max())*10 ])).title('Air concentration (mg/m^3)'),
                        y= alt.Y('naics_2022_sector_title:N', 
                                sort='-x',
                                axis=alt.Axis(title='NAICS sector',
                                titleX=-370))).configure_axis(labelLimit=1000)                                                                 
                    

                self.reported_container.altair_chart(box_and_whisker_mg_m3, use_container_width = True)
                # st.markdown('##### Figure note: Chart type is box-and-whisker, but for some substances only outliers can be seen')
            else:   
                self.reported_container.markdown('### No occupational inhalation data available')

        return osha