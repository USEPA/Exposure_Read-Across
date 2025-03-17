# Effort has been made to conform code to PEP8
# conventions, including the 79-character line limit.
# However, in rare cases, it was not possible to
# break lines prior to their length exceeding 79 characters
# while maintaining readability
import streamlit as st
from ctxpy import Exposure
import pandas as pd
import altair as alt


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
        chemexpo = "https://comptox.epa.gov/chemexpo"
        self.reported_container.subheader('Reported Information')
        self.reported_container.markdown(
                   'Reported consumer use information is'
                   ' obtained from the '
                   'EPA Chemicals and Products Database (CPDat) v4.0.'
                   'using CTX. To interactively '
                   'explore this database, visit the '
                   f'[Chemical Exposure Knowledgebase (ChemExpo)]({chemexpo}).')
        self.reported_container.info('Hover over any datapoint or bar on a '
                                     'figure to see its full set of labels')
            
        self.reported_container.divider()

        # Pulling reported functional use
        function_info = expo.search_cpdat(vocab_name="fc", 
                                          dtxsid=self.structure_dtxsid)
        function_info = pd.DataFrame(function_info)
           
        if function_info.empty:
            self.reported_container.error('No Reported Function '
                                          'information available')
        else:
            # Dropping id, dtxsid, docid columns 
            function_info.drop(columns=['id', 'dtxsid', 'docid'], inplace=True)
            
            #Dropping rows with no function information 
            function_info=function_info[function_info['reportedfunction'].notna()]

            function_info.rename(columns={'datatype':'Data Type',
                                          'doctitle':'Document Title',
                                          'docdate':'Document Date',
                                          'reportedfunction':'Reported Function',
                                          'functioncategory':'Function Category'},
                                 inplace=True)
            
            chart_function_info = (function_info
                                   .drop_duplicates(subset=['Document Title', 
                                                            'Reported Function',
                                                            'Function Category'],
                                                    ignore_index=True))
            
            grouped_cfi = (chart_function_info
                         .groupby(by=['Reported Function',
                                      'Function Category'])
                         .size()
                         .reset_index(name='Docs per function'))
            # Experimented with making this plot vertical,
            # but decided that it looked
            # worse, due to the bars being much shorter
            # and running into the legend
            scale = alt.Scale(scheme='category20')
            docs_per_function_cat = (alt
                                     .Chart(grouped_cfi)
                                     .mark_bar()
                                     .encode(x="Docs per function",
                                             y='Function Category',
                                             color=(alt
                                                    .Color('Reported Function',
                                                           scale=scale)))
                                     .configure_axis(labelLimit=1000))
            chemexpo = "https://comptox.epa.gov/chemexpo/functional_use_categories/"
            self.reported_container.markdown('### Functional Use Information')
            (self
             .reported_container
             .markdown("Functional uses facilitates the "
                       "understanding of the different"
                       " ways the substance may be used "
                       "both in commerce as well as in "
                       "consumer products and industrial processes. For more "
                       "information about Function Categories "
                       "and how they are used in "
                       f"CPDat, visit [the Function Categories page]({chemexpo}) on "
                       "ChemExpo. Data used to generated "
                       "the figure below can be "
                       "obtained as a CSV file by "
                       "clicking the `Download Data` button "
                       "below the figure caption."))
            
            self.reported_container.altair_chart(docs_per_function_cat,
                                                 use_container_width=True)
            
            function_info_csv = function_info.to_csv(index=False)
            
            (self
             .reported_container
             .caption(f"The number of documents containing {self.structure_dtxsid} in "
                      "CPDat. Documents are grouped by the unique Function Categories "
                      "(FCs) the chemical is reported to serve. Within each FC, the "
                      "count of documents having "
                      "different reported functions are "
                      "shown by the differnt colored "
                      "sections of the bar. Documents "
                      "that had no functional use information are excluded."))

            self.reported_container.download_button(
                            label="Download Data",
                            data=function_info_csv,
                            file_name=("reported_functional_use_"
                                       f"{self.structure_dtxsid}.csv")
                            )

    # Now pulling list presence information
    def list_presence(self):
        expo = Exposure()
        chemexpo = "https://comptox.epa.gov/chemexpo/list_presence_tags/"
        self.reported_container.markdown('### List Presence Information')
        (self
         .reported_container
         .markdown("List Presence Keywords inform "
                   "the broad use of chemicals across "
                   "different sectors, products,"
                   " geographic location, and even use in "
                   "regulatory or government documetns. "
                   "For more information about "
                   "List Presence Keywords and how "
                   "they are used in CPDat, visit the "
                   f"[List Presence Keywords page]({chemexpo}) on ChemExpo."))
        list_presence = expo.search_cpdat(vocab_name='lpk',
                                          dtxsid=self.structure_dtxsid)
        list_presence = pd.DataFrame(list_presence)
        if list_presence.empty:
            self.reported_container.error('No List Presence '
                                          'information available')
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
            (self
             .reported_container
             .caption("The List Presence Keywords "
                      "associated with documents that "
                      f"contain the chemical {self.structure_dtxsid} " 
                      "on which the "
                      "entered substance is present."))
       
    # Pulling product-inclusion data 
    def product_inclusion(self):
        expo = Exposure()
    
        self.reported_container.markdown('### Products Reporting '
                                         'This Substance')
        puc_url = 'https://comptox.epa.gov/chemexpo/pucs/'
        (self
         .reported_container
         .markdown("Product use categories provide "
                   "information on how a consumer "
                   "product may be used. These categories "
                   "have specifically been "
                   "developed for use in exposure modeling "
                   "and assessment. For more "
                   "information about Product Use Categories and how they are "
                   f"used in CPDat, visit the [Product Use Categories page]({puc_url}) "
                   "on ChemExpo. Using the drop-down "
                   "menu below, you can select the "
                   "PUC-level to use for aggregation. "
                   "Checking the `Show all PUCs` box "
                   "below will limit the number of "
                   "categories for display (if there is "
                   "a large number). Data used to "
                   "generated this figure can be "
                   "obtained from the `Download Data` button below the figure "
                   "caption."))
        pucs = expo.search_cpdat(vocab_name='puc', 
                                 dtxsid=self.structure_dtxsid)
        pucs_available = pd.DataFrame(pucs)
        if pucs_available.empty:
            self.reported_container.error('No PUC Information Available')
        else:
            try:
                (pucs_available
                 .rename(columns={'doctitle':'Document Title', 
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
                                  'component':'Component',},
                         inplace=True))
                
                pucs_available_csv = pucs_available.to_csv(index=False)
                
                pucs_available.drop(columns=['id', 'dtxsid', 'docid',], 
                                    inplace=True)
                puc_show_opts = (self
                                 .reported_container
                                 .selectbox(("Which level of PUC results "
                                            "would you like to see?"),
                                            ("General Category",
                                             ('PUC general category + '
                                              'PUC product family'))))

                if puc_show_opts == 'PUC general category + PUC product family':
                    # Duplicates are dropped because 
                    # the metric of interest is the 
                    # number of distinct products, not the number of records 
                    pucs = ['Product Name',
                            'General Category',
                            'Product Family']
                    no_dups_pucs_available = (pucs_available
                                              .drop_duplicates(subset=pucs,
                                                               ignore_index=True))
                    grouped_pucs = (no_dups_pucs_available
                                    .groupby(by=['General Category', 
                                                 'Product Family'])
                                    .size()
                                    .reset_index(name='Number of product names'))
                    #String concatenation to get the combined names 
                    col = 'PUC general category + PUC product family'
                    grouped_pucs[col] = (grouped_pucs['General Category'] +
                                         " + " +
                                         grouped_pucs['Product Family'])
                    #Dropping empty cells
                    grouped_pucs = grouped_pucs[grouped_pucs['Product Family'] != '']
                    grouped_pucs = grouped_pucs[grouped_pucs['General Category'] != '']     
                
                if puc_show_opts == "General Category":
                    # Duplicates are dropped because
                    # the metric of interest is the
                    # number of distinct products, not the number of records 
                    pucs = ['Product Name', 'General Category',]
                    no_dups_pucs_available = (pucs_available
                                              .drop_duplicates(subset=pucs,
                                                               ignore_index=True))
                    grouped_pucs = (no_dups_pucs_available
                                    .groupby(by=['General Category'])
                                    .size()
                                    .reset_index(name='Number of product names'))
                    
                    #Dropping empty cells 
                    grouped_pucs = grouped_pucs[grouped_pucs['General Category'] != '']

                show_all_check = self.reported_container.checkbox("Show all PUCs")
                (self
                 .reported_container
                 .markdown("PUC categories represented by "
                           "only a small number of "
                           "products have been hidden. "
                           "Check the `Show all PUCs` box"
                           "above to show them."))
                if not show_all_check:
                    (grouped_pucs
                     .sort_values(by='Number of product names', 
                                  ascending=False, 
                                  ignore_index=True, 
                                  inplace=True))
                    topval = grouped_pucs['Number of product names'][0]
                    idx = grouped_pucs['Number of product names'] > (topval*.025)
                    grouped_pucs = grouped_pucs[idx]

                prods_per_catfam = (alt
                                    .Chart(grouped_pucs)
                                    .mark_bar()
                                    .encode(x=alt.X('Number of product names:Q'),
                                            y=alt.Y(puc_show_opts + ':N',
                                                    sort='-x',
                                                    axis=alt.Axis(title=puc_show_opts,
                                                                  titleX=-350)))
                                    .configure_axis(labelLimit=1000))
                self.reported_container.altair_chart(prods_per_catfam,
                                                     use_container_width=True)
                (self
                 .reported_container
                 .caption("Number of unique CPDat products that contain "
                          f"{self.structure_dtxsid} and have "
                          "been assigned a Product Use Category (PUC)."))
                (self
                 .reported_container
                 .download_button(label="Download Data",
                                  data=pucs_available_csv,
                                  file_name=("product_use_data_"
                                             f"{self.structure_dtxsid}.csv")))      
                
            # Bare "except" used due to inform user
            # that error was encountered without program halting
            except:
                (self
                 .reported_container
                 .error("Error encountered when "
                        "retrieving product use information"))

    # Underscores used in method argument names to
    # tell streamlit not to hash the arguments
    # (as recommended by error message)
    @st.cache_data
    def import_osha(_self):
        
        """
        Reads in the OSHA CEHD dataset produced by Jeff Minucci's model 

        Returns:
            Dataframe: The DataFrame of osha data
        """
        osha = pd.read_feather(_self.script_location / 
                               "data" / 
                               "osha_processed_full.feather")
        osha.rename(columns={'conc_mgm3': 'exposure_level',
                             'sample_type': 'measure_unit_id',
                             'subsector_name': 'naics_2022_subsector_title',
                             'sector_name': 'naics_2022_sector_title'
                             }, inplace=True)
        return osha

    def osha_info(self):
        (self
         .reported_container
         .markdown('### OSHA Occupational Monitoring Information'))
        (self
         .reported_container
         .markdown("OSHA's Chemical Exposure Health Data "
                   "(CEHD) provides industrial-"
                   "hygiene samples of chemicals "
                   "measured in a workplace by OSHA "
                   "compliance officers from 1984 to "
                   "2022. While CEHD provides dermal "
                   "wipes and air measurements, "
                   "only the air concentrations are "
                   "provided by the Exposure Data "
                   "Viewer. These air concentration "
                   "measurements are taken by sampling "
                   "either a worker's personal "
                   "breathing space or the shared working area."))

        osha = self.import_osha()
        
        self.reported_container.container
        osha_of_interest = osha[osha['dtxsid'] == self.structure_dtxsid].copy()

        if osha_of_interest.empty:
            self.reported_container.error('No exposure data available')
        else:
            #Include a note that all zero-values have been dropped
            #osha_of_interest.drop('index', axis=1, inplace=True)
            
            csv_osha_of_interest = osha_of_interest.to_csv(index=False)
            self.reported_container.markdown('### OSHA Workplace '
                                             'Air Monitoring Data')
            (self
             .reported_container
             .download_button(label="Download Data",
                              data=csv_osha_of_interest,
                              file_name=("osha_cehd_data_"
                                         f"{self.structure_dtxsid}.csv")))
            # This dataframe will be displayed as a box and whisker plot
            # Dropping entries with a value of zero, 
            # as these do not contain useful 
            # information
            idx = (osha_of_interest['exposure_level'] == '0' |
                   osha_of_interest['exposure_level'] == '0.0')
            osha_of_interest = osha_of_interest[~idx]

            #Dropping "Fibers per cubic centimeter" and empty rows.
            #Is not necessary to have a separate clause for null entries, 
            # as those are classified as not being in the list.
            idx = osha_of_interest['measure_unit_id'].isin(['M', 'P'])
            osha_of_interest = osha_of_interest[idx]
            
            # Conversion of units to mg/m^3 is not necessary,
            # as it has already been done  
            
            #Change the labeled unit type when the conversion is done
            osha_of_interest['measure_unit_id'] = 'M'
            osha_of_interest['measure_unit_name'] = 'Milligrams per cubic meter'
        
            #The units that are convertable to mg/m^3:
                #P - mg/m3, #X- micrograms, Y - milligrams,
            #Not convertable:
            # F - fibers/cc, % - percentage of bulk material    

            #Dropping any rows where the ppm_to_mg function returns a value of zero, 
            # as these cannot be placed on a log scale 
            
            osha_of_interest = osha_of_interest[osha_of_interest['exposure_level'] != 0]

            naics_sector_or_sub = (self
                                   .reported_container
                                   .selectbox(("Which level of NAICS results "
                                               "would you "
                                               "like to see?"),
                                              ("NAICS Sector",
                                               "NAICS Subsector")))
            if not osha_of_interest.empty:
                (self
                 .reported_container
                 .markdown("#### Occupational inhalation "
                           "exposure data on the entered "
                           "substance from filtered version "
                           "of the dataset CEHD, "
                           "converted to units of mg/m^3"))
                if naics_sector_or_sub == "NAICS Subsector":
                    (self
                     .reported_container
                     .markdown("Data is separated by NAICS subsector. "
                               "Any entries with units that "
                               "cannot be converted to "
                               "$mg/m^3$ or which have a "
                               "non-inhalation sample type "
                               "are not displayed."))
                    expo_range = [(osha_of_interest['exposure_level'].min())/10,
                                  (osha_of_interest['exposure_level'].max())*10]
                    box_plot = (alt
                                .Chart(osha_of_interest)
                                .mark_boxplot()
                                .encode(x=(alt.X('exposure_level:Q',
                                                 scale=alt.Scale(type="log",
                                                                 domain=expo_range))
                                            .title('Air concentration (mg/m^3)')),
                                        y= alt.Y('naics_2022_subsector_title:N',
                                                 sort='-x',
                                                 axis=alt.Axis(title='NAICS Subsector',
                                                               titleX=-370)))
                                .configure_axis(labelLimit=1000))
                    
                if naics_sector_or_sub == "NAICS Sector":
                    idx = osha_of_interest['naics_2022_sector_title'].notna()
                    osha_of_interest_nona = osha_of_interest[idx]
                    (self
                     .reported_container
                     .markdown("Data is separated by NAICS sector. "
                               "Any entries with "
                               "units that cannot be converted to "
                               "mg/m^3 or which have "
                               "a non-inhalation sample type "
                               "are not displayed."))
                    expo_range = [(osha_of_interest_nona['exposure_level'].min())/10,
                                  (osha_of_interest_nona['exposure_level'].max())*10]
                    box_plot = (alt.Chart(osha_of_interest_nona)
                                .mark_boxplot()
                                .encode(x=(alt.X('exposure_level:Q',
                                                  scale=alt.Scale(type="log",
                                                                  domain=expo_range))
                                            .title('Air concentration (mg/m^3)')),
                                        y=alt.Y('naics_2022_sector_title:N',
                                                 sort='-x',
                                                 axis=alt.Axis(title='NAICS sector',
                                                               titleX=-370)))
                                .configure_axis(labelLimit=1000))

                self.reported_container.altair_chart(box_plot,
                                                     use_container_width=True)
            else:   
                (self
                 .reported_container
                 .error('No occupational inhalation data available'))

        return osha
