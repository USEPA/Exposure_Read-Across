import streamlit as st
from ctxpy import Exposure
import pandas as pd
import altair as alt

def cpdat_displays(structure_dtxsid):
    expo = Exposure(x_api_key='aaa69edc-d6d6-4d60-83d1-d9bd8e82f12f')
    with st.container(border=True):
                #Pulls reported functional use, predicted functional use, and presence in consumer/industrial formulations or articles 
                #Dataset displays are sequential due to being too wide to display as side-by-side columns 
                
                st.markdown('## Usage Information ')
                st.markdown('##### Mouse-over any entries whose labels are cut off to see the full labels')
                st.markdown('##### Data pulled from CPDAT by DTXSID ')
                #Cannot find a way to control the size of header text;
                #none of the options desrcribed in the documentation work
                st.header('', divider=True)

                #Pulling reported functional use
                function_info = expo.search_cpdat(vocab_name="fc", dtxsid=structure_dtxsid)
                function_info = pd.DataFrame(function_info)
                function_info_csv = function_info.to_csv(index=False)
                st.download_button(
                                    label="Download all substance functional use data as CSV",
                                    data=function_info_csv,
                                    file_name="functional_use_data.csv"
                                    )
                    
                
                if function_info.empty:
                    st.markdown('### No Reported Function information available')
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
                    st.markdown('### Functional Use Information')
                    st.altair_chart(docs_per_function_cat, use_container_width = True)
                    st.markdown("##### Figure note: Collected documents containing entered substance grouped by reported function within a function category. " 
                                "Only data with a reported function has been included. "
                                "This data facilitates the understanding of the uses of the entered substance by displaying the substance's reported functions."
                                )  
                
                #Now pulling list presence information 
                st.markdown('### List Presence Information')
                st.markdown('##### The lists on which the entered substance is present. ' 
                            'These provide additional information about the contexts in which the substance is used.')
                list_presence = expo.search_cpdat(vocab_name='lpk', dtxsid=structure_dtxsid)
                list_presence = pd.DataFrame(list_presence)
                if list_presence.empty:
                    st.markdown('### No List Presence information available')
                else:
                    list_presence.drop(columns=['id','dtxsid','docid',], inplace=True)
                    list_presence.rename(columns={'doctitle':'Document Title', 
                                                'docdate':'Document Date',
                                                'reportedfunction':'Reported Function',
                                                'functioncategory':'Function Category',
                                                'keywordset':'Key Word Set',
                                                'docsubtitle':'Document Subtitle',
                                                'organization':'Organization',
                                                }, inplace=True)
                    st.dataframe(list_presence) 
                
                #Pulling product-inclusion data 
                st.markdown('### Products Reporting This Substance')
                pucs = expo.search_cpdat(vocab_name='puc', dtxsid=structure_dtxsid)
                pucs_available = pd.DataFrame(pucs)
                if pucs_available.empty:
                    st.markdown('### No PUC Information Available')
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
                        st.download_button(
                                    label="Download product use data as CSV",
                                    data=pucs_available_csv,
                                    file_name="functional_use_data.csv"
                                    )
                        pucs_available.drop(columns=['id','dtxsid','docid',], inplace=True)
                        puc_show_opts = st.selectbox("Which level of PUC results would you like to see?", ("General Category", 'PUC general category + PUC product family'))

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
                        


                        show_all_check = st.checkbox("Show all PUC categories")
                        st.write("PUC categories represented by only a small number of products have been hidden. Click the checkbox above to show them") 
                        if not show_all_check:
                            grouped_pucs.sort_values(by='Number of product names', ascending=False, ignore_index=True, inplace=True)
                            topval = grouped_pucs['Number of product names'][0]
                            grouped_pucs = grouped_pucs[grouped_pucs['Number of product names']>(topval*.025)]                            

                        prods_per_catfam = alt.Chart(grouped_pucs).mark_bar().encode(
                            x=alt.X('Number of product names:Q'),
                            y=alt.Y(puc_show_opts + ':N',  
                                    sort='-x',  
                                    axis=alt.Axis(title=puc_show_opts,
                                    titleX=-350),
                            )).configure_axis(labelLimit=1000) 
                        st.altair_chart(prods_per_catfam, use_container_width = True)
                       # st.markdown('##### Figure note: Number of products associated with each combination of PUC general category and product family')
                    except:
                        st.error("Error encountered when retrieving product use information")