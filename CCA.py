import pandas as pd
import os

def get_std_table(layer, site):
    """This function creates an OTU coverage dataframe and it's corresponding ordered normalized standard deviations. The OTU coverage table is saved as a csv file.
    Input:
    - layer: The kingdom specified, could either be Bacteria, Eukaryota, or Archaea. String format.
    - site: The station of interest. Could either be S, P, or CD. String format
    Output:
    - sub_result: The coverage file of the station given the kingdom
    - std_df: The normalized standard deviations dataframe. 
    """
    #Initiate an empty dataframe
    sub_result = pd.DataFrame()
    #Iterate over the files within the given directory
    for file in os.listdir('16s_18s/files'):
        #Check if the file ends with csv
        if file.endswith('.csv'):
            #And starts with Phyla_of_.. The specified layer and the given site
            if file.startswith('Phyla_of_' + layer + '_' + site):
                #Print the file name
                print(file)
                #Read the file
                data = pd.read_csv('16s_18s/files/' + file)
                #Change the column names
                data.columns = ['Phyla', layer + file[len(layer) + 9: -4]]
                #Set the phyla as index
                data = data.set_index('Phyla')
                #Add the dataframe to the previous initiated dataframe
                sub_result = pd.concat([sub_result, data], axis=1)
    #Fill the missing values with 0
    sub_result.fillna(0, inplace = True)
    #Print the top rows
    print(sub_result.head(5))


    #filter the column names containing the specified layer within the lineage and the site name
    filter = [col for col in sub_result if col.startswith(layer + '_' + site)]

    #Calculate the mean and std of the counts
    std_df = sub_result[filter].std(axis = 1)
    mean_df = sub_result[filter].mean(axis = 1)
    #Normalize the standard deviations by dividing by the mean
    std_df = std_df.divide(mean_df)
    #Sort the standard deviation values
    std_df.sort_values(inplace = True, ascending = False)
    #Filter the subresults dataframe
    sub_result = sub_result[filter]
    #Save the generated results to csv file
    sub_result.to_csv('16s_18s/files/OTU_table_%s_%s.csv' % (layer, site))
    #Return the dataframe containin the counts and the ordered standard deviation
    return sub_result, std_df
