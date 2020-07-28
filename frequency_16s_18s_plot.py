import pandas as pd
import numpy as np
from collections import Counter
from prettytable import PrettyTable
import matplotlib.pyplot as plt
import argparse
import operator
from CCA import get_std_table
import os


def get_percentage(mseq, layer):
    """This function extracts the lineage of the given mseq file and creates a dictionary of the relative occurence of a specified level (Kingdom or Phyla) in percentages.
    Input:
    - mseq: The mseq file
    - layer: The level of the lineage of interest. Could either be Kingdom or Phyla. String format.
    Output:
    - result: The raw occurence of the specified level
    - norm_kingdom: The percentages of the occurence of the specified level
    - organisms: The whole lineage per hit. """
    #Initate 2 empty lists
    organisms = []
    hits = []
    #Iterate over the lines of the given mseq file
    for line in mseq:
        #read the line
        hit = line.decode("utf-8")
        #Check if 'NS' is present in the selected line
        if 'NS' in hit:
            #Split the line on the tab
            hit = hit.split('\t')
            #Save the first value of the line to the list
            hits.append(hit[0])
            #Save the lineage
            layers = hit[13].split(';')[:2]
            #Check the length of the lineage
            #If it's only 1..
            if len(layers) == 1:
                #Add only one time unknown to the lineage
                layers.append('unknown')
            #If it's only 0..
            elif len(layers) == 0:
                #Add 2 times unknown to the lineage
                layers.append(['unknown', 'unknown'])
            #Save the created lineage to the organisms list
            organisms.append(layers)
    #Check the specified level
    #If kingdom is specified..
    if layer == 'Kingdom':
        #Count the occurence of the kingdoms within the lineages
        result = Counter(x[0] for x in organisms)
    #If Phyla is speciefied..
    elif layer == 'Phyla':
        #Count the occurence of the phyla within the lineages
        result = Counter(x[1].split(' ', 1)[0] for x in organisms)
    #The number of 16S/18S rRNA hits
    total = len(mseq)

    # Normalize the occurence dictionary to percentages
    norm_kingdom = {k: v / total * 100 for k, v in result.items()}
    #Return The occurence, the Normalized occurence in percentages, the list containing the lineages
    return result, norm_kingdom, organisms

def create_table(field_name, result):
    """This function prints a table containing the frequency given a dictionary.
    Input:
    - field_name: The name of the column besides Frequency
    - result: dictionary containing the raw counts
    Output:
    - x: the created table containing the frequencies. """
    #Initiate the PrettyTable function
    x = PrettyTable()
    # Create the table names
    x.field_names = [field_name, "Frequency"]
    #Iterate over the dictionary
    for key, value in result.items():
        #Check if the level contains phyla
        if key.endswith('[phylum]'):
            #Change the name of the key
            key = key[:-8]
        #Add the frequency and updated key to the created table
        x.add_row([key, value])
    #Print the table
    print(x)
    #Return the table
    return x


ap = argparse.ArgumentParser(description = 'This script creates multiple stacked plots and csv files containing the raw counts of the 16S/18S annotation at phylum and kingdom level.')

ap.add_argument('-n', '--name', required = True, help = 'Give the base name of how the created graphs should be saved.')
ap.add_argument('-m', '--mseq', required = True, nargs = '+', default = [], help = 'Mseq file(s) of interest')
ap.add_argument('-s', '--sample_list', required = True, nargs='+', default = [], help = 'Give a list of the names of the samples in the same order as the samples are given.')
ap.add_argument('-top', '--top', required = False, default = 'all', help ='Give the top number of phyla to be plotted ')
ap.add_argument('-ex', '--excluding_phylum', required = False, nargs = '+', default = [], help = 'Give the names of the Phyla to be excluded from the graphs if necessary')

args, leftovers = ap.parse_known_args()

#Create 3 empty lists
bar_chart_list = []
bar_chart_names = []
mseq_file_list = []
#Iterate over the given mseq files
for mseq_file in args.mseq:
    #Open the mseq file in read modus
    mseq = open(mseq_file, 'rb')
    #Read the mseq file per line
    mseq = mseq.read().splitlines()
    #Print the top of the mseq file
    print(mseq[:5])
    #Append the mseq file to the list
    mseq_file_list.append(mseq)
    #Initiate 2 empty lists
    sub_bar_chart_list = []
    sub_list_names = []
    #Get the needed dictionaries given the mseq file on Kingdom level
    result, norm_kingdom, organisms = get_percentage(mseq, 'Kingdom')
    #Convert raw counts dictionary to dataframe
    data = pd.DataFrame.from_dict(result, orient='index').reset_index()
    #Rename the columns of the dataframe
    data.columns = ['Kingdom', 'Frequency']
    #Save the dataframe as csv file
    data.to_csv('Kingdom_of_' + i + '_' + args.name)
    #Print the frequency table on Kingdom level
    create_table('Kingdom', result)
    #Iterate over the keys of the raw counts dictionary on Kingdom level
    for i in result.keys():
        #Check if the key is empty, if so, continue
        if i == '':
            continue
        #Initiate 2 empty lists
        domain_count = []
        sub_names = []
        #Iterate over the lineages
        for org in organisms:
            #If the given Phylum is present within the lineage..
            if i in org:
                #Append the Phylum name in the right format to the list
                if ' ' in org[1]:
                    if(org[1][0].isupper()):
                        domain_count.append(org[1].split(' ', 1)[0])
                else:
                    domain_count.append(org[1])
        #Check the occurence of the phyla within the list
        result_p = Counter(domain_count)
        #Convert raw counts dictionary to dataframe
        data = pd.DataFrame.from_dict(result_p, orient='index').reset_index()
        #Rename the columns of the dataframe
        data.columns = ['Phyla', 'Frequency']
        #Save the dataframe as csv file
        data.to_csv('Phyla_of_' + i + '_' + args.name)

        #Sum all the the values
        total = sum(Counter(domain_count).values())
        #Normalize the counts by returning the percentages of occurence
        norm_p = {k: v / total * 100 for k, v in result_p.items()}
        #Print the frequency table of the raw counts
        create_table('Phyla of ' + i, result_p)

#Initiate an empty list
bar_chart_names = []
#Iterate over the normalized kingdom values
for key, value in norm_kingdom.items():
    #Save the keys to the list
    bar_chart_names.append(key)

#Remove duplicates from the reads
keys = list(set(bar_chart_names))
#Create a list of zeros
old_value = [0] * len(args.sample_list)
#Initiate 2 empty lists
key_list = []
df_list = []

#Iterate over the unique keys
for i in keys:
    #Check if empty, if so, continue
    if i == '':
        continue
    #Get a dataframe containing the frequency per phylum and their corresponding standard deviation (ordered)
    otu, std = get_std_table(i, args.name)
    #Combine them
    combi = pd.concat([otu, std], axis = 1)
    #Change last column name
    combi.columns.values[-1] = "std"

    #Sort the DataFrame again on standard deviation
    combi.sort_values(by = 'std', inplace = True, ascending = False)

    #Check if a top value is given
    #If so..
    if args.top != 'all':
        #Select the top rows
        combi = combi.head(int(args.top))
    #Initiate 3 empty lists
    value_list = []
    sub_values_chart = []
    sub_names_chart = []
    #Start a counter
    count = 0
    #Create a list of zeros
    sub_old_value = [0] * len(args.sample_list)
    #Initiate an empty dictionary
    total_dict = {}
    #append the key to a list
    key_list.append(i)
    #Iterate over all hits in all mseq files
    for mseq in mseq_file_list:
        #Get the occurence, the Normalized occurence in percentages, the list containing the lineages, at kingdom level
        result, norm, orgi = get_percentage(mseq, 'Kingdom')
        #append the percentage of a given key
        try:
            value_list.append(norm[i])
        #If the key is not present in the dictionary, append 0
        except:
            value_list.append(0)
            #Print the missing key
            print(args.sample_list[count], i)
        #Iterate over the lineages
        for org in orgi:
            #Check if the key is present within the lineage
            #If so..
            if i in org:

                #Check if the phylum is present in the standard deviation created dataframe
                if org[1] in combi.index:
                    #Convert phylum to the correct format and append it to a list
                    if ' ' in org[1]:
                        if(org[1][0].isupper()):
                            sub_names_chart.append(org[1].split(' ', 1)[0])
                    elif org[1] == 'unknown':
                        continue
                    #Check if the phylum should not be excluded
                    elif org[1] in args.excluding_phylum:
                        continue
                    else:
                        sub_names_chart.append(org[1])
        #Create frequency dictionary
        result_p = Counter(sub_names_chart)
        #Convert it to dataframe
        data = pd.DataFrame.from_dict(result_p, orient='index').reset_index()
        #Change the name of the columns
        try:
            data.columns = ['Phyla', 'Frequency']
        except:
            data = pd.DataFrame(columns = ['Phyla', 'Frequency'])
            #Add the lineage to the phylum column and add zeros to the frequency column
            data['Phyla'] = orgi
            data['Frequency'] = np.zeros(len(orgi))

        #Get the total number of counts
        total = sum(Counter(sub_names_chart).values())
        #Normalize them to percentages
        norm_p = {k: v / total * 100 for k, v in result_p.items()}
        #save the normalized dictionary to another dictionary
        total_dict[args.sample_list[count]] = norm_p
        #add 1 to the counter
        count += 1
    #append the converted and transposed total DataFrame generated from the dictionary of dictionaries
    df_list.append(pd.DataFrame(total_dict).T)
    #plot a stacked bar plot of the percentages at kingdom level and save it
    plt.bar(np.arange(len(args.sample_list)), tuple(value_list), width = 0.35, bottom = old_value)
    old_value = list(map(operator.add, old_value,value_list))
plt.xticks(np.arange(len(args.sample_list)), tuple(args.sample_list))
plt.legend(key_list, loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
#If directory does not exist, create it
if os.path.isdir('16s_18s/frequency_plots') == False:
    os.mkdir('16s_18s/frequency_plots')
plt.savefig('16s_18s/frequency_plots/' + args.name + '_kingdom.png', bbox = 'tight')
plt.show()
#Initiate a counter
counter = 0
#Iterate over the dataframes per mseq file
for df in df_list:
    #Plot a stacked barplot of the percentages at phylum level
    df.plot(kind="bar", stacked=True).legend(fontsize=6,ncol=2, bbox_to_anchor=(1, 1))
    plt.tight_layout()
    #Create a name how to save the plot
    save_name = '16s_18s/frequency_plots/Phyla_of_' + key_list[counter] + '_' + args.name
    if len(args.excluding_phylum) != 0:
        for exc in args.excluding_phylum:
            save_name += '_no_' + exc
    if  args.top != 'all':
        save_name += '_top' + str(args.top)
    save_name += '.png'
    #Save and show the plot
    plt.savefig(save_name, bbox='tight')
    plt.show()
    #Add 1 to the counter
    counter  += 1
