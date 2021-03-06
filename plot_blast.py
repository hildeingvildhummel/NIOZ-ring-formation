import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import re
import pandas as pd
from scipy import stats

ap = argparse.ArgumentParser(description='Get the lineage of the protein functions as annotated by KEGG')
ap.add_argument('-f', '--file', nargs='+', required = True, help = 'The output files of BLAST generated by the same reference database')
ap.add_argument('-s', '--samples', nargs='+', required = True, help = 'Sample names')
ap.add_argument('-refs', '--ref_sample', required = True, help = 'Name of the Reference sample used during BLAST')
ap.add_argument('-o', '--output', required = True, help = 'Output name')
ap.add_argument('-t', '--table', required = True, help = 'Give the PROKKA tbl file of the reference sample')


args, leftovers = ap.parse_known_args()
#Initiate an empty list
sample_list = []
#Iterate over the sample names
for i in args.samples:
    #Save the name to the list
    sample_list.append(i)
#Append the name of the referece sample to the list
sample_list.append(args.ref_sample)
#Get the base name of the samples
samples_sub = [x[:-1] for x in sample_list]
samples_sub = list(set(samples_sub))
#Initiate a counter
counter = 0
#Create an empty dictionary
dict = {}
#Iterate over the given BLAST output files
for f in args.file:
    #Open and read the file per line
    file = open(f, 'r')
    file = file.read().splitlines()
    #Iterate over the lines
    for i in file:
        #Split it on tab
        i = i.split('\t')
        #Select the name of the gene
        key = i[1]
        #Select the percent identity
        pident = i[9]
        #Select the coverage
        pcovs = i[10]
        #If the Coverage is below 80\% continue
        if int(pcovs) < 80:
            continue

        #If the counter is equal to 0...
        if counter == 0:
            #Check if the key is present within the dictionary, if so continue
            if key in dict:
                continue
            #If not...
            else:
                #Save the percent identity to the given gene name
                dict[key] = [pident]
        #If the counter is unequal to 0...
        else:
            #Check if the genename is already in the dictionary
            if key in dict:
                #Check the lenght of the value of the given key
                if len(dict[key]) == counter:
                    #If it is equal to the counter append the pident to the value of the given key
                    dict[key].append(pident)
                #Otherwise continue
                else:
                    continue
            #If the key is not present within the dictionary..
            else:
                #Create a list of zeros
                list_1 = [0] * counter
                #Append the percent identity to this list
                list_1.append(pident)
                #Save the list with its corresponding gene name
                dict[key] = list_1

    #Iterate over the dictionary keys
    for k in dict.keys():
        #Check if all the values have a correct length
        if len(dict[k]) != counter + 1:
            #Otherwise, append 0
            dict[k].append(0)
    #Add 1 to the counter
    counter += 1
#Print the number of genes
print('key length: ', len(dict.keys()))

#Open the prokka table file
annotation_file = open(args.table, 'r')
#Get the base name of the genes
d = list(dict.keys())[0].split('_')[0]
#Print the base name of the genes
print(d)
#Split the table file based on the base name
annotation_file = [d+e for e in annotation_file.read().split(d) if e]

#Create an empty list
annotation = []
#Iterate over the dictionary keys
for k in dict.keys():
    #Select the parts of the table file containing the key
    interest = [s for s in annotation_file if k+'\n' in s]
    #Save the annotation of prokka
    annotation.append([i.split('\t')[4].split('\n')[0] for i in interest])
#Print the number of annotations
print('length annotation: ', len(annotation))

#Flat the list of annotations
flat_list = [item for sublist in annotation for item in sublist]
#save the dictionary values to a list
data = list(dict.values())
#Convert the list to an array
an_array = np.array(data, dtype=np.float64)

#Set the array in a dataframe with sample names as columns and annotations as index
df = pd.DataFrame(data=an_array, index=flat_list, columns=args.samples)
#Add a reference sample column to the dataframe only containing 100%
df[args.ref_sample] = [100]*len(df)
#Print the length of the dataframe
print('length dataframe: ', len(df))
#Save the dataframe as csv file
df.to_csv(args.output[:-4] + '.csv')

#Create heatmap of the dataframe
g = sns.clustermap(df, cmap="vlag")
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
#Save and show the heatmap
plt.savefig(args.output, bbox_inches = "tight")
plt.show()

#Check the name of the reference sample
#Based on the name of the reference sample, extract the rows where the replicates of the reference are equal to 100 and all the other samples are unequal to 100.
if args.ref_sample == 'PO1':
    new_df = df[(df[args.ref_sample] == 100) & (df[args.ref_sample[:-1] + '2'] == 100) & (df[args.ref_sample[:-1] + '3'] == 100) & (df['PR1'] != 100) & (df['PR2'] != 100) & (df['PR3'] != 100)]
elif args.ref_sample == 'PR1':
    new_df = df[(df[args.ref_sample] == 100) & (df[args.ref_sample[:-1] + '2'] == 100) & (df[args.ref_sample[:-1] + '3'] == 100) & (df['PO1'] != 100) & (df['PO2'] != 100) & (df['PO3'] != 100)]
elif args.ref_sample == 'SR1':
    new_df = df[(df[args.ref_sample] == 100) & (df[args.ref_sample[:-1] + '2'] == 100) & ((df['SO1'] != 100) & (df['SO1'] > 60)) & ((df['SO2'] != 100) & (df['SO2'] > 60))]
elif args.ref_sample == 'SO1':
    new_df = df[(df[args.ref_sample] == 100) & (df[args.ref_sample[:-1] + '2'] == 100) & (df['SR1'] != 100) & (df['SR2'] != 100)]
elif args.ref_sample == 'CDI1':
    new_df = df[(df[args.ref_sample] == 100) & (df[args.ref_sample[:-1] + '2'] == 100) & (df[args.ref_sample[:-1] + '3'] == 100) & (df['CDO1'] != 100) & (df['CDO2'] != 100) & (df['CDO3'] != 100) & (df['CDR2'] != 100) & (df['CDR3'] != 100)]
elif args.ref_sample == 'CDR2':
    new_df = df[(df[args.ref_sample] == 100) & (df[args.ref_sample[:-1] + '3'] == 100) & (df['CDO1'] != 100) & (df['CDO2'] != 100) & (df['CDO3'] != 100) & (df['CDI1'] != 100) & (df['CDI2'] != 100) & (df['CDI3'] != 100)]
elif args.ref_sample == 'CDO1':
    new_df = df[(df[args.ref_sample] == 100) & (df[args.ref_sample[:-1] + '2'] == 100) & (df[args.ref_sample[:-1] + '3'] == 100) & (df['CDI1'] != 100) & (df['CDI2'] != 100) & (df['CDI3'] != 100) & (df['CDR2'] != 100) & (df['CDR3'] != 100)]
#Save the selected rows as a csv file
new_df.to_csv(args.output[:-4] + '_difference.csv')
#Print the number of selected rows.
print('new df length: ', len(new_df))
