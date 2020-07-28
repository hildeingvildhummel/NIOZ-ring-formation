import pandas as pd
import json
import matplotlib.pyplot as plt

"""This script creates coverage file on genus level per sample based on the station name. Furthermore, it creates a stacked bar plot of the most abundant phyla as annotated by DIAMOND"""

#Initiate the sation name
sample = 'P'
#Open the gene annotation file
with open('KEGG/lineage_{}_tags.txt'.format(sample)) as f:
    data = json.load(f)
#Check the sample name and get the corresponding sample names and the names of the gene coverage files per sample
if sample == 'S':
    samples = ['SO1', 'SO2', 'SR1', 'SR2']
    files = ['KEGG/SO1_samples.count', 'KEGG/SO2_samples.count', 'KEGG/SR1_samples.count', 'KEGG/SR2_samples.count']
elif sample == 'P':
    samples = ['PO1', 'PO2', 'PO3', 'PR1', 'PR2', 'PR3']
    files = ['KEGG/PO1_samples.count', 'KEGG/PO2_samples.count', 'KEGG/PO3_samples.count', 'KEGG/PR1_samples.count', 'KEGG/PR2_samples.count', 'KEGG/PR3_samples.count']
elif sample == 'CD':
    samples = ['CDO1', 'CDO2', 'CDO3', 'CDI1', 'CDI2', 'CDI3', 'CDR2', 'CDR3']
    files = ['KEGG/CDO1_samples.count', 'KEGG/CDO2_samples.count', 'KEGG/CDO3_samples.count', 'KEGG/CDI1_samples.count', 'KEGG/CDI2_samples.count', 'KEGG/CDI3_samples.count', 'KEGG/CDR2_samples.count', 'KEGG/CDR3_samples.count']

#Initiate 2 empty dictionaries
dict = {}
phylum_dict = {}
#Iterate over the gene annotations
for k,lineage in data.items():
    #Save the annotation at genus level, otherwise at the highest known level.
    if len(lineage) < 6 or not lineage[6]:
        if len(lineage) < 5 or not lineage[5]:
            if len(lineage) < 4 or not lineage[4]:
                if len(lineage) < 3 or not lineage[3]:
                    if len(lineage) < 2 or not lineage[2]:
                        if len(lineage) < 0  or not lineage[0]:
                            dict[k] = 'root'
                        else:
                            dict[k] = lineage[0]
                    else:
                        dict[k] = lineage[2]
                else:
                    dict[k] = lineage[3]
            else:
                dict[k] = lineage[4]
        else:

            dict[k] = lineage[5]
    else:

        dict[k] = lineage[6]
    #Save the annotation at phylum level if known, otherwise save 'unknown'
    if len(lineage) < 2 or not lineage[2]:
        phylum_dict[k] = 'unknown'
    else:
        phylum_dict[k] = lineage[2]
#Iniitate a counter
counter = 0
#Iterate over the coverage file names
for i in files:
    #Print the counter
    print(counter)
    #Open and read the file
    f = open(i, 'r')
    f = f.read().split('\n')

    #Create 2 empty lists
    counts = []
    genes = []
    #Iterate over the lines within the coverage file
    for j in f[:-5]:
        #Split on tab
        j = j.split('\t')
        #Append the coverage to the list
        counts.append(j[1])
        #Append the gene name to the list
        genes.append(j[0])
    #Check if counter is 0...
    if counter == 0:
        #If so, create a DataFrame with the genenames as index
        super_df = pd.DataFrame(index = genes)
        #Select the sample name as column name and save the coverage to that column
        super_df[samples[counter]] = list(map(int, counts))
        #Print the length of the created dataframe
        print(len(super_df))
    #if the counter is unequal to 0..
    else:
        #Create a dataframe with the gene names as index
        df = pd.DataFrame(index = genes)
        #Select the sample name as column name and save the coverage to that column
        df[samples[counter]] = list(map(int, counts))
        #Append the created DataFrame to the previously created dataframe
        super_df = pd.concat([super_df, df], axis=1, sort=False)
        #Print the length of the dataframe
        print(len(df))
    #Add 1 to the counter
    counter += 1

#Save the gene names
s = super_df.index.to_series()
#Map the genus annotations to the gene names and save the annotations as index
super_df.index = s.map(dict)
#Initiate a new column with the index as values
super_df['column'] = super_df.index
#Select the columns of interest based on the station name
if sample == 'CD':
    df_g = super_df[['column', 'CDO1', 'CDO2', 'CDO3', 'CDI1', 'CDI2', 'CDI3', 'CDR2', 'CDR3']]
elif sample == 'S':
    df_g = super_df[['column', 'SO1', 'SO2', 'SR1','SR2']]
elif sample == 'P':
    df_g = super_df[['column', 'PO1', 'PO2', 'PO3', 'PR1','PR2', 'PR3']]
#Sum the values if they are annotated the same
df_g = df_g.groupby(['column']).sum()
#Print the top rows of the dataframe
print(df_g.head(5))
#Save the dataframe as csv file
df_g.to_csv('DIAMOND/{}_counts_genus.csv'.format(sample))

#Map the phyla annotations to the gene names and save the annotations as index
super_df.index = s.map(phylum_dict)
#Initiate a new column with the index as values
super_df['column'] = super_df.index
#Select the columns of interest based on the station name
if sample == 'CD':
    df_p = super_df[['column', 'CDO1', 'CDO2', 'CDO3', 'CDI1', 'CDI2', 'CDI3', 'CDR2', 'CDR3']]
elif sample == 'S':
    df_p = super_df[['column', 'SO1', 'SO2', 'SR1','SR2']]
elif sample == 'P':
    df_p = super_df[['column', 'PO1', 'PO2', 'PO3', 'PR1','PR2', 'PR3']]
#Sum the values if they are annotated the same
df_p = df_p.groupby(['column']).sum()
#Sort the values based on the sum of the rows
df_p = df_p.assign(tmp=df_p.sum(axis=1)).sort_values('tmp', ascending=False).drop('tmp', 1)
#Convert the coverages to percentages
df_p = df_p/df_p.sum(axis = 0)*100
#Transpose the dataframe
df_p = df_p.T
#Initiate a figure
fig, ax = plt.subplots()
#Select the top 10 rows
abundant = df_p.iloc[:, :10]
#Create a stacked bar plot of the top rows
abundant.plot(kind='bar', stacked=True, ax = ax)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set(ylabel='Percentage')
#Save and show the figure
plt.savefig('DIAMOND/{}_phyla_abundant_10.png'.format(sample), bbox_inches='tight')
plt.show()
