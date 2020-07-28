import argparse
import pandas as pd
import json
import numpy as np
from sklearn.decomposition import PCA
import random
import itertools
import matplotlib.pyplot as plt

def stacked_bar(nitrogen_df, title, save, n_col = 1):
    """Create a stacked bar plot of a dataframe.
    Input:
    - nitrogen_df: The dataframe
    - title: Title of the bar plot
    - save: The save name of the plot
    - n_col = the number of columns of the legend. Default is 1"""
    #Create a stacked bar plot
    fig, ax = plt.subplots()
    nitrogen_df.plot(kind='bar', stacked=True, ax = ax)
    #Set the title
    ax.set_title(title)
    #Set the legend
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #Set the y-label
    ax.set(ylabel='Percentage')
    #Save and show the plot
    plt.savefig(save, bbox_inches='tight')
    plt.show()


def plot_stacked_bar(super_df, pathway, save):
    """This function creates stacked bar graphs showing the percentage abundance of genes within a specified pathway.
    Input:
    - super_df: A dataframe containing the information to plot
    - pathway: Name of the patway of interest. If no pathway of interest is given, set to True for distribution of all the pathways and set to False to get phylum distribution of all the pathways.
    - save: The name of the plot to be saved.
    Output:
    - stacked bar plot"""
    #If a pathway is not specified
    if pathway == True:
        #Extract pathway information from the dataframe and sum the coverage if they are annotated the same
        nitrogen_df = super_df.groupby('Pathway').sum()
        #Print the created dataframe
        print(nitrogen_df)
    #If a pathway is specified..
    elif pathway != False:
        #Extract the pathway of interest from the dataframe
        nitrogen_df = super_df.loc[super_df['Pathway'] == pathway]
        #Group by module and the sum the coverages if they are annotated the same
        nitrogen_df = nitrogen_df.groupby('Module').sum()
    #If nothing is specified..
    else:
        #Extract phylum information from the dataframe and sum the coverages if they are annotated the same
        nitrogen_df = super_df.groupby("Phylum").sum()
        #Sort the values based on the total abundance
        nitrogen_df = nitrogen_df.assign(tmp=nitrogen_df.sum(axis=1)).sort_values('tmp', ascending=False).drop('tmp', 1)
    #Scale the coverages and transform it to percentages
    nitrogen_df = nitrogen_df/nitrogen_df.sum(axis = 0)*100
    #Transpose the dataframe
    nitrogen_df = nitrogen_df.T
    # If a pathway is specified
    if pathway != False:
        stacked_bar(nitrogen_df, 'Module distribution', save)
    else:
        stacked_bar(nitrogen_df, 'Phylum distribution', save, ncol=2)
        abundant = nitrogen_df.iloc[:,:10]
        stacked_bar(abundant, 'Most abundant phylum distribution', save[:-4] + '_abundant_10.png')
        rare = nitrogen_df[nitrogen_df.columns[(nitrogen_df<=0.05).any()]]
        stacked_bar(rare, 'Rare phylum distribution', save[:-4] + '_rare.png')

ap = argparse.ArgumentParser(description='Get the lineage of the function annotated by KEGG. Get stacked bar plots of the distributions of the pathways of interest and there corresponding most abundant organisms. Finally, a total image of the annotation is created by PCA. The abundance is saved as a csv file, if a gene is assigned to more than 1 K-number, this corrected by dividing the abundance by the number of assinged K-numbers.')
ap.add_argument('-f', '--file', nargs='+', required = True, help = 'The output files of htseq, seperate the files by space')
ap.add_argument('-s', '--samples', nargs='+', required = True, help = 'Sample names')
ap.add_argument('-K', '--kegg', required = True, help = 'Kegg file')
ap.add_argument('-db', '--database', required = True, help = 'The KEGG database KO')
ap.add_argument('-o', '--output', required = True, help = 'Output name')
ap.add_argument('-a', '--gene_annotation', required = True, help = 'the _tags.txt output of get_lineage.py')


args, leftovers = ap.parse_known_args()

#Try to read the csv file
try:
    super_df = pd.read_csv(args.output)
#If the file does not exist yet, create the csv file
except:
    #Initiate a counter
    counter = 0
    #Open the Kegg file and read it
    kegg_f = open(args.kegg, 'r')
    kegg_file = kegg_f.read()
    #Open the KO database
    with open(args.database) as json_file:
        data = json.load(json_file)
    #Initiate an empty dictionary
    dict = {}
    #Iterate over the levels within the database
    for item in data['children']:
        for i in item['children']:
            for j in i['children']:
                try:

                    for k in j['children']:
                        #extract the annotation
                        string = ' '.join(k['name'].split(' ')[2:])
                        #create dictionary with the K numbers and its function
                        dict[k['name'].split(' ')[0]] = string.split(';')[1]
                #If no K-number is present within this level..
                except:
                    #Print the level
                    print(j['name'])

    #Create the delimiter
    d = '\tM'
    #Split the kegg file on module
    kegg = ['M'+e for e in kegg_file.split(d) if e]
    #Split the kegg file on pathway
    kegg_path = kegg_file.split('\n\n')
    #Open the gene annotation as created by get_lineage.py
    with open(args.gene_annotation) as json_file:
        gene_lineage = json.load(json_file)
    #Iterate over htseq files containing the coverage per gene per sample.
    for i in args.file:
        #Print the counter
        print(counter)
        #Open and read the coverage file in the correct format
        f = open(i, 'r')
        f = f.read().split('\n')
        #Initiate 10 empty lists
        counts = []
        genes = []
        pathways = []
        modules = []
        phyla = []
        classes = []
        orders = []
        families = []
        genuses = []
        functions = []
        #Iterate over the lines within the coverage file containing the coverage information
        for j in f[:-3]:
            #Split the lines on the tab
            j = j.split('\t')
            #Get the indexes of the modules containing the gene
            indices = [i for i, x in enumerate(kegg) if j[0] in x]
            #Check the number of hits
            #If there are any hits found..
            if len(indices) != 0:
                #Iterate over the hits
                for i in indices:
                    #Extract the hits in the right format
                    res = [i for i in kegg[i].split('\n\t\t')]
                    #Extract the final information
                    indices_res = [i for i, x in enumerate(res) if j[0] in x]
                    #Iterate over the indexes of interest
                    for index in indices_res:
                        #Check the format and correct for it, extracting the name of the gene
                        if '\n\n\n' in res[index]:
                            gene_names = res[index].split('\n\n\n')[0]
                            gene_names = gene_names.split(' ')[1: ]
                        else:
                            gene_names = res[index].split(' ')[1: ]
                        #Split the genenames
                        gene_names = [i.split(',') for i in gene_names]
                        #Convert list to the correct format
                        gene_names = list(itertools.chain(*gene_names))
                        #Extract the K-number
                        k = res[index].split(' ')[0]
                        #Check whether an '\n' is present within the gene name, if so remove it.
                        try:
                            gene_names[-1] = gene_names[-1].split('\n')[0]
                        except:
                            gene_names[-1] = gene_names[-1]
                        #Set present to No
                        present = 'No'

                        #Iterate over the gene names
                        for g in gene_names:
                            #Check if it is an exact match
                            if g == j[0]:
                                #If so, set present to yes
                                present = 'Yes'
                        #If present is No, continue
                        if present == 'No':
                            continue
                        #If present is Yes
                        else:
                            #Append the function of the K-number to a list
                            functions.append(dict[k])
                            #Append the gene name to a list
                            genes.append(j[0])
                            #Append the coverage (normalized for the number of appearances) to a list
                            counts.append(int(j[1])/len(indices))
                            #Check the counter, is it is 0..
                            if counter == 0:
                                #Extract the Module from the Kegg file
                                Module = ' '.join(kegg[i].split('\n\t\t')[0].split(' ')[1:])
                                #Save the module to the list
                                modules.append(Module)
                                #Try to append the phylum to a list, if it is known.
                                try:
                                    phyla.append(gene_lineage[j[0]][2])
                                except:
                                    phyla.append('unknown')
                                #Try to append the class to a list, if it is known.
                                try:
                                    classes.append(gene_lineage[j[0]][3])
                                except:
                                    classes.append('unknown')
                                #Try to append the order to a list, if it is known.
                                try:
                                    orders.append(gene_lineage[j[0]][4])
                                except:
                                    orders.append('unknown')
                                #Try to append the family to a list, if it is known.
                                try:
                                    families.append(gene_lineage[j[0]][5])
                                except:
                                    families.append('unknown')
                                #Try to append the genus to a list, if it is known.
                                try:
                                    genuses.append(gene_lineage[j[0]][6])
                                except:
                                    genuses.append('unknown')
                                #Get the index of interest from the the Pathway splitted kegg file
                                indices_path = [i for i, x in enumerate(kegg_path) if Module in x]
                                #Extract the pathway the Module belongs to
                                for p in indices_path:
                                    path = kegg_path[p]
                                #Convert the pathway to the correct format
                                str_list = list(filter(None, path.split('\n')))
                                pathway = str_list[0]
                                #Save the pathway to a list
                                pathways.append(pathway)
        #If the counter is equal to 0 ..
        if counter == 0:
            #Initiate a dataframe with the gene names as index
            super_df = pd.DataFrame(index = genes)
            #Create a column called Pathway and save the found pathways to it
            super_df['Pathway'] = pathways
            #Create a column called Module and save the found modules to it
            super_df['Module'] = modules
            #Create a column called Function and save the found functions to it
            super_df['Function'] = functions
            #Create a column called Phylum and save the found phyla to it
            super_df['Phylum'] = phyla
            #Create a column called Class and save the found classes to it
            super_df['Class'] = classes
            #Create a column called Family and save the found families to it
            super_df['Family'] = families
            #Create a column called Genus and save the found genuses to it
            super_df['Genus'] = genuses
            #Create a column named after the sample of interest and save the coverage to it
            super_df[args.samples[counter]] = counts
            #Print the length of the dataframe
            print(len(super_df))
        #If the counter is unequal to 0..
        else:
            #Create a dataframe with the genes as index
            df = pd.DataFrame(index = genes)
            #Create a column named after the sample of interest and save to coverage to it
            df[args.samples[counter]] = counts
            #Combine the DataFrame to the DataFrame created by counter == 0 by index
            super_df = pd.concat([super_df, df], axis=1, sort=False)
            #Print the length of the dataframe
            print(len(df))
        #Add 1 to the counter
        counter += 1
    #Save the concatenated dataframe as csv file
    super_df.to_csv(args.output)

#Initiate a counter
count = 0
#Create an empty list
pathway_plot = []
#Initialize a figure
fig, axs = plt.subplots(1, 4)

#Iterate over the Pathway column of the dataframe
for selected_pathway, df_pathway in super_df.groupby('Pathway'):
    #append the name of the pathway to the list
    pathway_plot.append(selected_pathway)
    #Group the pathways by phylum
    gk = df_pathway.groupby('Phylum')
    #Get the total coverage
    gk = gk.sum()
    #Scale the coverage by their total coverage
    gk = gk/gk.sum(axis = 0)
    #Transpose the created dataframe
    gk = gk.T
    #Select only the phyla with an relative abundance above 0.05
    gk = gk[gk.columns[(gk>0.05).any()]]
    #Create an empty list
    SO1_list = []
    #Check if the counter is 0..
    if count == 0:
        #Create a dictionary with the column values as keys and random values
        col_dict = {k:np.random.rand(3,) for k in gk.columns.unique()}
        #Convert the keys to a list
        labels = list(col_dict.keys())
        #Create the box for the legend of the plot
        handles = [plt.Rectangle((0,0),1,1, color=col_dict[label]) for label in labels]
    #If the counter is unequal to 0...
    else:
        #Iterate over the unique column values of the dataframe
        for key in gk.columns.unique():
            #Check if the key is not present within the random value dictionary..
            if key not in col_dict.keys():
                #If not present, add the key with a random value
                col_dict[key] = np.random.rand(3,)
            #Overwrite the labels list and the box of the legend with the newly created dictionary
            labels = list(col_dict.keys())
            handles = [plt.Rectangle((0,0),1,1, color=col_dict[label]) for label in labels]
    #Iterate over the column values of the dataframe
    for column in gk.columns:
        #Append a bar of the barplot to a list of a single sample with its corresponding random value of the dictionary as color.
        SO1_list.append(axs[count].bar(args.samples, gk[column], align='center', width= 0.4, color=[col_dict[column]]))
        #Set the name of the pathway as the title
        axs[count].set_title(selected_pathway)
        #Set the sample names as x labels.
        axs[count].set_xticklabels(args.samples, rotation=90, fontsize=7)
        #If count is 0..
        if count == 0:
            #Set the y label to Ratio
            axs[count].set(ylabel='Ratio')
    #Add 1 to the counter
    count += 1
#make sure the complete plot is visible
fig.tight_layout()
#Add the legend
plt.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
#Save and show the plot
plt.savefig(args.output[:-4] + '.png', bbox_inches='tight')
plt.show()

#get the stacked bar plot of the Phyla
plot_stacked_bar(super_df, pathway = False, save = args.output[:-4] + '_Phyla.png')

#Get the stacked bar plot of the distribution of the pathways
plot_stacked_bar(super_df, pathway = True, save = args.output[:-4] + '_Pathways.png')

#Get the stacked bar plot of the distribution of the Nitrogen metabolism
plot_stacked_bar(super_df, pathway = 'Nitrogen metabolism', save = args.output[:-4] + '_Nitrogen.png')

#Get the stacked bar plot of the distribution of the Sulfur metabolism
plot_stacked_bar(super_df, pathway = 'Sulfur metabolism', save = args.output[:-4] + '_Sulfur.png')

#Get the stacked bar plot of the distribution of the Carbon fixation
plot_stacked_bar(super_df, pathway = 'Carbon fixation', save = args.output[:-4] + '_Carbon.png')

#Get the stacked bar plot of the distribution of the Photosynthesis
plot_stacked_bar(super_df, pathway = 'Photosynthesis', save = args.output[:-4] + '_Photosynthesis.png')

#Scale the coverage
total = super_df[args.samples]/super_df[args.samples].sum(axis = 0)
#Get the coverage of the samples within the given pathways and transpose it
total = total.T
#Get PCA function with 2 PCAs
pca = PCA(n_components=2)
#Fit the PCA transform to the KEGG annoted genes
PCA_proteins = pca.fit_transform(total)
#Create an empty list
explained_variance = []
#Iterate over the explained variances
for i in pca.explained_variance_ratio_:
    #Transform it to percentages
    a = float(i) * 100
    #Save it to the list
    explained_variance.append("%.2f" % a)
#Create a PCA dataframe
principal_df = pd.DataFrame(data = PCA_proteins, columns = ['principal component 1', 'principal component 2'])
#Initiate a figure
plt.figure(figsize=(10,10))
plt.xticks(fontsize=12)
plt.yticks(fontsize=14)
#Set the explained variances percentages as x and y labels
plt.xlabel('Principal Component - 1 ({}%)'.format(explained_variance[0]),fontsize=20)
plt.ylabel('Principal Component - 2 ({}%)'.format(explained_variance[1]),fontsize=20)
#Set the title
plt.title("Principal Component Analysis",fontsize=20)

#iteratate over the samples
for target in args.samples:
    #Keep the information about the samples
    indicesToKeep = total.index == target
    #Create the scatter plot of the PCA dataframe
    plt.scatter(principal_df.loc[indicesToKeep, 'principal component 1']
               , principal_df.loc[indicesToKeep, 'principal component 2'], c = np.random.rand(3,), s = 50)
#Iterate over the samples again
for i, txt in enumerate(args.samples):
    #Assign the correct label to the points within the plot
    plt.annotate(txt, (principal_df.loc[i, 'principal component 1'], principal_df.loc[i, 'principal component 2']))

#Save and show the figure
plt.savefig(args.output[:-4] + '_PCA.png', bbox_inches='tight')
plt.show()
