import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
from bin_contigs import get_bin_csv
from os import listdir
from os.path import isfile, join

"""This script creates heatmaps of the given created bins per station of interest.
Input:
- sample: The station of interest. Could either be S, P, or C. String format
- path: The path where the txt files of the bins are saved
Output:
- Heatmap of the bins of interest within the given path of the station of interest """

#Specify the name of the site. Could either be P, S or C. String format
sample = 'C'
#Specify the path to the txt files containing the contigs per bin
path = 'Bins/concoct/'
# path = 'Bins/maxbin/'
# path = 'Bins/metabat/'

#Initiate an empty list
bins = []
#Extract the files of interest from the path
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
#Iterate over these files
for i in onlyfiles:
    #Check whether these files belong to the station of interest and if it is a txt file
    if sample in i and i.endswith('.txt'):
        #Extract the name of the bin
        number = i.split('_')[1].split(sample)[0]
        #Save it to the list
        bins.append(number)
#Create a list of station of the number of bins
sites = [sample]*len(bins)
#Print the bins
print(bins)
#Iterate over the list of stations and the bin names
for b, s in zip(bins, sites):
    #If annotation is set to True..
    annotation = True
    #Try to open the csv file
    try:
        file = pd.read_csv(path + 'bin_%s%s.csv' % (s, b), index_col = 0)
    #If it does not exist, create it
    except:
        file = get_bin_csv(b, s, path)
    #Open the coverage file
    total = pd.read_csv('bowtie/{}_bowtie.csv'.format(s), index_col=0)
    #Get the total coverage per sample
    total_sum = total.sum(axis = 0, skipna = True)
    #Perform total sum scaling
    norm_file = file/total_sum
    #If annotation is set to True..
    if annotation == True:
        #Open the annotation file
        with open('lineage_{}.txt'.format(s)) as json_file:
            data = json.load(json_file)
        #Initiate an empty list
        annotation = []
        #Iterate over the contigs within the bin
        for x in norm_file.index:
            #Try to annotate it
            try:
                annotation.append(data[x])
            #If no annotation is known, remove the contig
            except:
                norm_file.drop(x, inplace = True)
        #Add the annotation to the dataframe
        norm_file['annotated'] = annotation
        #Sum the coverage of the contigs that are annotated the same
        norm_file = norm_file.groupby("annotated").sum()
    #Print the top rows of the created dataframe
    print(norm_file.head(5))
    #Try to create an heatmap of the dataframe
    try:
        g = sns.clustermap(norm_file, cmap="vlag")
    except:
        plt.close()
        try:
            norm_file = file/total_sum
            g = sns.clustermap(norm_file, cmap="vlag")
        #If not possible, only a single contig is present within the bin..
        except:
            print('single sample bin')
            plt.close()
            continue
    #Rotate the labels
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    #save and close the figure
    plt.savefig(path + '%s_%s_contigs_summed.png' % (s, b), bbox_inches = "tight")
    plt.close()
