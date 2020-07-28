import pandas as pd
import numpy as np
import os
import itertools
import matplotlib.pyplot as plt
from scipy import stats



def bowtie_df_idxstats(site, path, save_name):
    """This function creates a csv file containing a dataframe with the coverage per sample per created contig by the co-assembly.
    Input:
    - site: The name of the station, could either be P, S, or CD. String format
    - path: The path name containg the txt files with the coverage per sample
    - save_name: The base name of the output csv file
    Output:
    - total_df: The dataframe containing the coverage of the contigs per sample """
    #Initiate an empty dataframe
    total_df = pd.DataFrame()
    #Iterate over all the files given the path name
    for f in os.listdir(path):
        #Check if the file ends with contigs.txt and starts with the name of the Site (S, P or CD)
        if  f.endswith('contigs.txt') and f.startswith(site):
            #Split on . and select the first part
            name = f.split('.')[0]
            #Open the file in read modus
            file = open(path + f, 'r')
            #Open the file per line
            file = file.read().splitlines()
            #Split the file by tab
            file =[x.split('\t') for x in file]
            #flatten the list
            new_list = [sublist for sublist in file]
            #Convert the list to an array
            new_list = np.array(new_list)
            #Print the file name
            print(name)
            #Create a dataframe with 2 columns, Contig and the name of the file
            df = pd.DataFrame(columns =['Contig', name])
            #Add the created array to the Contig column
            df['Contig'] = new_list[:, 0]
            #Add the counts within the created array as integers to the column with the file name
            df[name] = list(map(int,new_list[:, 2]))
            #Set the contigs as the index
            df = df.set_index('Contig')
            #Save the counts within the dataframe
            values = df.values
            #Convert the counts to a list
            merged = list(itertools.chain.from_iterable(values))
            #Convert the values within the list to integers
            merged = [int(x) for x in merged]
            #Concat the created dataframe to the previous initiated dataframe
            total_df = pd.concat([total_df, df], axis=1)
            #Print the information of the created dataframe and the input
            print('Original: ', len(file))
            print('Kept: ', len(df))
            print('Sum of the reads: ', sum(merged))

            print('Percentage kept: ', float("{:.2f}".format(len(df)/len(file)*100)), '%')
    #Fill the missing values with 0
    total_df.fillna(0, inplace=True)
    #Print the top values of the total dataframe
    print(total_df.head(5))
    #Print the lenght of the total dataframe
    print(len(total_df))
    #Save the total dataframe as a csv file
    total_df.to_csv(path + save_name + '.csv')
    return total_df


def KW_test_diversity(array1, array2, array3=None):
    """This function performs the Kruskal-Wallis test given at least 2 array.
    Input:
    - array1: The first numpy array
    - array2: The second numpy array
    - array3: Optional, the third numpy array)
    Output:
    - Print the statistical measure with its corresponding P-value"""
    #If a third array is given
    if array3 != None:
        #Perform Kruskal-Wallis test
        print(stats.kruskal(array1, array2, array3))
    #If only 2 arrays are given
    else:
        #Perform the Kruskal-Wallis test
        print(stats.kruskal(array1, array2))

"""16S/18S Eukaryotes diversity arrays"""
# chao_CDO = [58.75, 63.2, 68.6]
# chao_CDI = [25.000000, 33.5, 39]
# chao_CDR = [46, 53.2]
# shan_CDO = [1.977174, 2.029174, 2.089919]
# shan_CDI = [1.700087, 1.792492, 1.956412]
# shan_CDR = [1.910991, 1.914247]
#
# chao_PO = [59, 70, 86]
# chao_PR = [89.5, 98, 107.5]
# shan_PO = [2.529132, 2.492969, 2.568359]
# shan_PR = [2.583616, 2.450368, 2.367848]
#
# chao_SO = [40.5, 58.5]
# chao_SR = [72.2, 78.5]
# shan_SO = [2.060866, 2.212016]
# shan_SR = [1.930972, 1.693287]

"""16S/18S Bacteria diversity arrays"""
# chao_CDO = [36, 36, 36]
# chao_CDI = [29,32,34]
# chao_CDR = [35,35]
# shan_CDO = [2.1, 2.11,2.12]
# shan_CDI = [2.02,2.03,2.03]
# shan_CDR = [2.07,2.09]
#
# chao_PO = [25,31,32]
# chao_PR = [33,33,33]
# shan_PO = [1.85,1.85,1.84]
# shan_PR = [1.85,1.85,1.85]
#
# chao_SO = [30,32]
# chao_SR = [34,36]
# shan_SO = [1.8,1.86]
# shan_SR = [1.9,1.89]

"""16S/18S Archaea diversity arrays"""
# chao_CDO = [6,6,6]
# chao_CDI = [5,6,6]
# chao_CDR = [6,6]
# shan_CDO = [0.66, 0.74, 0.75]
# shan_CDI = [0.7,0.75,0.7]
# shan_CDR = [0.69,0.68]
#
# chao_PO = [3,5,6]
# chao_PR = [6,6,6]
# shan_PO = [0.38,0.56,0.92]
# shan_PR = [0.94,0.95,0.92]
#
# chao_SO = [3,4]
# chao_SR = [5,5]
# shan_SO = [0.62,0.76]
# shan_SR = [0.7,0.67]

"""Genes annotated by DIAMOND diversity arrays"""
chao_CDO = [4004.60000, 3967.53191, 3934.88535]
chao_CDI = [4007.51429, 3918.33929, 3873.94118]
chao_CDR = [3901.36683, 3927.44385]
shan_CDO = [5.191524, 5.309860, 5.000179]
shan_CDI = [5.530389, 5.466342, 5.500696]
shan_CDR = [4.966367, 5.058536]
#
chao_PO = [3141.76404, 3140.90576, 3155.87429]
chao_PR = [3218.82258, 3164.89796, 3300.44262]
shan_PO = [4.816715, 4.660854, 4.914905]
shan_PR = [4.896646, 4.694754, 5.014400]

chao_SO = [3074.48551, 2991.48734]
chao_SR = [3018.53191, 2996.55152]
shan_SO = [5.103714, 4.237788]
shan_SR = [4.486575, 4.925228]

#Print the sample name, the diversity measure, the results of the KW test and the corresponding mean and standard deviation
print('P')
print('Chao')
KW_test_diversity(chao_PO, chao_PR)
print('PO: ', np.mean(chao_PO), np.std(chao_PO))
print('PR: ', np.mean(chao_PR), np.std(chao_PR))
print('Shannon')
KW_test_diversity(shan_PO, shan_PR)
print('PO: ', np.mean(shan_PO), np.std(shan_PO))
print('PR: ', np.mean(shan_PR), np.std(shan_PR))
print('S')
print('Chao')
KW_test_diversity(chao_SO, chao_SR)
print('SO: ', np.mean(chao_SO), np.std(chao_SO))
print('SR: ', np.mean(chao_SR), np.std(chao_SR))
print('Shannon')
KW_test_diversity(shan_SO, shan_SR)
print('SO: ', np.mean(shan_SO), np.std(shan_SO))
print('SR: ', np.mean(shan_SR), np.std(shan_SR))

print('C')
print('Chao')
KW_test_diversity(chao_CDI, chao_CDO, chao_CDR)
print('CDI: ', np.mean(chao_CDI), np.std(chao_CDI))
print('CDO: ', np.mean(chao_CDO), np.std(chao_CDO))
print('CDR: ', np.mean(chao_CDR), np.std(chao_CDR))
print('Shannon')
KW_test_diversity(shan_CDI, shan_CDO, shan_CDR)
print('CDI: ', np.mean(shan_CDI), np.std(shan_CDI))
print('CDO: ', np.mean(shan_CDO), np.std(shan_CDO))
print('CDR: ', np.mean(shan_CDR), np.std(shan_CDR))

#Create contig coverage dataframe
path = 'bowtie/'
# bowtie_df_idxstats('S', path, 'S_bowtie')
# bowtie_df_idxstats('P', path, 'P_bowtie')
# bowtie_df_idxstats('C', path, 'C_bowtie')
