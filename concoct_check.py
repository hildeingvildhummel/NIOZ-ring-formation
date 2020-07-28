import pandas as pd
import argparse

ap = argparse.ArgumentParser(description = 'This script converts the csv coverage file of the Concoct output to the correct format txt file containing a single txt file per bin containing at least 5 significant contigs, with the contigs within the bin.')
ap.add_argument('-c', '--concoct', required = True, help = 'The concoct csv file')
ap.add_argument('-f', '--textfile', required = True, help = 'Text file containing the significant contigs')
ap.add_argument('-o', '--output', required = True, help = 'The output basename')

args, leftovers = ap.parse_known_args()

#Open the Concoct coverage file
file = pd.read_csv(args.concoct, header = None)
#Inititae the column names
file.columns = ['Contig', 'Bin']
#Set the contig names as index
file_2 = file.set_index('Contig')
#Open the txt file with the contigs of interest
significant = open(args.textfile, 'r')
#Read the txt file per line
significant = significant.read().splitlines()
#Select the bin numbers of the contigs of interest
sub = file_2.loc[significant].dropna()
#Reset the index
sub = sub.reset_index()
#Group by the unique bin number and sum the occurence of the bin number
counts = sub.groupby('Bin')['Contig'].nunique()
#Print the created dataframe
pd.set_option('display.max_rows', None)
print(counts)

#Selected the bins if the occurence is above 5
interesting_bins = counts.where(counts > 5).dropna().index.astype(int).to_list()
#Iterate over the bins of interest
for i in interesting_bins:
    #Select the contigs of the interesting bins
    contigs = file.loc[file['Bin'] == i]['Contig'].to_list()
    #Save the contigs of this specific bin to a txt file in the correct format
    with open('bin_' + str(i) + args.output + '.txt', 'w') as f:
        for item in contigs:
            f.write("%s\n" % item)
