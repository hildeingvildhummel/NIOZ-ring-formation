import json
import argparse
import subprocess
import numpy as np

ap = argparse.ArgumentParser('This function extracts the gene names if they are annotated as a indicated organism. It returns a list file containing a single gene name by line.')
ap.add_argument('-f', '--json_file', required = True, help = 'The json file containing the annotations per Prokka tag')
ap.add_argument('-org', '--organism', required = True, help = 'Name of the organism to search for')
args, leftovers = ap.parse_known_args()

def search(values, searchFor):
    """This function looks for an item within the values of a given dictionary.
    Input:
    - values: The dictionary to look into
    - searchFor: The item to look for
    Output:
    - key_list: List of keys where the values contains the item"""
    #Initiate an empty list
    key_list = []
    #Iterate over the dictionary
    for k, v in values.items():
        #If the item is in the value..
        if searchFor in v:
            #Save the key to the list
            key_list.append(k)
    #Return the list with the keys
    return key_list

try:
#Try to open the created tags file by the get_lineage.py
    with open(args.json_file) as f:
        data = json.load(f)
except:
    #If this file does not exist, print this message
    print('The required json file does not excist yet. Run the get_lineage.py script first')

#Look in the tags dictionary for the given organism
keys = search(data, args.organism)

#Open a txt file in write modus
with open(args.organism + '_name.list', 'w') as f:
    #Save the list of gene names as txt file 
    for item in keys:
        f.write("%s\n" % item)
