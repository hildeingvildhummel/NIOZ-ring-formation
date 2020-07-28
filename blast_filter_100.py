import argparse

ap = argparse.ArgumentParser(description = 'Remove the genes if all the samples have a 100% identity to the reference database.')
ap.add_argument('-b', '--blast', nargs='+', required = True, help = 'The file(s) containing the blast output')
ap.add_argument('-index', '--index', required = True, help = 'The index of the column containing the pIDENT, starting the count with zero. So, the first column is 0, the second column is 1, etc.')
ap.add_argument('-o', '--output', required = True, help = 'The output name of the list containing the IDs to keep. The output is saved as list file')

args, leftovers = ap.parse_known_args()

#Create 2 empty dictionaries
dict = {}
dict_list = {}
#Iterate over the blast output files
for f in args.blast:
    #Print the name of the blast output file
    print(f)
    #Create an empty list
    key_list = []
    #Open the blast output file in read modus
    blast = open(f, 'r')
    #Read the Blast output file per line, creating a list
    blast = blast.read().splitlines()
    #Select the column within the blast file containing the pIDENT
    blast_ident = [s.split('\t')[int(args.index)] for s in blast]
    #Select the gene ID of PROKKA of the reference database
    ID = [s.split('\t')[1] for s in blast]
    #Select the gene ID of PROKKA of the query
    ID_u = [s.split('\t')[0] for s in blast]
    #Iterate over both gene names (reference and query) and the pIDENT
    for i, j, k in zip(ID, blast_ident, ID_u):
        #Check if the reference ID is already in the empty list
        if i in key_list:
            #If so, append the ID of the query to the dictionary with the given reference ID as reference
            dict_list[i].append(k)
            continue
        #If the refrence ID is already noted within the second dictionary..
        if i in dict:
            #Add the given pIDENT to the already existing value of the key given the reference ID
            dict[i] = float(dict[i]) + float(j)
            #Append the query ID to the dictionary given the reference ID
            dict_list[i].append(k)
            #Append the reference ID to the list
            key_list.append(i)

        else:
            #If the reference ID is not saved in either the list or the first dictionary..
            # Add the pIDENT to the dictionary with the reference ID as key
            dict[i] = float(j)
            # Add the query ID as a list to the dictionary with the reference ID as key
            dict_list[i] = [k]
            #Save the reference ID to the list
            key_list.append(i)
#Create a new dictionary, only containing the values and keys if the value is unequal to the number of samples * 100
d = { k : v for k,v in dict.items() if v != len(args.blast) * 100}

#Save the query ID values if their corresponding key is present within the filtered dictionary d
ID_list = [dict_list[x] for x in list(d.keys())]
#Flatten the list
flatten_list = [val for sublist in ID_list for val in sublist]

#Save the output
with open(args.output + '.list', 'w') as f:
    for item in list(set(flatten_list)):
        f.write("%s\n" % item)
