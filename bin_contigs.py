import pandas as pd

def get_bin_csv(b, s, path):
    """This function creates an abundance dataframe of the contigs, given a bin. The abundance dataframe is saved as a csv file.
    Input:
    - b: The bin number in string format
    - s: The sample, could either be S, P, or C. In string format
    - path: The path to the text files containing the contig names per bin
    Output:
    - selected: The dataframe with the abundance of the given contigs within the selected bin."""
    # Print the Site and the Bin number
    print(s, b)
    # Open the text file containing the Contig names of the specified bin
    file = open(path + 'bin_%s%s.txt' % (b, s), 'r')
    #Read the file per line, creating a list
    file = file.read().splitlines()

    #Read the dataframe containing the abundance per contig, given the Site
    counts = pd.read_csv('bowtie/{}_bowtie.csv'.format(s), index_col = 0, header = 0)
    # Select the rows of the dataframe that contain the Contig names within the text file
    try:
        selected = counts.loc[file]
    except:
        #if error occurs, split the contig name prior to selecting the rows of the dataframe
        file = [i.split(' ')[0] for i in file]
        selected = counts.loc[file]
    #Print the top 5 rows of the created dataframe
    print(selected.head(5))
    #Save the created dataframe as csv file
    selected.to_csv(path + 'bin_%s%s.csv' % (s, b))
    #Return the created dataframe 
    return selected
