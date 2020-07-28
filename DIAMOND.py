import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np

def plot_PCA(csv_file, save_name):
    """This function creates a PCA plot based on the annotation by KEGG
    Input
    - csv_file: The annotation csv file of the KEGG annotations
    - save_name: The name how to save the plot
    Output:
    A save PCA figure."""
    #Open and read the csv file
    total = pd.read_csv(csv_file, index_col=0)
    #Select all rows and from the seventh column
    total = total.iloc[:, 7:]
    #Initiate a PCA with 2 components
    pca = PCA(n_components=2)
    #Perform the PCA on the data
    PCA_proteins = pca.fit_transform(total.T)
    #Create an empty list
    explained_variance = []
    #Iterate over the explained variance of the PCA
    for i in pca.explained_variance_ratio_:
        #Convert it to a percentage
        a = float(i) * 100
        #Save the percentage to the list with 2 decimals
        explained_variance.append("%.2f" % a)
    #Create a PCA dataframe
    principal_df = pd.DataFrame(data = PCA_proteins, columns = ['principal component 1', 'principal component 2'])
    #Craete a figure
    plt.figure(figsize=(10,10))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=14)
    #Set the explained variance as the x and y labels
    plt.xlabel('Principal Component - 1 ({}%)'.format(explained_variance[0]),fontsize=20)
    plt.ylabel('Principal Component - 2 ({}%)'.format(explained_variance[1]),fontsize=20)
    #Get the columns
    targets = list(total.columns)
    #Print them
    print(targets)
    #Iterate over the column names
    for target in targets:
        #Convert it to indices
        indicesToKeep = total.columns == target
        #Plot the scatter plot
        plt.scatter(principal_df.loc[indicesToKeep, 'principal component 1']
                   , principal_df.loc[indicesToKeep, 'principal component 2'], c = np.random.rand(3,), s = 50)
    #Annotate the points in the PCA scatter plot
    for i, txt in enumerate(targets):
        plt.annotate(txt, (principal_df.loc[i, 'principal component 1'], principal_df.loc[i, 'principal component 2']))
    #Save and show the figure
    plt.savefig(save_name, bbox_inches='tight')
    plt.show()

def print_most_abundant_organisms(lineage_file, only_out = False):
    """This function prints the most abundant organisms within the annotation file
    Input:
    - lineage_file: lineage csv file
    - only_out: if set to True, only the control mat will be taken into account """
    #Read the lineage file
    df = pd.read_csv(lineage_file, index_col = 0)
    #Sum per sample
    total_sum = df.sum(axis = 0, skipna = True)
    #Normalize by TSS
    df = df/total_sum
    #Check if only_out is set to True..
    if only_out == True:
        #Filter columns on 'O'
        df = df[df.filter(like='O').columns]
    #Get the mean and sort by it
    df['mean'] = df.mean(axis = 1)
    df = df.sort_values(by=['mean'], ascending = False)
    #Iterate over the normalized coverages
    for i in range(len(df)):
        #If the total sum is below 0.6..
        if df.iloc[:i + 1, -1].sum() >= 0.6:
            #Print the total normalized coverage of the selected samples
            print('The total coverage of these organism(s) is/are : {:.2f}%'.format(df.iloc[:i + 1, -1].sum() * 100))
            #Remove the mean column
            df.drop(['mean'], axis=1, inplace = True)
            #Print the selected rows
            print(df.iloc[:i+1, :])
            break

def scatter_plot_significance(bowtie_csv, top_value_csv, save_name, top_value_csv_2 = None, top_value_csv_3 = None):
    """This function plots all contigs by PCA and colors the significant contigs red.
    Input
    - bowtie_csv: The coverage file per contig
    - top_value_csv: The results of the ALDEx2 test
    - top_value_csv_2: If more than one result is generated..
    - top_value_csv_3: If more than one result is generated..
    Output:
    Show and save the generated PCA plots """
    #Read the coverage file
    df = pd.read_csv(bowtie_csv, header=0, index_col = 0)
    #Save the index to list
    all_contigs = df.index.to_list()
    #Read the file containing the results of ALDEx2
    sign_df = pd.read_csv(top_value_csv, header=0, index_col =0)
    #Save the index to list
    contigs = sign_df.index.to_list()
    #Check if more than 1 result is given..
    if top_value_csv_2 != None:
        #Read the file containing the results of ALDEx2
        sign_df_2 = pd.read_csv(top_value_csv_2, header=0, index_col =0)
        #Save the index to list
        contigs_2 = sign_df.index.to_list()
        #Read the file containing the results of ALDEx2
        sign_df_3 = pd.read_csv(top_value_csv_3, header=0, index_col =0)
        #Save the index to list
        contigs_3 = sign_df.index.to_list()
    #Initate an empty list
    labels = []
    #If more than 1 result is given...
    if top_value_csv_2 != None:
        #Iterate over all the contigs
        for i in all_contigs:
            #If present in one of the results..
            #Save 'Significant' to the list
            if i in contigs:
                labels.append('Significant')
            elif i in contigs_2:
                labels.append('Significant')
            elif i in contigs_3:
                labels.append('Significant')
            #Otherwise..
            else:
                #Save 'Non-significant' to the list
                labels.append('Non-significant')
    #If only 1 result is given..
    else:
        #Iterate over all the contigs
        for i in all_contigs:
            #Check if it is in the significant contig file..
            if i in contigs :
                #If so, save 'Significant' to the list
                labels.append('Significant')
            #Otherwise..
            else:
                #Save 'Non-significant' to the list.
                labels.append('Non-significant')
    #Append the list to the dataframe
    df['target'] = labels
    #Set the labels as index
    df.set_index('target', inplace = True)
    #Print the top rows of the generated DataFrame
    print(df.head(5))

    #Initate a PCA with 2 components
    pca = PCA(n_components=2)
    #Perform the PCA
    PCA_proteins = pca.fit_transform(df)
    #Initiate an empty list
    explained_variance = []
    #Iterate over the explained variance
    for i in pca.explained_variance_ratio_:
        #Convert the explained variance to percentages
        a = float(i) * 100
        #Save it to the list with 2 decimals
        explained_variance.append("%.2f" % a)
    #Create a PCA dataframe
    principal_df = pd.DataFrame(data = PCA_proteins, columns = ['principal component 1', 'principal component 2'])
    #Initate a figure
    plt.figure(figsize=(10,10))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=14)
    #Set the explained variance as x and y label of the plot
    plt.xlabel('Principal Component - 1 ({}%)'.format(explained_variance[0]),fontsize=20)
    plt.ylabel('Principal Component - 2 ({}%)'.format(explained_variance[1]),fontsize=20)
    #Set Non-significant as black and significant as red
    targets = ['Non-significant', 'Significant']
    colors = ['k', 'r']
    #Iterate over the significance..
    for target, color in zip(targets, colors):
        #Select the corresponding contigs
        indicesToKeep = df.index == target
        #Create the scatter plot
        plt.scatter(principal_df.loc[indicesToKeep, 'principal component 1']
                , principal_df.loc[indicesToKeep, 'principal component 2'], c = color)
    #Save and show the figure
    plt.savefig(save_name[:-4] + '_pertype.png', bbox_inches='tight')
    plt.show()


print_most_abundant_organisms('lineage_S.csv', only_out = False)
plot_PCA('KEGG/CD_KEGG_annotated.csv', 'KEGG/CD_KEGG_PCA.png')
scatter_plot_significance('bowtie/C_bowtie_clr.csv', 'bowtie/Top_values_C_I_O.csv', 'bowtie/C_PCA_significance.png', 'bowtie/Top_values_C_I_R.csv', 'bowtie/Top_values_C_O_R.csv')
