import subprocess
import argparse

ap = argparse.ArgumentParser(description = 'This script converts the Concoct output to fasta files ')
ap.add_argument('-f', '--file', required = True, help = 'The Concoct output coverage file')
ap.add_argument('-o', '--fasta', required = True, help = 'The directory containing the empty fasta files to save the sequences to')

args, leftovers = ap.parse_known_args()
#Open the concoct coverage file
file = open(args.file, 'r')
#Read per line
file = file.read().splitlines()
#Iterate over the lines
for f in file:
    #Split on tab
    f = f.split('\t')
    #Open a fasta file in append modus
    with open(args.fasta + f[1] + ".fasta", "a+") as fasta:
        #Append the fasta sequence to the created fasta file 
        subprocess.call(["grep", f[0] + " ", "-A", "1", "final.contigs.fa"], stdout=fasta)
