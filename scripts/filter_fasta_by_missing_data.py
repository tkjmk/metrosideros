#!/usr/bin/python3
"""Program that filters fastas by the amount of missing data in the sequence and the number of individuals that have more than a set amount of missing data. """ 
__author__ = 'Tane Kafle (tk916@ic.ac.uk)'
__version__ = '0.0.3'
__date__ = 'May 2018'


"""
Changelog:
----------
v0.0.3 - May 2025
    - Refactored to include main() function and argparse handling

v0.0.2 - May 2018
    - Added option to output files in phylip format,  it can be controlled using the <output_format> flag
    - Log file is now output optionally with flag <log_file>

v0.0.1 - May 2018
    - Initial version
"""


import argparse, os, sys
from os import listdir
from os.path import isfile, join
from Bio import AlignIO 
from Bio.Align import MultipleSeqAlignment
import numpy as np

def restricted_float(x):
    '''
    Defining a limit for floats between 0 and 1 inclusive. To ensure that threshold input for <minimum_nonN_seq_threshold> in argparse is between 0 and 1.
    '''
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

def parse_arguments():
    parser = argparse.ArgumentParser(usage = 'python3 filter_fasta_by_missing_data.py -i STR -o STR -p1 STR -p2 STR [-m INT -s FLOAT[0:1] -f STR -l BOOL -h]\n\nPython 3 required!!\n\nThis program takes in fasta files from <input_directory>, and if they pass the filter, they are output into <output_directory> as a fasta or phylip file depending on <output_format>. The main filter is the <minimum_nonN_seq_threshold>, which is the proportion (given as a float between 0 and 1) of the sequence you want to be actual data (as opposed to missing data [N]). The second part of the filter is the <minimum_indv_per_population_threshold>, where an integer is given and at least that many individuals from each population have to pass the <minimum_nonN_seq_threshold> for the fasta to pass the filter. Populations are defined using the <pop1> and <pop2> flags, listing individuals names from fasta seperated by commas. When entering input directory, a <file_extension> can also be stated, so that only files with that extension are input into program. A log file called "filter_fasta_by_missing_data.log" will also be output in the output directory if <log_file> is set to True. \n\nPlease ensure the output directory has been created! \n\nPlease ensure that <input_directory> and <output_directory> pathways end with a forward slash.\n\nUse -h to show full help screen.\n\nAuthor: Tane Kafle\n\n')
    parser.add_argument("-i", "--input_directory", required=True, help="Input directory that contains the fasta files. Please end with '/' at the end of pathway.")
    parser.add_argument("-o", "--output_directory", required=True, help="The directory where you would like to output the fasta files that pass the filter. A log file will also be output here. Please end with '/' at the end of pathway.")
    parser.add_argument("-f", "--file_extension", help="Put file ending e.g. '.fas' to ensure only fasta files are selected for analysis.", default=".fas", type=str)
    parser.add_argument("-s", "--minimum_nonN_seq_threshold", help="Please enter the minimum threshold of sequence that is required (i.e. the amount of sequence that is not Ns). For example 0.6 means that sequences with more than 40% Ns in sequences will not pass.", default = 0.8, type=restricted_float)
    parser.add_argument("-p1", "--pop1", required=True, help="A comma seperated list (no spaces) of individuals from population 1, as written in FASTA.", type=str)
    parser.add_argument("-p2", "--pop2", required=True, help="A comma seperated list (no spaces) of individuals from population 2, as written in FASTA.", type=str)
    parser.add_argument("-m", "--minimum_indv_per_population_threshold", help="Argument to copy (c) or move (m) the files from input directory to output directory.", default = 2 , type=int)
    parser.add_argument("-f", "--output_format",  help="Output format, either fasta (default) or phylip.", required=False, default="fasta", type=str)
    parser.add_argument("-l", "--log_file",  help="Log file output or not in output directory. Please input boolean value (True or False).", required=False, default=True , choices = [True, False], type=bool)
    return parser.parse_args()

def main():
    args = parse_arguments()

    inpath = args.input_directory
    outpath = args.output_directory

    logf = args.log_file
    outfmt = args.output_format

    # Building a list of filenames within the directory (with or without specific file extension depending on <file_extension> flag).

    fastafiles = [f for f in listdir(inpath) if f.endswith(args.file_extension) and isfile(join(inpath, f))]
    if len(fastafiles) == 0:
        raise ValueError("There are no files with a " + args.file_extension + " extension in this directory, please select <input_directory> that contains files with that extension or <file_extension> found within this directory.")
    if logf == True:
        log = open(outpath + "filter_fasta_by_missing_data.log",'w') #Creating log file.
        log.write("#MIN_SEQ_THRESHOLD = " + str(args.minimum_nonN_seq_threshold)  +"\n#MIN_INDV_PER_POP_THRESHOLD = " + str(args.minimum_indv_per_population_threshold) + "\nfasta, result, pop1_passed, pop2_passed\n") #Inputting first bit of information.


    Nthreshold = 1 - args.minimum_nonN_seq_threshold #The maximum proportion of Ns allowed in a sequence.
    pop1 = args.pop1 #Converting the individuals names given into a list.
    pop1 = [x.strip() for x in pop1.split(',')] 
    pop2 = args.pop2
    pop2 = [x.strip() for x in pop2.split(',')]
    min_pop_thresh = args.minimum_indv_per_population_threshold


    #These are dictionaries for:
    aln = {} #the read in fasta.
    #seqLen = {} #length of sequences
    #nSeq = {} #number of sequences in fasta
    #filt_aln = {} 
    #aln_a = {} #these are the sequences.
    #aln_r = {} #this is length of the sequence stored in a range.
    #aln_id = {} #this is to store the IDs of each individual in sequence.
    aln_dict = {} #this will be a dictionary of dictionaries where the keys for each fasta will be the individual name, and the item will be its sequence as a list.
    for i in fastafiles: #Reading in fasta files and other stats into dictionary.
        aln[i] = AlignIO.read(inpath + i, "fasta")
    # seqLen[i] = aln[i].get_alignment_length() #Length of aligned sequences in FASTA
    # nSeq[i] = len(aln[i]) #Number of sequences from FASTA
    # filt_aln[i] = MultipleSeqAlignment([])
    # aln_a[i] = np.char.array([list(rec) for rec in aln[i]]) #Array with alignments
    # aln_r[i] = range(seqLen[i]) #Range of alignments 0:seqLen-1
    # aln_id[i] = np.char.array([rec.id for rec in aln[i]])
        aln_dict[i] = {rec.id:list(rec) for rec in aln[i]} #dictionary comprehension storing the name of sequence with the sequence.

    fas_passed = [] #List to store the names of fasta's that passed filter.
    for fasta in aln_dict: #Looping through fastas.
        p1_counter = 0 #These are counters that go up everytime an individual passes the missing data filter.
        p2_counter = 0
        for seq in aln_dict[fasta]: #Looping through sequences within the fasta.
            number_of_Ns = aln_dict[fasta][seq].count("N") #Counting the number of Ns in the sequence.
            fraction_of_Ns = number_of_Ns / len(aln_dict[fasta][seq]) #Dividing this by the length of the sequence to get the proportion of missing data.
            if fraction_of_Ns < Nthreshold: #If the amount of missing data is less than the threshold than the sequence has passed, and depending on the population, the counter goes up.
                if seq in pop1:
                    p1_counter += 1
                elif seq in pop2:
                    p2_counter += 1
                else:
                    continue
        if p1_counter >= min_pop_thresh and p2_counter >= min_pop_thresh: #This then ensures that the minimum number of individuals from each population have passed the missing data filter.
            fas_passed.append(fasta)
            if logf == True:
                log.write(fasta + ", Passed, " + str(p1_counter) + ", " + str(p2_counter) + "\n") #Entering information to log file.
            if outfmt == "fasta":
                outfile = open(outpath + fasta, 'w') #Writing them to the output directory, if they have passed.
                AlignIO.write(aln[fasta], outfile, "fasta")
                outfile.close()
            elif outfmt == "phylip":
                fastasplit = fasta.rsplit('.', 1)
                phylip = fastasplit[0] + ".phy"
                outfile = open(outpath + phylip, 'w') #Writing them to the output directory, if they have passed.
                AlignIO.write(aln[fasta], outfile, "phylip-sequential")
                outfile.close()
            else:  
                print("Output format must be either fasta or phylip\n")
                sys.exit()
        else:
            if logf == True:
                log.write(fasta + ", Failed, " + str(p1_counter) + ", " + str(p2_counter) + "\n") #Entering information to log file.
            continue

    if logf == True:
        log.close()

    print("Finished, " + str(len(fas_passed)) + " of " + str(len(fastafiles)) + " passed the filters. These have been output to " + outpath + "\n") #Printing a nice output


if __name__ == "__main__":
    main()
