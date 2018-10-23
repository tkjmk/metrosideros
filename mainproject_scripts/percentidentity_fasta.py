#!/usr/bin/python3

"""Program that replaces mismatched sequences with Ns based on the number of matches in a specified amount of sequence.""" 
__author__ = 'Tane Kafle (tk916@ic.ac.uk)'
__version__ = '0.0.1'
__date__ = 'March 2018'

from Bio import AlignIO 
from Bio.Align import MultipleSeqAlignment
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(usage = "python3 percentidentity_fasta.py -i STR -p STR\n\nPython 3 required!\n\nThis program takes in fasta files for each contig that have individuals (from 2 populations) in the same order in each fasta and each fasta is of the same length. It will compare each sequence for matching bases (that are not Ns; note that gaps are not recognised, please convert to Ns) and print the percentage identity as a decimal to terminal. The script outputs the percent identity within the 2 populations and between them. \n\nPlease ensure: \n- that the directory should only contain the fasta files (can have directories, but no other files)\n- all sequences within a fasta files should be the same length, without gaps (-) and with individuals in the same order (and the same number of individuals; add a sequence of Ns if you do not have sequence data for that individual).\n\nAuthor: Tane Kafle\n\n")
parser.add_argument("-i", "--input_directory", help="input directory. Directory should only contain the fasta files to be read in.", required=True)
parser.add_argument("-p", "--population_of_individuals", help="Enter a list seperate by commas of the population of each individual in the exact order the individuals appear in the fasta", type=str, required=True)



args = parser.parse_args()

def findfiles(directory):
    objects = os.listdir(directory)  # find all objects in a dir

    files = []
    for i in objects:  # check if very object in the folder ...
        if isFile(directory + i):  # ... is a file.
            files.append(i)  # if yes, append it.
    return files

def isFile(object):
    try:
        os.listdir(object)  # tries to get the objects inside of this object
        return False  # if it worked, it's a folder
    except Exception:  # if not, it's a file
        return True

fastadir = args.input_directory
#fastadir = "/home/tane/Project_CMEE/Sandbox/testing/"
fastafiles = findfiles(fastadir) #Creating a list with the fastas.

#These are dictionaries for:
aln = {} #the read in fasta.
seqLen = {} #length of sequences
nSeq = {} #number of sequences in fasta
filt_aln = {} 
aln_a = {} #these are the sequences.
aln_r = {} #this is length of the sequence stored in a range.

for i in fastafiles: #Reading in fasta files and other stats into dictionary.
    aln[i] = AlignIO.read(fastadir + i, "fasta")
    seqLen[i] = aln[i].get_alignment_length() #Length of aligned sequences in FASTA
    nSeq[i] = len(aln) #Number of sequences from FASTA
    filt_aln[i] = MultipleSeqAlignment([])
    aln_a[i] = np.char.array([list(rec) for rec in aln[i]]) #Array with alignments
    aln_r[i] = range(seqLen[i]) #Range of alignments 0:seqLen-1

def get_unique_list(lst):
    '''
    Function that will return a list with the unique values from that list.
    '''
    if isinstance(lst,list):
        return list(set(lst))

pop = args.population_of_individuals
pop = [x.strip() for x in pop.split(',')]
pop_unq = get_unique_list(pop)


def triangular_number(n):
    '''
    Returns triangular number i.e. the factorial of a number, but addition rather than multipliplication). Used in function consecutive_pairs.
    '''
    return n * (n + 1) // 2

def consecutive_pairs(int):
    '''
    Returns consecutive pairs of numbers i.e. for comparison, will give it to use as indices for use in python e.g. if you give the number 3, you will get returned [[0,1], [0,2], [1,2]]
    ''' 
    int = int - 1
    pairs = []
    firstitem = 0
    seconditem = 0
    for i in range(triangular_number(int)):
        seconditem += 1
        pairs.append([firstitem, seconditem])
        if seconditem == int:
            firstitem += 1
            seconditem = firstitem
    return pairs


#Lists that store all the records of the probability of similar base per site.
withinpop1 = []
withinpop2 = []
betweenpop = []

for fasta in aln_a: #Goes through each fasta from the directory
    comparisons = consecutive_pairs(len(aln_a[fasta])) #Generate index pairs to compare each of the sequences from this fasta.
    for pair in comparisons:
        #Assigning the sequences to variables
        seq1 = aln_a[fasta][pair[0]]
        seq2 = aln_a[fasta][pair[1]]
        # Variables to count the:
        seq1N = 0 #number of Ns in sequence 1.
        seq2N = 0 #number of Ns in sequence 2.
        samebase = 0 #number of bases that are the same, that are not an N.
        comparisonsmade = 0 #number of total bases that were compared i.e. sites where both sequences did not have an N.
        for base in range(len(seq1)): #Looping through sequences.
            if seq1[base] == 'N': #Counting Ns in sequence 1.
                seq1N += 1
            if seq2[base] == 'N': #Counting Ns in sequence 2.
                seq2N += 1
            if seq1[base] != 'N' and seq2[base] != 'N': #Counting where bases are the same and neither sequence has an N.
                comparisonsmade += 1
                if seq1[base] == seq2[base]:
                    samebase += 1
        
    
        
        if seq1N/len(seq1) < 0.5 and seq2N/len(seq2) < 0.5: #Ensuring that their are less than 50% N in both sequences and then adding the percent identity to the appropriate list: between populations and within populations 1 or 2.
            if pop[pair[0]] != pop[pair[1]]:
                betweenpop.append(samebase/comparisonsmade)
            else:
                if pop[pair[0]] == pop_unq[0]:
                    withinpop1.append(samebase/comparisonsmade)
                if pop[pair[0]] == pop_unq[1]:
                    withinpop2.append(samebase/comparisonsmade)
#Calculating mean percentage identity of bases (given as decimal)
between_pop_mean = sum(betweenpop)/len(betweenpop)
within_pop1_mean = sum(withinpop1)/len(withinpop1)
within_pop2_mean = sum(withinpop2)/len(withinpop2)

#Printing output to terminal.
print("\n\nWithin Population 1 (" + pop_unq[0] + ") percent identity: " + '%.10f' % within_pop1_mean + " (difference: " + str('{:.3e}'.format(1 - within_pop1_mean)) + ")\n\nWithin Population 2 (" +  pop_unq[1] + ") percent identity: " + '%.10f' % within_pop2_mean + " (difference: " + str('{:.3e}'.format(1 - within_pop2_mean)) + ")\n\nBetween Populations percent identity: " + '%.10f' % between_pop_mean + " (difference: " + str('{:.3e}'.format(1 - between_pop_mean)) + ")\n\n")