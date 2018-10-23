#!/usr/bin/python2

"""Program that replaces mismatched sequences with Ns based on the number of matches in a specified amount of sequence.""" 
__author__ = 'Tane Kafle (tk916@ic.ac.uk)'
__version__ = '0.0.1'
__date__ = 'Feb 2018'

from Bio import AlignIO 
#import glob
from Bio.Align import MultipleSeqAlignment
#from Bio.SeqRecord import SeqRecord
import numpy as np
import argparse

parser = argparse.ArgumentParser(usage = "python2 remove_mistmatched_seq.py -i STR [-o STR -O STR -w INT -t INT]\n\nPython 2 required!\n\nThis program filters pairwise fasta alignments for areas that have mismatches. It reads in both aligned sequences from a single fasta file and goes through the sequence by the <window_size> and if there are less matches in this window than the <matches_in_window> threshold, then that window will be replaced by 'N's. The final output will be a fasta file with the two aligned sequences and any 'mismatched' areas replaced by 'N's.\n\nNote that gaps (dashes) in either sequence are considered as a match. The two sequences must also be of same length.\n\nAuthor: Tane Kafle\n\n")
parser.add_argument("-i", "--input_file_name", help="input file name. Can be a pathway.", required=True)
parser.add_argument("-o", "--output_directory", help="path to output folder (must exist)")
parser.add_argument("-O", "--output_file_name", help="output file name, please also put desired file extension")
parser.add_argument("-w", "--window_size", help="number of bases to be analysed at a time for matches", type=int, default=10)
parser.add_argument("-t", "--matches_in_window", help="threshold value, in terms of number of bases that must match between the two sequences for the window to be output in the filtered alignments.", type=int, default=7)

#OUTPUT FILE NAME IF NONE PROVIDED
args = parser.parse_args()

if args.output_file_name is not None:
    outputfilename = args.output_file_name
else:
    inputfile = args.input_file_name
    inpname = list(reversed(inputfile.rsplit('/',1)))[0]
    inpname = inpname.split(".")
    fileext = inpname.pop((len(inpname)-1)) #remove and stores file extension name
    inpname.append("mismatch.filt")
    inpname.append(fileext)
    outputfilename = ".".join(inpname)


#FASTAPATH
fastafile = args.input_file_name
aln = AlignIO.read(fastafile, "fasta") #Read in FASTA
seqLen = aln.get_alignment_length() #Length of aligned sequences in FASTA
nSeq = len(aln) #Number of sequences from FASTA
if nSeq != 2:
    raise ValueError('Please input fasta file with exactly 2 aligned sequences')
if len(aln[0]) != len(aln[1]):
    raise ValueError('Please ensure the two sequences in the fasta file are of teh same length')
filt_aln = MultipleSeqAlignment([])
aln_a = np.array([list(rec) for rec in aln], np.character) #Array with alignments
aln_r = range(seqLen) #Range of alignments 0:seqLen-1
step = args.window_size
threshold = args.matches_in_window
if threshold > step:
    raise ValueError('Please input <matches_in_window> integer less than or equal to <window_size> integer')

if step > seqLen:
    raise ValueError('Please ensure <window_size> is less than the length of the sequence')

aln_r_st = range(0, seqLen, int(step)) #Step range
cons_ar = np.array([[],[]]) #Concensus Array
bases = ["a","t","c","g"]

for i in aln_r_st: #Represents every 10 bases of sequence
    chunk = aln_a[:,i:i+step] #Even if last chunk is less than 10, Python won't count indices that aren't there.
    counter = 0 #Initialising counter, bases that are the same in the chunk or if either include a dash, lead to a 1 being added to the counter. Base mismatches will not add anything.
    for base in range(len(chunk[0])):
        if chunk[0,base] == '-' or chunk[1,base] == '-':
            counter += 1
        elif chunk[0,base] == chunk[1,base]:
            counter += 1
    if (counter * (float(step)/len(chunk[0]))) >= threshold: #This division of ste/len(chunk[0]) is to ensure the final alignment (which for example could be 5 bases and therefore even a perfect match will not create a counter greater than 7, assuming that is the threshold value)
        cons_ar = np.concatenate((cons_ar,chunk), axis=1)
    else:
        N_string = ["N" for x in range(len(chunk[0]))]
        N_array = np.array([N_string, N_string])
        cons_ar = np.concatenate((cons_ar, N_array), axis=1)

cons_list = list()
cons_list = cons_ar.tolist()

for i in range(len(aln)):
    filt_aln.add_sequence(descriptor = aln[i].id,sequence = str(cons_list[i]).translate(None, ",' []"))

if args.output_directory is not None:
    outputdirectory = args.output_directory
    outputfile = outputdirectory + outputfilename
else:
    outputfile = outputfilename

count = AlignIO.write(filt_aln, outputfile, "fasta")

# base = base
# base = -
# base != base

# +1
# +1
# +0



#   -   -   -   -   -   -   A   -   -   -                               #sequence 1
#   A   T   G   C   A   C   A   G   T   -                               #sequence 2
#   A   T   G   C   A   C   A   G   T   -   #Score here would be 10     #concensus sequence

#   A   T   G   C   A   C   A   G   T   A
#   G   A   C   A   T   C   G   A   C   G
#   N   N   N   N   N   C   N   N   N   N   #Score here would be 1

#   -   -   -   -   -   -   A   -   -   -
#   -   -   -   G   -   -   A   C   -   -  
#   -   -   -   G   -   -   A   C   -   -   #Score here would be 10

#   -   T   -   -   C   G   A   -   -   -
#   -   C   -   G   A   T   A   -   C   -  
#   -   N   -   G   N   N   A   -   C   -   #Score here would be 7   