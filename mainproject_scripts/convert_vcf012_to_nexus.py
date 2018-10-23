#!/usr/bin/python3 -W ignore
__author__ = 'Tane Kafle (tk916@ic.ac.uk)'
__version__ = '0.0.1'
__date__ = 'May 2018'
import argparse
parser = argparse.ArgumentParser(usage = "python3 convert_vcf012_to_nexus.py -i STR -o STR [-h]\n\nPython 3 required!!\n\nThis program filters converts an 012 file generated from a VCF using vcftools --012 flag to a nexus file. The 012 flag provides a matrix for each individual from each SNP, represented by a 0, 1 or 2; the number represents the number of alternative alleles (as opposed to reference) present in that individual for that SNP. Please just input the file ending with .012 into the <input_012> argument to get your output nexus file. The .012.indv file is also required to be present in the same directory (the individual file should have the exact same file name, just with the .indv extension added). You can optionally also change the location and name of the output nexus file using the <output_nexus> argument.\n\nAuthor: Tane Kafle\n\n")
parser.add_argument("-i", "--input_012", required=True, help="Please enter the location of the .012 file outputted from the vcftools --012 flag, and ensure the .012.indv file is also present in this directory.")
parser.add_argument("-o", "--output_nexus", help="Please enter the output name for the nexus file. If no argument is passed, the file name will be that of the 012 file with .nex at the end in the same directory.", default=None)

args = parser.parse_args()

file012 = args.input_012

in012 = open(file012, "r" ) #File with the 012 matrix generated from a VCF using VCFtools.
in012names = open(file012 + ".indv", "r") #File with individual names from VCF.

individuals = [] #Reading in individual names.
for indiv in in012names:
    individuals.append(indiv)

individuals = [i.strip("\n") for i in individuals] #Removing new line characters from individual names.


sequences = [] #Reading in the sequences to a list.
for lines in in012:
    sequences.append(lines)


sequences = [i.split("\t") for i in sequences] #Splitting out 012 SNP sequence characters into individual items of a list.

for i in range(0,len(sequences)):
    sequences[i][len(sequences[i]) - 1] = sequences[i][len(sequences[i]) - 1].strip("\n") #Removing the new line character of the last item of each list.
    sequences[i] = [j.replace("-1", "-") for j in sequences[i]] #Replacing "-1" with "-" to represent missing data.
    sequences[i].pop(0) #Removing the first item in the list, which is not part of the 012 sequence (it is the individual name).

numbofseq = len(sequences) #Number of sequences from the 012 file.
seqlen = len(sequences[0]) #Length of sequences from the 012 file.

#Writing out all the information into the nexus file.
if args.output_nexus == None:
    outnex = open(file012.replace(".012", "") + ".nex", "w")
else:
    outnex = open(args.output_nexus, "w")

outnex.write("#NEXUS\n\nBegin data;\n\tDimensions ntax=" + str(numbofseq) + " nchar=" + str(seqlen) + ';\n\tFormat datatype=integerdata symbols="012" gap=-;\n\tMatrix\n')

for i in range(0,len(sequences)):
    outnex.write(individuals[i] + "  \t\t" + "".join(sequences[i]) + "\n")

outnex.write("\t;\nEnd;\n")
outnex.close()
in012names.close()
in012.close()