#!/usr/bin/python2 -W ignore
__author__ = 'Tane Kafle (tk916@ic.ac.uk)'
__version__ = '0.0.2'
__date__ = 'May 2018'


"""
Changelog:
----------
v0.0.2 - May 2025
    - Refactored to include main() function and argparse handling

v0.0.1 - May 2018
    - Initial version
"""


import os, sys, argparse

def parse_arguments():
	parser = argparse.ArgumentParser(usage = "python2 leastmissingdataperSNP_twopop.py -i STR -o STR -p1 INT -p2 INT [-p1n STR -p2n STR -t INT -h]\n\nPython 2 required!!\n\nThis program filters VCF files with SNP data. It first ensures that there are at least <threshold> (default = 2) individuals with data for that SNP from each population. If there are, it then selects the SNP with the least missing data over all individuals from both populations to be kept. If two SNPs have the same amount of missing data, the first one will remain. A tab seperated .stat output will be put with contig and the amount of data from each population for that SNP. Please ensure pop1 samples are the first samples in the VCF, followed by the pop2 samples. This script used Adam Ciezarek's LeastMissingDataPerSNP.py script as a template.\n\nAuthor: Tane Kafle\n\n")
	parser.add_argument("-i", "--input_vcf", required=True, help="input VCF file pathway with desired VCF name at end.")
	parser.add_argument("-o", "--output_vcf", required=True, help="output VCF file pathway with desired VCF name at end.")
	parser.add_argument("-p1", "--pop1", required=True, help="Number of individuals/samples in the VCF from population 1.", type=int)
	parser.add_argument("-p2", "--pop2", required=True, help="Number of individuals/samples in the VCF from population 2.", type=int)
	parser.add_argument("-p1n", "--pop1_name", help="Population 1 name/identifier.", default="pop1", type=str)
	parser.add_argument("-p2n", "--pop2_name", help="Population 2 name/identifier.", default="pop2", type=str)
	parser.add_argument("-t", "--threshold", help = "The minimum number of individuals from each population", default=2, type=int)
	return parser.parse_args()

def main():
	args = parse_arguments()

	#Getting the number of individuals per population, so we can get the location via indices of what sample data on VCF is from which population.
	pop1len = args.pop1
	pop2len = args.pop2

	p1vl = 9 #population 1 VCF indices lower bound.
	p1vh = 9 + pop1len #population 1 VCF indices high i.e. what index is that of the last sample of Pop 1.
	p2vl = 9 + pop1len #population 2 VCF indices low
	p2vh = 9 + pop1len + pop1len #population 1 VCF indices high

	vcfin = open(args.input_vcf, 'r')
	vcfout = open(args.output_vcf, 'w')

	headers = {}

	# For loop to add the SNPs with the least missing data per contig to the headers dictionary.
	for line in vcfin:
		if line.startswith('#'):
			vcfout.write(line)
		else:
			#Counting for missing data per sample first and if they both don't have at least 2, then this line is ignored.
			pop1samples = line.split('\t')[p1vl:p1vh]
			pop2samples = line.split('\t')[p2vl:p2vh]
			pop1missingdatacount = pop1samples.count("./.:.:.:.")+ pop1samples.count("./.:.:.:.\n")
			pop2missingdatacount = pop2samples.count("./.:.:.:.")+ pop2samples.count("./.:.:.:.\n")
			if (pop1len - pop1missingdatacount < args.threshold) or (pop2len - pop2missingdatacount < args.threshold):
				continue
			else:
				if line.split('\t')[0] not in headers: #Checks if the contig has been looked at and if it hasn't, adds it to headers and calculates the amount of missing data.
					headers[line.split('\t')[0]] = line+"____"+str(line.count("./.:.:.:."))
				else: #If the contig has previously been looked at, it counts the missing data and then adds the one with less missing data to header, then the cycle starts again.
					a = line.count("./.:.:.:.")
					b = headers[line.split('\t')[0]].split('____')[1]
					if b < a:
						headers[line.split('\t')[0]] = line+"____"+str(line.count("./.:.:.:.")) #Headers has the vcfline followed by some underscores and then the amount of missing data in that line.				
					
	for i in headers:
		vcfout.write(str(headers[i].split('____')[0]))

	## This is to create the first part of information required for gppfst, the chromosome/loci + the amount of missing data per species (specific to my species ME_NE and ME_SC) - would work for someone else if they have 10 species of both following each other.
	#This part is editing the VCF input name to create a VCF output name with .stats instead of .vcf.
	vcfstatoutname = args.output_vcf.split('.')
	vcfstatoutname.pop() #removing the ending of the file name i.e. .vcf.
	vcfstatoutname = '.'.join(vcfstatoutname)
	vcfstatoutname = vcfstatoutname + ".stats"

	vcfstatout = open(vcfstatoutname, 'w')
	vcfstatout.write("Contig\tPos\t" + args.pop1_name + ".samples\t" + args.pop2_name + ".samples\n") #Creating header of stats file.

	for i in headers:
		VCFline = headers[i].split('____')[0] #Getting the VCF input line.
		VCFline = VCFline.split('\t') #Splitting individual line into individual sections.
		POP1 = VCFline[9:(9+pop1len)] #These are the sample indices for the 10 POP1 samples, next line counts the missing data for these samples. The same is then done fro ME_SC.
		POP1data = str(pop1len - (POP1.count("./.:.:.:.")+ POP1.count("./.:.:.:.\n")))
		POP2 = VCFline[(9+pop1len):(9+pop1len+pop2len)]
		POP2data = str(pop2len - (POP2.count("./.:.:.:.") + POP2.count("./.:.:.:.\n")))
		position = VCFline[1]
		statslines = [i, position, POP1data, POP2data] #Creating list with the stats for that contig: which contig it is and missing data.
		vcfstatout.write('\t'.join(statslines) + '\n') #Writing out to file.
		
		
	print("\nNumber of SNPs: " + str(len(headers))) # Counting number of SNPs in VCF (there are one per contig, so also number of contigs)

if __name__ == "__main__":
    main()