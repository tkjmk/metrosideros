#!/usr/bin/python3
__author__ = 'Tane Kafle (tk916@ic.ac.uk)'
__version__ = '0.0.2'
__date__ = 'June 2018'
import os, sys, argparse
parser = argparse.ArgumentParser(usage = "python3 estimate_dxy.py -i STR -o STR -p1 INT -p2 INT [-p1n STR -p2n STR -h]\n\nPython 3 required!!\n\nThis program calculates Dxy as according to Irwin et al. (2016), Molecular Ecology 26(18), pp4488-4507. This between-group nucleotide differentiation (Dxy)  is determed as p1(1-p2) + p2(1-p1), where p1 is the grequecy of a given allele in the first population and p2 is the frequency of that allele in the second population and ranges from 0 to 1. This was originally made for both gppfst and minotaur. It also calculates pi for both populations and the difference between them. A .dxy output file is also printed to a desired directory.\n\nAuthor: Tane Kafle\n\n")
parser.add_argument("-i", "--input_vcf", required=True, help="input VCF file pathway.")
parser.add_argument("-o", "--output_dxy", required=True, help="output DXY directory pathway. Please include final forward slash.")
parser.add_argument("-p1", "--pop1", required=True, help="Number of individuals/samples in the VCF from population 1.", type=int)
parser.add_argument("-p2", "--pop2", required=True, help="Number of individuals/samples in the VCF from population 2.", type=int)
parser.add_argument("-p1n", "--pop1_name", help="Population 1 name/identifier.", default="pop1", type=str)
parser.add_argument("-p2n", "--pop2_name", help="Population 2 name/identifier.", default="pop2", type=str)

#Version 0.0.2 on 11th June 2018 has pi calculations added.

args = parser.parse_args()

#Getting the number of individuals per population, so we can get the location via indices of what sample data on VCF is from which population.
pop1len = args.pop1
pop2len = args.pop2

p1vl = 9 #population 1 VCF indices lower bound.
p1vh = 9 + pop1len #population 1 VCF indices high i.e. what index is that of the last sample of Pop 1.
p2vl = 9 + pop1len #population 2 VCF indices low
p2vh = 9 + pop1len + pop1len #population 1 VCF indices high

vcfin = open(args.input_vcf, 'r')

dxyoutname = args.input_vcf.split('/')[-1] #extracting file name
dxyoutname = dxyoutname.split('.')
dxyoutname.pop() #removing the ending of the file name i.e. .vcf.
dxyoutname = '.'.join(dxyoutname)
dxyoutname = args.output_dxy + dxyoutname + ".dxy"

dxyout = open(dxyoutname, 'w')
dxyout.write("Contig\tReference\tAlternative\tDxy\t" + args.pop1_name + "_Indv\t" + args.pop2_name + "_Indv\t" + args.pop1_name + "_RefAllele_Count\t" + args.pop2_name + "_RefAllele_Count\t" + args.pop1_name + "_AltAllele_Count\t" + args.pop2_name + "_AltAlle_Count\t" + "Pi_" + args.pop1_name + "\tPi_" + args.pop2_name + "\tDelta_Pi_" + args.pop2_name + args.pop1_name) #Creating header of stats file.

for line in vcfin:
    if line.startswith('#'):
        continue
    else:
        contig_nam = line.split('\t')[0]
        ref_allele = line.split('\t')[3]
        alt_allele = line.split('\t')[4]
     
        #Extracting allele information for population 1
        pop1samples = line.split('\t')[p1vl:p1vh] #These are the pop1 sample information, where each sample is an item in this list.
        pop1sampleinfo = [component for components in pop1samples for component in components.split(':')] #Splitting the sample information into individual components.
        pop1alleles = [genos for genos in pop1sampleinfo if '/' in genos] #Only taking the allele information (if it is reference or alternative i.e. 0 or 1)
        pop1alleles_sep = [allele for alleles in pop1alleles for allele in alleles.split('/')] #Putting all alleles into a single pool to be counted.
        #Extracting allele information for population 2
        pop2samples = line.split('\t')[p2vl:p2vh]
        pop2sampleinfo = [component for components in pop2samples for component in components.split(':')]
        pop2alleles = [genos for genos in pop2sampleinfo if '/' in genos]
        pop2alleles_sep = [allele for alleles in pop2alleles for allele in alleles.split('/')]

        #Allele frequency counts.
        pop1_ref_freq = pop1alleles_sep.count('0') / (pop1alleles_sep.count('0') + pop1alleles_sep.count('1'))
        pop1_alt_freq = pop1alleles_sep.count('1') / (pop1alleles_sep.count('0') + pop1alleles_sep.count('1'))
        pop2_ref_freq = pop2alleles_sep.count('0') / (pop2alleles_sep.count('0') + pop2alleles_sep.count('1'))
        pop2_alt_freq = pop2alleles_sep.count('1') / (pop2alleles_sep.count('0') + pop2alleles_sep.count('1'))

        dxy_ref = (pop1_ref_freq * (1 - pop2_ref_freq)) + (pop2_ref_freq * (1 - pop1_ref_freq ))
        dxy_alt = (pop1_alt_freq * (1 - pop2_alt_freq)) + (pop2_alt_freq * (1 - pop1_alt_freq ))

        #Count number of individuals without data
        pop1_inv = pop1len - (pop1alleles_sep.count('.') / 2)
        pop2_inv = pop2len - (pop2alleles_sep.count('.') / 2)

        #Estimating pi
        pop1_n = pop1alleles_sep.count("0") + pop1alleles_sep.count("1") #total_alleles_count
        pop1_k = pop1alleles_sep.count("0") #Count of reference allele (you get same results if you use alternative allele - as all are biallelic)
        pop1_pi = (2 *pop1_k*(pop1_n - pop1_k))/(pop1_n*(pop1_n -1))

        pop2_n = pop2alleles_sep.count("0") + pop2alleles_sep.count("1") #total_alleles_count
        pop2_k = pop2alleles_sep.count("0")
        pop2_pi = (2 *pop2_k*(pop2_n - pop2_k))/(pop2_n*(pop2_n -1))

        delt_pop2pop1 = pop2_pi - pop1_pi

        if abs(round(dxy_ref,3) - round(dxy_alt,3)) >= 0.1: #The dxy using ref or alt allele should give the same result if only biallelic SNPs are used. Rounded to 3sf and difference greater than 0.1, because sometimes there are slight different results from calculated and rounded values.
            print("There are multiallelic SNPs present. Please remove them for this analysis.")
            sys.exit()
        outputline = [contig_nam, ref_allele, alt_allele, str(dxy_ref), str(int(pop1_inv)), str(int(pop2_inv)), str(pop1alleles_sep.count('0')), str(pop2alleles_sep.count('0')), str(pop1alleles_sep.count('1')), str(pop2alleles_sep.count('1')), str(pop1_pi), str(pop2_pi), str(delt_pop2pop1)]
        dxyout.write("\n" + "\t".join(outputline))
