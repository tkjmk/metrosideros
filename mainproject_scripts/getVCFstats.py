#!/usr/bin/python3
__author__ = 'Tane Kafle (tk916@ic.ac.uk)'
__version__ = '0.0.1'
__date__ = 'August 2018'
import os, sys, argparse
parser = argparse.ArgumentParser(usage = "python3 getVCFstats.py -i STR -o STR -p1 INT -p2 INT [-p1n STR -p2n STR -h]\n\nPython 3 required!!\n\nThis program calculates VCF stats. It works for VCFs containing two populations and can give stats on fixed polymorphisms between the two, shared polymorphisms, private to population 1 or population 2, and if they are monomorphic. The stat code is:\n\t0 = fixed\n\t1 = private to pop1\n\t2 = private to pop2\n\t3 = shared\n\t4 = monomorphic\n\t5+ = error. A .vcs output file containing the variant call stat code for each SNP is also printed to a desired directory.\n\nAuthor: Tane Kafle\n\n")
parser.add_argument("-i", "--input_vcf", required=True, help="input VCF file pathway.")
parser.add_argument("-o", "--output_vcs", required=True, help="output variant call statistics directory pathway. Please include final forward slash.")
parser.add_argument("-p1", "--pop1", required=True, help="Number of individuals/samples in the VCF from population 1.", type=int)
parser.add_argument("-p2", "--pop2", required=True, help="Number of individuals/samples in the VCF from population 2.", type=int)
parser.add_argument("-p1n", "--pop1_name", help="Population 1 name/identifier.", default="pop1", type=str)
parser.add_argument("-p2n", "--pop2_name", help="Population 2 name/identifier.", default="pop2", type=str)

# I want shared polymorphisms. I want fixed differences. I want those different in one and those different in the other.
# Counting reference allele
# 1 0 fixed difference, 0 1 fixed difference.
# 0.7 0.5 shared polymoprhism, 0.5 0.7.
# 1 0.3 only in one, 0.3 1 only in the other
# 0 means fixed, 1 means found only in pop1, 2 means found only in pop 2 and 3 means found in both, 4 means monomorphic and 5+ means error.


args = parser.parse_args()

#Getting the number of individuals per population, so we can get the location via indices of what sample data on VCF is from which population.
pop1len = args.pop1
pop2len = args.pop2

p1vl = 9 #population 1 VCF indices lower bound.
p1vh = 9 + pop1len #population 1 VCF indices high i.e. what index is that of the last sample of Pop 1.
p2vl = 9 + pop1len #population 2 VCF indices low
p2vh = 9 + pop1len + pop1len #population 1 VCF indices high

vcfin = open(args.input_vcf, 'r')

vcsoutname = args.input_vcf.split('/')[-1] #extracting file name
vcsoutname = vcsoutname.split('.')
vcsoutname.pop() #removing the ending of the file name i.e. .vcf.
vcsoutname = '.'.join(vcsoutname)
vcsoutname = args.output_vcs + vcsoutname + ".vcs"

vcsout = open(vcsoutname, 'w')
vcsout.write("Contig\tReference\tAlternative\t" + args.pop1_name + "_Indv\t" + args.pop2_name + "_Indv\t" + args.pop1_name + "_RefAllele_Count\t" + args.pop2_name + "_RefAllele_Count\tStat_Code") #Creating header of stats file.

# Variables to store the number of appearances of SNPs falling under these groups.
fixed = 0
share = 0
popu1 = 0
popu2 = 0
monom = 0
pfix1 = 0 # These are all potential ones, where there is complete missing data in some SNPs.
pfix2 = 0
ppop1 = 0
ppop2 = 0
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
        if (pop1alleles_sep.count('0') + pop1alleles_sep.count('1')) == 0:
            pop1_ref_freq = -1
            pop2_ref_freq = pop2alleles_sep.count('0') / (pop2alleles_sep.count('0') + pop2alleles_sep.count('1'))
        elif (pop2alleles_sep.count('0') + pop2alleles_sep.count('1')) == 0:
            pop2_ref_freq = -1
            pop1_ref_freq = pop1alleles_sep.count('0') / (pop1alleles_sep.count('0') + pop1alleles_sep.count('1'))
        else:
            pop1_ref_freq = pop1alleles_sep.count('0') / (pop1alleles_sep.count('0') + pop1alleles_sep.count('1'))
#            pop1_alt_freq = pop1alleles_sep.count('1') / (pop1alleles_sep.count('0') + pop1alleles_sep.count('1'))
            pop2_ref_freq = pop2alleles_sep.count('0') / (pop2alleles_sep.count('0') + pop2alleles_sep.count('1'))
#            pop2_alt_freq = pop2alleles_sep.count('1') / (pop2alleles_sep.count('0') + pop2alleles_sep.count('1'))
        
        # Doing checks of reference allele frequencies to determine if SNP is fixed, shared, private to one of the populations, monomorphic, or throws up an error.
        
        if pop1_ref_freq == -1:
            if pop2_ref_freq == -1:
                print("There are SNPs with no data for either population. Please remove.")
                sys.exit()
            elif pop2_ref_freq == 0 or pop2_ref_freq == 1:
                stat = 6.0 # Fixed in pop2 but missing in pop1.
                pfix2 +=1
            elif pop2_ref_freq > 0 and pop2_ref_freq < 1:
                stat = 6.2 # Fixed difference
                ppop2 +=1
        elif pop2_ref_freq == -1:
            if pop1_ref_freq == 0 or pop1_ref_freq == 1:
                stat = 6.0 # Fixed in pop1 but missing in pop1.
                pfix1 +=1
            elif pop1_ref_freq > 0 and pop1_ref_freq < 1:
                stat = 6.1 #Fixed difference
                ppop1 +=1
        elif pop1_ref_freq == 1:
            if pop2_ref_freq == 0:
                stat = 0
                fixed +=1
            elif pop2_ref_freq == 1:
                stat = 4
                monom +=1
            elif pop2_ref_freq > 0 and pop2_ref_freq < 1:
                stat = 2
                popu2 +=1
            else:
                stat = 5.1
        elif pop1_ref_freq == 0:
            if pop2_ref_freq == 1:
                stat = 0
                fixed +=1
            elif pop2_ref_freq == 0:
                stat = 4
                monom +=1
            elif pop2_ref_freq > 0 and pop2_ref_freq < 1:
                stat = 2
                popu2 +=1
            else:
                stat = 5.2
        elif pop1_ref_freq > 0 and pop1_ref_freq < 1:
            if pop2_ref_freq > 0 and pop2_ref_freq < 1:
                stat = 3
                share +=1
            elif pop2_ref_freq == 0 or pop2_ref_freq == 1:
                stat = 1
                popu1 += 1
            else:
                stat = 5.3
        else:
            stat = 5.0
                
            
        #Count number of individuals without data
        pop1_inv = pop1len - (pop1alleles_sep.count('.') / 2)
        pop2_inv = pop2len - (pop2alleles_sep.count('.') / 2)

        outputline = [contig_nam, ref_allele, alt_allele, str(int(pop1_inv)), str(int(pop2_inv)), str(pop1alleles_sep.count('0')), str(pop2alleles_sep.count('0')), str(stat)]
        vcsout.write("\n" + "\t".join(outputline))

print("Stats\nNumber of FIXED: " + str(fixed) + "\nNumber of SHARED: " + str(share) + "\nNumber of PRIVATE (" + args.pop1_name + ") : " + str(popu1) + "\nNumber of PRIVATE (" + args.pop2_name + ") : " + str(popu2) + "\nNumber of MONOMORPHIC: " + str(monom) + "\n") # Writing stats to terminal.

print("MISSING DATA IN ONE POPULATION:\nPotentially fixed (" + args.pop1_name + " & " + args.pop2_name + ") : " + str(pfix1) + " + " + str(pfix2) + "\nPotentially private (" + args.pop1_name + ") : " + str(ppop1) + "\nPotentially private (" + args.pop2_name + ") : " + str(ppop2))