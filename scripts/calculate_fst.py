#!/usr/bin/python3

__author__ = 'Tane Kafle (tk916@ic.ac.uk)'
__version__ = '0.0.2'
__date__ = 'August 2018'

"""
Changelog:
----------
v0.0.2 - May 2025
    - Refactored to include main() function and argparse handling

v0.0.1 - August 2018
    - Initial version
"""

import os, sys, argparse



def FST_Wright_pairwise (col): #According to Chen
    NpopA = float(col[0]) #Number of chromosomes there are for samples in each population
    NpopB = float(col[2])
    
    popAcount= int(col[1]) #Number of particular allele
    popBcount= int(col[3])
    
    pop1freq = popAcount / float(NpopA )
    pop2freq = popBcount / float(NpopB )

    return (((pop1freq-pop2freq)**2) / ((pop1freq+pop2freq) * (2-(pop1freq+pop2freq))))


def FST_W_pairwise (col): #According to https://github.com/pontussk/popstats/blob/master/popstats.py
    NpopA = float(col[0]) #Number of chromosomes there are for samples in each population
    NpopB = float(col[2])
    
    popAcount= int(col[1]) #Number of particular allele
    popBcount= int(col[3])
    

    npops= 2.0
    #nsamples = float(NpopA + NpopB)
    n_bar= (NpopA / npops) + (NpopB / npops)
    samplefreq = ( (popAcount+popBcount) / (NpopA + NpopB) )
    pop1freq = popAcount / float(NpopA )
    pop2freq = popBcount / float(NpopB )
    Npop1 = NpopA
    Npop2 = NpopB
    S2A= (1/ ( (npops-1.0) * n_bar) ) * ( ( (Npop1)* ((pop1freq-samplefreq)**2) ) + ( (Npop2)*((pop2freq-samplefreq)**2) ) )
    nc = 1.0/(npops-1.0) * ( (Npop1+Npop2) - (((Npop1**2)+(Npop2**2)) / (Npop1+Npop2)) )
    T_1 = S2A -( ( 1/(n_bar-1) ) * ( (samplefreq * (1-samplefreq)) -  ((npops-1)/npops)* S2A ) )
    T_2 = (( (nc-1) / (n_bar-1) ) * samplefreq *(1-samplefreq) )   +  (1.0 +   (((npops-1)*(n_bar-nc))  / (n_bar-1)))       * (S2A/npops)
    FST2 = T_1/T_2
    return (T_1,T_2,FST2)


def FST_W_pairwise_TK (col): #According to Chen
    NpopA = float(col[0]) #Number of chromosomes there are for samples in each population
    NpopB = float(col[2])
    
    popAcount= int(col[1]) #Number of particular allele
    popBcount= int(col[3])
    

    npops= 2.0
    #nsamples = float(NpopA + NpopB)
    n_bar= (NpopA / npops) + (NpopB / npops)
    samplefreq = ( (popAcount+popBcount) / (NpopA + NpopB) )
    pop1freq = popAcount / float(NpopA )
    pop2freq = popBcount / float(NpopB )
    Npop1 = NpopA
    Npop2 = NpopB
    S2A= (1/ ( (npops-1.0) * n_bar) ) * ( ( (Npop1)* ((pop1freq-samplefreq)**2) ) + ( (Npop2)*((pop2freq-samplefreq)**2) ) ) 
    nc = 1.0/(npops-1.0) * ( (Npop1+Npop2) - (((Npop1**2)+(Npop2**2)) / (Npop1+Npop2)) )
    T_1 = S2A -( ( 1/((2*n_bar)-1) ) * ( (samplefreq * (1-samplefreq)) -  ((npops-1)/npops)* S2A ) ) #the n_bar is multiplied by 2, as it is said in the Chen paper.
    T_2 = (( (nc-1) / (n_bar-1) ) * samplefreq *(1-samplefreq) )   +  (1.0 +   (((npops-1)*(n_bar-nc))  / (n_bar-1)))       * (S2A/npops)
    FST2 = T_1/T_2
    return (T_1,T_2,FST2)


def FST_H_pairwise (col): #According to https://github.com/pontussk/popstats/blob/master/popstats.py
    NpopA = float(col[0])
    NpopB = float(col[2])
    popAcount= int(col[1])
    popBcount= int(col[3])
    pop1freq = popAcount / float(NpopA )
    pop2freq = popBcount / float(NpopB )
    Npop1 = NpopA
    Npop2 = NpopB
    T_1=(pop1freq-pop2freq)**2 - ((pop1freq*(1.0-pop1freq))/(Npop1-1)) - ((pop2freq*(1.0-pop2freq))/(Npop2-1))
    T_2=(pop1freq*(1.0-pop2freq)) + (pop2freq*(1.0-pop1freq))
    FST3 = T_1/T_2
    #T_1=T_2
    # #T_2=1.0
    return (T_1,T_2, FST3)


def FST_H_pairwise_TK (col): #According to Chen
    NpopA = float(col[0])
    NpopB = float(col[2])
    popAcount= int(col[1])
    popBcount= int(col[3])
    pop1freq = popAcount / float(NpopA )
    pop2freq = popBcount / float(NpopB )
    Npop1 = NpopA
    Npop2 = NpopB
    npops = 2.0
    #samplefreq = ( (popAcount+popBcount) / (NpopA + NpopB) )

    HW=(1/npops) * (((2*pop1freq) * (1-pop1freq)) + ((2*pop2freq) * (1-pop2freq)))
    HB=(1/(npops*(npops-1))) * (((2*pop1freq)*(1-pop2freq)) + (2*pop2freq)*(1-pop1freq))

    T_1=(pop1freq-pop2freq)**2 - ((pop1freq*(1.0-pop1freq))/(Npop1-1)) - ((pop2freq*(1.0-pop2freq))/(Npop2-1))
    T_2=(pop1freq*(1.0-pop2freq)) + (pop1freq*(1.0-pop2freq))
    #FST3 = T_1/T_2
    FST3 = 1 - (HW/HB)
    #T_1=T_2
    # #T_2=1.0
    return (T_1,T_2, FST3)

def FST_H_pairwise_TK_unbiased (col): #According to Chen
    NpopA = float(col[0])
    NpopB = float(col[2])
    popAcount= int(col[1])
    popBcount= int(col[3])
    pop1freq = popAcount / float(NpopA )
    pop2freq = popBcount / float(NpopB )
    Npop1 = NpopA
    Npop2 = NpopB
    npops = 2.0
    #samplefreq = ( (popAcount+popBcount) / (NpopA + NpopB) )

    HW=(npops-1) * (((((2 * Npop1)/((2*Npop1)-1))*pop1freq) * (1-pop1freq)) + (((((2 * Npop2)/((2*Npop2)-1))*pop2freq) * (1-pop2freq))))
    HB=((2*pop1freq)*(1-pop2freq)) + ((2*pop2freq)*(1-pop1freq))

    T_1=(pop1freq-pop2freq)**2 - ((pop1freq*(1.0-pop1freq))/(Npop1-1)) - ((pop2freq*(1.0-pop2freq))/(Npop2-1))
    T_2=(pop1freq*(1.0-pop2freq)) + (pop1freq*(1.0-pop2freq))
    #FST3 = T_1/T_2
    FST3 = 1 - (HW/HB)
    #T_1=T_2
    # #T_2=1.0
    return (HW,HB, FST3)

# FST_W_pairwise([12,4,14,2])
# FST_W_pairwise_TK([12,4,14,2])
# FST_H_pairwise([12,4,14,2])
# FST_Wright_pairwise([12,4,14,2])

def parse_arguments():
    parser = argparse.ArgumentParser(usage = "python3 calculate_fst.py -i STR -o STR -p1 INT -p2 INT [-h]\n\nPython 3 required!!\n\nThis program calculates Fst as according to Chen et al. (2015), https://doi.org/10.1371/journal.pone.0135368. It calculates Wright's, Weir-Cockerham and Hudson. A .fst output file is also printed to a desired directory.\n\nAuthor: Tane Kafle\n\n")
    parser.add_argument("-i", "--input_vcf", required=True, help="input VCF file pathway.")
    parser.add_argument("-o", "--output_fst", required=True, help="output FST directory pathway. Please include final forward slash.")
    parser.add_argument("-p1", "--pop1", required=True, help="Number of individuals/samples in the VCF from population 1.", type=int)
    parser.add_argument("-p2", "--pop2", required=True, help="Number of individuals/samples in the VCF from population 2.", type=int)
    return parser.parse_args()


def main():
    args = parse_arguments() #parser.parse_args()

    #Getting the number of individuals per population, so we can get the location via indices of what sample data on VCF is from which population.
    pop1len = args.pop1
    pop2len = args.pop2

    p1vl = 9 #population 1 VCF indices lower bound.
    p1vh = 9 + pop1len #population 1 VCF indices high i.e. what index is that of the last sample of Pop 1.
    p2vl = 9 + pop1len #population 2 VCF indices low
    p2vh = 9 + pop1len + pop1len #population 1 VCF indices high

    vcfin = open(args.input_vcf, 'r')

    fstoutname = args.input_vcf.split('/')[-1] #extracting file name
    fstoutname = fstoutname.split('.')
    fstoutname.pop() #removing the ending of the file name i.e. .vcf.
    fstoutname = '.'.join(fstoutname)
    fstoutname = args.output_fst + fstoutname + ".fst"

    fstout = open(fstoutname, 'w')
    fstout.write("Contig\tReference\tAlternative\tPop1_indiv\tPop1_refcount\tPop2_indiv\tPop2_refcount\tWright_FST\tW_FST\tW_FST_TK\tH_FST\tH_FST_TK\tH_FST_TK_ub") #Creating header of stats file.

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
            pop1_ref_count = pop1alleles_sep.count('0')
            pop1_alt_count = pop1alleles_sep.count('1')
            pop2_ref_count = pop2alleles_sep.count('0')
            pop2_alt_count = pop2alleles_sep.count('1')

            #Count number of individuals without data
            pop1_inv = pop1len - (pop1alleles_sep.count('.') / 2)
            pop2_inv = pop2len - (pop2alleles_sep.count('.') / 2)

            stats_ref = [pop1_inv*2, pop1_ref_count, pop2_inv*2, pop2_ref_count] #To input into FST things
            stats_alt = [pop1_inv*2, pop1_alt_count, pop2_inv*2, pop2_alt_count]
            #print(contig_nam + " and " + str(pop1_inv*2) + " " + str(pop1_ref_count) + " " + str(pop2_inv*2) + " " + str(pop2_ref_count)) #Hash for no output.
            #print("ALT: " + str(pop1_alt_count) + " " + str(pop2_alt_count)) #Hash for no output.

            if pop2_inv*2 == pop2_ref_count:# If all individuals have the same alleles from pop 2...
                if pop1_inv*2 == pop1_ref_count: #...and pop1, the FST is 1 for all 3.
                    Wr_FST_ref, WC_FST_ref, WC_FST_TK_ref, Hu_FST_ref, Hu_FST_TK_ref,Hu_FST_TK_ref_ub = [0.0,0.0,0.0,0.0,0.0,0.0]
                else: #Swapping them around (and hoping errors don't pop up)
                    stats_ref = [pop2_inv*2, pop2_ref_count, pop1_inv*2, pop1_ref_count]
                    Wr_FST_ref = FST_Wright_pairwise(stats_ref)
                    Wr_FST_alt = FST_Wright_pairwise(stats_alt)
                    WC_FST_ref = FST_W_pairwise(stats_ref)[2]
                    WC_FST_TK_ref = FST_W_pairwise_TK(stats_ref)[2]
                    Hu_FST_ref = FST_H_pairwise(stats_ref)[2]
                    Hu_FST_TK_ref = FST_H_pairwise_TK(stats_ref)[2]
                    Hu_FST_TK_ref_ub = FST_H_pairwise_TK_unbiased(stats_ref)[2]
            elif pop1_ref_count == 0 and pop2_ref_count == 0: #If they are all 0, FST is 1.
                Wr_FST_ref, WC_FST_ref, WC_FST_TK_ref, Hu_FST_ref, Hu_FST_TK_ref,Hu_FST_TK_ref_ub = [0.0,0.0,0.0,0.0,0.0,0.0]
            else: #Do the calculations.
                Wr_FST_ref = FST_Wright_pairwise(stats_ref)
                Wr_FST_alt = FST_Wright_pairwise(stats_alt)
                WC_FST_ref = FST_W_pairwise(stats_ref)[2]
                WC_FST_TK_ref = FST_W_pairwise_TK(stats_ref)[2]
                Hu_FST_ref = FST_H_pairwise(stats_ref)[2]
                Hu_FST_TK_ref = FST_H_pairwise_TK(stats_ref)[2]
                Hu_FST_TK_ref_ub = FST_H_pairwise_TK_unbiased(stats_ref)[2]
            if abs(round(Wr_FST_ref,3) - round(Wr_FST_alt,3)) >= 0.1: #The fst using ref or alt allele should give the same result if only biallelic SNPs are used. Rounded to 3sf and difference greater than 0.1, because sometimes there are slight different results from calculated and rounded values.
                print("There are multiallelic SNPs present. Please remove them for this analysis.")
                sys.exit()
            outputline = [contig_nam, ref_allele, alt_allele, str(pop1_inv), str(pop1_ref_count), str(pop2_inv), str(pop2_ref_count), str(Wr_FST_ref), str(WC_FST_ref), str(WC_FST_TK_ref), str(Hu_FST_ref), str(Hu_FST_TK_ref), str(Hu_FST_TK_ref_ub)]
            fstout.write("\n" + "\t".join(outputline))


    fstout.close()


if __name__ == "__main__":
    main()
