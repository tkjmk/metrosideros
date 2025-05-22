These scripts were written for my analysis on the divergence of Lord Howe Island *Metrosideros*, and therefore were largely written for analysis comparing the two populations. Analysis on SNPs were based upon VCF files with information for two populations.

# calculate_fst.py
This script does calculates FST using various estimators by taking in a VCF file with SNP data for two populations.
It was used to calculate the FST values that were input into MINOTAUR and GppFst.
# leastmissingdataperSNP_twopop.py
This script filters a VCF file so that only one SNP remains per contig, based upon missing data. You can optionally also choose a minimum threshold of non-missing data required in both species.
It was used to remove the effects of linkage between SNPs for analyses, and 2 SNPs were required in both species to calculate dXY and FST.
# convert_vcf012_to_nexus.py
This script takes a .012 file generated from a VCF by VCFTOOLS to a nexus file.
It was used as input into BEAUti to generate an XML file for SNAPP.
# percentidentity_fasta.py
This script compares pairwise every individual from two populations and calculate the within species and between species percent identity in sequence.
It was used to estimate a prior for SNAPP.
# esimate_dxy.py
This script estimates dXY as p1(1 - p2) + p2(1 - p1) where p1 and p2 are the frequencies of a particular alelle within population 1 and population 2, respectively. It also calculates nucleotide diversity, pi.
This was to get locus-specific dXY from each SNP in VCF and input into MINOTAUR.
# remove_mismatched_seq.py
This filters pairwise FASTA sequence so that if there are more than a certain threshold of mismatches in a window, all missing data is converted to Ns.
This was used before the final steps of reference transcriptome assembly to prevent erroneous SNP calls.
# filter_fasta_by_missing_data.py
This script goes through FASTA sequence with contigs for all individuals, it will look through each one and see how much missing data there is for each individual, and if a certain threshold of sequence is not met the whole contig will be removed.
# getVCFstats.py
This script takes a VCF of two populations and returns information on the number of shared SNPs and private SNPs between two populations.
This was used to extract stats of VCFs before and after filtering.			
# select_random_files.py
This script will copy or move a selected number of files randomly from a directory to another directory.
This was written when trying to sample FASTAs for some long analyses.
# gppfst_experimentalparams.py
This script combines various stats to be input into GppFst, including a FST file generated using VCFtools for two populations; a STATS file that has the name of the contig, position of SNP, number of samples for population 1 and population 2 from a VCF (generated from my leastmissingdataperSNP_twopop.py script), dXY generated from my estimate_dxy.py script;  and finally a CDS-only BED file with the contig name, and start and end positions of the CDS.
This was used to generate a file that I ended up putting into GppFst and MINOTAUR.

Running any scripts with the -h flag will give more information on their usage.
