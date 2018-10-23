#!/usr/bin/python3
__author__ = 'Tane Kafle (tk916@ic.ac.uk)'
__version__ = '0.0.2'
__date__ = 'May 2018'
#Update to version 0.0.2 on 8th June 2018 to also include Dxy on output, negative Fst values are coverted to 0 and changed it to Python3.

import pandas as pd
import os, sys, argparse
parser = argparse.ArgumentParser(usage = 'python2 gppfst_experimentalparams.py -f STR -s STR -b STR -o STR [-h]\n\nPython 3 required!!\n\nThis program creates a tab seperated format file that can be used as experimental parameter inputs for gppfst on R. It requires a FST file generated using VCFtools for two populations; a STATS file that has the name of the contig, position of SNP, number of samples for population 1 and population 2 from a VCF (generated from my leastmissingdataperSNP_twopop.py script), Dxy generated from my estimate_dxy.py script;  and finally a CDS-only BED file with the contig name, and start and end positions of the CDS. \n\nAuthor: Tane Kafle\n\n')
parser.add_argument("-f", "--fst_file", required=True, help="fst file pathway.")
parser.add_argument("-s", "--stats_file", required=True, help="stats file pathway.")
parser.add_argument("-b", "--bed_file", required=True, help="bed file pathway.")
parser.add_argument("-d", "--dxy_file", required=True, help="dxy file pathway.")
parser.add_argument("-o", "--output_file", required=True, help="output file pathway.")

args = parser.parse_args()

vcf_fst = pd.read_csv(args.fst_file, sep="\t")
vcf_dxy = pd.read_csv(args.dxy_file, sep = "\t", header = 0, names = ["Contig", "Drop1", "Drop2", "Dxy", "Drop3", "Drop4", "ME_NE.ref_freq", "ME_SC.ref_freq", "ME_NE.alt_freq", "ME_SC.alt_freq", "Pi.ME_NE", "Pi.ME_SC", "DeltaPi.ME_SC-ME_NE"])
stats = pd.read_csv(args.stats_file, sep="\t")
bed = pd.read_csv(args.bed_file, sep="\t", header = None, names = ["Contig_bed", "Start", "End", "Info1", "Info2", "Info3", "Info4", "Info5", "Info6", "Info7"])
bed = bed.drop(["Info1", "Info2", "Info3", "Info4", "Info5", "Info6", "Info7"], axis=1)

stats_fst = vcf_fst.merge(stats, left_on="CHROM", right_on='Contig', how='left')
stats_fst_dxy = pd.merge(stats_fst, vcf_dxy[['Contig','Dxy', 'Pi.ME_NE', 'Pi.ME_SC', 'DeltaPi.ME_SC-ME_NE']], on='Contig', how='left')
sequencelength = bed[["End"]].sub(bed["Start"], axis='index')

bedwithlength = bed.assign(seqlen=sequencelength)

final = stats_fst_dxy.merge(bedwithlength, left_on="CHROM", right_on="Contig_bed", how='left')
final = final.drop(["Contig", "Pos", "Contig_bed"], axis=1)
final = final.rename(str.upper, axis='columns')
final = final.rename(index=str, columns={"CHROM": "CONTIG"})
final[["WEIR_AND_COCKERHAM_FST"]] = final[["WEIR_AND_COCKERHAM_FST"]].clip_lower(0)

final.to_csv(args.output_file, sep="\t", index = False)

