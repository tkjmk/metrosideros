#!/usr/bin/python3 -W ignore
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
import shutil as sh
import random
from os import listdir
from os.path import isfile, join


def parse_arguments():
    parser = argparse.ArgumentParser(usage = 'python3 select_random_files.py -i STR -o STR [-f STR -m STR -n FLOAT -a STR -h]\n\nPython 3 required!!\n\nThis program selects a random (or pseudo-random, if you prefer) subset of files from <input_directory> and then, via the selection of the <action> flag, copies or moves them into the <output_directory>. The number of files to be subsetted can be determined using the <method_of_subsetting> flag, followed by the <number_of_files_to_subset> flag, where a value between 0 and 1 is used if p (percentage) is selected or an integer below or equal to the total number of files is entered if a (absolute) is selected. The <file_extension> flag can be used so that only files with a certain extension are copied/moved from the <input_directory>.\n\nPlease ensure the output directory is created and empty! If there are files with the same name, these will not be copied/moved into the directory.\n\nPlease ensure that <input_directory> and <output_directory> pathways end with a forward slash.\n\nUse -h to show full help screen.\n\nAuthor: Tane Kafle\n\n')
    parser.add_argument("-i", "--input_directory", required=True, help="Input directory that contains the files. Please end with '/' at the end of pathway.")
    parser.add_argument("-o", "--output_directory", required=True, help="The directory where you would like to output the subset. Please end with '/' at the end of pathway.")
    parser.add_argument("-f", "--file_extension", help="Put file ending e.g. '.fas' to only choose a random subset of these files", default=None, type=str)
    parser.add_argument("-m", "--method_of_subsetting", help="Choose random subset of files as absolute (a) or percentage (p) of total number of files.", default='p', choices=['a','p'], type=str)
    parser.add_argument("-n", "--number_of_files_to_subset", help="If percentage chosen in <method_of_subsetting> put a value between 0 and 1, if absolute choose a whole number below the number of total files to be chosen from.", default=1, type=float)
    parser.add_argument("-a", "--action", help = "Argument to copy (c) or move (m) the files from input directory to output directory.", default='c', choices = ['c','m'], type=str)
    return parser.parse_args()

def main():
    args = parse_arguments()

    inpath = args.input_directory
    outpath = args.output_directory

    # Building a list of filenames within the directory (with or without specific file extension depending on <file_extension> flag).
    if args.file_extension == None:
        files = [f for f in listdir(inpath) if isfile(join(inpath, f))]
        if len(files) == 0:
            raise ValueError("There are no files in this directory, please select <input_directory> that contains files.")
    else:
        files = [f for f in listdir(inpath) if f.endswith(args.file_extension) and isfile(join(inpath, f))]
        if len(files) == 0:
            raise ValueError("There are no files with a " + args.file_extension + " extension in this directory, please select <input_directory> that contains files with that extension or <file_extension> found within this dirrectory.")

    # This part determines the number of files that are going to be subsetted depending on the <method_of_subsetting> and <number_of_files_to_subset>.
    if args.method_of_subsetting == "a":
        filenum = int(args.number_of_files_to_subset) #If a float is input e.g. 12.5, this will always be rounded down to 12.
    elif args.method_of_subsetting == "p":
        if args.number_of_files_to_subset > 1 or args.number_of_files_to_subset < 0:
            raise ValueError("Please input a <number_of_files_to_subset> value between 0 and 1 inclusive.")
        filenum = int(args.number_of_files_to_subset * len(files))

    if filenum > len(files):
        raise ValueError(str(filenum) + " files attempted to be subsetted from " + str(len(files)) + " total files.\nPlease lower <number_of_files_to_subset> so that it is below the total number of files to be subsetted.\n")


    subset_files = random.sample(files, filenum) # This line selects a random subset of files 


    # The following part of script either copies or moves the subsetted files to the output directory and then prints a message on the number of files from the total number of files that have been copied/moved over.

    if args.action == 'c':
        for filename in subset_files:
            sh.copy2(inpath + filename, args.output_directory)
        print(str(filenum) + " files have been subsetted from " + str(len(files)) + " total files and copied into " + outpath)
    elif args.action == 'm':
        for filename in subset_files:
            sh.move(outpath + filename, args.output_directory)
        print(str(filenum) + " files have been subsetted from " + str(len(files)) + " total files and moved into " + outpath)


if __name__ == "__main__":
    main()