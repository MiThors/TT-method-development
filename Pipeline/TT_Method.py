'''Scrip for running the TT method without an outgroup. 
The user can tweak parameters defined at the start if they so wish.
Functions used are found in the functions.py script.
Please see the Manual (ReadMe.md) or use the command 'TT_Method.py --help' for help.

Milo Thordarson: anth2886@student.uu.se'''

# Module loading
import TT_functions
import argparse
import os

# Filtering parameters that can be changed by the user
# These represent what you consider low and high coverage for a genotype position, as an int
low_coverage = 10
high_coverage = 500
# This list contains all acceptable values in the FILTERS column of a vcf file, default set to be passing all filters or a non-entry
vcf_filters = ['PASS','.']


# Arguments from commandline using argparse
parser = argparse.ArgumentParser()
parser.add_argument("-1", "--ind1",
                    required = True, 
                    nargs = '+',
                    help = "genotype files for the first population")
parser.add_argument("-2", "--ind2",
                    required = True,
                    nargs = '+',
                    help = "genotype files for the second population")
parser.add_argument("-t", "--type", 
                    required = True, 
                    choices = ['vcf', 'tped', 'bam'], 
                    help = "type of genotype files input")
parser.add_argument("-a", "--ancestral",  
                    required = True,
                    nargs = '+',
                    help = "files containing ancestral states")
parser.add_argument("-k", "--keywords",  
                    nargs = 2, 
                    default = ["ind1", "ind2"],
                    help = "names of 2 populations, MUST be same order as -1 and -2 flag")
parser.add_argument("-o", "--out", 
                    default = "TT_out_ind1_ind2",
                    help = "optional flag to indicate that all chromosomes are in the same file")
parser.add_argument("--test", action = "store_true")
args = parser.parse_args()

# Turn args into informatively named variables
files_ind1 = args.ind1
files_ind2 = args.ind2
file_type = args.type
files_anc = args.ancestral
ind1_key, ind2_key = args.keywords
out_dir = args.out

# Check that filepaths are the same lengths
file_tot = len(files_ind1)
if any(len(lst) != file_tot for lst in [files_ind2, files_anc]):
    print("Error: Unequal amount of files for ind1, ind2 or ancestral.")
    print(f'Files given for ind1 are: {files_ind1}')
    print(f'Files given for ind2 are: {files_ind2}')
    print(f'Files given for ancestral states are: {files_anc}')
    exit(1)

# Make output dir, will complain if it already exists, which is why this is so early in the script
if not args.test: os.mkdir(out_dir)

if file_type == 'vcf': 
    for i in range(file_tot):
        TT_functions.get_counts_vcf(files_ind1[i], files_ind2[i], files_anc[i], low_coverage, high_coverage, vcf_filters)
        # And then you will need to combine them somewhere here
