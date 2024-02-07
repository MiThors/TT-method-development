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
parser.add_argument("-c", "--chromosomes", 
                    action = "store_true", 
                    help = "optional flag to indicate that all chromosomes are in the same file")
parser.add_argument("-o", "--out", 
                    default = "TT_out_ind1_ind2",
                    help = "optional flag to indicate that all chromosomes are in the same file")
args = parser.parse_args()

# Turn args into informatively named variables
files_ind1 = args.ind1
files_ind2 = args.ind2
file_type = args.type
files_anc = args.ancestral
ind1_key, ind2_key = args.keywords
out_dir = args.out
chr_in_one_file = args.chromosomes

# Check that filepaths are the same lengths
length = len(files_ind1)
if any(len(lst) != length for lst in [files_ind2, files_anc]):
    print("Error: Unequal amount of files for ind1, ind2 or ancestral. Double check you are including all intended files.")
    exit(1)

# Other variables
nucl=['A','C','G','T']
nt_set=set(nucl)

# Make output dir, will complain if it already exists, which is why this is so early in the script
#os.mkdir(out_dir)

#if file_type == 'vcf': TT_functions.get_counts_vcf(file1_path, file2_path, anc_path, ind1_key, ind2_key, out_dir, chr_in_one_file)



