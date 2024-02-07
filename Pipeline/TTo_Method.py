'''Script for running the TTo method. 
The user can tweak parameters defined at the start if they so wish.
Functions used are found in the functions.py script.
Please see the Manual (ReadMe.md) or use the command 'TT_Method.py --help' for help.

Milo Thordarson: anth2886@student.uu.se'''

# Module loading
import sys
import gzip
import TT_functions
import argparse
from zipfile import ZipFile

# Filtering parameters that can be changed by the user
# These represent what you consider low and high coverage for a genotype position, as an int
LOW_COVERAGE = 10
HIGH_COVERAGE = 500
# This list contains all acceptable values in the FILTERS column of a vcf file, default set to be passing all filters or a non-entry
VCF_FILTERS = ['PASS','.']


# Arguments from commandline using argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files",
                    required = True, 
                    nargs = 3,
                    help = "3 genotype files for populations to compare and outgroup, outgroup MUST be last")
parser.add_argument("-t", "--type", 
                    required = True, 
                    choices = ['vcf', 'tped', 'bam'], 
                    help = "type of genotype files input")
parser.add_argument("-a", "--ancestral",  
                    required = True, 
                    help = "file containing ancestral states")
parser.add_argument("-k", "--keywords",  
                    nargs = 3, 
                    default = ["ind1", "ind2", "outgroup"],
                    help = "names of 2 populations, MUST be same order as in --files")
parser.add_argument("-c", "--chromosomes", 
                    action = "store_true", 
                    help = "optional flag to indicate that all chromosomes are in the same file")
parser.add_argument("-o", "--out", 
                    default = "TT_out_ind1_ind2",
                    help = "optional flag to indicate that all chromosomes are in the same file")
parser.add_argument("-T", "--TT", 
                    action = "store_true",
                    help = "optional flag to indicate that you want to include TT calculation of statistics")
args = parser.parse_args()

# Turn args into informatively named variables
file1_path, file2_path, outgroup = args.files
file_type = args.type
anc_path = args.ancestral
ind1_key, ind2_key, outgroup_key = args.keywords
out_dir = args.out
chr_in_one_file = args.chromosomes
include_TT_estimates = args.TT