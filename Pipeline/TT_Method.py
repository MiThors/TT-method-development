'''Scrip for running the TT method without an outgroup. 
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

# Arguments from commandline using argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--in",
                    required = True, 
                    nargs = 2,
                    help = "2 genotype files for populations to compare")
parser.add_argument("-t", "--type", 
                    required = True, 
                    choices = ['vcf', 'tped', 'bam'], 
                    help = "type of genotype files input")
parser.add_argument("-a", "--ancestral",  
                    required = True, 
                    help = "file containing ancestral states")
parser.add_argument("-k", "--keywords",  
                    nargs = 2, 
                    default = ["ind1", "ind2"],
                    help = "names of 2 populations, MUST be same order as in --in")
parser.add_argument("-c", "--chromosomes", 
                    action = "store_true", 
                    help = "optional flag to indicate that all chromosomes are in the same file")
parser.add_argument("-o", "--out", 
                    default = "TT_out_ind1_ind2",
                    help = "optional flag to indicate that all chromosomes are in the same file")
args = parser.parse_args()

# Filtering parameters that can be changed by the user

