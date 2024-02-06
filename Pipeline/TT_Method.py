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
parser.add_argument("-i", "--in",  nargs=2, help="2 genotype files for populations to compare", required=True)
parser.add_argument("-t", "--type", choices=['vcf', 'tped', 'bam'], help="type of genotype files input", required=True)
parser.add_argument("-a", "--ancestral", help="file containing ancestral states", required=True)
parser.add_argument("-k", "--keywords",  nargs=2, help="names of 2 populations, MUST be same order as in --in")
args = parser.parse_args()

# Filtering parameters that can be changed by the user

