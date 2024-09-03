'''Script for running the TTo method. 
The user can tweak parameters defined at the start if they so wish.
Functions used are found in the functions.py script.
Please see the Manual (ReadMe.md) or use the command 'python3 TT_Method.py --help' for help.

Milo Thordarson: anth2886@student.uu.se'''

# Module loading
import functions # TTo specific functions from functions.py file
import argparse
import os
import multiprocessing
import sys

# Filtering parameters that can be changed by the user
# These represent what you consider low and high coverage for a genotype position, as an int
low_coverage = 10
high_coverage = 500
# All acceptable values in the FILTERS column of a vcf file, default set to be passing all filters or a non-entry
vcf_filters = ['PASS','.']


# Arguments from commandline using argparse
parser = argparse.ArgumentParser()
parser.add_argument("-1", "--pop1",
                    required = True, 
                    nargs = '+',
                    help = "genotype files for the first population")
parser.add_argument("-2", "--pop2",
                    required = True,
                    nargs = '+',
                    help = "genotype files for the second population")
parser.add_argument("-og", "--outgroup",
                    required = True,
                    nargs = '+',
                    help = "genotype files for the outgroup population")
parser.add_argument("-t", "--type", 
                    required = True, 
                    choices = ['vcf', 'tped', 'bam'], 
                    help = "type of genotype files input")
parser.add_argument("-a", "--ancestral",  
                    required = True,
                    nargs = '+',
                    help = "files containing ancestral states")
parser.add_argument("-k", "--keywords",  
                    nargs = 3, 
                    default = ["pop1", "pop2", "outgroup"],
                    help = "names of 2 populations, MUST be same order as -1 and -2 flag then the outgroup")
parser.add_argument("-o", "--out", 
                    default = "TTo_out_pop1_pop2",
                    help = "name of the output directory")
parser.add_argument("-c", "--counts",
                    action="store_true",
                    help = "output a file with all counts per chromosome per window")
parser.add_argument("-w", "--window",
                    default = "5000000",
                    help = "set the window size for calculating local parameters")
# REMOVE --test: for testing purposes only
parser.add_argument("--test", action = "store_true")
args = parser.parse_args()

# Turn args into informatively named variables
files_pop1 = args.pop1
files_pop2 = args.pop2
files_outgroup = args.outgroup
file_type = args.type
files_anc = args.ancestral
pop1_key = args.keywords[0] 
pop2_key = args.keywords[1]
outg_key = args.keywords[2]
out_dir = args.out
print_counts = args.counts
win_size = args.window

# Check that filepaths are the same lengths
file_tot = len(files_pop1)
if any(len(lst) != file_tot for lst in [files_pop2, files_outgroup, files_anc]):
    print("Error: Unequal amount of files for pop1, pop2 or ancestral.")
    print(f'Files given for pop1 are: {files_pop1}')
    print(f'Files given for pop2 are: {files_pop2}')
    print(f'Files given for outgroup are: {files_outgroup}')
    print(f'Files given for ancestral states are: {files_anc}')
    sys.exit(1)

# Make output dir, will complain if it already exists, which is why this is so early in the script
if not args.test: os.mkdir(out_dir)

# Initialise the directory with all the counts per chromosome
counts_dict = {}

# For vcf filetype
if file_type == 'vcf': 
    # Create iterable list with all input parameters for counting
    iterables = [[files_pop1[i], files_pop2[i], files_anc[i], low_coverage, high_coverage, vcf_filters, win_size] for i in range(file_tot)]
    iterables_outgroup = [[files_pop1[i], files_pop2[i], files_outgroup[i], files_anc[i], low_coverage, high_coverage, vcf_filters, win_size] for i in range(file_tot)]
    # To avoid infinite recursion
    if __name__ == '__main__':
        with multiprocessing.Pool() as pool:
            # Computes for files in parallel using CPU cores available to user
            results = pool.map(functions.get_counts_vcf_TT, iterables)
            results_outgroup = pool.map(functions.get_counts_vcf_TTo, iterables_outgroup)

# From results get a list of the counts and if user selected, print the counts by chromosome and window and output to file
counts = []
windows = {}
counts_outgroup = []
windows_outgroup = {}

if print_counts:
    count_file = open(out_dir + "/" + pop1_key + pop2_key + "_TT_Counts.txt", 'w')
    count_file_outgroup = open(out_dir + "/" + pop1_key + pop2_key + "_TTo_Counts.txt", 'w')

for comparison in results:
    for key in comparison[1]:
        if not key in windows:
            windows.update({key: []})
    counts.extend(comparison[0][key])
    windows[key].extend(comparison[1][key])
    if print_counts:
        count_file.write("#" + key + "\n")
        for i in range(len(windows[key])):
                count_file.write(str(windows[key][i]) + "\t" + str(counts[i]) + "\n")

for comparison in results_outgroup:
    for key in comparison[1]:
        if not key in windows_outgroup:
            windows_outgroup.update({key: []})
    counts_outgroup.extend(comparison[0][key])
    windows_outgroup[key].extend(comparison[1][key])
    if print_counts:
        count_file_outgroup.write("#" + key + "\n")
        for i in range(len(windows_outgroup[key])):
                count_file_outgroup.write(str(windows_outgroup[key][i]) + "\t" + str(counts_outgroup[i]) + "\n")

if print_counts: 
    count_file.close()
    count_file_outgroup.close()

