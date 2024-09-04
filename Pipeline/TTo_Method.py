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
parser.add_argument("-T", "--TTCounts",
                    nargs = 1,
                    default = False,
                    help = "optional file location of TT counts if already run")
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
file_TT = args.TTCounts
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
    if not file_TT: iterables = [[files_pop1[i], files_pop2[i], files_anc[i], low_coverage, high_coverage, vcf_filters, win_size] for i in range(file_tot)]
    iterables_outgroup = [[files_pop1[i], files_pop2[i], files_outgroup[i], files_anc[i], low_coverage, high_coverage, vcf_filters, win_size] for i in range(file_tot)]
    # To avoid infinite recursion
    if __name__ == '__main__':
        with multiprocessing.Pool() as pool:
            # Computes for files in parallel using CPU cores available to user
            if not file_TT: results = pool.map(functions.get_counts_vcf_TT, iterables)
            results_outgroup = pool.map(functions.get_counts_vcf_TTo, iterables_outgroup)
        pool.close()

# From TT results get a list of the counts and if user selected, print the counts by chromosome and window and output to file
counts = []
outgroup_counts = []
if print_counts:
    if not file_TT: count_file = open(out_dir + "/" + pop1_key + pop2_key + "_TT_Counts.txt", 'w')
    outgroup_count_file = open(out_dir + "/" + pop1_key + pop2_key + "_TTo_Counts.txt", 'w')
# Comparison being one group of files if multiple were submitted
if not file_TT:
    for comparison in results:
        for chrom in comparison:
            # First list for each chromosome are the counts
            counts.extend(comparison[chrom][0])
            if print_counts:
                count_file.write("#" + chrom + "\n")
                # The second list contains the window positions
                for i in range(len(comparison[chrom][0])):
                    count_file.write(str(comparison[chrom][1][i]) + "\t" + str(comparison[chrom][0][i]) + "\n")
    if print_counts: count_file.close()

# Same is done for the TTo counts and if TT counts, chromosomes are compared
if file_TT:
    with open('TT_out_pop1_pop2/DSub100NSub100_TT_Counts.txt','rt',encoding='utf-8') as TT_counts:
        l = TT_counts.readline()
        chrom = l[1:].strip()
        TT_dict = {chrom: [[],[]]}
        l = TT_counts.readline()
        while l:
            if l[0] == '#' and l[1:] not in TT_dict: 
                chrom = l[1:].strip()
                TT_dict.update({chrom: [[],[]]})
            elif l.strip(): 
                all_values = l.split()
                window = tuple([int(x.strip('(),')) for x in all_values[0:2]])
                TT_dict[chrom][1].append(window)
                count = [int(x.strip('[],')) for x in all_values[2:]]
                TT_dict[chrom][0].append(count)
                counts.append(count)
                
            l = TT_counts.readline()
    TT_counts.close()
  
for comparison in results_outgroup:
    for chrom in comparison:
        if file_TT and chrom not in TT_dict: 
            print("Error: a chromosome generated from this TTo run was not found in the file from the TT run. Please make sure that the input files used for this TTo are exactly the same as were used in TT.")
            sys.exit(1)
        if file_TT and TT_dict[chrom][1] != comparison[chrom][1]:
            print("Error: the windows from this TTo run and in the file from the TT run were not the same. Please make sure that the input files used for this TTo are exacctly the same as were used in TT.")
            sys.exit(1)
        # First list for each chromosome are the counts
        outgroup_counts.extend(comparison[chrom][0])
        if print_counts:
            outgroup_count_file.write("#" + chrom + "\n")
            # The second list contains the window positions
            for i in range(len(comparison[chrom][0])):
                outgroup_count_file.write(str(comparison[chrom][1][i]) + "\t" + str(comparison[chrom][0][i]) + "\n")
if print_counts: outgroup_count_file.close()