'''Scrip for running the TT method. 
The user can tweak parameters defined at the start if they so wish.
Functions used are found in the functions.py script.
Please see the Manual (ReadMe.md) or use the command 'TT_Method.py --help' for help.

Milo Thordarson: anth2886@student.uu.se'''

# Module loading
import functions
import argparse
import os
import time
import multiprocessing

# Filtering parameters that can be changed by the user
# These represent what you consider low and high coverage for a genotype position, as an int
low_coverage = 10
high_coverage = 500
# This list contains all acceptable values in the FILTERS column of a vcf file, default set to be passing all filters or a non-entry
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
                    default = ["pop1", "pop2"],
                    help = "names of 2 populations, MUST be same order as -1 and -2 flag")
parser.add_argument("-o", "--out", 
                    default = "TT_out_pop1_pop2",
                    help = "name of the output directory")
parser.add_argument("--test", action = "store_true")
args = parser.parse_args()

# Turn args into informatively named variables
files_pop1 = args.pop1
files_pop2 = args.pop2
file_type = args.type
files_anc = args.ancestral
pop1_key, pop2_key = args.keywords
out_dir = args.out

# Check that filepaths are the same lengths
file_tot = len(files_pop1)
if any(len(lst) != file_tot for lst in [files_pop2, files_anc]):
    print("Error: Unequal amount of files for pop1, pop2 or ancestral.")
    print(f'Files given for pop1 are: {files_pop1}')
    print(f'Files given for pop2 are: {files_pop2}')
    print(f'Files given for ancestral states are: {files_anc}')
    exit(1)

# Make output dir, will complain if it already exists, which is why this is so early in the script
if not args.test: os.mkdir(out_dir)

# Initialise the directory with all the counts per chromosome
counts_dict = {}

# For vcf filetype
if file_type == 'vcf': 
    t0 = time.time()
    iterables = [[files_pop1[i], files_pop2[i], files_anc[i], low_coverage, high_coverage, vcf_filters] for i in range(file_tot)]
    if __name__ == '__main__':
        with multiprocessing.Pool() as pool:
            results = pool.map(functions.get_counts_vcf_TT, iterables)
    
    for dict in results:
        for key in dict:
            if key in counts_dict:
                counts_dict[key] = [counts_dict[key][i] + dict[key][i] for i in range(9)]
            else: counts_dict.update({key: dict[key]})

    #for i in range(file_tot):
    #    # Get the counts for a vcf file from this function, returns a dictionary with {chrom: m0, m1, m2, m3, m4, m5, m6, m7, m8} structure
    #    add_dict = functions.get_counts_vcf_TT(files_pop1[i], files_pop2[i], files_anc[i], low_coverage, high_coverage, vcf_filters)
    #    for key in add_dict:
    #        if key in counts_dict:
    #            counts_dict[key] = [counts_dict[key][i] + add_dict[key][i] for i in range(9)]
    #        else: counts_dict.update({key: add_dict[key]})
t1 = time.time()
print(counts_dict)
print(f'Time elapsed for parallelized = {t1-t0}')

loop_counts = {}

if file_type == 'vcf': 
    t0 = time.time()
    for i in range(file_tot):
        # Get the counts for a vcf file from this function, returns a dictionary with {chrom: m0, m1, m2, m3, m4, m5, m6, m7, m8} structure
        add_dict = functions.get_counts_vcf_TT(iterables[i])
        for key in add_dict:
            if key in loop_counts:
                loop_counts[key] = [loop_counts[key][i] + add_dict[key][i] for i in range(9)]
            else: loop_counts.update({key: add_dict[key]})

t1 = time.time()
print()
print(loop_counts)
print(f'Time elapsed for looped = {t1-t0}')