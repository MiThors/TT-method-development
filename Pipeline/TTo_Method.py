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

[alfa1,alfa2,test1,test2,y,tau2_1,tau2_2,tau3_1,tau3_2,B1,B2,U1,U2,V1,V2,tau_test,T1,T2,J1,J2,m_counts] = functions.get_estimates_vcf_TTo(counts, outgroup_counts)

# Open all output files, print the header and estimates, close all output files
alfa1_out=open(out_dir+'/alfa1_cond.res','w')
alfa2_out=open(out_dir+'/alfa2_cond.res','w')
test1_out=open(out_dir+'/test1_cond.res','w')
test2_out=open(out_dir+'/test2_cond.res','w')
y_out=open(out_dir+'/y_cond.res','w')
tau2_1_out=open(out_dir+'/tau2_1_cond.res','w')
tau2_2_out=open(out_dir+'/tau2_2_cond.res','w')
tau3_1_out=open(out_dir+'/tau3_1_cond.res','w')
tau3_2_out=open(out_dir+'/tau3_2_cond.res','w')
B1_out=open(out_dir+'/B1_cond.res','w')
B2_out=open(out_dir+'/B2_cond.res','w')
U1_out=open(out_dir+'/U1_cond.res','w')
U2_out=open(out_dir+'/U2_cond.res','w')
V1_out=open(out_dir+'/V1_cond.res','w')
V2_out=open(out_dir+'/V2_cond.res','w')
tau_test_out=open(out_dir+'/tau_test_cond.res','w')
T1_out=open(out_dir+'/T1_cond.res','w')
T2_out=open(out_dir+'/T2_cond.res','w')
J1_out=open(out_dir+'/J1_cond.res','w')
J2_out=open(out_dir+'/J2_cond.res','w')
m_counts_out=open(out_dir+'/m_counts_cond.res','w')

header= '\t'.join(['pop1','pop2','obs_mean','wbj_mean','wbj_var'])+'\n'

alfa1_out.write(header + '\t'.join([pop1_key,pop2_key,alfa1[0],alfa1[1],alfa1[2]])+'\n')
alfa2_out.write(header + '\t'.join([pop1_key,pop2_key,alfa2[0],alfa2[1],alfa2[2]])+'\n')
test1_out.write(header + '\t'.join([pop1_key,pop2_key,test1[0],test1[1],test1[2]])+'\n')
test2_out.write(header + '\t'.join([pop1_key,pop2_key,test2[0],test2[1],test2[2]])+'\n')
y_out.write(header + '\t'.join([pop1_key,pop2_key,y[0],y[1],y[2]])+'\n')
tau2_1_out.write(header + '\t'.join([pop1_key,pop2_key,tau2_1[0],tau2_1[1],tau2_1[2]])+'\n')
tau2_2_out.write(header + '\t'.join([pop1_key,pop2_key,tau2_2[0],tau2_2[1],tau2_2[2]])+'\n')
tau3_1_out.write(header + '\t'.join([pop1_key,pop2_key,tau3_1[0],tau3_1[1],tau3_1[2]])+'\n')
tau3_2_out.write(header + '\t'.join([pop1_key,pop2_key,tau3_2[0],tau3_2[1],tau3_2[2]])+'\n')
B1_out.write(header + '\t'.join([pop1_key,pop2_key,B1[0],B1[1],B1[2]])+'\n')
B2_out.write(header + '\t'.join([pop1_key,pop2_key,B2[0],B2[1],B2[2]])+'\n')
U1_out.write(header + '\t'.join([pop1_key,pop2_key,U1[0],U1[1],U1[2]])+'\n')
U2_out.write(header + '\t'.join([pop1_key,pop2_key,U2[0],U2[1],U2[2]])+'\n')
V1_out.write(header + '\t'.join([pop1_key,pop2_key,V1[0],V1[1],V1[2]])+'\n')
V2_out.write(header + '\t'.join([pop1_key,pop2_key,V2[0],V2[1],V2[2]])+'\n')
T1_out.write(header + '\t'.join([pop1_key,pop2_key,T1[0],T1[1],T1[2]])+'\n')
T2_out.write(header + '\t'.join([pop1_key,pop2_key,T2[0],T2[1],T2[2]])+'\n')
J1_out.write(header + '\t'.join([pop1_key,pop2_key,J1[0],J1[1],J1[2]])+'\n')
J2_out.write(header + '\t'.join([pop1_key,pop2_key,J2[0],J2[1],J2[2]])+'\n')
m_counts_out.write(header + '\t'.join([pop1_key,pop2_key,m_counts[0],m_counts[1],m_counts[2]])+'\n')