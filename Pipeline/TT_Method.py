'''Script for running the TT method. 
The user can tweak parameters defined at the start if they so wish.
Functions used are found in the functions.py script.
Please see the Manual (ReadMe.md) or use the command 'python3 TT_Method.py --help' for help.

Milo Thordarson: anth2886@student.uu.se'''

# Module loading
import functions # TT specific funtions from functions.py file
import argparse
import os
import multiprocessing
import sys

# Filtering parameters that can be changed by the user
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
                    help = "name of the output directory, defaults to creating one named 'TT_out_pop1_pop2' in current directory")
parser.add_argument("-c", "--counts", 
                    action="store_true",
                    help = "output a file with all counts per chromosome per window")
parser.add_argument("-d", "--depth_thresholds",
                    required = True,
                    nargs = 2,
                    type = int,
                    help = "user should define what is too low and too high depth for their given genomes, as it varies greatly sample to sample")
parser.add_argument("-w", "--window", 
                    default = 5000000,
                    type = int,
                    help = "set the window size for calculating local parameters, default is 5000000 which correspends to about 5 cM")

args = parser.parse_args()

# Turn args into informatively named variables
files_pop1 = args.pop1
files_pop2 = args.pop2
files_anc = args.ancestral
pop1_key = args.keywords[0]
pop2_key = args.keywords[1]
if args.out: out_dir = args.out
else: out_dir = "TT_out_" + pop1_key + "_" + pop2_key
print_counts = args.counts
low_coverage, high_coverage = args.depth_thresholds
win_size = args.window

# Check that filepaths are the same lengths
file_tot = len(files_pop1)
if any(len(lst) != file_tot for lst in [files_pop2, files_anc]):
    print("Error: Unequal amount of files for pop1, pop2 or ancestral.")
    print(f'Files given for pop1 are: {files_pop1}')
    print(f'Files given for pop2 are: {files_pop2}')
    print(f'Files given for ancestral states are: {files_anc}')
    sys.exit(1)

# Make output dir, will complain if it already exists, which is why this is so early in the script
os.mkdir(out_dir)

# Create iterable list for each set of files with all input parameters for multicore counting
iterables = [[files_pop1[i], files_pop2[i], files_anc[i], low_coverage, high_coverage, vcf_filters, win_size]for i in range(file_tot)]

# To avoid infinite recursion
if __name__ == '__main__':
    with multiprocessing.Pool() as pool:
        # Computes for files in parallel using CPU cores available to user
        results = pool.map(functions.get_counts_TT, iterables)
    pool.close()

# Combine results into a single dictionary of the counts and if user selected, print the counts by chromosome and window and output to file
counts = []
if print_counts:
    count_file = open(out_dir + "/" + pop1_key + pop2_key + "_TT_Counts.txt", 'w')
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

# Obtain estiamtes of the model parameters using the counts. For each parameter is the observed mean from all counts, the wbj mean and the wbj variance
[alfa1,alfa2,thetaA,mu_t1,mu_t2,mu_nu1,mu_nu2,mu_diff_t1_t2,drift1,drift2,theta1,theta2,W1ratio,W2ratio,D1,D2,P1,P2,P1_time,P2_time,Fst] = functions.get_estimates_TT(counts)

# Open all output files, print the header and estimates, close all output files
alfa1_out=open(out_dir+'/alfa1.res','w')
alfa2_out=open(out_dir+'/alfa2.res','w')
thetaA_out=open(out_dir+'/thetaA.res','w')
mu_t1_out=open(out_dir+'/mu_t1.res','w')
mu_t2_out=open(out_dir+'/mu_t2.res','w')
mu_nu1_out=open(out_dir+'/mu_nu1.res','w')
mu_nu2_out=open(out_dir+'/mu_nu2.res','w')
mu_diff_t1_t2_out=open(out_dir+'/mu_diff_t1_t2.res','w')
drift1_out=open(out_dir+'/drift1.res','w')
drift2_out=open(out_dir+'/drift2.res','w')
theta1_out=open(out_dir+'/theta1.res','w')
theta2_out=open(out_dir+'/theta2.res','w')
W1ratio_out=open(out_dir+'/W1ratio.res','w')
W2ratio_out=open(out_dir+'/W2ratio.res','w')
D1_out=open(out_dir+'/D1.res','w')
D2_out=open(out_dir+'/D2.res','w')
P1_out=open(out_dir+'/P1.res','w')
P2_out=open(out_dir+'/P2.res','w')
P1_time_out=open(out_dir+'/P1_time.res','w')
P2_time_out=open(out_dir+'/P2_time.res','w')
Fst_out=open(out_dir+'/Fst.res','w')

header= '\t'.join(['pop1','pop2','obs_mean','wbj_mean','wbj_var'])+'\n'

alfa1_out.write(header + '\t'.join([pop1_key,pop2_key,alfa1[0],alfa1[1],alfa1[2]])+'\n')
alfa2_out.write(header + '\t'.join([pop1_key,pop2_key,alfa2[0],alfa2[1],alfa2[2]])+'\n')
thetaA_out.write(header + '\t'.join([pop1_key,pop2_key,thetaA[0],thetaA[1],thetaA[2]])+'\n')
mu_t1_out.write(header + '\t'.join([pop1_key,pop2_key,mu_t1[0],mu_t1[1],mu_t1[2]])+'\n')
mu_t2_out.write(header + '\t'.join([pop1_key,pop2_key,mu_t2[0],mu_t2[1],mu_t2[2]])+'\n')
mu_nu1_out.write(header + '\t'.join([pop1_key,pop2_key,mu_nu1[0],mu_nu1[1],mu_nu1[2]])+'\n')
mu_nu2_out.write(header + '\t'.join([pop1_key,pop2_key,mu_nu2[0],mu_nu2[1],mu_nu2[2]])+'\n')
mu_diff_t1_t2_out.write(header + '\t'.join([pop1_key,pop2_key,mu_diff_t1_t2[0],mu_diff_t1_t2[1],mu_diff_t1_t2[2]])+'\n')
drift1_out.write(header + '\t'.join([pop1_key,pop2_key,drift1[0],drift1[1],drift1[2]])+'\n')
drift2_out.write(header + '\t'.join([pop1_key,pop2_key,drift2[0],drift2[1],drift2[2]])+'\n')
theta1_out.write(header + '\t'.join([pop1_key,pop2_key,theta1[0],theta1[1],theta1[2]])+'\n')
theta2_out.write(header + '\t'.join([pop1_key,pop2_key,theta2[0],theta2[1],theta2[2]])+'\n')
W1ratio_out.write(header + '\t'.join([pop1_key,pop2_key,W1ratio[0],W1ratio[1],W1ratio[2]])+'\n')
W2ratio_out.write(header + '\t'.join([pop1_key,pop2_key,W2ratio[0],W2ratio[1],W2ratio[2]])+'\n')
D1_out.write(header + '\t'.join([pop1_key,pop2_key,D1[0],D1[1],D1[2]])+'\n')
D2_out.write(header + '\t'.join([pop1_key,pop2_key,D2[0],D2[1],D2[2]])+'\n')
P1_out.write(header + '\t'.join([pop1_key,pop2_key,P1[0],P1[1],P1[2]])+'\n')
P2_out.write(header + '\t'.join([pop1_key,pop2_key,P2[0],P2[1],P2[2]])+'\n')
P1_time_out.write(header + '\t'.join([pop1_key,pop2_key,P1_time[0],P1_time[1],P1_time[2]])+'\n')
P2_time_out.write(header + '\t'.join([pop1_key,pop2_key,P2_time[0],P2_time[1],P2_time[2]])+'\n')
Fst_out.write(header + '\t'.join([pop1_key,pop2_key,Fst[0],Fst[1],Fst[2]])+'\n')

alfa1_out.close()
alfa2_out.close()
thetaA_out.close()
mu_t1_out.close()
mu_t2_out.close()
mu_nu1_out.close()
mu_nu2_out.close()
mu_diff_t1_t2_out.close()
drift1_out.close()
drift2_out.close()
theta1_out.close()
theta2_out.close()
W1ratio_out.close()
W2ratio_out.close()
D1_out.close()
D2_out.close()
P1_out.close()
P2_out.close()
P1_time_out.close()
P2_time_out.close()
Fst_out.close()
