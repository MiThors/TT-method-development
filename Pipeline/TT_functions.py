'''A script contianing all the defined functions for using both TT and TTo methods.

Milo Thordarson: anth2886@student.uu.se'''

import gzip

def get_var_form():
    return

def check_if_missingness():
    return

def parse_var_genotypes():
    return

def get_genotype():
    
    return

def orient_and_get_count():
    return

def check_if_ok_and_get_var_form_TTo():
    return

def get_sample_conf():
    return

def check_if_pass_coverage():
    return

def make_out_str():
    return

def bad_coverage(depth, low, high):
    '''Function to check whether or not the coverage of a position in a genome's vcf is between the set thresholds.
    Input: the contents of a vcf line as a list, the int index of which column the FORMAT information is stored.
    Output: True if the coverage is bad, False if it passed the check.'''
    if depth == '.' : 
        return True # Coverage is bad if undefined
    else: depth = int(depth)
    if depth > low and depth < high:
        return False
    else: return True

def get_configuration_index(nucl_A, genotype_1, genotype_2, ref_1, ref_2, alt_1, alt_2):
    if alt_1 == '.' and alt_2 == '.':
        if ref_1 == ref_2 == nucl_A: return 0 # m0 = O_{0,0}
        elif ref_1 != nucl_A: return 1 # m1 = O_{1, 0}
        elif ref_2 != nucl_A: return 2 # m2 = O_{0, 1}
        

    

def get_counts_vcf(pop1, pop2, anc, low_cov, high_cov, filters):
    '''Function for getting counts from a vcf file. Opens files, checks formatting for pop1 and pop2 columns is correct, aligns positions in all files, ignores lines that do not pass filters, then adds counts to the appropriate situation.
    pop1, pop2, anc = filepaths for all files, list of one or more
    only_one_file = true or false if all chromosomes are in one file
    low_cov, high_cov = coverage thresholds for filtering, set up in main method file
    filters = list of values considered acceptable for FILTER field of vcf file
    Returns all eight count scenarios in a dictionary, keys are chromosomes, one list of counts per chromosome'''
    # Variable initialisation
    nucl = ['A','C','G','T']
    nt_set = set(nucl)
    out_dict = {}
    # Opening the files
    with gzip.open(anc,'rt',encoding='utf-8') as ancestral:
        with gzip.open(pop1, 'rt', encoding='utf-8') as file_1:
            with gzip.open(pop2, 'rt', encoding='utf-8') as file_2:
                # This is a little loop to skip past the VCF file headers for both vcf files
                l1 = file_1.readline()
                l2 = file_2.readline()
                la = ancestral.readline()
                while l1[0:2] == '##':
                    l1 = file_1.readline()
                while l2[0:2] == '##':
                    l2 = file_2.readline()
                
                # The line that should be left is the names of all the columns, and so we can get what column the POS, QUAL and FILTER, etc. are at
                pop1_columns = l1.strip().split()
                pop2_columns = l2.strip().split()
                anc_columns = la.strip().split()
                try: 
                    # Indexing for VCF format so there are no magic numbers, and as a test of correct file format given
                    chrom_1_ind, chrom_2_ind = pop1_columns.index("#CHROM"), pop2_columns.index("#CHROM")
                    pos_1_ind, pos_2_ind = pop1_columns.index("POS"), pop2_columns.index("POS")
                    ref_1_ind, ref_2_ind = pop1_columns.index("REF"), pop2_columns.index("REF")
                    alt_1_ind, alt_2_ind = pop1_columns.index("ALT"), pop2_columns.index("ALT")
                    qual_1_ind, qual_2_ind = pop1_columns.index("QUAL"), pop2_columns.index("QUAL")
                    filter_1_ind, filter_2_ind = pop1_columns.index("FILTER"), pop2_columns.index("FILTER")
                    format_1_ind, format_2_ind = pop1_columns.index("FORMAT"), pop2_columns.index("FORMAT")
                    pos_A_ind = anc_columns.index("POS")
                    nucl_A_ind = anc_columns.index("NUCL")
                    print()
                except ValueError:
                    print(f"Could not find all columns in in vcf files {pop1} or {pop2}, or all columns in ancestral file. Please check that formatting is correct.")
                    exit(1)
                
                # While loop for to align positions if needed
                while l1 and l2 and la:
                    l1 = file_1.readline().strip().split()
                    l2 = file_2.readline().strip().split()
                    la = ancestral.readline().strip().split()
                    if not l1 or not l2 or not la : break
                    pos_1 = int(l1[pos_1_ind])
                    pos_2 = int(l2[pos_2_ind])
                    pos_A = int(la[pos_A_ind])
                    while not pos_1 == pos_2 == pos_A:
                        if pos_1 == min(pos_1, pos_2, pos_A):
                            l1 = file_1.readline().strip().split()
                            if l1: pos_1 = int(l1[pos_1_ind])
                            else: break
                        elif pos_2 == min(pos_1, pos_2, pos_A):
                            l2 = file_2.readline().strip().split()
                            if l2: pos_2 = int(l2[pos_2_ind])
                            else: break
                        elif pos_A == min(pos_1, pos_2, pos_A):
                            la = ancestral.readline().strip().split()
                            if la: pos_A = int(la[pos_A_ind])
                            else: break
                    # Check to make sure the positions are all the same. 
                    if not pos_1 == pos_2 == pos_A: 
                        print(f"Error: Files never managed to be reach at the same position, {anc} ended at {pos_A}, {pop1} at {pos_1}, and {pop2} at {pos_2}. Please check that correct files are being compared, or file formatting.")
                        exit(1)
                    
                    # Define variables
                    nucl_A = la[nucl_A_ind]
                    ref_1, ref_2 = l1[ref_1_ind], l2[ref_2_ind]
                    alt_1, alt_2 = l1[alt_1_ind], l2[alt_2_ind]
                    chrom_1, chrom_2 = l1[chrom_1_ind], l2[chrom_2_ind]
                    
                    # Series of quality and assumption checks to make sure that we can keep going
                    if chrom_1 != chrom_2: continue 
                    elif l1[qual_1_ind] == '.' or l2[qual_2_ind] == '.' : continue
                    elif l1[filter_1_ind] not in filters or l2[filter_2_ind] not in filters : continue
                    elif nucl_A not in nucl: continue # If the ancient nucleotide is not resolved, we skip
                    elif len(alt_1) > 1 or len(alt_2) > 1: continue # We skip multiallelic sites as they violate assumptions
                    elif len(set([nucl_A,ref_1,ref_2,alt_1,alt_2]).difference('.')) > 2: continue # Another check for multiallelic site
                    elif not set([nucl_A,ref_1,ref_2,alt_1,alt_2]).difference('.').issubset(nt_set): continue # All nucleotides should be A, T, C or G.
                    
                    # If it passes these checks, we get the genotype and coverage and check if they are also appropriate
                    l1_format_info = l1[format_1_ind].split(':')
                    l2_format_info = l2[format_2_ind].split(':')
                    try:
                        genotype_1_ind, genotype_2_ind = l1_format_info.index("GT"), l2_format_info.index("GT")
                        coverage_1_ind, coverage_2_ind = l1_format_info.index("DP"), l2_format_info.index("DP")
                    except ValueError:
                        continue
                    
                    l1_genotype_info, l2_genotype_info,  = l1[format_1_ind + 1].split(':'), l2[format_2_ind + 1].split(':')
                    coverage_1, coverage_2 = l1_genotype_info[coverage_1_ind], l2_genotype_info[coverage_2_ind]
                    genotype_1, genotype_2 = l1_genotype_info[genotype_1_ind], l2_genotype_info[genotype_2_ind]
                    
                    if bad_coverage(coverage_1, low_cov, high_cov) or bad_coverage(coverage_2, low_cov, high_cov) : continue # Check the coverage is within acceptable thresholds
                    elif genotype_1.count('.') > 0 or genotype_2.count('.') > 0: continue # Check if genotypes are undefined
                    # Check if current chromosome exists in the dict already, if not add another key for that
                    if chrom_1 not in out_dict: out_dict.update({chrom_1 : [0, 0, 0, 0, 0, 0, 0, 0, 0]})
                    # Get the type of sample configuration, represented as the index of m0, m1, ... m8
                    configuration_index = get_configuration_index(nucl_A, genotype_1, genotype_2, ref_1, ref_2, alt_1, alt_2)
                    # Add one count to the relevant chromosome and configuration count
                    out_dict[chrom_1][configuration_index] += 1
    if out_dict:
        return out_dict
    else:
        print("Error: It seems that every position in files {pop1} and {pop2} failed all checks and no counts were generated for these files. Please check file formatting or whether all positions truly violate assumptions.")
        exit(1)
