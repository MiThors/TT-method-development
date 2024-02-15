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

def check_if_ok_and_get_var_form_TT():
    return

def check_if_ok_and_get_var_form_TTo():
    return

def get_sample_conf():
    return

def check_if_pass_coverage():
    return

def make_out_str():
    return

def bad_coverage(line_list, FORMAT_index, low, high):
    '''Function to check whether or not the coverage of a position in a genome's vcf is between the set thresholds.
    Input: the contents of a vcf line as a list, the int index of which column the FORMAT information is stored.
    Output: True if the coverage is bad, False if it passed the check.'''
    genotype_info = line_list[FORMAT_index].split(':')
    if 'GT' in genotype_info and 'DP' in genotype_info:
        depth_index = genotype_info.index('DP')
        genotype_values = line_list[FORMAT_index + 1].split(':')
        depth = int(genotype_values[depth_index])
        if depth > low and depth < high:
            return False
        else: return True
    else: return True

def get_counts_vcf(ind1, ind2, anc, low_cov, high_cov, filters):
    '''Function for getting counts from a vcf file. Opens files, checks formatting for ind1 and ind2 columns is correct, aligns positions in all files, ignores lines that do not pass filters, then adds counts to the appropriate situation.
    ind1, ind2, anc = filepaths for all files, list of one or more
    only_one_file = true or false if all chromosomes are in one file
    low_cov, high_cov = coverage thresholds for filtering, set up in main method file
    filters = list of values considered acceptable for FILTER field of vcf file
    Returns all eight count scenarios in a dictionary, keys are chromosomes, one list of counts per chromosome'''
    with gzip.open(anc,'rt',encoding='utf-8') as ancestral:
        with gzip.open(ind1, 'rt', encoding='utf-8') as file_1:
            with gzip.open(ind2, 'rt', encoding='utf-8') as file_2:
                # This is a little loop to skip past the VCF file headers for both vcf files
                l1=file_1.readline()
                l2=file_2.readline()
                la=ancestral.readline()
                while l1[0:2]=='##':
                    l1=file_1.readline()
                while l2[0:2]=='##':
                    l2=file_2.readline()
                # The line that should be left is the names of all the columns, and so we can get what column the POS, QUAL and FILTER, etc. are at
                ind1_columns = l1.strip().split()
                ind2_columns = l2.strip().split()
                anc_columns = la.strip().split()
                try: ind1_POS = ind1_columns.index("POS")
                except ValueError:
                    print(f"Could not find 'POS' column in vcf file {ind1}. Please check that formatting is correct.")
                    exit(1)
                try: ind2_POS = ind2_columns.index("POS")
                except ValueError:
                    print(f"Could not find 'POS' column in vcf file {ind2}. Please check that formatting is correct.")
                    exit(1)
                try: ind1_QUAL = ind1_columns.index("QUAL")
                except ValueError:
                    print(f"Could not find 'QUAL' column in vcf file {ind1}. Please check that formatting is correct.")
                    exit(1)
                try: ind2_QUAL = ind2_columns.index("QUAL")
                except ValueError:
                    print(f"Could not find 'QUAL' column in vcf file {ind2}. Please check that formatting is correct.")
                    exit(1)
                try: ind1_FILTER = ind1_columns.index("FILTER")
                except ValueError:
                    print(f"Could not find 'FILTER' column in vcf file {ind1}. Please check that formatting is correct.")
                    exit(1)
                try: ind2_FILTER = ind2_columns.index("FILTER")
                except ValueError:
                    print(f"Could not find 'FILTER' column in vcf file {ind2}. Please check that formatting is correct.")
                    exit(1)
                try: ind1_FORMAT = ind1_columns.index("FORMAT")
                except ValueError:
                    print(f"Could not find 'FORMAT' column in vcf file {ind1}. Please check that formatting is correct.")
                    exit(1)
                try: ind2_FORMAT = ind2_columns.index("FORMAT")
                except ValueError:
                    print(f"Could not find 'FORMAT' column in vcf file {ind2}. Please check that formatting is correct.")
                    exit(1)
                # This stuff is HARD CODED so be careful of that
                try: anc_POS = anc_columns.index("POS")
                except ValueError:
                    print(f"Could not find 'POS' column in vcf file {anc}. Please check that formatting is correct.")
                    exit(1)
                try: anc_SUPPORT = anc_columns.index("SUPPORT")
                except ValueError:
                    print(f"Could not find 'SUPPORT' column in vcf file {anc}. Please check that formatting is correct.")
                    exit(1)
                try: anc_NUCL = anc_columns.index("NUCL")
                except ValueError:
                    print(f"Could not find 'NUCL' column in vcf file {anc}. Please check that formatting is correct.")
                    exit(1)
                # While loop for while the files exist and less than 100 sampels have been taken
                while l1 and l2 and la:
                    # Open the next lines which should have actual data
                    l1=file_1.readline().strip().split()
                    l2=file_2.readline().strip().split()
                    la=ancestral.readline().strip().split()
                    # Get positions and loop to align positions if needed. 
                    POS_1 = int(l1[ind1_POS])
                    POS_2 = int(l2[ind2_POS])
                    POS_A = int(la[anc_POS])
                    while not POS_1 == POS_2 == POS_A:
                        if POS_1 == min(POS_1, POS_2, POS_A):
                            l1=file_1.readline().strip().split()
                            POS_1 = int(l1[ind1_POS])
                        elif POS_2 == min(POS_1, POS_2, POS_A):
                            l2=file_2.readline().strip().split()
                            POS_2 = int(l2[ind2_POS])
                        elif POS_A == min(POS_1, POS_2, POS_A):
                            la=ancestral.readline().strip().split()
                            POS_A = int(la[anc_POS])
                        else: break
                    # Check to make sure the positions are all the same. 
                    if not POS_1 == POS_2 == POS_A: 
                        print("Error: Files never managed to be set at the same position.")
                        exit(1)
                    # Series of checks to make sure that we can keep going
                    if la[anc_SUPPORT] != '3': continue
                    elif l1[ind1_QUAL] == '.' or l2[ind2_QUAL] == '.': continue
                    elif l1[ind1_FILTER] not in filters or l2[ind2_FILTER] not in filters: continue
                    elif bad_coverage(l1, ind1_FORMAT, low_cov, high_cov) or bad_coverage(l2, ind2_FORMAT, low_cov, high_cov): continue
    return