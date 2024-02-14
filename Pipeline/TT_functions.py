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
                while l1[0:2]=='##':
                    l1=file_1.readline()
                while l2[0:2]=='##':
                    l2=file_2.readline()
                # The line that should be left is the names of all the columns, and so we can get what column the POS, QUAL and FILTER, etc. are at
                ind1_columns = l1.strip().split()
                ind2_columns = l1.strip().split()
                try:
                    ind1_POS = ind1_columns.index("POS")
                except ValueError:
                    print(f"Could not find 'POS' column in vcf file {ind1}. Please check that formatting is correct.")
                    exit(1)
                ind2_POS = ind2_columns.index("POS")
                ind1_QUAL = ind1_columns.index("QUAL")
                ind2_QUAL = ind2_columns.index("QUAL")
                ind1_FILTER = ind1_columns.index("FILTER")
                ind2_FILTER = ind2_columns.index("FILTER")
                ind1_FORMAT = ind1_columns.index("FORMAT")
                ind2_FORMAT = ind2_columns.index("FORMAT")
                # This stuff is HARD CODED so be careful of that
                ancestral_POS = 0 # HARD CODED
                ancestral_SUPPORT = 2
                sample = 0
                count = 0
                la = 'Init'
                # While loop for while the files exist and less than 100 sampels have been taken
                while l1 and l2 and la and (sample < 100):
                    count += 1
                    # Open the next lines which should have actual data
                    l1=ind_1.readline().strip().split()
                    l2=ind_2.readline().strip().split()
                    la=ancestral.readline().strip().split()
                    # Make sure that the ancestral that had to be initialised is actually existing.
                    if la == 'Init': 
                        print(f"Ancestral file can't be read, throw exception.")
                        break
                    # Get positions and loop to align positions if needed. 
                    POS_1 = int(l1[ind1_POS])
                    POS_2 = int(l2[ind2_POS])
                    POS_A = int(la[ancestral_POS])
                    while not POS_1 == POS_2 == POS_A:
                        if POS_1 == min(POS_1, POS_2, POS_A):
                            l1=ind_1.readline().strip().split()
                            POS_1 = int(l1[ind1_POS])
                        elif POS_2 == min(POS_1, POS_2, POS_A):
                            l2=ind_2.readline().strip().split()
                            POS_2 = int(l2[ind2_POS])
                        elif POS_A == min(POS_1, POS_2, POS_A):
                            la=ancestral.readline().strip().split()
                            POS_A = int(la[ancestral_POS])
                        else: break
                    # Check to make sure the positions are all the same. 
                    if not POS_1 == POS_2 == POS_A: 
                        print("Error: Files never managed to be set at the same position.")
                        break
                    # Series of checks to make sure that we can keep going
                    if la[ancestral_SUPPORT] != '3': continue
                    elif l1[ind1_QUAL] == '.' or l2[ind2_QUAL] == '.': continue
                    elif l1[ind1_FILTER] not in FILTERS or l2[ind2_FILTER] not in FILTERS: continue
                    elif bad_coverage(l1, ind1_FORMAT) or bad_coverage(l2, ind2_FORMAT): continue
                    
    return