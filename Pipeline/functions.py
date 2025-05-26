'''A script contianing all the defined functions for using both TT and TTo methods.

Milo Thordarson: anth2886@student.uu.se'''

import gzip
import sys
from math import log
import traceback

def get_indexes(columns_list):
    '''Function to get the indexes for the relevant columns in a vcf file. This avoids magic numbers and makes code easier to understand, and additionally helps to catch formatting issues. 
    In: the header line of a vcf, split on whitespace into a list
    Out: all 7 indexes as ints'''
    chrom_ind = columns_list.index("#CHROM")
    pos_ind = columns_list.index("POS")
    ref_ind = columns_list.index("REF")
    alt_ind = columns_list.index("ALT")
    qual_ind = columns_list.index("QUAL")
    filter_ind = columns_list.index("FILTER")
    format_ind = columns_list.index("FORMAT")
    return chrom_ind, pos_ind, ref_ind, alt_ind, qual_ind, filter_ind, format_ind

def quality_and_filter_check(quality,filter,acceptable_filters):
    '''Function simply to combine checking quality and filters columns in vcf as this needs to be repeated but is never changed much.
    Input: the quality column of a site, the filters column of a site, acceptable filters defined by user in TT or TTo method python script.
    Output: False if checks are all failed, True if it passed all checks.'''
    if '.' == quality: return False 
    elif filter not in acceptable_filters: return False
    else: return True

def bad_coverage(depth, low, high):
    '''Function to check whether or not the coverage of a position in a genome's vcf is between the set thresholds.
    Input: the depth of the genotype as a stirng obtained from the vcf, and the low and high thresholds for coverage as ints
    Output: True if the coverage is bad, False if it passed the check.'''
    if depth == '.' : 
        return True # Coverage is bad if undefined
    else: depth = int(depth)
    if depth > low and depth < high:
        return False
    else: return True

def passed_conditioning(outgroup_genotype, ref, alt, nucl_A):
    '''Function to check if a derived allele is found in the outgroup for TTo conditional counting.
    Input: outgroup genotype as a string, reference nucleotide for outgroup, alternate nucleotide for outgroup and ancestral nucleotide all as strings
    Output: True if no derived alleles are found in the outgroup, False if there is at least one derived allele'''
    # No derived
    if alt == '.' and ref == nucl_A: derived = 0
    # Both derived
    elif alt == '.' and ref != nucl_A: derived = 2
    # One or more derived, counting alternate nucleotide
    elif ref == nucl_A: derived = outgroup_genotype.count("1")
    # Reference is the derived nucleotide, one or more derived
    else: derived = outgroup_genotype.count("0")
    if derived > 0 : return True
    else: return False

def get_configuration_index(nucl_A, genotype_1, genotype_2, ref_1, ref_2, alt_1, alt_2):
    '''Function to get what configuration the genotypes of pop1 and pop2 are in, based on the numbers of derived and ancestral alleles. 
    Input: Ancestral nucleotide, genotypes of both populations, reference and altnerative nucleotides for both populations, all in strings
    Output: int values corresponding to the configuration'''
    # No derived in pop1
    if alt_1 == '.' and ref_1 == nucl_A: pop1_derived = 0
    # Both derived in pop1
    elif alt_1 == '.' and ref_1 != nucl_A: pop1_derived = 2
    # One or more derived, counting alternate nucleotide in pop1
    elif ref_1 == nucl_A: pop1_derived = genotype_1.count("1")
    # Reference is the derived nucleotide, one or more derived, counting reference nucleotide in pop1
    else: pop1_derived = genotype_1.count("0")
    # Same as above for pop 2
    if alt_2 == '.' and ref_2 == nucl_A: pop2_derived = 0
    elif alt_2 == '.' and ref_2 != nucl_A: pop2_derived = 2
    elif ref_2 == nucl_A: pop2_derived = genotype_2.count("1")
    else: pop2_derived = genotype_2.count("0")

    if pop1_derived == 0 and pop2_derived == 0: return 0 # m0 = O{0,0}
    elif pop1_derived == 1 and pop2_derived == 0: return 1 # m1 = O{1,0}
    elif pop1_derived == 0 and pop2_derived == 1: return 2 # m2 = O{0,1}
    elif pop1_derived == 2 and pop2_derived == 0: return 3 # m3 = O{2,0}
    elif pop1_derived == 0 and pop2_derived == 2: return 4 # m4 = O_{0,2}
    elif pop1_derived == 1 and pop2_derived == 1: return 5 # m5 = O{1,1}
    elif pop1_derived == 2 and pop2_derived == 1: return 6 # m6 = O{2,1}
    elif pop1_derived == 1 and pop2_derived == 2: return 7 # m7 = O{1,2}
    elif pop1_derived == 2 and pop2_derived == 2: return 8 # m8 = O{2,2}
    else: 
        print("Error: it seems genotype counting was able to obtain values other than (0, 1, 2) for one or both populations. According to our checks, that should not be possible, please check genotype information in vcfs.")
        sys.exit(1)

def get_TTo_configuration_index(nucl_A, genotype_1, genotype_2, ref_1, ref_2, alt_1, alt_2, ref_og, alt_og):
    '''Function to get what configuration the genotypes of pop1 and pop2 are in, based on the numbers of derived and ancestral alleles. 
    Input: Ancestral nucleotide, genotypes of both populations, reference and altnerative nucleotides for both populations, all in strings
    Output: int values corresponding to the configuration'''
    # No derived in pop1
    nucl_set = set([nucl_A, ref_1, ref_2, alt_1, alt_2, ref_og, alt_og]).difference('.')
    if len(nucl_set) == 1: #Case where using the outgroup, there is no variation
        if ref_1 == nucl_A: return 0
        else: return 8
    if len(nucl_set) == 2: 
        if alt_1=='.' and alt_2=='.' and alt_og=='.': #Case where when using the outgroup, there is variation, but none of the populations have a derived
            if ref_1 == nucl_A: return 0
            else: return 8
    if alt_1 == '.' and ref_1 == nucl_A: pop1_derived = 0
    # Both derived in pop1
    elif alt_1 == '.' and ref_1 != nucl_A: pop1_derived = 2
    # One or more derived, counting alternate nucleotide in pop1
    elif ref_1 == nucl_A: pop1_derived = genotype_1.count("1")
    # Reference is the derived nucleotide, one or more derived, counting reference nucleotide in pop1
    else: pop1_derived = genotype_1.count("0")
    # Same as above for pop 2
    if alt_2 == '.' and ref_2 == nucl_A: pop2_derived = 0
    elif alt_2 == '.' and ref_2 != nucl_A: pop2_derived = 2
    elif ref_2 == nucl_A: pop2_derived = genotype_2.count("1")
    else: pop2_derived = genotype_2.count("0")

    if pop1_derived == 0 and pop2_derived == 0: return 0 # m0 = O{0,0}
    elif pop1_derived == 1 and pop2_derived == 0: return 1 # m1 = O{1,0}
    elif pop1_derived == 0 and pop2_derived == 1: return 2 # m2 = O{0,1}
    elif pop1_derived == 2 and pop2_derived == 0: return 3 # m3 = O{2,0}
    elif pop1_derived == 0 and pop2_derived == 2: return 4 # m4 = O_{0,2}
    elif pop1_derived == 1 and pop2_derived == 1: return 5 # m5 = O{1,1}
    elif pop1_derived == 2 and pop2_derived == 1: return 6 # m6 = O{2,1}
    elif pop1_derived == 1 and pop2_derived == 2: return 7 # m7 = O{1,2}
    elif pop1_derived == 2 and pop2_derived == 2: return 8 # m8 = O{2,2}
    else: 
        print("Error: it seems genotype counting was able to obtain values other than (0, 1, 2) for one or both populations. According to our checks, that should not be possible, please check genotype information in vcfs.")
        sys.exit(1)

def get_counts_TT(iterable):
    '''Function for getting counts from a vcf file. Opens files, checks formatting for pop1 and pop2 columns is correct, aligns positions in all files, ignores lines that do not pass filters, then adds counts to the appropriate situation.
    Input:
    pop1, pop2, anc = filepaths for all files, list of one or more
    low_cov, high_cov = coverage thresholds for filtering, set up in main method file
    filters = list of values considered acceptable for FILTER field of vcf file
    Output: 
    out_dict = dictionary, keys are chromosomes, values are list of count lists per window and tuple of start and end position of each window
    '''
    # Variable initialisation
    pop1 = iterable[0]
    pop2 = iterable[1]
    anc = iterable[2]
    low_cov = int(iterable[3])
    high_cov = int(iterable[4])
    filters = iterable[5]
    window_size = int(iterable[6])
    nucl = ['A','C','G','T']
    nt_set = set(nucl)
    window_step = 0
    out_dict = {}
    local_count = []
    win_start = 0
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
                    chrom_1_ind, pos_1_ind, ref_1_ind, alt_1_ind, qual_1_ind, filter_1_ind, format_1_ind = get_indexes(pop1_columns)
                    chrom_2_ind, pos_2_ind, ref_2_ind, alt_2_ind, qual_2_ind, filter_2_ind, format_2_ind = get_indexes(pop2_columns)
                    pos_A_ind = anc_columns.index("POS")
                    nucl_A_ind = anc_columns.index("NUCL")
                except ValueError:
                    print(f"Error: Could not find all columns in in vcf files {pop1} or {pop2}, or all columns in ancestral file. Please check that formatting is correct.")
                    sys.exit(1)
                # Then loop for while file exists to do the counting
                while l1 and l2 and la:
                    l1 = file_1.readline().strip().split()
                    l2 = file_2.readline().strip().split()
                    la = ancestral.readline().strip().split()
                    # Case for the end of the file, update dictionaries and break
                    if not l1 or not l2 or not la: 
                        if local_count: out_dict[current_chrom][0].append(local_count)
                        out_dict[current_chrom][1].append((win_start, current_pos))
                        break
                    pos_1 = int(l1[pos_1_ind])
                    pos_2 = int(l2[pos_2_ind])
                    pos_A = int(la[pos_A_ind])
                    # If need to align positions
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
                    # Check for end of file, then to make sure the positions are all the same and throw and error that they never managed to align
                    if not l1 or not l2 or not la: 
                        if local_count: out_dict[current_chrom][0].append(local_count)
                        out_dict[current_chrom][1].append((win_start, current_pos))
                        break
                    if not pos_1 == pos_2 == pos_A: 
                        print(f"Error: Files never managed to be reach at the same position, {anc} ended at {pos_A}, {pop1} at {pos_1}, and {pop2} at {pos_2}. Please check that correct files are being compared, or file formatting.")
                        sys.exit(1)
                    # Once those checks are passed, start getting genotype and chromosome information
                    nucl_A = la[nucl_A_ind]
                    ref_1, ref_2 = l1[ref_1_ind], l2[ref_2_ind]
                    alt_1, alt_2 = l1[alt_1_ind], l2[alt_2_ind]
                    chrom_1, chrom_2 = l1[chrom_1_ind], l2[chrom_2_ind]
                    # Check to make sure the chromosomes are the same
                    if chrom_1 != chrom_2: 
                        print(f'Error: Files at same positions, but chromosomes in {pop1} and {pop2} are not the same, please check file formatting.')
                        sys.exit(1)
                    # Check if current chromosome exists in the dict already, if not add another key for that and if counting was done on a previous chromosome, update dictionaries
                    if chrom_1 not in out_dict: 
                        out_dict.update({chrom_1 : [[],[]]})
                        if local_count: 
                            out_dict[current_chrom][0].append(local_count)
                            out_dict[current_chrom][1].append((win_start, current_pos))
                        local_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                        window_step = 0
                        win_start = 0
                    # If a position in the file is skipped, this while loop will increse the steps of the window to account for them in the genome, and update the window if it reaches the end
                    while pos_1 - win_start != window_step:
                        window_step += 1
                        if window_step >= window_size:
                            if not local_count: local_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                            out_dict[chrom_1][0].append(local_count)
                            out_dict[chrom_1][1].append((win_start, win_start+window_size-1))
                            local_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                            window_step = 0
                            win_start += window_size
                    # If at the end of a window, update output directories and start a new window
                    if window_step >= window_size:
                        out_dict[chrom_1][0].append(local_count)
                        out_dict[chrom_1][1].append((win_start, current_pos))
                        local_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                        window_step = 0
                        win_start = pos_1
                    # To keep track of the previous chromosome and position for updating the directories correctly
                    current_chrom = chrom_1
                    window_step += 1
                    current_pos = pos_1
                    # Series of quality and assumption checks to make sure that we can keep going 
                    #if quality_and_filter_check(l1[qual_1_ind], l1[filter_1_ind], filters) == False or quality_and_filter_check(l2[qual_2_ind], l2[filter_2_ind], filters) == False: continue
                    if not '.' in [l1[qual_1_ind], l2[qual_2_ind]]:
                                if (l1[filter_1_ind] in ['PASS','.']) and (l2[filter_2_ind] in ['PASS','.']):
                                    if not (low_cov, high_cov) == (0,0):
                                        [coverage1,genotype1]=get_genotype(l1[9:])
                                        [coverage2,genotype2]=get_genotype(l2[9:])
                                        if not check_if_pass_coverage(coverage1,low_cov,high_cov) or not check_if_pass_coverage(coverage2,low_cov,high_cov): continue
                                    if nucl_A in nucl:
                                        var_form=check_if_ok_and_get_var_form_TT(nucl_A,ref_1,ref_2,alt_1,alt_2)
                                        if not var_form=='':
                                            if var_form=='OK_NO_VARIATION':
                                                if ref_1==nucl_A:
                                                    local_count[0] += 1
                                                else:
                                                    local_count[8] += 1
                                            elif var_form=='OK_POLY':
                                                der_count1=orient_and_get_count(genotype1,ref_1,alt_1,nucl_A)
                                                der_count2=orient_and_get_count(genotype2,ref_1,alt_2,nucl_A)
                                                local_count[get_sample_conf(der_count1,der_count2)] += 1
    # Check if the output exists, and if so return it, else there has been an error       
    if out_dict:
        return out_dict
    else:
        print(f"Error: It seems that every position in files {pop1} and {pop2} failed all checks and no counts were generated for these files. Please check file formatting or whether all positions truly violate assumptions.")
        sys.exit(1)

def get_genotype(a_list):
    b_geno=''
    coverage=0
    for x in a_list:
        d=x.split(':')
        if d[0]=='./.':
            return [0,'']
        b_geno=d[0]
        if len(d)>2:
            if d[2]=='.':
                return [0,'']
            coverage += int(d[2])
        else:
            coverage+=int(d[1])
    return [coverage,b_geno]

def check_if_pass_coverage(a_coverage,LOW_COV_THRESH,HIGH_COV_THRESH):
    if (a_coverage>LOW_COV_THRESH and a_coverage<HIGH_COV_THRESH):
        return 1
    else:
        return 0

def check_if_ok_and_get_var_form(anc_nt,ref_nt1,ref_nt2,ref_nt_ogrp,alt_nt1,alt_nt2,alt_nt_ogrp):
    set1=set([anc_nt,ref_nt1,ref_nt2,ref_nt_ogrp,alt_nt1,alt_nt2,alt_nt_ogrp]).difference('.')
    NUCL=['A','C','G','T']
    nt_set=set(NUCL)
    if set1.issubset(nt_set):
        if len(set1)==1:
            return 'OK_NO_VARIATION'
        if len(set1)==2:
            if alt_nt1=='.' and alt_nt2=='.' and alt_nt_ogrp=='.':
                return 'OK_NO_VARIATION'
            else:
                return 'OK_POLY'
    return ''

def orient_and_get_count(genotype,ref_nt,alt_nt,anc_nt):
    if genotype.count('2')>0:
        print('Genotype contains 2 derived.')
    if alt_nt=='.':
        if ref_nt==anc_nt:
            return 0
        else:
            return 2
    else:
        if ref_nt==anc_nt:
            return genotype.count('1')
        else:
            return genotype.count('0')

def get_sample_conf(der_count1,der_count2):
    if der_count1+der_count2==0:
        return 0
    if der_count1+der_count2==4:
        return 0
    if der_count1+der_count2==1:
        if der_count1==1:
            return 1
        else:
            return 2
    if der_count1+der_count2==2:
        if der_count1==2:
            return 3
        if der_count2==2:
            return 4
        else:
            return 5
    if der_count1+der_count2==3:
        if der_count1==2:
            return 6
        else:
            return 7
    print('THIS SHOULD NEVER HAPPEN')
    return 100

def check_if_ok_and_get_var_form_TT(anc_nt,ref_nt1,ref_nt2,alt_nt1,alt_nt2):
    set1=set([anc_nt,ref_nt1,ref_nt2,alt_nt1,alt_nt2]).difference('.')
    NUCL=['A','C','G','T']
    nt_set=set(NUCL)
    if set1.issubset(nt_set):
        if len(set1)==1:
            return 'OK_NO_VARIATION'
        if len(set1)==2:
            if alt_nt1=='.' and alt_nt2=='.':
                return 'OK_NO_VARIATION'
            else:
                return 'OK_POLY'
    return ''

def get_counts_TTo(iterable):
    '''Function for getting counts from a vcf file for the TTo method, TT counts are included if no TT count file was given. Opens files, checks formatting for pop1 and pop2 columns is correct, aligns positions in all files, ignores lines that do not pass filters, then adds counts to the appropriate situation.
    pop1, pop2, outgroup, anc = filepaths for all files, list of one or more
    low_cov, high_cov = coverage thresholds for filtering, set up in main method file
    filters = list of values considered acceptable for FILTER field of vcf file
    Returns all eight count scenarios in a dictionary, keys are chromosomes, one list of counts per chromosome'''
    # Variable initialisation
    pop1 = iterable[0]
    pop2 = iterable[1]
    outgroup = iterable[2]
    anc = iterable[3]
    low_cov = int(iterable[4])
    high_cov = int(iterable[5])
    filters = iterable[6]
    window_size = int(iterable[7])
    if iterable[8]: count_TT = False
    else: count_TT = True
    nucl = ['A','C','G','T']
    nt_set = set(nucl)
    out_dict = {}
    local_count = []
    out_dict_TT = {}
    local_count_TT = []
    win_start = 0
    # Opening the files
    with gzip.open(anc,'rt') as ancestral:
        with gzip.open(pop1, 'rt') as file_1:
            with gzip.open(pop2, 'rt') as file_2:
                with gzip.open(outgroup, 'rt') as file_og:
                    # This is a little loop to skip past the VCF file headers for both vcf files
                    l1 = file_1.readline()
                    l2 = file_2.readline()
                    lo = file_og.readline()
                    la = ancestral.readline()
                    while l1[0:2] == '##':
                        l1 = file_1.readline()
                    while l2[0:2] == '##':
                        l2 = file_2.readline()
                    while lo[0:2] == '##':
                        lo = file_og.readline()
                    # The line that should be left is the names of all the columns, and so we can get what column the POS, QUAL and FILTER, etc. are at
                    pop1_columns = l1.strip().split()
                    pop2_columns = l2.strip().split()
                    outgroup_columns = lo.strip().split()
                    anc_columns = la.strip().split()
                    try: 
                        # Indexing for VCF format so there are no magic numbers, and as a test of correct file format given
                        chrom_1_ind, pos_1_ind, ref_1_ind, alt_1_ind, qual_1_ind, filter_1_ind, format_1_ind = get_indexes(pop1_columns)
                        chrom_2_ind, pos_2_ind, ref_2_ind, alt_2_ind, qual_2_ind, filter_2_ind, format_2_ind = get_indexes(pop2_columns)
                        chrom_OG_ind, pos_OG_ind, ref_OG_ind, alt_OG_ind, qual_OG_ind, filter_OG_ind, format_OG_ind = get_indexes(outgroup_columns)
                        pos_A_ind = anc_columns.index("POS")
                        nucl_A_ind = anc_columns.index("NUCL")
                    except ValueError:
                        print(f"Error: Could not find all columns in in one or multiple vcf files {pop1}, {pop2} or {file_og}, or all columns in ancestral file. Please check that formatting is correct.")
                        sys.exit(1)
                    while l1 and l2 and la and lo:
                        l1 = file_1.readline().strip().split()
                        l2 = file_2.readline().strip().split()
                        lo = file_og.readline().strip().split()
                        la = ancestral.readline().strip().split()
                        if not l1 or not l2 or not la or not lo: 
                            if local_count:
                                out_dict[current_chrom][0].append(local_count)
                                out_dict[current_chrom][1].append((win_start, current_pos))
                                if count_TT:
                                    out_dict_TT[current_chrom][0].append(local_count_TT)
                                    out_dict_TT[current_chrom][1].append((win_start, current_pos))
                            break
                        pos_1 = int(l1[pos_1_ind])
                        pos_2 = int(l2[pos_2_ind])
                        pos_OG = int(lo[pos_OG_ind])
                        pos_A = int(la[pos_A_ind])
                        while not pos_1 == pos_2 == pos_A == pos_OG:
                            if pos_1 == min(pos_1, pos_2, pos_OG, pos_A):
                                l1 = file_1.readline().strip().split()
                                if l1: pos_1 = int(l1[pos_1_ind])
                                else: break
                            elif pos_2 == min(pos_1, pos_2, pos_OG, pos_A):
                                l2 = file_2.readline().strip().split()
                                if l2: pos_2 = int(l2[pos_2_ind])
                                else: break
                            elif pos_OG == min(pos_1, pos_2, pos_OG, pos_A):
                                lo = file_og.readline().strip().split()
                                if lo: pos_OG = int(lo[pos_OG_ind])
                                else: break
                            elif pos_A == min(pos_1, pos_2, pos_OG, pos_A):
                                la = ancestral.readline().strip().split()
                                if la: pos_A = int(la[pos_A_ind])
                                else: break
                        # Check for end of file
                        if not l1 or not l2 or not la or not lo: 
                            if local_count:
                                out_dict[current_chrom][0].append(local_count)
                                out_dict[current_chrom][1].append((win_start, current_pos))
                                if count_TT:
                                    out_dict_TT[current_chrom][0].append(local_count_TT)
                                    out_dict_TT[current_chrom][1].append((win_start, current_pos))
                            break
                        # Check to make sure the positions are all the same. 
                        if not pos_1 == pos_2 == pos_A == pos_OG: 
                            print(f"Error: Files never managed to be reach at the same position, {anc} ended at {pos_A}, {pop1} at {pos_1}, {pop2} at {pos_2} and {outgroup} at {pos_OG}. Please check that correct files are  being compared, or file formatting.")
                            sys.exit(1)
                        # Define variables
                        nucl_A = la[nucl_A_ind]
                        ref_1, ref_2, ref_OG = l1[ref_1_ind], l2[ref_2_ind], lo[ref_OG_ind]
                        alt_1, alt_2, alt_OG = l1[alt_1_ind], l2[alt_2_ind], lo[alt_OG_ind]
                        chrom_1, chrom_2, chrom_OG = l1[chrom_1_ind], l2[chrom_2_ind], lo[chrom_OG_ind]

                        if not chrom_1 == chrom_2 == chrom_OG:
                            print(f'Error: Files at same positions, but chromosomes in {pop1}, {pop2} and {file_og} are not the same, please check file formatting.')
                            sys.exit(1)
                        
                        # Check if current chromosome exists in the dict already, if not add another key for that
                        if chrom_1 not in out_dict: 
                            out_dict.update({chrom_1: [[],[]]})
                            out_dict_TT.update({chrom_1: [[],[]]})
                            if local_count: 
                                out_dict[chrom_1][0].append(local_count)
                                out_dict[chrom_1][1].append((win_start, win_start + window_size))
                                if count_TT:
                                    out_dict_TT[chrom_1][0].append(local_count_TT)
                                    out_dict_TT[chrom_1][1].append((win_start, win_start + window_size))
                            local_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                            if count_TT: local_count_TT = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                            win_start = 0
                        # If a position in the file is skipped, this while loop will increse the steps of the window to account for them in the genome, and update the window if it reaches the end
                        while pos_1 > (win_start + window_size):
                            out_dict[chrom_1][0].append(local_count)
                            out_dict[chrom_1][1].append((win_start, win_start + window_size))
                            local_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                            if count_TT:
                                out_dict_TT[chrom_1][0].append(local_count_TT)
                                out_dict_TT[chrom_1][1].append((win_start, win_start + window_size))
                                local_count_TT = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                            win_start += window_size
                        
                        current_chrom = chrom_1
                        current_pos = pos_1

                        if count_TT:
                            # Series of quality and assumption checks to make sure that we can keep going 
                            if not '.' in [l1[qual_1_ind], l2[qual_2_ind]]:
                                if (l1[filter_1_ind] in ['PASS','.']) and (l2[filter_2_ind] in ['PASS','.']):
                                    if not (low_cov, high_cov) == (0,0):
                                        [coverage1,genotype1]=get_genotype(l1[9:])
                                        [coverage2,genotype2]=get_genotype(l2[9:])
                                        if not check_if_pass_coverage(coverage1,low_cov,high_cov) or not check_if_pass_coverage(coverage2,low_cov,high_cov): continue
                                    if nucl_A in nucl:
                                        var_form=check_if_ok_and_get_var_form_TT(nucl_A,ref_1,ref_2,alt_1,alt_2)
                                        if not var_form=='':
                                            if var_form=='OK_NO_VARIATION':
                                                if ref_1==nucl_A:
                                                    local_count_TT[0] += 1
                                                else:
                                                    local_count_TT[8] += 1
                                            elif var_form=='OK_POLY':
                                                der_count1=orient_and_get_count(genotype1,ref_1,alt_1,nucl_A)
                                                der_count2=orient_and_get_count(genotype2,ref_1,alt_2,nucl_A)
                                                local_count_TT[get_sample_conf(der_count1,der_count2)] += 1
                        # Series of quality and assumption checks to make sure that we can keep going
                        if not '.' in [l1[qual_1_ind], l2[qual_2_ind], lo[qual_OG_ind]]:
                                if (l1[filter_1_ind] in ['PASS','.']) and (l2[filter_2_ind] in ['PASS','.'] and (lo[filter_OG_ind] in ['PASS','.'])):
                                    if not (low_cov, high_cov) == (0,0):
                                        [coverage1,genotype1]=get_genotype(l1[9:])
                                        [coverage2,genotype2]=get_genotype(l2[9:])
                                        [coverage_ogrp,genotype_ogrp]=get_genotype(lo[9:])
                                        if not check_if_pass_coverage(coverage1,low_cov,high_cov) or not check_if_pass_coverage(coverage2,low_cov,high_cov) or not check_if_pass_coverage(coverage_ogrp,low_cov,high_cov): continue
                                    if nucl_A in nucl:
                                        var_form=check_if_ok_and_get_var_form(nucl_A,ref_1,ref_2,ref_OG,alt_1,alt_2,alt_OG)
                                        if not var_form=='':
                                            der_count_ogrp=orient_and_get_count(genotype_ogrp,ref_OG,alt_OG,nucl_A)
                                            if der_count_ogrp>0:
                                                if var_form=='OK_NO_VARIATION':
                                                    if ref_1==nucl_A:
                                                        local_count[0] += 1
                                                    else:
                                                        local_count[8] += 1
                                                elif var_form=='OK_POLY':
                                                    der_count1=orient_and_get_count(genotype1,ref_1,alt_1,nucl_A)
                                                    der_count2=orient_and_get_count(genotype2,ref_1,alt_2,nucl_A)
                                                    local_count[get_sample_conf(der_count1,der_count2)] += 1
                        
    if out_dict and not count_TT:
        return out_dict, False
    elif out_dict and out_dict_TT:
        return out_dict, out_dict_TT
    else:
        print(f"Error: It seems that every position in files {pop1}, {pop2} and {outgroup} failed all checks and no counts were generated for these files. Please check file formatting or whether all positions truly violate assumptions.")
        sys.exit(1)

def estimate_param_TT(counts):
    '''Function for estimating the parameters of the model for TT counts. Methods and math are not mine (Milo), see the original TT paper for more information (on the github page).
    Input: list of counts for each of the eight scenarios for either the entire genome (observed) or just one window (local)
    Output: list of all the different parameters'''
    # The first and last cases are not used except in the totals, therefore only the other 7 need to be initialised
    n1, n2, n3, n4, n5, n6, n7 = counts[1:8]
    n_tot=1.0*sum(counts)
    alfa1=2.0*n5/(n5+2.0*n6)
    alfa2=2.0*n5/(n5+2.0*n7)

    mu_t1_t2_diff=1.0*(n3-n4+0.5*(n1+n6-n2-n7))/n_tot

    thetaA=(3.0/n_tot)*(n5+2.0*n6)*(n5+2.0*n7)/(8.0*n5)

    mu_t1_p1=1.0*n3+0.5*n1
    mu_t1_p2=(n5+2.0*n6)*(n5+6.0*n7)/(8.0*n5)
    mu_t1=(mu_t1_p1-mu_t1_p2)/n_tot 
    mu_t2_p1=1.0*n4+0.5*n2
    mu_t2_p2=(n5+6.0*n6)*(n5+2.0*n7)/(8.0*n5)
    mu_t2=(mu_t2_p1-mu_t2_p2)/n_tot 
    if n5<2*n6:
        drift1=-1.0*log(alfa1)
        theta1=mu_t1/drift1
        the_nom=2.0*n6*(n1+n7)-1.0*n5*(n1+4.0*n3-n7)
        the_den=2.0*(2.0*n6-1.0*n5)*n_tot
        mu_nu1=the_nom/the_den
        the_nom=(2.0*n6+n5)*(n5*(8.0*n3+2.0*n7+n5)-2.0*n6*(6.0*n7+n5))
        the_den=(2.0*n6-n5)*(n5*(4.0*n1+8.0*n3-6.0*n7-n5)-2.0*n6*(6.0*n7+n5))
        W1ratio=the_nom/the_den
        logpt=1.0/(log(2.0*n5)-log(2.0*n6+n5))
        pt1=(2.0*n7+n1)*(2.0*n6+n5)/(2.0*n6-n5)
        pt2=(2.0*n3-n1)/logpt
        pt3=(2.0*n6+n5)*(6.0*n7+n5)/(4.0*n5)
        pt4=4.0*n5/(2.0*n6-n5)
        D1=(pt1+pt2-(pt3*(pt4+logpt)))/(2.0*n_tot)  
    else:
        drift1='NaN'
        theta1='NaN'
        mu_nu1='NaN'
        W1ratio='NaN'
        D1='NaN'    
    if n5<2*n7:
        drift2=-1.0*log(alfa2)
        theta2=mu_t2/drift2
        the_nom=2.0*n7*(n2+n6)-1.0*n5*(n2+4.0*n4-n6)
        the_den=2.0*(2.0*n7-1.0*n5)*n_tot
        mu_nu2=the_nom/the_den  
        the_nom=(2.0*n7+n5)*(n5*(8.0*n4+2.0*n6+n5)-2.0*n7*(6.0*n6+n5))
        the_den=(2.0*n7-n5)*(n5*(4.0*n2+8.0*n4-6.0*n6-n5)-2.0*n7*(6.0*n6+n5))
        W2ratio=the_nom/the_den 
        logpt=1.0/(log(2.0*n5)-log(2.0*n7+n5))
        pt1=(2.0*n6+n2)*(2.0*n7+n5)/(2.0*n7-n5)
        pt2=(2.0*n4-n2)/logpt
        pt3=(6.0*n6+n5)*(2.0*n7+n5)/(4.0*n5)
        pt4=4.0*n5/(2.0*n7-n5)
        D2=(pt1+pt2-(pt3*(pt4+logpt)))/(2.0*n_tot)  
    else:
        drift2='NaN'
        theta2='NaN'
        mu_nu2='NaN'
        W2ratio='NaN'
        D2='NaN'
    ## Method from Schlebusch et al 2012
    Ccount1=1.0*n3+0.5*n6
    D1orD2count1=0.5*n5+1.0*n7
    Ccount2=1.0*n4+0.5*n7
    D1orD2count2=0.5*n5+1.0*n6
    P1=-log((3.0/2.0)*D1orD2count1/(D1orD2count1+Ccount1))
    P2=-log((3.0/2.0)*D1orD2count2/(D1orD2count2+Ccount2))
    P1_time=P1*thetaA
    P2_time=P2*thetaA
    Fst=(2.0*n3+2.0*n4-1.0*n5)/(1.0*n1+1.0*n2+2.0*n3+2.0*n4+1.0*n5+1.0*n6+1.0*n7)
    ##
    return [alfa1,alfa2,thetaA,mu_t1,mu_t2,mu_nu1,mu_nu2,mu_t1_t2_diff,drift1,drift2,theta1,theta2,W1ratio,W2ratio,D1,D2,P1,P2,P1_time,P2_time,Fst]

def filter_list(g,n,a_list,site_counts):
    b_list=[]
    b_snps=[]
    num_rem_wins=0
    num_rem_snps=0
    for i in range(g):
        a_count=site_counts[i]
        try:
            b_list.append(1.0*float(a_list[i]))
            b_snps.append(a_count)
        except:
            num_rem_wins+=1
            num_rem_snps+=a_count
    return [g-num_rem_wins,n-num_rem_snps,b_list,b_snps]

def get_WBJ_mean_var(g,n,obs_mean,a_list,site_counts):
    a_mean=0
    if obs_mean=='NA':
        return ('NA','NA','NA')
    else:
        obs_mean=1.0*float(obs_mean)
        #print  'BEFORE',[g,n,a_list,site_counts]
        [g,n,a_list,site_counts]=filter_list(g,n,a_list,site_counts)
        #print  'AFTER',[g,n,a_list,site_counts]
        for i in range(g):
            a_count=1.0*site_counts[i]
            a_term=1.0-1.0*(a_count/n)
            a_mean+=a_term*a_list[i]
        a_mean=g*obs_mean-a_mean
        a_var=0
        for i in range(g):
            a_count=1.0*site_counts[i]
            hj=1.0*n/a_count
            a_term=hj*obs_mean-(hj-1.0)*a_list[i]-a_mean
            a_var+=a_term*a_term/(hj-1.0)
        a_var=a_var/g
        return (str(obs_mean),str(a_mean),str(a_var))

def get_estimates_TT(count_list):
    '''Function to obtain all parameter estimates and wbj statistics for TT method. Observed parameters for whole genome is calculated, then local parameters per window are calcualted, and both used to obtian wbj estimates. 
    Input: List of count lists per window, will be in order of chromosome but that information is not needed
    Output: a list containing lists for each parameter of observed parameter values, wbj mean and wbj variance.'''
    l_alfa1 = []
    l_alfa2 = []
    l_thetaA = []
    l_mu_t1 = []
    l_mu_t2 = []
    l_mu_nu1 = []
    l_mu_nu2 = []
    l_mu_diff_t1_t2 = []
    l_drift1 = []
    l_drift2 = []
    l_theta1 = []
    l_theta2 = []
    l_W1ratio = []
    l_W2ratio = []
    l_D1 = []
    l_D2 = []
    l_P1 = []
    l_P2 = []
    l_P1_time = []
    l_P2_time = []
    l_Fst = []
    num_sites = []
    g=0
    n=0
    obs_d = [0 for i in range(9)]
    # Using list comprehension and zipping, sum up the list of lists to obtain the observed data
    obs_d = [sum(count_window) for count_window in zip(*count_list)]
    try:
        [obs_alfa1,obs_alfa2,obs_thetaA,obs_mu_t1,obs_mu_t2,obs_mu_nu1,obs_mu_nu2,obs_mu_diff_t1_t2,obs_drift1,obs_drift2,obs_theta1,obs_theta2,obs_W1ratio,obs_W2ratio,obs_D1,obs_D2,obs_P1,obs_P2,obs_P1_time,obs_P2_time,obs_Fst]=estimate_param_TT(obs_d)
    except ZeroDivisionError as zeros:
        traceback.print_exc()
        print()
        print("Help: This error occurs when one of the genotype situation counts in the offending line of code is 0. If counts are all 0, check file formatting, the script can't read your genomes. If only some cases are 0, make sure the windows are wide enough or reconsider the comparison between the populations.")
        sys.exit(1)
    # Then go through and estimate parameters per window
    for count_window in count_list:
        if sum(count_window)>0:
            g+=1
            n+=sum(count_window)
            local_count=[obs_d[0]-count_window[0],obs_d[1]-count_window[1],obs_d[2]-count_window[2],obs_d[3]-count_window[3],obs_d[4]-count_window[4],obs_d[5]-count_window[5],obs_d[6]-count_window[6],obs_d[7]-count_window[7],obs_d[8]-count_window[8]]
            try:
                [alfa1,alfa2,thetaA,mu_t1,mu_t2,mu_nu1,mu_nu2,mu_diff_t1_t2,drift1,drift2,theta1,theta2,W1ratio,W2ratio,D1,D2,P1,P2,P1_time,P2_time,Fst]=estimate_param_TT(local_count)
            except ZeroDivisionError as zeros:
                traceback.print_exc()
                print()
                print("Help: This error occurs when one of the genotype situation counts in the offending line of code is 0. If counts are all 0, check file formatting, the script can't read your genomes. If only some cases are 0, make sure the windows are wide enough or reconsider the comparison between the populations.")
                sys.exit(1)
            # Add the local parameters to the lists of all local parameters
            l_alfa1.append(alfa1)
            l_alfa2.append(alfa2)
            l_thetaA.append(thetaA)
            l_mu_t1.append(mu_t1)
            l_mu_t2.append(mu_t2)
            l_mu_nu1.append(mu_nu1)
            l_mu_nu2.append(mu_nu2)
            l_mu_diff_t1_t2.append(mu_diff_t1_t2)
            l_drift1.append(drift1)
            l_drift2.append(drift2)
            l_theta1.append(theta1)
            l_theta2.append(theta2)
            l_W1ratio.append(W1ratio)
            l_W2ratio.append(W2ratio)
            l_D1.append(D1)
            l_D2.append(D2)
            l_P1.append(P1)
            l_P2.append(P2)
            l_P1_time.append(P1_time)
            l_P2_time.append(P2_time)
            l_Fst.append(Fst)
            num_sites.append(sum(count_window))
    # Obtain wbj results for each parameter
    res=[get_WBJ_mean_var(g,n,obs_alfa1,l_alfa1,num_sites)]
    res.append(get_WBJ_mean_var(g,n,obs_alfa2,l_alfa2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_thetaA,l_thetaA,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_mu_t1,l_mu_t1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_mu_t2,l_mu_t2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_mu_nu1,l_mu_nu1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_mu_nu2,l_mu_nu2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_mu_diff_t1_t2,l_mu_diff_t1_t2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_drift1,l_drift1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_drift2,l_drift2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_theta1,l_theta1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_theta2,l_theta2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_W1ratio,l_W1ratio,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_W2ratio,l_W2ratio,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_D1,l_D1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_D2,l_D2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_P1,l_P1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_P2,l_P2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_P1_time,l_P1_time,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_P2_time,l_P2_time,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_Fst,l_Fst,num_sites))
    return res

def get_cond_estimates(counts):
    n1, n2, n3, n4, n5, n6, n7 = counts[1:8]
    n_tot = 1.0*sum(counts)
    alfa1=1.0*(n1+n7+n5)/(n1+2.0*n3+n6+0.5*n6)
    alfa2=1.0*(n2+n6+n5)/(n2+2.0*n4+n7+0.5*n5)
    test1=(2.0*n1+n5)/(2.0*n2+n5)-(2.0*n7+n5)/(2.0*n6+n5)
    test2=((n1-n2)+2.0*(n3-n4)+(n6-n7))/n_tot

    return [alfa1,alfa2,test1,test2]

def estimate_param_TTo(counts, outgroup_counts):
    '''Function for estimating the parameters of the model for TTo counts. Methods and math are not mine (Milo), see the original TT paper for more information (on the github page).
    Input: list of counts and outgroup confirmed counts for each of the eight scenarios for either the entire genome (observed) or just one window (local)
    Output: list of all the different parameters'''
    [alfa1,alfa2,test1,test2]=get_cond_estimates(outgroup_counts)
    n1, n2, n3, n4, n5, n6, n7 = counts[1:8]
    n_tot=1.0*sum(counts)
    if alfa1*alfa2>0:
        y=(9.0*n5)/(n_tot*2.0*alfa1*alfa2)
        tau2_1=(3.0* (2.0*alfa1*n6-(1.0-alfa1)*n5) ) /(n_tot*2.0*alfa1*alfa2)
        tau2_2=(3.0* (2.0*alfa2*n7-(1.0-alfa2)*n5) ) /(n_tot*2.0*alfa1*alfa2)
        tau3_1=((5.0-2.0*alfa1)*n5-4.0*alfa1*n6 ) /(n_tot*2.0*alfa1*alfa2)
        tau3_2=((5.0-2.0*alfa2)*n5-4.0*alfa2*n7 ) /(n_tot*2.0*alfa1*alfa2)
        B1=(0.5*n1+n3+0.5*n6+0.25*n5-(5.0*n5/(4.0*alfa1*alfa2)) )/n_tot
        B2=(0.5*n2+n4+0.5*n7+0.25*n5-(5.0*n5/(4.0*alfa1*alfa2)) )/n_tot
        mean_tau_3=0.5*(tau3_1+tau3_2)
        mean_tau_2=0.5*(tau2_1+tau2_2)
        T1=B1-0.5*mean_tau_3
        T2=B2-0.5*mean_tau_3
        J1=B1-(3.0*mean_tau_3)*(3.0*mean_tau_3)/(6.0*mean_tau_2)
        J2=B2-(3.0*mean_tau_3)*(3.0*mean_tau_3)/(6.0*mean_tau_2)
    else:
        y='NaN'
        tau2_1='NaN'
        tau2_2='NaN'
        tau3_1='NaN'
        tau3_2='NaN'
        B1='NaN'
        B2='NaN'
        T1='NaN'
        T2='NaN'
        J1='NaN'
        J2='NaN'
    U1=(0.5*(1.0-alfa1)*n1-alfa1*n3+0.5*(1.0-alfa2)*n7 +0.25*(2.0-alfa2)*n5 )/n_tot
    U2=(0.5*(1.0-alfa2)*n2-alfa2*n4+0.5*(1.0-alfa1)*n6 +0.25*(2.0-alfa1)*n5 )/n_tot
    if alfa1<1:
        V1=(0.5*n1-  (alfa1/(1.0-alfa1))*n3  +0.5*n7/(1.0-alfa1)  )/n_tot
    else:
        V1='NaN'
    if alfa2<1:
        V2=(0.5*n2-(alfa2/(1.0-alfa2))*n4+0.5*n6/(1.0-alfa2))/n_tot

    else:
        V2='NaN'
    tau_test=y-1.5*(tau2_1+tau2_2)

    return [alfa1,alfa2,test1,test2,y,tau2_1,tau2_2,tau3_1,tau3_2,B1,B2,U1,U2,V1,V2,tau_test,T1,T2,J1,J2]

def get_estimates_TTo(count_list, outgroup_count_list):
    '''Function to obtain all parameter estimates and wbj statistics for TTo method. How does this work?
    Input: List of count lists per window, will be in order of chromosome but that information is not needed, list of outgroup count lists, the same order
    Output: a list containing lists for each parameter of observed parameter values, wbj mean and wbj variance.'''
    l_alfa1 = []
    l_alfa2 = []
    l_test1 = []
    l_test2 = []
    l_y = []
    l_tau2_1 = []
    l_tau2_2 = []
    l_tau3_1 = []
    l_tau3_2 = []
    l_B1 = []
    l_B2 = []
    l_U1 = []
    l_U2 = []
    l_V1 = []
    l_V2 = []
    l_tau_test = []
    l_T1 = []
    l_T2 = []
    l_J1 = []
    l_J2 = []
    num_sites = []
    g=0
    n=0
    obs_d = [0 for i in range(9)]
    obs_cond_d = [0 for i in range(9)]
    obs_d = [sum(count_window) for count_window in zip(*count_list)]
    obs_cond_d = [sum(count_window) for count_window in zip(*outgroup_count_list)]
    print(obs_d, obs_cond_d)
    try:
        [obs_alfa1,obs_alfa2,obs_test1,obs_test2,obs_y,obs_tau2_1,obs_tau2_2,obs_tau3_1,obs_tau3_2,obs_B1,obs_B2,obs_U1,obs_U2,obs_V1,obs_V2,obs_tau_test,obs_T1,obs_T2,obs_J1,obs_J2]=estimate_param_TTo(obs_d,obs_cond_d)
    except ZeroDivisionError as zeros:
        traceback.print_exc()
        print()
        print("Help: This error occurs when one of the genotype situation counts in the offending line of code is 0. If counts are all 0, check file formatting, the script can't read your genomes. If only some cases are 0, make sure the windows are wide enough or reconsider the comparison between the populations.")
        sys.exit(1)
    for i in range(len(outgroup_count_list)):
        count = count_list[i]
        outgroup_count = outgroup_count_list[i]
        if sum(count)>0:
            g+=1
            n+=sum(count)
            local_count=[obs_d[0]-count[0],obs_d[1]-count[1],obs_d[2]-count[2],obs_d[3]-count[3],obs_d[4]-count[4],obs_d[5]-count[5],obs_d[6]-count[6],obs_d[7]-count[7],obs_d[8]-count[8]]
            local_outgroup_count=[obs_cond_d[0]-outgroup_count[0],obs_cond_d[1]-outgroup_count[1],obs_cond_d[2]-outgroup_count[2],obs_cond_d[3]-outgroup_count[3],obs_cond_d[4]-outgroup_count[4],obs_cond_d[5]-outgroup_count[5],obs_cond_d[6]-outgroup_count[6],obs_cond_d[7]-outgroup_count[7],obs_cond_d[8]-outgroup_count[8]]
            try:
                [alfa1,alfa2,test1,test2,y,tau2_1,tau2_2,tau3_1,tau3_2,B1,B2,U1,U2,V1,V2,tau_test,T1,T2,J1,J2]=estimate_param_TTo(local_count,local_outgroup_count)
            except ZeroDivisionError as zeros:
                traceback.print_exc()
                print()
                print("Help: This error occurs when one of the genotype situation counts in the offending line of code is 0. If counts are all 0, check file formatting, the script can't read your genomes. If only some cases are 0, make sure the windows are wide enough or reconsider the comparison between the populations.")
                sys.exit(1)
            l_alfa1.append(alfa1)
            l_alfa2.append(alfa2)
            l_test1.append(test1)
            l_test2.append(test2)
            l_y.append(y)
            l_tau2_1.append(tau2_1)
            l_tau2_2.append(tau2_2)
            l_tau3_1.append(tau3_1)
            l_tau3_2.append(tau3_2)
            l_B1.append(B1)
            l_B2.append(B2)
            l_U1.append(U1)
            l_U2.append(U2)
            l_V1.append(V1)
            l_V2.append(V2)
            l_tau_test.append(tau_test)
            l_T1.append(T1)
            l_T2.append(T2)
            l_J1.append(J1)
            l_J2.append(J2)

            num_sites.append(sum(count))
    res=[get_WBJ_mean_var(g,n,obs_alfa1,l_alfa1,num_sites)]
    res.append(get_WBJ_mean_var(g,n,obs_alfa2,l_alfa2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_test1,l_test1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_test2,l_test2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_y,l_y,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_tau2_1,l_tau2_1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_tau2_2,l_tau2_2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_tau3_1,l_tau3_1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_tau3_2,l_tau3_2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_B1,l_B1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_B2,l_B2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_U1,l_U1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_U2,l_U2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_V1,l_V1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_V2,l_V2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_tau_test,l_tau_test,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_T1,l_T1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_T2,l_T2,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_J1,l_J1,num_sites))
    res.append(get_WBJ_mean_var(g,n,obs_J2,l_J2,num_sites))
    res.append([str(x) for x in obs_d])
    return res

