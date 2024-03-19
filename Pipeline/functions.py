'''A script contianing all the defined functions for using both TT and TTo methods.

Milo Thordarson: anth2886@student.uu.se'''

import gzip
import sys


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


def derived_not_in_outgroup(outgroup_genotype, ref, alt, nucl_A):
    '''Function to check if the derived allele '''
    if alt == '.' and ref == nucl_A: derived = 0
    elif alt == '.' and ref != nucl_A: derived = 2
    elif ref == nucl_A: derived = outgroup_genotype.count("1")
    else: derived = outgroup_genotype.count("0")
    if derived > 0 : return False
    else: return True


def get_configuration_index(nucl_A, genotype_1, genotype_2, ref_1, ref_2, alt_1, alt_2):
    if alt_1 == '.' and ref_1 == nucl_A: pop1_derived = 0
    elif alt_1 == '.' and ref_1 != nucl_A: pop1_derived = 2
    elif ref_1 == nucl_A: pop1_derived = genotype_1.count("1")
    else: pop1_derived = genotype_1.count("0")
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
        

def get_counts_vcf_TT(iterable):
    '''Function for getting counts from a vcf file. Opens files, checks formatting for pop1 and pop2 columns is correct, aligns positions in all files, ignores lines that do not pass filters, then adds counts to the appropriate situation.
    pop1, pop2, anc = filepaths for all files, list of one or more
    only_one_file = true or false if all chromosomes are in one file
    low_cov, high_cov = coverage thresholds for filtering, set up in main method file
    filters = list of values considered acceptable for FILTER field of vcf file
    Returns all eight count scenarios in a dictionary, keys are chromosomes, one list of counts per chromosome'''
    pop1 = iterable[0]
    pop2 = iterable[1]
    anc = iterable[2]
    low_cov = iterable[3]
    high_cov = iterable[4]
    filters = iterable[5]
    window_size = int(iterable[6])
    # Variable initialisation
    nucl = ['A','C','G','T']
    nt_set = set(nucl)
    window_step = 0
    out_dict = {}
    local_count = []
    win_pos = {}
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
                
                while l1 and l2 and la:
                    l1 = file_1.readline().strip().split()
                    l2 = file_2.readline().strip().split()
                    la = ancestral.readline().strip().split()
                    if not l1 or not l2 or not la: 
                        if local_count: out_dict[current_chrom].append(local_count)
                        win_pos[current_chrom].append((win_start, current_pos))
                        break
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
                        sys.exit(1)

                    # Define variables
                    nucl_A = la[nucl_A_ind]
                    ref_1, ref_2 = l1[ref_1_ind], l2[ref_2_ind]
                    alt_1, alt_2 = l1[alt_1_ind], l2[alt_2_ind]
                    chrom_1, chrom_2 = l1[chrom_1_ind], l2[chrom_2_ind]

                    if chrom_1 != chrom_2: 
                        print(f'Error: Files at same positions, but chromosomes in {pop1} and {pop2} are not the same, please check file formatting.')
                        sys.exit(1)

                    # Check if current chromosome exists in the dict already, if not add another key for that
                    if chrom_1 not in out_dict: 
                        out_dict.update({chrom_1 : []})
                        win_pos.update({chrom_1: []})
                        if local_count:
                            out_dict[current_chrom].append(local_count)
                        local_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                    
                    if window_step == 0:
                        win_start = pos_1

                    if window_step >= window_size:
                        out_dict[chrom_1].append(local_count)
                        local_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                        window_step = 0
                        win_pos[chrom_1].append((win_start, current_pos))
                        win_start = pos_1

                    current_chrom = chrom_1
                    window_step += 1
                    current_pos = pos_1
                    
                    # Series of quality and assumption checks to make sure that we can keep going 
                    if '.' in [l1[qual_1_ind], l2[qual_2_ind]] : continue
                    elif l1[filter_1_ind] not in filters or l2[filter_2_ind] not in filters : continue
                    elif nucl_A not in nucl: continue # If the ancient nucleotide is not resolved, we skip
                    elif len(set([nucl_A,ref_1,ref_2,alt_1,alt_2]).difference('.')) > 2: continue # Check for multiallelic site
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
                    elif '.' in [genotype_1, genotype_2] : continue # Check if genotypes are undefined
                    elif "2" in [genotype_1, genotype_2] : continue # Check for multiallelic
                    # Get the type of sample configuration, represented as the index of m0, m1, ... m8
                    configuration_index = get_configuration_index(nucl_A, genotype_1, genotype_2, ref_1, ref_2, alt_1, alt_2)
                    # Add one count to the relevant chromosome and configuration count
                    local_count[configuration_index] += 1
                    
    if out_dict and win_pos:
        return out_dict, win_pos
    else:
        print(f"Error: It seems that every position in files {pop1} and {pop2} failed all checks and no counts were generated for these files. Please check file formatting or whether all positions truly violate assumptions.")
        sys.exit(1)


def get_counts_vcf_TTo(iterable):
    '''Function for getting counts from a vcf file. Opens files, checks formatting for pop1 and pop2 columns is correct, aligns positions in all files, ignores lines that do not pass filters, then adds counts to the appropriate situation.
    pop1, pop2, anc = filepaths for all files, list of one or more
    only_one_file = true or false if all chromosomes are in one file
    low_cov, high_cov = coverage thresholds for filtering, set up in main method file
    filters = list of values considered acceptable for FILTER field of vcf file
    Returns all eight count scenarios in a dictionary, keys are chromosomes, one list of counts per chromosome'''
    # Variable initialisation
    pop1 = iterable[0]
    pop2 = iterable[1]
    outgroup = iterable[2]
    anc = iterable[3]
    low_cov = iterable[4]
    high_cov = iterable[5]
    filters = iterable[6]
    window_size = int(iterable[7])
    nucl = ['A','C','G','T']
    nt_set = set(nucl)
    window_step = 0
    out_dict = {}
    local_count = []
    win_pos = {}
    win_start = {}
    # Opening the files
    with gzip.open(anc,'rt',encoding='utf-8') as ancestral:
        with gzip.open(pop1, 'rt', encoding='utf-8') as file_1:
            with gzip.open(pop2, 'rt', encoding='utf-8') as file_2:
                with gzip.open(outgroup, 'rt', encoding='utf-8') as file_og:
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

                    while l1 and l2 and la:
                        l1 = file_1.readline().strip().split()
                        l2 = file_2.readline().strip().split()
                        lo = file_og.readline()
                        la = ancestral.readline().strip().split()
                        if not l1 or not l2 or not la : 
                            if local_count:
                                out_dict[current_chrom].append(local_count)
                                win_pos[current_chrom].append((win_start, current_pos))
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
                            out_dict.update({chrom_1 : []})
                            win_pos.update({chrom_1: []})
                        if local_count:
                            out_dict[current_chrom].append(local_count)
                        local_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                    
                        if window_step == 0:
                            win_start = pos_1

                        if window_step >= window_size:
                            out_dict[chrom_1].append(local_count)
                            local_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                            window_step = 0
                            win_pos[chrom_1].append((win_start, current_pos))
                            win_start = pos_1

                        current_chrom = chrom_1
                        window_step += 1
                        current_pos = pos_1

                        # Series of quality and assumption checks to make sure that we can keep going
                        if '.' in [l1[qual_1_ind], l2[qual_2_ind], l2[qual_OG_ind]] : continue
                        elif l1[filter_1_ind] not in filters or l2[filter_2_ind] not in filters or lo[filter_OG_ind] not in filters : continue
                        elif nucl_A not in nucl: continue # If the ancient nucleotide is not resolved, we skip
                        elif len(set([nucl_A,ref_1,ref_2,ref_OG,alt_1,alt_2,alt_OG]).difference('.')) > 2: continue # Check for multiallelic sites
                        elif not set([nucl_A,ref_1,ref_2,ref_OG,alt_1,alt_2,alt_OG]).difference('.').issubset(nt_set): continue # All nucleotides should be A, T, C or G.

                        # If it passes these checks, we get the genotype and coverage and check if they are also appropriate
                        l1_format_info = l1[format_1_ind].split(':')
                        l2_format_info = l2[format_2_ind].split(':')
                        lo_format_info = lo[format_OG_ind].split(':')
                        try:
                            genotype_1_ind, genotype_2_ind, genotype_OG_ind = l1_format_info.index("GT"), l2_format_info.index("GT"), lo_format_info.index("GT")
                            coverage_1_ind, coverage_2_ind, coverage_OG_ind = l1_format_info.index("DP"), l2_format_info.index("DP"), lo_format_info.index("DP")
                        except ValueError:
                            continue

                        l1_genotype_info, l2_genotype_info, lo_genotype_info = l1[format_1_ind + 1].split(':'), l2[format_2_ind + 1].split(':'), lo[format_OG_ind + 1].split(':')
                        coverage_1, coverage_2, coverage_OG = l1_genotype_info[coverage_1_ind], l2_genotype_info[coverage_2_ind], lo_genotype_info[coverage_OG_ind]
                        genotype_1, genotype_2, genotype_OG = l1_genotype_info[genotype_1_ind], l2_genotype_info[genotype_2_ind], lo_genotype_info[genotype_OG_ind]

                        if bad_coverage(coverage_1, low_cov, high_cov) or bad_coverage(coverage_2, low_cov, high_cov) or bad_coverage(coverage_OG, low_cov, high_cov): continue # Check the coverage is within acceptable thresholds
                        elif '.' in [genotype_1, genotype_2, genotype_OG] : continue # Check if genotypes are undefined
                        elif "2" in [genotype_1, genotype_2, genotype_OG] : continue # Check for multiallelic
                        elif derived_not_in_outgroup(genotype_OG, ref_OG, alt_OG, nucl_A): continue
                        # Get the type of sample configuration, represented as the index of m0, m1, ... m8
                        configuration_index = get_configuration_index(nucl_A, genotype_1, genotype_2, ref_1, ref_2, alt_1, alt_2)
                        # Add one count to the relevant chromosome and configuration count
                        out_dict[chrom_1][configuration_index] += 1
    if out_dict and win_pos:
        return out_dict, win_pos
    else:
        print(f"Error: It seems that every position in files {pop1}, {pop2} and {outgroup} failed all checks and no counts were generated for these files. Please check file formatting or whether all positions truly violate assumptions.")
        sys.exit(1)
