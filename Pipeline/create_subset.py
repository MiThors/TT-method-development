import sys
import gzip
import get_file_name
from zipfile import ZipFile

ind1_vcf = " ../00_Example_Run/00_Data/Denisova_chr10.vcf.gz" # The path to the file
ind2_vcf = "../00_Example_Run/00_Data/Neandertal-Altai_chr10.vcf.gz" 
ancestral_txt = "../00_Example_Run/00_Data/Ancestral_states/chr10.txt.gz"

with gzip.open(ancestral_txt,'rt',encoding='utf-8') as ancestral:
    with gzip.open(ind1_vcf, 'rt', encoding='uft-8') as ind_1:
        with gzip.open(ind2_vcf, 'rt', encoding='uft-8') as ind_2:
            # This is a little loop to skip past the VCF file headers for both vcf files
            l1='##'
            while l1[0]=='##':
                l1=ind_1.readline()
            l2='##'
            while l2[0]=='##':
                l2=ind_2.readline()
            la=ancestral.readline()
            print(f'l1 = {l1}')
            print(f'l2 = {l2}')
            print(f'la = {la}')