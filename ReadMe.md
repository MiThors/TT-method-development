---------------------------------------------------------
# The TT Method Manual
---------------------------------------------------------

This repository contains guidelines and scripts for running the TT method for estimating divergence times between populations. This quick and transparent method requires two haploid genomes (or a single diploid genome) from each of the two populations, and the ancestral state for their common ancestor. Optionally, an additional two haploid genomes (or a single diploid genome) can be included, one that shares the same ancestral states but diverged earlier than the split between the study populations. This additional genome is referred to as the outgroup, and the method as TT-Outgroup, or TTO. Please see the published paper for full details (reference at the bottom). 

---------------------------------------------------------
## Setup and Data Formatting

The directory 'Pipeline' contains all the scripts necessary to run both the TT and the TTO method. Users should download all scripts to a suitable location. Please see the next section on running the scripts. 

Genomes for the two populations, and the outgroup for TTO, should be in VCF format. The method requires genotype data and depth data for each position and will disregard any sites missing those fields in the VCF. The thresholds to consider "good depth" can be easily changed by the user. Errors due to formatting assumptions have been written to be as informative as possible, but note that most issues arise from formatting. All VCF files should be zipped.

Both TT and TTO methods require the ancestral states at all positions in the genome. The scripts expect a zipped .txt file with two columns: POS for the position on the genome and NUCL for the ancestral nucleotide. It is up to the user to determine what is appropriate as an estimation for the ancestral nucleotide. 
An example of such files created for the human genome based on consensus among three species of apes can be found at the zenodo DOI below. We consider a particular position informative only when the ancestral state has consensus support among all three of Gorilla, Chimpanzee and Orangutan, hence there is an additional column for that information. If intending to estimate divergence times for human populations, the files in 'Ancestral_states.zip' can be downloaded and used, though please note that the method no longer assumes this final column (so other users using for example 2 ancestral species do not need to overly change their formatting), so please re-format to remove those rows with <3 consensus. 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4441887.svg)](https://doi.org/10.5281/zenodo.4441887)

---------------------------------------------------------

Brief description of included scripts:

'get_file_name.py' - this links keywords to full vcf file names and paths for ease of implementation. The User should edit to include vcf file paths and a relevant keyword for each individual's vcfs.

'count_sample_confs_per_ind_TT.py' & 'count_sample_confs_per_ind_TTO.py'.<br/>
These scripts take all-sites vcfs as input, and return counts of sample configurations in 5MB blocks of the genome. Users should edit each script to include file paths for ancestral state and vcf files. Resulting counts are outputted to 'DIR_counts_per_5cm_TT/' and 'DIR_counts_per_5cm_TTO/' respectively.

'get_counts_TT.sh' and 'get_counts_TTO.sh'.<br/>
These are example SLURM submission scripts that can be used to implement the above scripts for the 22 autosomes of the human genome, and for as many pairwise individual comparisons as desired. Users should edit to include relevant vcf keywords and SLURM commands.

'get_estimates_TT.py'.<br/>
This uses the sample configuration counts previously obtained to estimate parameters including divergence times. Results are outputted to 'DIR_estimates_TT/'.

'get_estimates_TTO.py'.<br/>
This script uses the sample configuration counts previously obtained to estimate parameters including divergence times. Results are outputted to 'DIR_estimates_TTO/'.

'plot_TT.R'.<br/>
This R script will create plots of divergence time estimates present in 'DIR_estimates_TT/', and output plots to 'DIR_plots/'.

'plot_TTO.R'.<br/>
This R script will create plots of divergence time estimates present in 'DIR_estimates_TTO/', and output plots to 'DIR_plots/'.

'wbj.py'.<br/>
This script contains functions used by 'get_estimates_TT.py' and 'get_estimates_TTO.py' to perform weighted bloack jack-knife estimation of parameters.

---------------------------------------------------------

Implementation:

The User should create the following directories at the same location as scripts:<br/>
DIR_counts_per_5cm_TT<br/>
DIR_counts_per_5cm_TTO_$OUTGROUP<br/>
DIR_error_TT<br/>
DIR_error_TTO<br/>
DIR_estimates_TT<br/>
DIR_estimates_TTO/$OUTGROUP_res/<br/> 
DIR_plots/TTO_$OUTGROUP<br/>
(where $OUTGROUP is the keyword of an outgroup individual from 'get_file_name.py')

Once the necessary script have been edited to include vcf locations, the TT & TTO methods can be implemented simply by using:<br/>
bash get_counts_TT.sh<br/>
python get_estimates_TT.py<br/>
Rscript TT_plot.R

bash get_counts_TTO.sh<br/>
python get_estimates_TTO.py $OUTGROUP<br/> 
Rscript TTO_plot.R $OUTGROUP<br/> 
(where OUTGROUP is the keyword of an outgroup individual from 'get_file_name.py')


---------------------------------------------------------

vcf file names:

Both the TT and TTO method are set up to take compressed all-sites vcfs as input files, with one vcf per chromosome (22 vcfs per individual). To be useable, vcf file names should follow the naming convention format "...chr1.vcf.gz" for chromosome 1, and "...chr2.vcf.gz" for chromosome 2 etc. For example, in the 'get_file_names.py' script provided, the keyword 'Neanderthal' points to a general vcf file name of 'AltaiNea.hg19_1000g.dq.bqual.RG.realn-snpAD_chr.vcf.gz'. The actual vcf file names should be 'AltaiNea.hg19_1000g.dq.bqual.RG.realn-snpAD_chr1.vcf.gz' for chromosome 1, and so on.

---------------------------------------------------------

For reference:<br/>
Estimating divergence times from DNA sequences.<br/>
Per Sjödin, James McKenna, Mattias Jakobsson.<br/>
https://doi.org/10.1093/genetics/iyab008








