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
## TT Arguments
All arguments can also be found by typing "python3 TT_Method.py --help" on the commandline. 

-1 --pop1: intakes one or more strings, the file path to the vcf(s) for the diploid individual from population 1. Required.

-2 --pop2: intakes one or more strings, the file path to the vcf(s) for the diploid individual from population 2. Required. 

-a --ancestral: intakes one or more strings, the file path to the ancestral states file. Information on how to generate these is found [below](#creating-ancestral-state-files). Required.

-k --keywords: intakes two strings, keywords for the individuals, and should be in the same order as the input parameters. Optional, and defaults to "pop1" and "pop2". 

-o --out: intakes one string, the name for the output directory. Optional, defaults "TT_out_pop1_pop2", using the keywords given.

-w, --window: intakes one int, the window size for weighted-block jackknife. Optional, default is 2,500,000 which corresponds to about 2.5cM in humans. 

-c --counts: flag to print counts per chromosome per window to a counts file in the output directory. Optional, can be used for inspecting when errors during calculations happen, or when the results of calculations look strange.

-d --depth_thresholds: required, two ints in the order of Low_Coverage then High_Coverage, the user should define what is too low or too high coverage for their particular dataset

---------------------------------------------------------
## TTo Arguments
All arguments can also be found by typing "python3 TTo_Method.py --help" on the commandline. 

-1 --pop1: intakes one or more strings, the file path to the vcf(s) for the diploid individual from population 1. Required.

-2 --pop2: intakes one or more strings, the file path to the vcf(s) for the diploid individual from population 2. Required. 

-og --outgroup: intakes one or more strings, the file path to the vcf(s) for the diploid individual from the outgroup population. Required. 

-a --ancestral: intakes 1 string, ancestral states file (with path).

-k --keywords: optional, intakes 3 strings, keywords for the individuals and the outgroup (MUST be in the same order as in -f), default = "ind1" "ind2" "outgroup".

-o --out: optional, intakes 1 string, a user-given name for the output directory. Default "TTo_out_ind1_ind2_outgroup".

-w, --window: intakes one int, the window size for weighted-block jackknife. Optional, default is 5,000,000 which corresponds to about 5cm. 

-c --counts: flag to print counts per chromosome per window to a counts file in the output directory. Optional, can be used for inspecting when errors during calculations happens, or when results of calculations look strange.

-T -TT: optional flag to include TT estimates of statistics (i.e. without using an outgroup) for comparison.

-d --depth_thresholds: required, two ints in the order of Low_Coverage then High_Coverage, the user should define what is too low or too high coverage for their particular dataset

---------------------------------------------------------
## File Input
Vcfs for populations and outgroups, a tab seperated text file of ancestral states with the headers "POS" and "NUCL" for two columns, the position on the genome and what nucleotide it is.


---------------------------------------------------------
## Creating Ancestral State Files
To Do Advice


---------------------------------------------------------
## Choosing an Outgroup (TTo)
Is a different outgroup to what I say in ancestral, this needs to be an in-group with the other genomes comapred to the ancestral states, but should be an outgorup to the compared genomes. 

---------------------------------------------------------
## References
Original paper:<br/>
Estimating divergence times from DNA sequences.<br/>
Per Sjödin, James McKenna, Mattias Jakobsson.<br/>
https://doi.org/10.1093/genetics/iyab008








