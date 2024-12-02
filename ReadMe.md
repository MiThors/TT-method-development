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

-d --depth_thresholds: required, two ints in the order of Low_Coverage then High_Coverage, the user should define what is too low or too high coverage for their particular dataset.

---------------------------------------------------------
## TTo Arguments
All arguments can also be found by typing "python3 TTo_Method.py --help" on the commandline. 

-1 --pop1: intakes one or more strings, the file path to the vcf(s) for the diploid individual from population 1. Required.

-2 --pop2: intakes one or more strings, the file path to the vcf(s) for the diploid individual from population 2. Required. 

-og --outgroup: intakes one or more strings, the file path to the vcf(s) for the diploid individual from the outgroup population. Required. 

-a --ancestral: intakes one or more strings, the file path to the ancestral states file. Information on how to generate these is found [below](#creating-ancestral-state-files). Required.

-k --keywords: intakes three strings, keywords for the individuals, and should be in the same order as the input parameters. Optional, and defaults to "pop1", "pop2" and "outgroup". 

-o --out: intakes one string, the name for the output directory. Optional, defaults "TT_out_pop1_pop2_outgroup", using the keywords given.

-w, --window: intakes one int, the window size for weighted-block jackknife. Optional, default is 2,500,000 which corresponds to about 2.5cM in humans. 

-c --counts: flag to print counts per chromosome per window to a counts file in the output directory. Optional, can be used for inspecting when errors during calculations happen, or when the results of calculations look strange.

-T -TT: optionally can provide the path to a TT counts file that was generated for the same two populations, to speed up the process as TT counts will not need to be re-counted. 

-d --depth_thresholds: required, two ints in the order of Low_Coverage then High_Coverage, the user should define what is too low or too high coverage for their particular dataset.

---------------------------------------------------------
## File Input
Vcfs for populations and outgroups, a tab-separated text file of ancestral states with the headers "POS" and "NUCL" for two columns, the position on the genome and what nucleotide it is. All chromosomes can be given in the same file, however, the scripts are parallelized based on providing files separately for each chromosome for each population. 

---------------------------------------------------------
## Creating Ancestral State Files
As mentioned, the human ancestral state files were generated using a consensus between the three great apes, which are the closest relatives of humans. From the appendix of the original paper: "We considered that the ancestral state of a site could be reliably inferred if the site had the same allele (A, C, T, or G) in the three apes and that there were at most two alleles among the three apes and the individuals being analyzed." This consensus of three is not necessarily required for any analysis. Users could use fewer or possibly more (to a point, else they would be too far away) species to do the analysis depending on what information they have available. But the distance should ideally be on a species level, rather than population.


---------------------------------------------------------
## Choosing an Outgroup (TTo)
Is a different outgroup to what I say in ancestral, this needs to be an in-group *population* with the other genomes compared to the ancestral states, but should simply have diverged earlier than the split between the two study populations. The point of the outgroup is to ascertain that derived alleles arose in the ancestral population before the divergence of the sub-populations. For example, a comparison between two non-African populations can use an African individual as this outgroup, as the out-of-Africa divergence will be older than any of the splits in the rest of the world. Therefore, using an outgroup that is another species entirely may be too distant to represent a useful ascertainment as derived alleles may be different. Users can run TT first, then try TTo with several different outgroups providing the TT file to assess the affect of outgroup choice. 

---------------------------------------------------------
## References
Original paper:<br/>
Estimating divergence times from DNA sequences.<br/>
Per Sjödin, James McKenna, Mattias Jakobsson.<br/>
https://doi.org/10.1093/genetics/iyab008








