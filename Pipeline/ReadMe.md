# The TT and TTo Methods
This repository contains guidelines and scripts for running the TT method for estimating divergence times between populations. This quick and transparent method requires two haploid genomes (or a single diploid genome) from each of the two populations, and the ancestral state for their common ancestor. Optionally, an additional two haploid genomes (or a single diploid genome) can be included, one that shares the same ancestral states but diverged earlier than the split between the study populations. This additional genome is referred to as the outgroup, and the method as TT-Outgroup, or TTO. Please see the published paper for full details (reference at the bottom). 

## TT Arguments
-1 --pop1: intakes one or more strings, the file path to the vcf(s) for the individual from population 1. 

-2 --pop2: intakes one or more stirngs, the file path to the vcf(s) for the individual from population 2.  

-a --ancestral: intakes 1 string, ancestral states file (with path).

-k --keywords: optional, intakes 2 strings, keywords for the individuals (MUST be in the same order as in -f), default = "ind1" "ind2".

-o --out: optional, intakes 1 string, a user-given name for the output directory. Default "TT_out_ind1_ind2".

## TTo Arguments
-i --in: intakes 3 strings, genotype file (with path) for the two individuals and the outgroup (MUST be in that order). Assumes files are separated by chromosome, see -C if not. 

-t --type: intakes 1 int, options for input filetypes, 1 = vcf, 2 = tped, 3 = bam. 

-a --ancestral: intakes 1 string, ancestral states file (with path).

-k --keywords: optional, intakes 3 strings, keywords for the individuals and the outgroup (MUST be in the same order as in -f), default = "ind1" "ind2" "outgroup".

-o --out: optional, intakes 1 string, a user-given name for the output directory. Default "TTo_out_ind1_ind2_outgroup".

-T -TT: optional flag to include TT estimates of statistics (i.e. without using an outgroup) for comparison.

## File Input
- Something or other

## Creating Ancestral State Files
- Tbh not like I have a clue. Just kinda vibe. Distance should be not too far, but should be faily distant. 

## Choosing an Outgroup (TTo)
- Is a different outgroup to what I say in ancestral, this needs to be an in-group with the other genomes comapred to the ancestral states, but should be an outgorup to the compared genomes. 

## References
Original TT and TTo Paper:<br/>
Estimating divergence times from DNA sequences.<br/>
Per Sjödin, James McKenna, Mattias Jakobsson.<br/>
https://doi.org/10.1093/genetics/iyab008
