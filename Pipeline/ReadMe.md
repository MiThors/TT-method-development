# The TT and TTo Methods
- Cite paper and stuff

## TT Arguments
-i --in: intakes 2 strings, genotype file (with path) for the two individuals. Assumes files are separated by chromosome, see -C if not. 

-t --type: intakes 1 int, options for input filetypes, 1 = vcf, 2 = tped, 3 = bam. 

-a --ancestral: intakes 1 string, ancestral states file (with path).

-k --keywords: optional, intakes 2 strings, keywords for the individuals (MUST be in the same order as in -f), default = "ind1" "ind2".

-C -Chromosomes: optional flag to indicate that all chromosomes are in the same file. Default false. 

-o --out: optional, intakes 1 string, a user-given name for the output directory. Default "TT_out_ind1_ind2".

## TTo Arguments
-i --in: intakes 3 strings, genotype file (with path) for the two individuals and the outgroup (MUST be in that order). Assumes files are separated by chromosome, see -C if not. 

-t --type: intakes 1 int, options for input filetypes, 1 = vcf, 2 = tped, 3 = bam. 

-a --ancestral: intakes 1 string, ancestral states file (with path).

-k --keywords: optional, intakes 3 strings, keywords for the individuals and the outgroup (MUST be in the same order as in -f), default = "ind1" "ind2" "outgroup".

-C -Chromosomes: optional flag to indicate that all chromosomes are in the same file. Default false. 

-o --out: optional, intakes 1 string, a user-given name for the output directory. Default "TTo_out_ind1_ind2_outgroup".

-T -TT: optional flag to include TT estimates of statistics (i.e. without using an outgroup) for comparison.

## File input
- Something or other

## Creating ancestral state files
- Tbh not like I have a clue. Just kinda vibe. Distance should be not too far, but should be faily distant. 

## Choosing an outgroup (TTo)
- Is a different outgroup to what I say in ancestral, this needs to be an in-group with the other genomes comapred to the ancestral states, but should be an outgorup to the compared genomes. 

