# PolyGEE
PolyGEE: A generalized estimating equation approach to the efficient and robust estimation of polygenic effects in large-scale association studies


This is a beta version of the software PolyGEE for the efficient and robust estimation of polygenic effects in large-scale association studies.


The command to run the tool is

./estimate input_file nsums scal bivar ld_score_dim file_lengths file_haps

nsums: # SNPs in the input file
scal: scale parameter to scale the estimated LD matrices
bivar: set this to 0 at the moment
ld_score_dim: set this to 1 at the moment
file_lengths: a file with 23 rows, 0 in the first 0, describing the number of SNPs for the haplotypes for each chromosome. An example file corresponding to the 1000 Genomes data is in the figshare repository. It is called "lengths_hap".

file_haps: a file with 22 rows, describing the file location for each chromosome-haplotype file. An example file is "files_hap" in the figshare repository.

format of the input file:


SNP | CHR | BP | z-score | LD-Score | index_1000Genomes | cM | cluster_id


The information for the last 3 columns can be extracted from the "ld_database" file in the figshare reppsitory

https://figshare.com/account/home

The input file should have more than 600k SNPs.

IMPORTANT: the executable requires the GSL library and the 1000 Genomes data
