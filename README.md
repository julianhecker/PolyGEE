# PolyGEE
PolyGEE: A generalized estimating equation approach to the efficient and robust estimation of polygenic effects in large-scale association studies


This is a beta version of the software PolyGEE for the efficient and robust estimation of polygenic effects in large-scale association studies.


The command to run the tool is

./estimate input_file nsums scal bivar ld_score_dim

nsums: # SNPs in the input file
scal: scale parameter to scale the estimated LD matrices
bivar: set this to 0 at the moment
ld_score_dim: set this to 1 at the moment

format of the input file:


SNP | CHR | BP | z-score | LD-Score | index_1000Genomes | cM | cluster_id


The information for the last 3 columns can be extracted from the "ld_panel" file in this repository.
The input file should have more than 600k SNPs.
