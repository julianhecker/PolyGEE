# PolyGEE
PolyGEE: A generalized estimating equation approach to the efficient and robust estimation of polygenic effects in large-scale association studies


This is a test version of the software "PolyGEE" for the efficient and robust estimation of polygenic effects in large-scale association studies.


### DOWNLOAD OF EXTERNAL DATA AND PREPARATION OF INPUT

1.) Download of 1000 Genomes data

To run the PolyGEE tool you need to download the 1000Genomes haplotype data from the figshare repository

https://figshare.com/s/9bc34e7c25e152297819

The files are called 1kghapX.impute.hap, where X denotes the chromosome. In addition, you need to prepare two files to provide the information about the 1000 Genomes data to the software. The first file should have 23 lines, a 0 in the first line and then the number of SNPs for all 22 chromosomes. The example file for this 1000Genomes data in the figshare repository is in the same directory and called "lenghts_hap". The second file lists the paths to the 22 haplotype files. The corresponding example file in the figshare repository is called "files_hap". Note: the haplotype data is for the european population and it is assumed that your study population matches this criteria.

2.) Preparation of the summary statistics input file

We assume that you have a qc-ed summary statistics file for a large-scale GWAS or meta-analysis (>5k samples). You can use for example the tool "ldsc" for LD Score regression to perform this step. 
Then, you should merge your data with the file "ld_database" in the figshare repository to create a file with the following columns:

SNP | CHR | BP | z-score | LD-Score | index_1000Genomes | cM | cluster_id

The last 4 columns come from the ld_database file. Your file should have at least 500k SNPs.


### RUN THE ANALYSIS

The command to run the tool is

./estimate input_file nsums scal bivar ld_score_dim file_lengths file_haps



nsums: # SNPs in the input file

scal: scale parameter to scale the estimated LD matrices

bivar: set this to 0 at the moment

ld_score_dim: set this to 1 at the moment

file_lengths: a file with 23 rows, 0 in the first 0, describing the number of SNPs for the haplotypes for each chromosome.

file_haps: a file with 22 rows, describing the file location for each chromosome-haplotype file. 



IMPORTANT: the executable requires the GSL library and the 1000 Genomes data. You can compile the tool using the make_file.

Output: The tool provides the estimates for 1+Na and for h^2 * C_study with standard error in brackets.

### REPRODUCING THE REAL DATA RESULTS

To reproduce the real data results for the population-based studies, please download the corresponding summary statistics from

http://www.med.unc.edu/pgc/results-and-downloads

and prepare the input according to the description above.
