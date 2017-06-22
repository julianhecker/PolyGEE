#ifndef _DATA_HPP_
#define _DATA_HPP_
#include "setup.hpp"
#include "estimator.hpp"
double** get_matrix(const char* filename,int row,int col);
double* get_vector(const char* filename,int n);
int* get_vector_i(const char* filename,int n);
summary_stat* get_sum(const char* filename,int n, int ld_score_dim,int bivar);
summary_stat* get_sum_new(const char* filename,int n, int ld_score_dim,int bivar);
block* read_in(summary_stat* sum,int nsums,int* num_blocks);
block* explore_block_structure(summary_stat* sum,int nsums,int* num_blocks);
void get_cm_map_summary_stat(summary_stat* sum,int n_sum,int* num_blocks);
void get_matrixdims(const char* filename,int* num, int* num_cum, int n);
void get_block_objects(block* blocks,int* blocksizes,int n_blocks,summary_stat* sum,int n_sum);
void get_cluster_ids(summary_stat* sum,int n_sum,int* num_blocks);
void get_block_sizes(int* blocksizes,summary_stat* sum,int n_sum);
estimator* prepare(summary_stat* sum,int nsums,block* blocks,int num_blocks,double scal,int bivar,const char* , const char* );
char** get_hap_files(const char* filename);
#endif /* _DATA_HPP_ */
