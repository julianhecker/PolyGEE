#ifndef _ESTIMATOR_HPP_
#define _ESTIMATOR_HPP_
#include <stddef.h>
#include "cmatrix.hpp"
#include "cluster.hpp"
#include <math.h>
#include <iostream>
#include "util.hpp"
using namespace std;
class estimator
{

   public:
   int n_cluster;
   int parameter_dim;
   int n_snp;
   
   cluster* cl;
   
   double* current_parameter;
   double* new_parameter;
   
   double* score;
   double** obs_fisher;
   double** obs_fisher_inv;
   
   double** cov;
   double** cov_tot;
   
   double** H_sum;
   double** H_sum_inv;
   
   double dispersion;
   double cm_width;
   
   
   estimator()
   {
      
	  cl=NULL;
   }
   ~estimator()
   {
      
	  if(cl!=NULL)delete[] cl;
	  if(score!=NULL) d_destroyVector(score);
	  if(obs_fisher!=NULL) d_destroyMatrix(obs_fisher,parameter_dim);
	  if(obs_fisher_inv!=NULL) d_destroyMatrix(obs_fisher_inv,parameter_dim);
	  if(current_parameter!=NULL) d_destroyVector(current_parameter);
	  if(new_parameter!=NULL) d_destroyVector(new_parameter);
	  if(cov!=NULL) d_destroyMatrix(cov,parameter_dim);
	  if(cov_tot!=NULL) d_destroyMatrix(cov_tot,parameter_dim);
	  if(H_sum!=NULL) d_destroyMatrix(H_sum,parameter_dim);
	  if(H_sum_inv!=NULL) d_destroyMatrix(H_sum_inv,parameter_dim);
	  
   }
   void init(int para_dim);
   void get_score_sum();
   void get_H_sum();
   void get_H_sum2();
   int update_step();
   void estimate(double* startpara,int len);
   void covar();
   //double get_corr_var(double h2_1,double h2_2, double h12, double var1, double var2, double v12);
   double get_corr_var(int ind1,int ind2,int ind3);
   double get_corr_var2(int ind1,int ind2,int ind3);


};


#endif /* _ESTIMATOR_HPP_ */