#include <iostream>

#include <algorithm>
#include <time.h>
#include <vector>
#include <math.h>
#include <string.h>
#include "util.hpp"


#include <stdio.h>
#include <stdlib.h>

#include "cluster.hpp"
#include "estimator.hpp"
#include "setup.hpp"
#include "data.hpp"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#define NCHR 22

#include <sys/time.h>
using namespace std;




int main(int argc,char *argv[])
{
  int i,j,z,k,it;
  const char* sumfile; // input file
  int nsums; // # SNPs in input file
  double scal; // scal of LD matrix
  int bivar; // bivar option for future extensions
  int ld_score_dim; // ld_score
  char* file_lengths;
  char* files_hap;
  if(argc==8)
  {
  sumfile=argv[1]; // input file
  nsums= atoi(argv[2]); // # SNPs in input file
  scal=atof(argv[3]); // scal of LD matrix
  bivar=atoi(argv[4]); // bivar option for future extensions
  ld_score_dim=atoi(argv[5]); // ld_score
  file_lengths=argv[6];
  files_hap=argv[7];
  }
  else
  {
     cout << "usage: estimate input_file nsums scal bivar ld_score_dim"<<endl;
	 exit(1);
  }
  summary_stat* sums=get_sum_new(sumfile,nsums,ld_score_dim,bivar);
  int nbl;

  block* blocks=explore_block_structure(sums,nsums,&nbl);

  estimator* est=prepare(sums,nsums,blocks,nbl,scal,bivar,file_lengths, files_hap);
  
  est->cm_width=1.0; // correlation restricted to 1 cM
  double* para=new double[3]; // 2 parameters: 1 for 1+Na and 1 for h^2*c_study
  para[1]=1.0;
  para[2]=0.001;// reasonable initial values
  for(i=1;i<=est->n_cluster;i++)
  {
       est->cl[i].dispersion_variance(1.0);
  }
  timeval start, end;
  gettimeofday(&start, 0);
  cout << "start estimation of parameters."<<endl;
  est->estimate(para,2);
  est->covar();
	
  gettimeofday(&end, 0);
  cout << "time in seconds: "<<end.tv_sec- start.tv_sec <<endl;
   
  cout << " 1+Na: "<<est->current_parameter[1]<<"("<<sqrt(est->cov_tot[1][1])<<")    h^2*c_study: "<< est->current_parameter[2]<<"("<<sqrt(est->cov_tot[2][2])<<")" <<endl;
  
 

 }











