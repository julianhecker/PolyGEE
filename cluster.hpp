#ifndef _CLUSTER_HPP_
#define _CLUSTER_HPP_
#include <stddef.h>
#include "cmatrix.hpp"
#include <math.h>
#include <iostream>
using namespace std;
class cluster
{

   public:
   int parameter_dim;
   int size;
   double** reg;
   double* y;
   cmatrix W;
   double* constants;
   double* eta;
   double* res;
   double* var;
   double** D;
   double** H;
   double* score;
   int group;
   double start_cm;
   double end_cm;
   double dispersion;
   int bivar;
   
   cluster ();
   cluster (int,int);
   ~cluster();
   void dispersion_variance(double disp);
   void get_eta(double* parameter,int len);
   void get_res();
   void get_var_gamma();
   void get_score(double* parameter,int len);
   void get_D_T_gamma();
   void get_H(double* parameter, int len);
   void init(int sizep,int parameter_dimp);
   void mult_res_with_inv_matrix(double* tmpVec,double* tmpRes);
};


#endif /* _CLUSTER_HPP_ */