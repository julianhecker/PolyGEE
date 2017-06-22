#ifndef _SETUP_HPP_
#define _SETUP_HPP_
#include <stddef.h>
#include "cmatrix.hpp"
#include "util.hpp"
#include <math.h>
#include <iostream>
using namespace std;
class block
{
   public:
   
   int size;
   int chr;
   int* indices;

};
class summary_stat
{

   public:
   
   char rs_id[80];
   unsigned int bp;
   double z1;
   double z2;
   int num_ld_scores;
   double* ld_score;
   int chr;
   int ref_index;
   double cm;
   int cluster_id;
   
   summary_stat (){
   ld_score=NULL;};
   summary_stat (int dim)
   {
     ld_score=d_createVector(dim);
   };
   ~summary_stat()
   {
     if(ld_score!=NULL) d_destroyVector(ld_score);
   };
   void init(int dim);

   
};

class summary
{

   public:
   
   char rs_id[80];
   unsigned int bp;
   double z_score;
   double z1;
   double z2;
   double info;
   double ld_score;
   int chr;
   int ref_index;
   double cm;
   
   
   int cluster_id;
   
   summary (){};
   ~summary(){};
  
   
};

class ref_SNP
{

   public:
   
   char rs_id[80];
   unsigned int bp;
   int chr;
   int include;
   double cm;
   double ld_score;
   double ld_w_score;
   
   ref_SNP (){};
   ~ref_SNP(){};
  
   
};




#endif /* _SETUP_HPP_ */