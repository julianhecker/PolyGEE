#include <math.h>
#include "data.hpp"
#include "estimator.hpp"

#include "setup.hpp"
#include "util.hpp"
#include <iostream>
#include <string.h>
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


#include <stdio.h>
#include <stdlib.h>
using namespace std;
#define NCHR 22
#define NINV 1006 //758


#define MAX_win 10000

double** get_matrix(const char* filename,int row,int col)
{
    FILE* fp;

    fp=fopen(filename,"r");
    double** mat=d_createMatrix(row,col);
    int i,j;

    double f;

    for(i=1;i<=row;i++)
    {
	  for(j=1;j<=col;j++)
	  {
        fscanf(fp,"%le",&f);
	    mat[i][j]=f;
		//cout<<  f << endl;

	  }
	  fscanf(fp,"\n");

    }
	fclose(fp);
	return mat;
}
char** get_hap_files(const char* filename)
{
    FILE* fp;

    fp=fopen(filename,"r");
    char** vec=new char*[NCHR+1];
    int i;
    for(i=1;i<=NCHR;i++)
    {
      vec[i]=new char[80];

    }
    for(i=1;i<=NCHR;i++)
    {

        fscanf(fp,"%s",vec[i]);

	    //cout << vec[i]<<endl;

    }
	fclose(fp);
	return vec;
}
double* get_vector(const char* filename,int n)
{
    FILE* fp;

    fp=fopen(filename,"r");
    double* vec=d_createVector(n);
    int i,j;

    double f;

    for(i=1;i<=n;i++)
    {

        fscanf(fp,"%le",&f);
	    vec[i]=f;
	    //cout << vec[i]<<endl;

    }
	fclose(fp);
	return vec;
}
void get_matrixdims(const char* filename,int* num, int* num_cum, int n)
{
    FILE* fp;

    fp=fopen(filename,"r");

    int i;
    int f;

    for(i=1;i<=n;i++)
    {

        fscanf(fp,"%i",&f);
	    num[i-1]=f;
	    //cout << vec[i]<<endl;

    }
	num_cum[0]=0;
	for(i=1;i<=n-1;i++)
	{
	  num_cum[i]=num_cum[i-1]+num[i];
	}
	fclose(fp);

}
int* get_vector_i(const char* filename,int n)
{
    FILE* fp;

    fp=fopen(filename,"r");
    int* vec=i_createVector(n);
    int i,j;

    int f;
    vec[0]=0;
    for(i=1;i<=n;i++)
    {

        fscanf(fp,"%i",&f);
	    vec[i]=f;
	    //cout << vec[i]<<endl;

    }
	fclose(fp);
	return vec;
}
summary_stat* get_sum_new(const char* filename,int n, int ld_score_dim,int bivar)
{
    FILE* fp;
    int i,j;
    fp=fopen(filename,"r");
    summary_stat* sum=new summary_stat[n+1];
    for(i=0;i<=n;i++)
    {
	  sum[i].init(ld_score_dim);
    }

    double z1,z2;
    z2=0.0;
    char rs_tmp[80];
    int chr;
    unsigned int bp;
    int ref_index;
    double ld;
    double cmf;
    int id;


    for(i=1;i<=n;i++)
    {

	    fscanf(fp,"%s",rs_tmp);
	    fscanf(fp,"%i",&chr);
	    fscanf(fp,"%i",&bp);
	    fscanf(fp,"%le",&z1);
	    if(bivar==1) fscanf(fp,"%le",&z2);
	    for(j=1;j<=ld_score_dim;j++)
	    {
		  fscanf(fp,"%le",&ld);
		  sum[i].ld_score[j]=ld;
	    }
	    fscanf(fp,"%i",&ref_index);
	    fscanf(fp,"%le",&cmf);
        fscanf(fp,"%i",&id);
	    //if(i%1000==0) cout << rs_tmp<<"    " << chr <<"    " <<bp<< "   "<<ld_score[1]<<"  "<< cmf<<endl;

	    strncpy ( sum[i].rs_id,rs_tmp,80);
	    sum[i].chr=chr;
	    sum[i].bp=bp;
	    sum[i].z1=z1;
	    sum[i].z2=z2;
	    sum[i].cm=cmf;
	    sum[i].ref_index=ref_index;

        sum[i].cluster_id=id;


    }
    fclose(fp);
    return sum;
}
summary_stat* get_sum(const char* filename,int n, int ld_score_dim,int bivar)
{
    FILE* fp;
    int i,j;
    fp=fopen(filename,"r");
    summary_stat* sum=new summary_stat[n+1];
    for(i=0;i<=n;i++)
	{
	  sum[i].init(ld_score_dim);
	}

    double z1,z2;
	z2=0.0;
	char rs_tmp[80];
	int chr;
	unsigned int bp;
	int ref_index;
	double ld;
    double cmf;
    for(i=1;i<=n;i++)
    {

	    fscanf(fp,"%s",rs_tmp);
		fscanf(fp,"%i",&chr);
		fscanf(fp,"%i",&bp);
		fscanf(fp,"%le",&z1);
		if(bivar==1) fscanf(fp,"%le",&z2);
		for(j=1;j<=ld_score_dim;j++)
		{
		  fscanf(fp,"%le",&ld);
		  sum[i].ld_score[j]=ld;
		}
		fscanf(fp,"%i",&ref_index);
		fscanf(fp,"%le",&cmf);
	    //if(i%1000==0) cout << rs_tmp<<"    " << chr <<"    " <<bp<< "   "<<ld_score[1]<<"  "<< cmf<<endl;

		strncpy ( sum[i].rs_id,rs_tmp,80);
		sum[i].chr=chr;
		sum[i].bp=bp;
		sum[i].z1=z1;
		sum[i].z2=z2;
		sum[i].cm=cmf;
		sum[i].ref_index=ref_index;


    }
	fclose(fp);
	return sum;
}
estimator* prepare(summary_stat* sum,int nsums,block* blocks,int num_blocks,double scal,int bivar,const char* file_lengths, const char* file_names)
{
   int i,j,k;
   int* num=i_createVector(NCHR);
   int* num_cum=i_createVector(NCHR);
   get_matrixdims(file_lengths,num,num_cum,NCHR+1);

   char** filenames=get_hap_files(file_names);


   int num_scores=sum[1].num_ld_scores;
   int nparameter=num_scores+1;

   estimator* est=new estimator;
   est->init(nparameter);
   est->n_cluster=num_blocks;
   est->cl=new cluster[est->n_cluster+1];
   for(i=1;i<=est->n_cluster;i++)
   {
     est->cl[i].bivar=bivar;
   }


   double** tmpmat;int* indices;
   int tmp_chr=1;
   tmpmat=get_matrix(filenames[1],num[1],NINV);

   for(i=1;i<=num_blocks;i++)
   {
        if(blocks[i].chr==tmp_chr+1){ d_destroyMatrix(tmpmat,num[tmp_chr]); tmpmat=get_matrix(filenames[blocks[i].chr],num[blocks[i].chr],NINV); tmp_chr=blocks[i].chr;}
        est->cl[i].init(blocks[i].size,nparameter);
		indices=i_createVector(blocks[i].size);

		est->cl[i].start_cm=sum[blocks[i].indices[1]].cm;
		est->cl[i].end_cm=sum[blocks[i].indices[blocks[i].size]].cm;
		est->cl[i].group=blocks[i].chr;

		for(j=1;j<=blocks[i].size;j++)
	    {
						est->cl[i].reg[j][1]=1.0;

						for(k=1;k<=num_scores;k++)
						{
						  est->cl[i].reg[j][k+1]=sum[blocks[i].indices[j]].ld_score[k];
						}
						est->cl[i].constants[j]=0.0;
						est->cl[i].y[j]=pow(sum[blocks[i].indices[j]].z1,2.0);
						indices[j]=sum[blocks[i].indices[j]].ref_index-num_cum[blocks[i].chr-1];


		}
		//cout << i<< "    "<<blocks[i].size <<"    " <<blocks[i].chr<<endl;
	    get_working_correlation(tmpmat,indices,num[est->cl[i].group],blocks[i].size,scal,est->cl[i].W.matrix);
		est->cl[i].W.invert();
		i_destroyVector(indices);


   }




  return est;
}
void get_cluster_ids(summary_stat* sum,int n_sum,int* num_blocks)
{

    int i;
    int id=1;
    int current_data_id;
    current_data_id=sum[1].cluster_id;
    sum[1].cluster_id=id;

    for(i=2;i<=n_sum;i++)
    {
      if(sum[i].cluster_id>current_data_id){ current_data_id=sum[i].cluster_id; id++;}
      sum[i].cluster_id=id;
    }
    *num_blocks=id;
}
block* explore_block_structure(summary_stat* sum,int nsums,int* num_blocks)
{
   int i;
   int num_cm_blocks;
   get_cluster_ids(sum,nsums,&num_cm_blocks);
   int* blocksizes=i_createVector(num_cm_blocks);
   block* blocks=new block[num_cm_blocks+1];
   for(i=1;i<=num_cm_blocks;i++){ blocksizes[i]=0;}
   get_block_sizes(blocksizes,sum,nsums);
   get_block_objects(blocks,blocksizes,num_cm_blocks,sum,nsums);
   i_destroyVector(blocksizes);
   *num_blocks=num_cm_blocks;
   return blocks;
}

block* read_in(summary_stat* sum,int nsums,int* num_blocks)
{
   int i;
   int num_cm_blocks;
   get_cm_map_summary_stat(sum,nsums,&num_cm_blocks);
   int* blocksizes=i_createVector(num_cm_blocks);
   block* blocks=new block[num_cm_blocks+1];
   for(i=1;i<=num_cm_blocks;i++){ blocksizes[i]=0;}
   get_block_sizes(blocksizes,sum,nsums);
   get_block_objects(blocks,blocksizes,num_cm_blocks,sum,nsums);
   i_destroyVector(blocksizes);
   *num_blocks=num_cm_blocks;
   return blocks;
}
void get_cm_map_summary_stat(summary_stat* sum,int n_sum,int* num_blocks)
{

    int i,j;
    j=1;
	int old_chr;
    int id;
    double tmp_start_cm;
	int num_tmp_block;
	int block_ctr;
	id=1;
	old_chr=1;
	tmp_start_cm=sum[1].cm;
	num_tmp_block=0;
	for(i=1;i<=n_sum;i++)
    {

	   if(num_tmp_block>MAX_win || fabs(tmp_start_cm-sum[i].cm)>=1.0 || sum[i].chr==old_chr+1){ id++; num_tmp_block=0; num_tmp_block=0; tmp_start_cm=sum[i].cm; old_chr=sum[i].chr;}
	   sum[i].cluster_id=id;
	   num_tmp_block++;
	   //cout << i <<  "    " <<sum[i].cluster_id<< "    " << sum[i].chr<<"    " <<sum[i].cm<<endl;

	}
	*num_blocks=id;
}
void get_block_sizes(int* blocksizes,summary_stat* sum,int n_sum)
{
   int i;

   for(i=1;i<=n_sum;i++)
   {
      blocksizes[sum[i].cluster_id]++;
   }
}
void get_block_objects(block* blocks,int* blocksizes,int n_blocks,summary_stat* sum,int n_sum)
{
   int i;
   int num_tot=0;
   int min_blocksize=100000;
   int max_blocksize=0;
   for(i=1;i<=n_blocks;i++)
   {
      //cout << i<<"    "<<blocksizes[i]<<endl;
      num_tot+=blocksizes[i];
	  if(blocksizes[i]>max_blocksize){ max_blocksize=blocksizes[i];}
	  if(blocksizes[i]<min_blocksize){ min_blocksize=blocksizes[i];}

      blocks[i].size=blocksizes[i];
	  blocks[i].indices=i_createVector(blocksizes[i]);
	  blocksizes[i]=1;

   }

   for(i=1;i<=n_sum;i++)
   {
      blocks[sum[i].cluster_id].indices[blocksizes[sum[i].cluster_id]]=i;
	  blocks[sum[i].cluster_id].chr=sum[i].chr;
	  blocksizes[sum[i].cluster_id]++;
	  //cout << i << "   " <<sum[i].chr<<"    " <<sum[i].cluster_id<<"    "<< sum[i].cm<<endl;
   }
   cout << "number of clusters: " << n_blocks<< " maximum cluster size: " << max_blocksize<< " minimum cluster size: " << min_blocksize<<endl;



}
