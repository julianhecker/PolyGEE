#include <math.h>
#include "util.hpp"
#include <iostream>
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
using namespace std;

#define TINY 1.0e-20

#define BIG 1.0e30
#define NINV 1006
#define Ne 11418
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
/*Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On
output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
a. nrot returns the number of Jacobi rotations that were required.*/
void jacobi(double **a, int n, double *d, double **v, int *nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=d_createVector(n);//vector(1,n);
	z=d_createVector(n);//vector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
		    //delete[] z;
			//delete[] b;
			d_destroyVector(z);
			d_destroyVector(b);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	
}
void multmat_trans2(double** res,double** a,double** b,int n1,int n2,int n3)
{
  int i,j,k;
  for(i=1;i<=n1;i++)
  {
    for(j=1;j<=n3;j++)
	{
	   res[i][j]=0.0;
	   for(k=1;k<=n2;k++)
	   {
	      res[i][j]+=a[i][k]*b[j][k];
	   }
	}
  }
}
void multmat_trans(double** res,double** a,double** b,int n1,int n2,int n3)
{
  int i,j,k;
  for(i=1;i<=n1;i++)
  {
    for(j=1;j<=n3;j++)
	{
	   res[i][j]=0.0;
	   for(k=1;k<=n2;k++)
	   {
	      res[i][j]+=a[k][i]*b[k][j];
	   }
	}
  }
}
void multmat_selftrans(double** res,double** a,int n1,int n2)
{
  int i,j,k;
  for(i=1;i<=n1;i++)
  {
    for(j=1;j<=n1;j++)
	{
	   res[i][j]=0.0;
	   for(k=1;k<=n2;k++)
	   {
	      res[i][j]+=a[i][k]*a[j][k];
	   }
	}
  }
}
void multmat_with_diag_right(double** res,double** a,double* b,int n1,int n2)
{
  int i,j;
  for(i=1;i<=n1;i++)
  {
    for(j=1;j<=n2;j++)
	{
	   res[i][j]=a[i][j]*b[j];
	}
  }
}
void multmat(double** res,double** a,double** b,int n1,int n2,int n3)
{
  int i,j,k;
  for(i=1;i<=n1;i++)
  {
    for(j=1;j<=n3;j++)
	{
	   res[i][j]=0.0;
	   for(k=1;k<=n2;k++)
	   {
	      res[i][j]+=a[i][k]*b[k][j];
	   }
	}
  }
}
void multmat_tr(double* tr,double** a,double** b,int n)
{
  int i,j,k;
  double diag=0.0;
  double tmp;
  for(i=1;i<=n;i++)
  {
    
	   tmp=0.0;
	   for(k=1;k<=n;k++)
	   {
	      tmp+=a[i][k]*b[k][i];
	   }
	   diag+=tmp;
	
  }
  *tr=diag;
}
void multmatvec(double* res,double** a,double* b,int col,int row)
{
  int i,j,k;
  for(i=1;i<=row;i++)
  {
    
	   res[i]=0.0;
	   for(k=1;k<=col;k++)
	   {
	      res[i]+=a[i][k]*b[k];
	   }
	
  }
}
void mat_plus_vecxvect(double** K,double* v1,double* v2,int n)
{
   int i,j;
   for(i=1;i<=n;i++)
   {
     for(j=1;j<=n;j++)
	 {
	    K[i][j]+=v1[i]*v2[j];
		
	 }
   }
}

void multmatvec_tr(double* res,double** a,double* b,int col,int row)
{
  int i,j,k;
  for(i=1;i<=col;i++)
  {
    
	   res[i]=0.0;
	   for(k=1;k<=row;k++)
	   {
	      res[i]+=a[k][i]*b[k];
	   }
	
  }
}


void eigen_symGSL(double** mat,int n,double* d,double** V,double* min, double* max)
{
  int dim=n;
  int i,j;
  double* tmp=d_createVector(n*n);// 1 too much
  for(i=1;i<=n;i++)
  {
    for(j=1;j<=n;j++)
	{
	   tmp[(i-1)*n+j-1]=mat[i][j];
	
	}
	
  }
  gsl_vector *eval = gsl_vector_alloc (dim);
  gsl_matrix *evec = gsl_matrix_alloc (dim, dim);

  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (dim);
  gsl_matrix_view m = gsl_matrix_view_array (tmp, dim, dim);
  
  gsl_eigen_symmv (&m.matrix, eval, evec, w);
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  double eval_i;
  double mineval=9999;
  double maxeval=-9999;
  int ctr=0;
  for(i=1;i<=n;i++)
  {
        
		eval_i = gsl_vector_get (eval, i-1); 
		//cout << eval_i<<endl;
		if(eval_i>=pow(10,-4))
		{
			
			d[i]=eval_i;
			for(j=1;j<=n;j++) V[j][i]=gsl_matrix_get (evec,j-1,i-1);
			ctr++;
			if(eval_i>maxeval){ maxeval=eval_i;}
			if(eval_i<mineval){ mineval=eval_i;}
		}
		
  }
  //*rank=ctr;
  *min=mineval;
  *max=maxeval;
  
  
  
  
 
  gsl_eigen_symmv_free (w);
  d_destroyVector(tmp);
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  

}
 
void GSL_QR_proj(double** mat,int n,int col,double** R)
{
 
  int i,j;
  
  gsl_matrix *m=gsl_matrix_calloc (n,col);
  gsl_matrix *A=gsl_matrix_calloc (n,col);
  for(i=1;i<=n;i++)
  {
    for(j=1;j<=col;j++)
	{
	    gsl_matrix_set (m, i-1, j-1,mat[i][j]);
	    gsl_matrix_set (A, i-1, j-1,mat[i][j]);
	}
	
  }
  gsl_vector* tau = gsl_vector_calloc (col);
  gsl_linalg_QR_decomp (m,tau);
  
  /*gsl_linalg_QR_QTmat (m, tau, A);
  for(i=1;i<=n;i++)
  {
    for(j=1;j<=col;j++)
	{
	    R[i][j]=gsl_matrix_get (A, i-1, j-1);
	   
	}
	
  }*/
  gsl_matrix_free (m);
  gsl_matrix_free (A);
  gsl_vector_free (tau);

}
void GSL_SVD(double** mat,int n,int col,int* rank,double** P)
{
 
  int i,j;
  
  gsl_matrix *m=gsl_matrix_calloc (n,col);
  gsl_matrix *A=gsl_matrix_calloc (n,col);
  gsl_matrix *X=gsl_matrix_calloc (col,col);
  gsl_matrix *V=gsl_matrix_calloc (col,col);
  gsl_vector* S = gsl_vector_calloc (col);
  gsl_vector* work = gsl_vector_calloc (col);
  for(i=1;i<=n;i++)
  {
    for(j=1;j<=col;j++)
	{
	    
	    gsl_matrix_set (A, i-1, j-1,mat[i][j]);
	}
	
  }
  
  //gsl_linalg_SV_decomp ( A, V,  S, work);
  int ctr=0;
  gsl_linalg_SV_decomp_mod (A, X, V, S, work);
  for(j=1;j<=col;j++)
  {
    if(gsl_vector_get(S,j-1)>=pow(10,-4))
	{
	     //cout << gsl_vector_get(S,j-1)<<endl;
		 for(i=1;i<=n;i++)
		 {
			
			 P[i][j]=gsl_matrix_get (A, i-1, j-1);
		 }
		 ctr++;
	}
  }
  *rank=ctr;
  
  gsl_matrix_free (m);
  gsl_matrix_free (A);
  gsl_matrix_free (V);
  gsl_matrix_free (X);
  gsl_vector_free (S);
  gsl_vector_free (work);

}
void inverse(double** mat,double** mat_inv,int dimension)
{
   int i,j;
   int s;
   //cout << "dim  " <<dimension<<endl;
   gsl_permutation * p = gsl_permutation_alloc (dimension);
   gsl_matrix * m = gsl_matrix_alloc (dimension,dimension);
   gsl_matrix * inverse = gsl_matrix_alloc (dimension,dimension);
   
   for(i=1;i<=dimension;i++)
   {
	     for(j=1;j<=dimension;j++)
		 {
		    gsl_matrix_set(m,i-1,j-1,mat[i][j]);
			//cout << S2[i][j]<<"   ";
		 }
		 //cout <<endl;
   }
   gsl_linalg_LU_decomp (m, p, &s);
   gsl_linalg_LU_invert (m, p, inverse);
   for(i=1;i<=dimension;i++)
   {
	     for(j=1;j<=dimension;j++)
		 {
		    mat_inv[i][j]=gsl_matrix_get(inverse,i-1,j-1);
		 }
   }
	
	  
   gsl_permutation_free (p);
   gsl_matrix_free (m);
   gsl_matrix_free (inverse);
}
void get_working_correlation(double** hap,int* inds,int nrows,int len,double scal,double** res)
{
   int i,j,k;
   double** tmpmat=d_createMatrix(len,NINV);
   double** corr=d_createMatrix(len,len);
   double* m=d_createVector(len);
   for(i=1;i<=len;i++)
   {
      for(j=1;j<=NINV;j++)
	  {
	    tmpmat[i][j]=hap[inds[i]][j];
	  }
	  //cout <<inds[i]<<endl;
   }
   for(i=1;i<=len;i++)
   {
      m[i]=0.0;
      for(j=1;j<=NINV;j++)
	  {
	    m[i]+=tmpmat[i][j];
	  }
	  m[i]/=(double)NINV;
	  for(j=1;j<=NINV;j++)
	  {
	    tmpmat[i][j]=(tmpmat[i][j]-m[i])/sqrt(m[i]*(1.0-m[i]));
	  }
	  
   }
   multmat_selftrans(corr,tmpmat,len,NINV);
   //cout.precision(2);
   
 
 
   for(i=1;i<=len;i++)
   {
      for(j=1;j<=len;j++)
	  {
	    corr[i][j]/=(double)NINV;
	    
	    corr[i][j]*=scal;
		if(i==j){ corr[i][j]+=(1.0-scal); }
	  }
	  
	  
   }
 
   for(i=1;i<=len;i++)
   {
      for(j=1;j<=len;j++)
	  {
	    res[i][j]=pow(corr[i][j],2.0);
		
	  }
	 
   }
   d_destroyMatrix(corr,len); d_destroyMatrix(tmpmat,len);
   d_destroyVector(m);
}
void get_working_correlation_recombshrink(double** hap,int* inds,int nrows,int len,double scal,double** res,double* cm)
{
   int i,j,k;
   double** tmpmat=d_createMatrix(len,NINV);
   double** corr=d_createMatrix(len,len);
   double* m=d_createVector(len);
   for(i=1;i<=len;i++)
   {
      for(j=1;j<=NINV;j++)
	  {
	    tmpmat[i][j]=hap[inds[i]][j];
	  }
	  //cout <<inds[i]<<endl;
   }
   for(i=1;i<=len;i++)
   {
      m[i]=0.0;
      for(j=1;j<=NINV;j++)
	  {
	    m[i]+=tmpmat[i][j];
	  }
	  m[i]/=(double)NINV;
	  for(j=1;j<=NINV;j++)
	  {
	    tmpmat[i][j]=(tmpmat[i][j]-m[i])/sqrt(m[i]*(1.0-m[i]));
	  }
	  
   }
   multmat_selftrans(corr,tmpmat,len,NINV);
   //cout.precision(2);
   for(i=1;i<=len;i++)
   {
      for(j=1;j<=len;j++)
	  {
	    corr[i][j]/=(double)NINV;
	    
	    corr[i][j]*=exp(-4*(double)Ne*fabs(cm[i]-cm[j])/((double)100*(double)NINV))*0.95;
		if(i==j){ corr[i][j]+=0.05; }
	  }
	  
	  
   }
 
   for(i=1;i<=len;i++)
   {
      for(j=1;j<=len;j++)
	  {
	    res[i][j]=pow(corr[i][j],2.0);
		
	  }
	 
   }
   d_destroyMatrix(corr,len); d_destroyMatrix(tmpmat,len);
}


void	axpy(double alpha, double* x, double* y, int len){
	int i;
	for(i=1;i<=len;i++){
		y[i]=alpha*x[i] + y[i];
	}
}
void	scal(double alpha, double* x, int len){
	int i;
	for(i=1;i<=len;i++){
		x[i]=alpha*x[i];
	}
}
void  copy(double* x,double* y, int len){
	int i;
	for(i=1;i<=len;i++){
		y[i]=x[i];
	}
}
double	dot(double* x, double* y, int len){
	int i;
	
	double Summe=0;
	for(i=1;i<=len;i++){
				Summe+=x[i]*y[i];
	}
	
	return Summe;
}

double	nrm2(double* x, int len){
	double A=sqrt(dot(x,x,len));
	return A;
}
double	nrm2_diff(double* x,double* y, int len){
    int i;
    double* tmp=d_createVector(len);
	for(i=1;i<=len;i++) tmp[i]=x[i]-y[i];
	double A=sqrt(dot(tmp,tmp,len));
	d_destroyVector(tmp);
	return A;
}







double*	d_createVector(int size)
{
    double *Array=0;
    
	Array=new double[size+1];
	
	
	return Array;
}
void	d_destroyVector(double* vector){
	delete[] vector;
}

int*	i_createVector(int size)
{
	int* Vector;
	Vector=new int[size+1];
	return Vector;
}
void	i_destroyVector(int* vector)
{
	delete[] vector;
}

double**	d_createMatrix(int sizeX, int sizeY){
	int i;
	double **Array=new double*[sizeX+1];
	for(i=0; i<=sizeX; i++){
	        Array[i]=new double[sizeY+1];
	}
	return Array;
}



void	d_destroyMatrix(double** matrix, int rows){
	    int i;
		for(i=0; i<=rows; i++){
			       delete[] matrix[i];
		}
		delete[] matrix;
}














