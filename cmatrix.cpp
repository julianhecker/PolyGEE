#include "cluster.hpp"
#include "util.hpp"
cmatrix::cmatrix()
{
  this->dimension=0;
  this->matrix=NULL;
  this->inv_matrix=NULL;
}
cmatrix::cmatrix(int dimensionp)
{
  this->dimension=dimensionp;
  
  this->matrix=d_createMatrix(dimension,dimension);
  this->inv_matrix=d_createMatrix(dimension,dimension);
  
}
cmatrix::~cmatrix ()
{
  if(this->matrix!=NULL)d_destroyMatrix(this->matrix,dimension);
  
  if(this->inv_matrix!=NULL)d_destroyMatrix(this->inv_matrix,dimension);
  
}
double cmatrix::smallest_eigenvalue(double* mine, double* maxe)
{
    double* tmp_d=d_createVector(dimension);
	double** tmpmat=d_createMatrix(dimension,dimension);
	double max,min;
    eigen_symGSL(matrix,dimension,tmp_d,tmpmat,&min,&max);
	*mine=min;
	*maxe=max;
	return 0;
}
void cmatrix::init(int dimensionp)
{
  int i,j;
  this->dimension=dimensionp;
  //cout << "dim init  "<<dimension<<endl;
  this->matrix=d_createMatrix(dimension,dimension);
  this->inv_matrix=d_createMatrix(dimension,dimension);
  for(i=1;i<=this->dimension;i++)
  { 
     for(j=1;j<=this->dimension;j++)
	 {
	   this->matrix[i][j]=0.0; this->inv_matrix[i][j]=0.0;
     }
  }
  
}
void cmatrix::invert()
{
  //cout << "dim invert  " << dimension<<endl; 
  inverse(matrix,inv_matrix,dimension);
  
}