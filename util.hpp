#ifndef _UTIL_HPP_
#define _UTIL_HPP_

void jacobi(double **a, int n, double *d, double **v, int *nrot);
void multmat_trans(double** res,double** a,double** b,int n1,int n2,int n3);
void multmat_selftrans(double** res,double** a,int n1,int n2);
void multmat_with_diag_right(double** res,double** a,double* b,int n1,int n2);
void multmat(double** res,double** a,double** b,int n1,int n2,int n3);
void multmat_tr(double* tr,double** a,double** b,int n);
void multmat_trans2(double** res,double** a,double** b,int n1,int n2,int n3);
void multmatvec(double* res,double** a,double* b,int col,int row);
void mat_plus_vecxvect(double** K,double* v1,double* v2,int n);

void multmatvec_tr(double* res,double** a,double* b,int col,int row);

//void eigen_symGSL(double** mat,int n,double* d,double** V,int* rank);
void eigen_symGSL(double** mat,int n,double* d,double** V,double* min, double* max);

void GSL_QR_proj(double** mat,int n,int col,double** R);
void GSL_SVD(double** mat,int n,int col,int* rank,double** P);
void inverse(double** mat,double** mat_inv,int dimension);

void	axpy(double alpha, double* x, double* y, int len);
void	scal(double alpha, double* x, int len);
void  copy(double* x,double* y, int len);
double	dot(double* x, double* y, int len);

double	nrm2(double* x, int len);
double	nrm2_diff(double* x,double* y, int len);

double*	d_createVector(int size);
void	d_destroyVector(double* vector);
int*	i_createVector(int size);
void	i_destroyVector(int* vector);
double**	d_createMatrix(int sizeX, int sizeY);
void	d_destroyMatrix(double** matrix, int rows);

void get_working_correlation_recombshrink(double** hap,int* inds,int nrows,int len,double scal,double** res,double* cm);
void get_working_correlation(double** hap,int* inds,int nrows,int len,double scal,double** res);

#endif /* _UTIL_HPP_ */
