#ifndef _CMATRIX_HPP_
#define _CMATRIX_HPP_
#include <stddef.h>
#include <iostream>
using namespace std;
class cmatrix
{

   public:
   int dimension;
   double** matrix;
   double** inv_matrix;
   
    cmatrix ();
    cmatrix (int);
	~cmatrix();
	void invert();
    void init(int dimensionp);
	double smallest_eigenvalue(double* mine, double* maxe);


};
#endif /* _CMATRIX_HPP_ */