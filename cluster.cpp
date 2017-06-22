#include "cluster.hpp"
#include "util.hpp"
cluster::cluster ()
{
  this->parameter_dim=0;
  this->size=0;
  this->reg=NULL;
  this->y=NULL;
  this->constants=NULL;
  this->eta=NULL;
  this->var=NULL;
  this->res=NULL;
  this->score=NULL;
  this->D=NULL;
  this->H=NULL;
  
}
cluster::cluster (int sizep,int parameter_dimp):W(sizep)
{
  this->parameter_dim=parameter_dimp;
  this->size=sizep;
  this->reg=d_createMatrix(size,parameter_dim);
  this->y=d_createVector(size);
  this->constants=d_createVector(size);
  this->res=d_createVector(size);
  this->eta=d_createVector(size);
  this->var=d_createVector(size);
  this->score=d_createVector(parameter_dim);
  this->D=d_createMatrix(parameter_dim,size);
  this->H=d_createMatrix(parameter_dim,parameter_dim);
 
  
}
cluster::~cluster ()
{
  if(this->reg!=NULL)d_destroyMatrix(this->reg,size);
  
  if(this->y!=NULL)d_destroyVector(this->y);
  if(this->constants!=NULL)d_destroyVector(this->constants);
  if(this->var!=NULL)d_destroyVector(this->var);
  if(this->eta!=NULL)d_destroyVector(this->eta);
  if(this->res!=NULL)d_destroyVector(this->res);
  if(this->score!=NULL)d_destroyVector(this->score);
  if(this->D!=NULL)d_destroyMatrix(this->D,parameter_dim);
  if(this->H!=NULL)d_destroyMatrix(this->H,parameter_dim);
  
}
void cluster::init(int sizep,int parameter_dimp)
{
  this->parameter_dim=parameter_dimp;
  this->size=sizep;
  this->reg=d_createMatrix(size,parameter_dim);
  this->y=d_createVector(size);
  this->constants=d_createVector(size);
  this->res=d_createVector(size);
  this->eta=d_createVector(size);
  this->var=d_createVector(size);
  this->score=d_createVector(parameter_dim);
  this->D=d_createMatrix(parameter_dim,size);
  this->H=d_createMatrix(parameter_dim,parameter_dim);
  if(bivar==0){this->W.init(sizep);}
  if(bivar==1){this->W.init(sizep/3);}
}
void cluster::get_eta(double* parameter,int len)
{
   int i,j;
   if(len==parameter_dim)
   {
	   for(i=1;i<=size;i++)
	   {
		 eta[i]=0.0;
		 for(j=1;j<=parameter_dim;j++)
		 {
		   eta[i]+=reg[i][j]*parameter[j];
		 }
		 eta[i]+=constants[i];
		 //cout << "beta  " << eta[i]<<endl;
	   }
	}
}
void cluster::get_res()
{
   int i,j;
   for(i=1;i<=size;i++)
   {
     res[i]=y[i]-eta[i];
	 //cout << "res  " << res[i]<<endl;
   }
}
void cluster::get_var_gamma()
{
   int i,j;
   for(i=1;i<=size;i++)
   {
     var[i]=pow(eta[i],2.0);
	 //cout << "var  " << var[i]<<endl;
   }
}
void cluster::get_D_T_gamma() //load D^T
{
   int i,j;
   double** tmpMat=d_createMatrix(parameter_dim,size);
   double* tmpVec=d_createVector(size);
   for(i=1;i<=parameter_dim;i++) //get X^T
   {
     for(j=1;j<=size;j++)
	 {
	    tmpMat[i][j]=reg[j][i];
	 }
	 
   }
   for(j=1;j<=size;j++) //u^' * A^0.5
   {
	    tmpVec[j]=pow(eta[j],-2.0)*sqrt(var[j]);
   }
   multmat_with_diag_right(D,tmpMat,tmpVec,parameter_dim,size);
   
   d_destroyVector(tmpVec);
   d_destroyMatrix(tmpMat,parameter_dim);
}
void cluster::get_score(double* parameter,int len)
{
   int i,j;
   double* tmpVec=d_createVector(size);
   double* tmpRes=d_createVector(size);
   
   
   //cout << "eta" <<endl;
   get_eta(parameter,len);
   //cout << "res" <<endl;
   get_res();
   //cout << "var "<<endl;
   get_var_gamma();
   //cout << "D" <<endl;
   get_D_T_gamma();
   
   for(j=1;j<=size;j++)
   {
	 if(fabs(var[j])<=pow(10,-6)){ cout << "error. variance equal to zero." <<endl; }		
	 tmpRes[j]=res[j]/sqrt(var[j]);
   }
   
   //multmatvec(tmpVec,W.inv_matrix,tmpRes,size,size);
   //cout << "inv" <<endl;
   mult_res_with_inv_matrix(tmpVec,tmpRes);  
   //cout << "mult" <<endl;
   multmatvec(score,D,tmpVec,size,parameter_dim);
	

   d_destroyVector(tmpVec);
   d_destroyVector(tmpRes);
}
void cluster::mult_res_with_inv_matrix(double* tmpVec,double* tmpRes)
{
  int i,j;
  int rsize;
  if(bivar==0){ multmatvec(tmpVec,W.inv_matrix,tmpRes,size,size);}
  if(bivar==1)
  {
      rsize=size/3;
      double* vec1=d_createVector(rsize);
	  double* res1=d_createVector(rsize);
	  double* vec2=d_createVector(rsize);
	  double* res2=d_createVector(rsize);
	  double* vec3=d_createVector(rsize);
	  double* res3=d_createVector(rsize);
      for(i=1;i<=rsize;i++)
	  { 
	     res1[i]=tmpRes[i];
		 res2[i]=tmpRes[i+rsize];
		 res3[i]=tmpRes[i+2*rsize];
	  }
	  multmatvec(vec1,W.inv_matrix,res1,rsize,rsize);
	  multmatvec(vec2,W.inv_matrix,res2,rsize,rsize);
	  multmatvec(vec3,W.inv_matrix,res3,rsize,rsize);
	  for(i=1;i<=rsize;i++)
	  { 
	    tmpVec[i]=vec1[i];//1.5*vec1[i]+0.5*vec2[i]-vec3[i];
	  }
	  
	  for(i=1;i<=rsize;i++)
	  { 
	    tmpVec[i+rsize]=vec2[i];//0.5*vec1[i]+1.5*vec2[i]-vec3[i];
	  }
	  for(i=1;i<=rsize;i++)
	  { 
	    tmpVec[i+2*rsize]=vec3[i];//-vec1[i]-vec2[i]+2.0*vec3[i];
	  }
	  
	  d_destroyVector(vec1);d_destroyVector(vec2);d_destroyVector(vec3);
	  d_destroyVector(res1);d_destroyVector(res2);d_destroyVector(res3);
	 
  }
}
void cluster::dispersion_variance(double disp)
{
   int i,j;
   this->dispersion=disp;
   for(i=1;i<=size;i++)
   {
     for(j=1;j<=size;j++)
	 {
	    W.inv_matrix[i][j]/=dispersion;
	 }
   }
}
void cluster::get_H(double* parameter, int len) //load D^T
{
   int i,j;
   get_eta(parameter,len);
   get_res();
   get_var_gamma();
   get_D_T_gamma();
   int rsize;
   double** tmpMat=d_createMatrix(size,parameter_dim);
   if(bivar==1)
   {
      rsize=size/3;
      double** tmpWinv=d_createMatrix(size,size);
	  for(i=1;i<=rsize;i++)
	  {
	    for(j=1;j<=rsize;j++)
		{
		   tmpWinv[i][j]=W.inv_matrix[i][j];//1.5*W.inv_matrix[i][j];
		   tmpWinv[rsize+i][rsize+j]=W.inv_matrix[i][j];//1.5*W.inv_matrix[i][j];
		   tmpWinv[2*rsize+i][2*rsize+j]=W.inv_matrix[i][j];//2.0*W.inv_matrix[i][j];
		   
		   tmpWinv[rsize+i][j]=0.0;//0.5*W.inv_matrix[i][j];
		   tmpWinv[i][rsize+j]=0.0;//0.5*W.inv_matrix[i][j];
		   tmpWinv[2*rsize+i][j]=0.0;//-1.0*W.inv_matrix[i][j];
		   tmpWinv[i][2*rsize+j]=0.0;//-1.0*W.inv_matrix[i][j];
		   tmpWinv[2*rsize+i][rsize+j]=0.0;//-1.0*W.inv_matrix[i][j];
		   tmpWinv[rsize+i][2*rsize+j]=0.0;//-1.0*W.inv_matrix[i][j];
		   
		}
	  }
	  multmat_trans2(tmpMat,tmpWinv,D,size,size,parameter_dim);
	  d_destroyMatrix(tmpWinv,size);
   }
   if(bivar==0) {multmat_trans2(tmpMat,W.inv_matrix,D,size,size,parameter_dim);}
   multmat(H,D,tmpMat,parameter_dim,size,parameter_dim);
   

   d_destroyMatrix(tmpMat,size);
}