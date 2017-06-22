#include "estimator.hpp"
#define EPS pow(10,-13)
#define CM 4
void estimator::covar()
{
   int i,j;
   int B1,B2;
   double** tmpMat=d_createMatrix(parameter_dim,parameter_dim);
   for(i=1;i<=parameter_dim;i++)
   { 
     for(j=1;j<=parameter_dim;j++)
	 {
	   cov[i][j]=0.0;
	 }
   }
   for(j=1;j<=n_cluster;j++)
   {
	  cl[j].get_score(current_parameter,parameter_dim);
   }
   for(i=1;i<=n_cluster;i++)
   { 
       mat_plus_vecxvect(cov,cl[i].score,cl[i].score,parameter_dim);
       
       
	   B1=1;
	   while(i-B1>=1 && cl[i].group==cl[i-B1].group && fabs(cl[i].start_cm-cl[i-B1].end_cm)<=cm_width)
	   {
	     mat_plus_vecxvect(cov,cl[i].score,cl[i-B1].score,parameter_dim);
		 B1++;
	   }
	   B2=1;
	   while(i+B2<=n_cluster && cl[i].group==cl[i+B2].group && fabs(cl[i].end_cm-cl[i+B2].start_cm)<=cm_width)
	   {
	     mat_plus_vecxvect(cov,cl[i].score,cl[i+B2].score,parameter_dim);
		 B2++;
	   }
	   //cout << i << "   " << B1-1<< "    " <<B2-1<< "    " <<cl[i].group<<"     " <<cl[i].start_cm <<"      " <<cl[i].end_cm<<endl;
	   /*for(j=1;j<=B;j++)
	   {
	      if(i-j>=1 && cl[i].group==cl[i-j].group) {mat_plus_vecxvect(cov,cl[i].score,cl[i-j].score,parameter_dim);}
		  if(i+j<=n_cluster && cl[i].group==cl[i+j].group) {mat_plus_vecxvect(cov,cl[i].score,cl[i+j].score,parameter_dim);}
	   }*/
	
	  
   }
   
   get_H_sum();
   
   inverse(H_sum,H_sum_inv,parameter_dim);
   multmat(tmpMat,cov,H_sum_inv,parameter_dim,parameter_dim,parameter_dim);
   multmat(cov_tot,H_sum_inv,tmpMat,parameter_dim,parameter_dim,parameter_dim);
   d_destroyMatrix(tmpMat,parameter_dim);
  
}

void estimator::init(int para_dim)
{
  this->parameter_dim=para_dim;
  this->score=d_createVector(parameter_dim);
  this->obs_fisher=d_createMatrix(parameter_dim,parameter_dim);
  this->obs_fisher_inv=d_createMatrix(parameter_dim,parameter_dim);
  this->cov=d_createMatrix(parameter_dim,parameter_dim);
  this->cov_tot=d_createMatrix(parameter_dim,parameter_dim);
  this->current_parameter=d_createVector(parameter_dim);
  this->new_parameter=d_createVector(parameter_dim);
  this->H_sum=d_createMatrix(parameter_dim,parameter_dim);
  this->H_sum_inv=d_createMatrix(parameter_dim,parameter_dim);

  
}

void estimator::get_score_sum()
{
  int i,j;
  for(i=1;i<=parameter_dim;i++)
  {
     this->score[i]=0.0;
  }
  for(j=1;j<=n_cluster;j++)
  {
	  cl[j].get_score(current_parameter,parameter_dim);
	  for(i=1;i<=parameter_dim;i++)
	  {
			this->score[i]+=cl[j].score[i];
			
	  }
	  //cout << cl[j].score[1]<< "    " <<cl[j].score[2]<<endl;
  }

  
  
  
}
void estimator::estimate(double* startpara,int len)
{  
   int i,j;
   
   for(i=1;i<=len;i++)
   {
     this->current_parameter[i]=startpara[i];
   }
   int conv=0;
   int it=0;
   while(conv==0 && it<=30)
   {
      conv=update_step();
	  //cout << "Iteration: " << it << "     " <<this->current_parameter[1]<< "    "<<this->current_parameter[2]<<endl;
	  it++;
   }
   //cout << "beta_1: " <<this->current_parameter[1]<< "  beta_2: "<<this->current_parameter[2]<<"  #iterations: "<<it<<endl;
    cout << "estimation done. #iterations computed: "<<it<<endl;
}
int estimator::update_step()
{
   int i,j;
   get_score_sum();
   get_H_sum();
   
   double* update_vec=d_createVector(parameter_dim);
   inverse(H_sum,H_sum_inv,parameter_dim);
   
   
   multmatvec(update_vec,H_sum_inv,score,parameter_dim,parameter_dim);
   for(i=1;i<=parameter_dim;i++){ new_parameter[i]=current_parameter[i]+update_vec[i];}
   double diff=nrm2_diff(new_parameter,current_parameter,parameter_dim);
   
   for(i=1;i<=parameter_dim;i++){ current_parameter[i]=new_parameter[i];}
  
   
   d_destroyVector(update_vec);
   if(diff<=EPS){ return 1;} else{ return 0;}
}

void estimator::get_H_sum()
{
  int i,j,k;
  for(i=1;i<=parameter_dim;i++)
  {
     for(j=1;j<=parameter_dim;j++)
	 {
	   this->H_sum[i][j]=0.0;
	 }
  }
  for(k=1;k<=n_cluster;k++)
  {
	  cl[k].get_H(current_parameter,parameter_dim);
	  //if(cl[k].group==6)
	  //{
	  for(i=1;i<=parameter_dim;i++)
	  {
			for(j=1;j<=parameter_dim;j++)
	        {
			   this->H_sum[i][j]+=cl[k].H[i][j];
	        }
	  }
	  //}
  }
  //cout << this->H_sum[1][1]<<"      "<< this->H_sum[2][2]<<endl;
  
}
void estimator::get_H_sum2()
{
  int i,j,k;
  for(i=1;i<=parameter_dim;i++)
  {
     for(j=1;j<=parameter_dim;j++)
	 {
	   this->H_sum[i][j]=0.0;
	 }
  }
  for(k=1;k<=n_cluster;k++)
  {
	  cl[k].get_H(current_parameter,parameter_dim);
	  if(cl[k].group==6)
	  {
	  for(i=1;i<=parameter_dim;i++)
	  {
			for(j=1;j<=parameter_dim;j++)
	        {
			   this->H_sum[i][j]+=cl[k].H[i][j];
	        }
	  }
	  }
  }
 
  
}

double estimator::get_corr_var(int ind1,int ind2,int ind3)
 {
    double sig1=current_parameter[ind1];
	double sig2=current_parameter[ind2];
	double sig12=current_parameter[ind3];
	
    double rho=pow(sig12,2.0)/(sig1*sig2);
  
    double vare=pow(rho,1.0)*(cov_tot[ind1][ind1]/(4*pow(sig1,2.0))+cov_tot[ind2][ind2]/(4*pow(sig2,2.0)));
    vare+=pow(rho,1.0)*cov_tot[ind1][ind2]/(2*sig1*sig2);
    vare+=1/(sig1*sig2)*(cov_tot[ind3][ind3]-cov_tot[ind1][ind3]*sig12/(sig1)-cov_tot[ind2][ind3]*sig12/(sig2));
    return vare; 
}
double estimator::get_corr_var2(int ind1,int ind2,int ind3)
 {
    double sig1=current_parameter[ind1];
	double sig2=current_parameter[ind2];
	double sig12=current_parameter[ind3];
	
    double rho=pow(sig12,2.0)/(sig1*sig2);
  
    double vare=pow(rho,1.0)*(cov_tot[ind1][ind1]/(4*pow(sig1,2.0))+cov_tot[ind2][ind2]/(4*pow(sig2,2.0)));
    //vare+=pow(rho,1.0)*cov_tot[ind1][ind2]/(2*sig1*sig2);
    vare+=1/(sig1*sig2)*(cov_tot[ind3][ind3]);//-cov_tot[ind1][ind3]*sig12/(sig1)-cov_tot[ind2][ind3]*sig12/(sig2));
    return vare; 
}

