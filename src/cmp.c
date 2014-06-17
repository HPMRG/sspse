/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "cmp.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

double cmp(int x, double llambda, double nu, double lzcmp, int give_log)
  {
     double dev;
     dev=x*llambda-nu*lgamma(x+1.0)-lzcmp;
     if(give_log){
      return( dev );
     }else{
      return( exp(dev) );
     }
  }

double zcmp(double lambda, double nu, double err, int K, int give_log)
  {
     double mj=0.0, z, aj;
     int j;
     z = 1.0 + lambda;
     aj = lambda;
     for (j = 2; j < K; j++){
       mj=lambda/pow((double)j,nu);
       aj*=mj;
       z+=aj;
//Rprintf("j %d\n", j);
     }
//Rprintf("cmp terms %d aj %f err %f\n", j, aj, err*(1-mj));
//   while (aj < err*(1.-mj) && j < 1000){
//   
//   Next is relative error which we ignore
     while (aj > err*z && j < 200){
//if(j > 180){
//       	Rprintf("cmp terms %d aj %f err*z %f z %f\n", j, aj, err*z, z);
//Rprintf("nu %f lambda %f\n", nu, lambda);
//}
       j++;
       mj=lambda/pow((double)j,nu);
       aj*=mj;
       z+=aj;
//Rprintf("j %d\n", j);
     }
//
//
//if(j==1000 | (aj > 1.)){
if(j>=200 && nu > 0.01){
       mj=pow(lambda,1.0/nu);
       aj=pow(mj,(1.0-nu)/2.0)*exp(nu*mj)/(pow(2.0*M_PI,(nu-1.0)/2.0)*sqrt(nu));
//Rprintf("nu %f lambda %f j %d ztilde %f z %f\n", nu, lambda, j, aj, z);
       if(aj > 1.0 + lambda){z=aj;}
//}else{
//Rprintf("nu %f lambda %f j %d z %f\n", nu, lambda, j, z);
//if((mj > 1.) | (aj > 1.)){ Rprintf("mj %f aj %f\n", mj, aj); }
//if((mj > 1.) | (aj > 1.)){ return(-200000.0); }
}

//   Add approx to remainder term
//   mj=lambda/pow((double)(j+1),nu);
//   aj*=mj;
//   z+=(aj/(1.-mj));
//if(j!=100){
//Rprintf("cmp terms %d aj %e err %e bdd %e z %e\n", j, aj, err*(1-mj), aj/(1.-mj),z);
//}
//
     if(give_log){
      return( log(z) );
     }else{
      return( z );
     }
  }

void dcmp (int *x, double *lambda, double *nu, int *n, double *err, int *give_log, double *val) {
  int i, give_log1=1, xmax;
  double lzcmp;
  xmax=x[(*n)-1];
  for (i = 0; i < (*n)-1; i++){
    if(x[i] > xmax){xmax=x[i];}
  }
  lzcmp = zcmp(*lambda, *nu, *err, 2*xmax, give_log1);
//Rprintf("lzcmp %f \n", lzcmp);
  for (i = 0; i < *n; i++){
    val[i] = cmp(x[i], log(*lambda), *nu, lzcmp, *give_log);
  }
}

void rcmp (int *x, double *lambda, double *nu, int *n, int *K, double *err) {
  int i, Ki, ni, popi, give_log0=0, give_log1=1;
  double gb, lzcmp, llambda;
  double *pi = (double *) malloc(sizeof(double) * (*K));
  llambda = log(*lambda);
  GetRNGstate();  /* R function enabling uniform RNG */
  Ki = (*K);
  ni = (*n);
  lzcmp = zcmp(*lambda, *nu, *err, 2*Ki, give_log1);
  for (i = 0; i < Ki; i++){
    pi[i] = cmp(i, llambda, *nu, lzcmp, give_log0);
  }
  for (i=1; i< Ki; i++){
    pi[i]=pi[i-1]+pi[i];
  }
  for (i = 0; i < ni; i++){
   gb = pi[Ki-1] * unif_rand();
   popi = 0;
   while(gb > pi[popi]){popi++;}
   x[i] = popi;
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(pi);
}
