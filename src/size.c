/****************************************************************/
/* Compute of log-likelihood */
/****************************************************************/

#include "size.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

double bnw_llikN(int *N, int *K, int *n, int *s, int *nk, int *Nk){
    int i, k, Nkk, Nc, ddd;
    double ll;
    ll=0.;
    Nc=0;
//  ddd=0;
//  Rprintf("K=%d n=%d nk[0]=%d\n",*K,*n,nk[0]);
//    for(k=0;k<(*K);k++){
//	  ddd+=Nk[k];
//    }
//  Rprintf("N= %d\n",ddd);
    for(k=0;k<(*K);k++){
	  Nkk=Nk[k];
	  Nc+=(k+1)*Nkk;
    Rprintf("k=%d Nkk=%d nk[k]=%d\n",k,Nkk,nk[k]);
	  if(Nkk<nk[k]) return(-10000.);
	  ll+=lgammafn(Nkk+1.)-lgammafn(Nkk-nk[k]+1.);
    }
//    Nc=*N;
//  Rprintf("K=%d n=%d ll=%f Nc=%d\n",*K,*n,ll,Nc);
    for(i=0;i<(*n);i++){
	  if(Nc<=0) return(-10000.);
	  ll+=log(1.*s[i]/(1.*Nc));
  Rprintf("i=%d s[i]=%d Nc=%d\n",i,s[i],Nc);
	  Nc-=s[i];
    }
    return(ll);
}

double dmultinorm(int *N, int *K, int *Nk, double *lprob){
    int k, Nkk;
    double ll;
    ll=lgammafn((*N)+1.);
    for(k=0;k<*K;k++){
      Nkk=Nk[k];
      if(Nkk>0){
        ll+=Nkk*lprob[k]-lgammafn(Nkk+1.);
      }
    }
    return(ll);
}

void bnw_NC(int *N, int *K, int *n, int *s, int *nk, int *Nk,
	    double *prob,
	    double *mu,
	    double *rho, int *M, double *unpos){
    int i, k, Nkk;
    double htn, htd, htnq, q, q1, htnq1;
    double htn2, htd2;
    double Nd, Nc, Mc, Kd, lM;
    double *lprob = (double *) malloc(sizeof(double) * (*K));
    double *qprob = (double *) malloc(sizeof(double) * (*K));
    double *lqprob = (double *) malloc(sizeof(double) * (*K));
//    Rprintf("K=%d n=%d\n",*K,*n);
    Nd=(double)(*N);
    Kd=(double)(*K);
    lM=log(*M);
//    Rprintf("N=%f ldwar=%f\n",Nkk,ldwar(&Nkk,mu,rho));
    Mc=lgammafn(Nd+Kd)-lgammafn(Kd)-lgammafn(Nd+1.);
//   Rprintf("K=%d Mc=%f\n",*K, Mc);
    Nc=0.;
    for(k=0;k<*K;k++){
      Nkk=k+1;
      lprob[k]=ldwarint(&Nkk,mu,rho);
      prob[k]=exp(lprob[k]);
      Nc+=prob[k];
//   Rprintf("k=%d prob=%f\n",k,prob[k]);
    }
    for(k=0;k<*K;k++){
      prob[k]=prob[k]/Nc;
      lprob[k]=log(prob[k]);
    }
    Nc=0.;
    for(k=0;k<*K;k++){
      qprob[k]=prob[k]*(k+1);
      Nc+=qprob[k];
//   Rprintf("k=%d prob=%f\n",k,prob[k]);
    }
    for(k=0;k<*K;k++){
      qprob[k]=qprob[k]/Nc;
      lqprob[k]=log(qprob[k]);
    }
/*  Generate a random draw from the population of sizes */
   Rprintf("N=%d\n",*N);
  for(k=0;k<*K;k++){
   Rprintf(" %f",prob[k]);
  }
   Rprintf("\n");

    GetRNGstate();  /* R function enabling uniform RNG */

    htn=0.;
    htd=0.;
    htn2=0.;
    htd2=0.;
    for(i=0;i<(*M);i++){
//   for(k=0;k<(*K);k++){
//    Nk[k]=0;
//   }
     rmultinom(*N, qprob, *K, Nk);
   Rprintf("N=%d\n",*N);
  for(k=0;k<*K;k++){
   Rprintf(" %d",Nk[k]);
  }
   Rprintf("\n");
     Nkk=0;
     for(k=0;k<(*K);k++){
      Nkk+=(Nk[k]*(k+1));
     }
//   Nc=bnw_llikN(N,K,n,s,nk,Nk);
     q=-dmultinorm(N,K,Nk,lqprob);
     Rprintf("q=%f\n",q);
     htnq=bnw_unposN(N,K,n,s,nk,Nk,lprob)+q;
   if(i == 0) {
    q1=q;
    htnq1=htnq;
   }
//   if(htnq > -9999.) htn+=exp(htnq-q-htnq1+q1);
//   htd+=exp(-q+q1);
//   htn+=exp(htnq-q-htnq1+q1);
     htn+=htnq;
     htn2+=htnq*htnq;
//   if(htnq > -9999.) htn+=exp(htnq-q);
//   htd+=exp(-q+q1);
     htd+=q;
     htd2+=(q)*(q);
//   Rprintf("i=%d Nk=%d ll=%f %f\n",i+1,Nk[0],ll,bnw_unposN(N,K,n,s,nk,Nk,lprob)-Mc);
    }
    Rprintf("htn=%f htd=%f\n",log(htn),log(htd));
//  htn=log(htn)+htnq1-q1-log(htd)+q1+Mc;
    q1=htn/(*M);
    htn=q1+0.5*(htn2/(*M)-q1*q1);
    q1=htd/(*M);
    htd=q1+0.5*(htd2/(*M)-q1*q1);
    htn=htn-htd+Mc;
//  htn=log(htn)-log(htd)+Mc;
//    cpos=log(cpos/(1.*(*M)));
    PutRNGstate();  /* Disable RNG before returning */
    free(lprob);
//    Rprintf("K=%d Mc=%f cpos=%f\n",*K, Mc,cpos);
   Rprintf("N=%d #pop ties=%d #mean ties %f\n",(*N),Nkk,(Nkk*1.)/(*N));
    *unpos=htn;
}

double ldwarint(int *N, double *mu, double *rho){
    double p, lpmf, r;
    r = (*rho);
    p = (*mu)*(r-2.) - r + 1.;
    lpmf = log(r - 1.) + lgammafn(*N + p) + lgammafn(r + p) - 
           lgammafn(1. + p) - lgammafn(*N + r + p);
    return(lpmf);
}

double bnw_unposN(int *N, int *K, int *n, int *s, int *nk, int *Nk,
		  double *lprob){
    int i, k, Nkk, Nc;
    double ll;
    Nc=0;
    ll=lgammafn((*N)+1.);
    for(k=0;k<*K;k++){
	  Nkk=Nk[k];
	  Nc+=(k+1)*Nkk;
//    Rprintf("k=%d Nkk=%d nk[k]=%d\n",k,Nkk,nk[k]);
	  if(Nkk<nk[k]) return(-10000.);
	  if(Nkk>0){
	    ll+=lgammafn(Nkk+1.)-lgammafn(Nkk-nk[k]+1.);
	    ll+=Nkk*lprob[k]-lgammafn(Nkk+1.);
	  }
    }
//  Nc=*N;
    for(i=0;i<(*n);i++){
//    Rprintf("i=%d Nc=%d s[i]=%d p=%f\n",i,Nc,s[i],log(1.*s[i]/(1.*Nc)));
          if(Nc<=0) return(-10000.);
	  ll+=log(1.*s[i]/(1.*Nc));
	  Nc-=s[i];
    }
//  Rprintf("N=%d ll=%f\n",*N,ll);
    return(ll);
}
