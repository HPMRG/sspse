/****************************************************************/
/* Compute of log-likelihood */
/****************************************************************/

#include "size.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void bnw_unpos(int *K, int *n, int *s, int *snk, double *Nk, double *mu,
	       double *rho, double *unpos){
    int i, k;
    double N=0., ll=0., Nc, Nkk;
//    Rprintf("K=%d n=%d\n",*K,*n);
	  Nkk=Nk[0];
//    Rprintf("N=%f ldwar=%f\n",Nkk,ldwar(&Nkk,mu,rho));
    for(k=0;k<*K;k++){
	  Nkk=Nk[k];
	  if(Nkk>0){
	    N+=(k+1)*Nkk;
	    ll+=lgammafn(Nkk+1.)-lgammafn(Nkk-snk[k]+1.);
	    ll+=Nkk*ldwar(&Nkk,mu,rho)-lgammafn(Nkk+1.);
    Rprintf("k=%d Nkk=%f ldwar=%f\n",k+1,Nkk,ldwar(&Nkk,mu,rho));
	  }
    }
    ll+=lgammafn(N+1.);
//  Rprintf("N=%f llik=%f\n",N,ll);
    Nc=N;
    for(i=0;i<(*n);i++){
//  Rprintf("ll=%f s=%d Nc=%f\n",ll,s[i],Nc);
	  ll+=log(s[i]/Nc);
	  Nc-=s[i];
    }
//  Rprintf("llik=%f",ll);
    *unpos=ll;
}

void dwarC(int *N, double *mu, double *rho, double *pmf){
    double p;
    p = 1./(*mu);
    *pmf = log(*rho - 1.) + lgammafn(*N + p) + lgammafn(*rho + p) - 
          lgammafn(1. + p) - lgammafn(*N + *rho + p);
}

double ldwar(double *N, double *mu, double *rho){
    double p, lpmf;
    p = 1./(*mu);
    lpmf = log(*rho - 1.) + lgammafn(*N + p) + lgammafn(*rho + p) - 
          lgammafn(1. + p) - lgammafn(*N + *rho + p);
    return(lpmf);
}

void bnw_s_llik(int *K, int *n, int *s, double *Nk, double *llik){
    int i, k;
    double Nc=0., ll=0.;
    int *snk = (int *) malloc(sizeof(int) * (*K));

    for(k=0;k<*K;k++) snk[k]=0.;
    for(i=0;i< *n;i++) snk[s[i] - 1]++;

//    Rprintf("K=%d n=%d\n",*K,*n);
    for(k=0;k<*K;k++){
	  Nc+=(k+1)*Nk[k];
	  ll+=lgammafn(Nk[k]+1.)-lgammafn(Nk[k]-snk[k]+1.);
    }
//    Rprintf("llik=%f\n",ll);
    for(i=0;i<(*n);i++){
//    Rprintf("ll=%f s=%d Nc=%f\n",ll,s[i],Nc);
	  ll+=log(s[i]/Nc);
	  Nc-=s[i];
    }
//    Rprintf("llik=%f",ll);
    free(snk);
    *llik=ll;
}

