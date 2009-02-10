/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "size.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void bnw_llik(int *K, int *n, int *s, int *nk, double *Nk, double *llik){
    int i, k;
    double Nc=0., ll=0.;
    for(k=0;k<*K;k++){
	  Nc+=(k+1)*Nk[k];
	  if(Nk[k]<nk[k]){*llik=-1000000.;return;}
	  ll+=lgammafn(Nk[k]+1.)-lgammafn(Nk[k]-nk[k]+1.);
    }
    for(i=0;i<(*n);i++){
	  if(Nc<=0.){*llik=-1000000.;return;}
	  ll+=log(s[i]/Nc);
	  Nc-=s[i];
    }
    *llik=ll;
}
double bnw_llikN(int *K, int *n, int *s, int *nk, int *Nk){
    int i, k, Nkk, Nc;
    double ll;
    ll=0.;
    Nc=0;
    for(k=0;k<(*K);k++){
	  Nkk=Nk[k];
	  Nc+=(k+1)*Nkk;
	  if(Nkk<nk[k]) return(-1000000.);
	  ll+=lgammafn(Nkk+1.)-lgammafn(Nkk-nk[k]+1.);
    }
    for(i=0;i<(*n);i++){
	  if(Nc<=0) return(-1000000.);
	  ll+=log(1.*s[i]/(1.*Nc));
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
    int i, k, Nkk, Mi;
    int nvalidhtn, nvalidhtd, nvalidthistime, minvalid;
    double mhtni, vhtni, mhtdi, vhtdi, htn, htd, q;
    double mhtn, vhtn, mhtd, vhtd;
    double Nd, Nc, lcardN, Kd, Md, lM;
    double *lprob = (double *) malloc(sizeof(double) * (*K));
    double *qprob = (double *) malloc(sizeof(double) * (*K));
    double *lqprob = (double *) malloc(sizeof(double) * (*K));
    double *htnv = (double *) malloc(sizeof(double) * (*M));
    double *htdv = (double *) malloc(sizeof(double) * (*M));

    Mi=(int)(*M);
    Nd=(double)(*N);
    Md=(double)(*M);
    Kd=(double)(*K);
    lM=log(*M);
//  This is the number of populations of size N with K classes
    lcardN=lgammafn(Nd+Kd)-lgammafn(Kd)-lgammafn(Nd+1.);

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
// Rprintf("N=%d\n",*N);
//for(k=0;k<*K;k++){
// Rprintf(" %f",prob[k]);
//}
// Rprintf("\n");

    GetRNGstate();  /* R function enabling uniform RNG */

    mhtni=0.;
    mhtdi=0.;
    vhtni=0.;
    vhtdi=0.;
    mhtn=0.;
    vhtn=0.;
    mhtd=0.;
    vhtd=0.;
    nvalidhtn=0;
    nvalidhtd=0;
    minvalid=trunc(0.05*Mi);
    if(minvalid < 100) minvalid=100;
    while(nvalidhtn < minvalid){
    for(i=0;i<Mi;i++){
//   for(k=0;k<(*K);k++){
//    Nk[k]=0;
//   }
     rmultinom(*N, qprob, *K, Nk);
// Rprintf("N=%d\n",*N);
//for(k=0;k<*K;k++){
// Rprintf(" %d",Nk[k]);
//}
// Rprintf("\n");
//   Nkk=0;
//   for(k=0;k<(*K);k++){
//    Nkk+=(Nk[k]*(k+1));
//   }
//   Nc=bnw_llikN(K,n,s,nk,Nk);
     q=-dmultinorm(N,K,Nk,lqprob);
     htdv[i]=q;
     htnv[i]=bnw_unposN(N,K,n,s,nk,Nk,lprob)+q;
    }
//  Rprintf("mhtni=%f mhtdi=%f\n",log(mhtni),log(mhtdi));
    mhtni=0.;
    mhtdi=0.;
    nvalidthistime=0;
    for(i=0;i<Mi;i++){
      if(htnv[i] > -900000.){
        mhtni+=htnv[i]/Md;
        nvalidthistime++;
      }
      mhtdi+=htdv[i]/Md;
    }
    mhtni=Md*mhtni/nvalidthistime;
    vhtni=0.;
    vhtdi=0.;
    for(i=0;i<Mi;i++){
      if(htnv[i] > -900000.){
        vhtni+=(htnv[i]-mhtni)*(htnv[i]-mhtni)/nvalidthistime;
      }
      vhtdi+=(htdv[i]-mhtdi)*(htdv[i]-mhtdi)/Md;
    }
//  mhtni=log(mhtni)+htnq1-q1-log(mhtdi)+q1+lcardN;
    mhtn=(nvalidhtn*mhtn+nvalidthistime*mhtni)/(nvalidhtn+nvalidthistime);
    vhtn=(nvalidhtn*vhtn+nvalidthistime*vhtni)/(nvalidhtn+nvalidthistime);
    mhtd=(nvalidhtd*mhtd+Mi*mhtdi)/(nvalidhtd+Mi);
    vhtd=(nvalidhtd*vhtd+Mi*vhtdi)/(nvalidhtd+Mi);
    nvalidhtn+=nvalidthistime;
    nvalidhtd+=Mi;
    Rprintf("N=%d minvalid=%d nvalidhtn=%d nvalidhtd=%d\n",
	     *N,minvalid,nvalidhtn, nvalidhtd);
    }
    htn=mhtn+0.5*vhtn+log(1.*nvalidhtn);
    htd=mhtd+0.5*vhtd+lM;
//  htn=mhtn+0.5*vhtn;
//  htd=mhtd+0.5*vhtd;
    *unpos=htn-htd+lcardN;
    Rprintf("N=%d Pct Valid=%f\n",(*N),(nvalidhtn*1.)/Md);
//  htn=htn-htd;
//  htn=log(htn)-log(htd)+lcardN;
//    cpos=log(cpos/(1.*(*M)));
    PutRNGstate();  /* Disable RNG before returning */
    free(lprob);
    free(qprob);
    free(lqprob);
//  free(htnv);
//  free(htdv);
//    Rprintf("K=%d lcardN=%f cpos=%f\n",*K, lcardN,cpos);
// Rprintf("N=%d #pop ties=%d #mean ties %f\n",(*N),Nkk,(Nkk*1.)/(*N));
}
void bnw_mp(int *N, int *lenN, int *K, int *n, int *s, int *nk,
	    double *Cval,
	    double *prob,
	    int *Nprior,
	    int *Nmle,
	    double *mu,
	    double *rho, int *M){
    int i, k, Nkk, Mi, Nlen;
    int Ntotprior, Ntotpriori;
    double fact, lCval, q;
    double Nc, Md, lM;
    double *dprob = (double *) malloc(sizeof(double) * (*K));

    Mi=(int)(*M);
    Md=(double)(*M);
    Nlen=(*lenN);
    lM=log(*M);
    lCval=log(*Cval);

    fact=1.;
    Nc=0.;
    for(k=0;k<*K;k++){
      Nkk=k+1;
      dprob[k]=exp(ldwarint(&Nkk,mu,rho));
      Nc+=dprob[k];
    }
    for(k=0;k<*K;k++){
      dprob[k]=dprob[k]/Nc;
    }
    for(i=0;i<Nlen;i++){
      prob[i]=0.;
    }

    GetRNGstate();  /* R function enabling uniform RNG */

    i=0;
    while(i < Mi){
/*  Generate a random draw from the size prior (here uniform) */
     Ntotpriori = trunc(Nlen*unif_rand()+1.);
     Ntotprior = N[Ntotpriori];
//   Rprintf("N= %d",Ntotprior);
/*  Generate a random draw from the population of sizes */
     rmultinom(Ntotprior, dprob, *K, Nprior);
// for(k=0;k<*K;k++){
//  Rprintf(" %d",Nprior[k]);
// }
//  Rprintf("\n");
//   Nkk=0;
//   for(k=0;k<(*K);k++){
//    Nkk+=(Nk[k]*(k+1));
//   }
//   
//   q=-dmultinorm(N,K,Nk,lprob);
//   
     q=bnw_llikN(K,n,s,nk,Nprior);
     if(q > lCval){
       Rprintf("Warning: Rejection sampling bound log(C)=%f exceeded\n",lCval);
       Rprintf("         by drawn value of %f.\n",q);
       Rprintf("         Resetting the value to 110 percent of the draw.\n");
       i=0;
       fact =1.1;
       lCval=q+log(fact);
       for(k=0;k<(*K);k++){
         Nmle[k]=Nprior[k];
       }
     }
// Rprintf("i=%d Ntotprior=%d lCval=%f q=%f\n",i,Ntotprior, lCval, q);
     if(lCval+log(unif_rand()) < q){
      prob[Ntotpriori]=prob[Ntotpriori]+1;
      i++;
      if( ((10*i) % Mi)==0 && Mi > 500){
       Rprintf("Sampled %d from %d\n", i, Mi);}
     }

    }
    for(k=0;k<Nlen;k++){
        prob[k]=prob[k]/Md;
    }
    *Cval=exp(lCval)/fact;
    PutRNGstate();  /* Disable RNG before returning */
    free(dprob);
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
	  if(Nkk<nk[k]) return(-1000000.);
	  if(Nkk>0){
	    ll+=lgammafn(Nkk+1.)-lgammafn(Nkk-nk[k]+1.);
	    ll+=Nkk*lprob[k]-lgammafn(Nkk+1.);
	  }
    }
    for(i=0;i<(*n);i++){
          if(Nc<=0) return(-1000000.);
	  ll+=log(1.*s[i]/(1.*Nc));
	  Nc-=s[i];
    }
    return(ll);
}
