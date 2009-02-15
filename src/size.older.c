/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "size.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void bnw_NCwar(int *N, int *K, int *n, int *s, int *nk, int *Nk,
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
	    double *lbound,
	    double *dprob,
	    double *prob,
	    double *NtotMLE,
	    int *Nprior,
	    int *Nmle,
	    int *M){
    int i, j, k, Mi, Nlen;
    int Ntotprior, Ntotpriori;
    double fact, bound, q;
    double Md, lM;

    Mi=(int)(*M);
    Md=(double)(*M);
    Nlen=(*lenN);
    lM=log(*M);
    bound=(*lbound);

    fact=1.;
    for(i=0;i<Nlen;i++){
      prob[i]=0.;
      NtotMLE[i]=-1000000.;
    }

    GetRNGstate();  /* R function enabling uniform RNG */

    i=0;
    j=0;
    while(i < Mi){
/*  Generate a random draw from the size prior (here uniform) */
     Ntotpriori = trunc(Nlen*unif_rand());
     Ntotprior = N[Ntotpriori];
/*   Generate a random draw from the population of sizes */
     rmultinom(Ntotprior, dprob, *K, Nprior);
     q=bnw_llikN(K,n,s,nk,Nprior);
     j++;
     if(q > bound){
       Rprintf("Warning: Rejection sampling bound log(C)=%f exceeded\n",bound);
       Rprintf("         by drawn value of %f.\n",q);
       Rprintf("         Resetting the value to 110 percent of the draw.\n");
       i=0;
       fact=1.1;
       bound=q+log(fact);
       for(k=0;k<(*K);k++){
         Nmle[k]=Nprior[k];
       }
     }
// Rprintf("i=%d Ntotprior=%d lbound=%f q=%f\n",i,Ntotprior, lbound, q);
      if( (10000*(j / 10000))==j){
       Rprintf("Sampled %d i=%d N=%d =%f bound=%f\n",j,i,Ntotprior,q,bound);}
     if(bound+log(unif_rand()) < q){
      prob[Ntotpriori]=prob[Ntotpriori]+1;
      i++;
      if( ((10*i) % Mi)==0 && Mi > 500){
       Rprintf("Sampled %d from %d\n", i, Mi);}
     }
     if(q > NtotMLE[Ntotpriori]){NtotMLE[Ntotpriori]=q;}
    }
    for(k=0;k<Nlen;k++){
        prob[k]=prob[k]/Md;
    }
//  *bound=exp(lbound)/fact;
    *lbound=bound-log(fact);
    PutRNGstate();  /* Disable RNG before returning */
}

void bnw_mpwar(int *N, int *lenN, int *K, int *n, int *s, int *nk,
	    double *lbound,
	    double *prob,
	    double *NtotMLE,
	    int *Nprior,
	    int *Nmle,
	    double *mu,
	    double *rho, int *M){
    int i, k, Nkk, Mi, Nlen;
    int Ntotprior, Ntotpriori;
    double fact, bound, q;
    double Nc, Md, lM;
    double *dprob = (double *) malloc(sizeof(double) * (*K));

    Mi=(int)(*M);
    Md=(double)(*M);
    Nlen=(*lenN);
    lM=log(*M);
    bound=(*lbound);

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
      NtotMLE[i]=-1000000.;
    }

    GetRNGstate();  /* R function enabling uniform RNG */

    i=0;
    while(i < Mi){
/*  Generate a random draw from the size prior (here uniform) */
     Ntotpriori = trunc(Nlen*unif_rand());
     Ntotprior = N[Ntotpriori];
/*  Generate a random draw from the population of sizes */
     rmultinom(Ntotprior, dprob, *K, Nprior);
// for(k=0;k<*K;k++){
//  Rprintf(" %d",Nprior[k]);
// }
//  Rprintf("\n");
//   Rprintf("N= %d",Ntotprior);
//   Nkk=0;
//   for(k=0;k<(*K);k++){
//    Nkk+=(Nk[k]*(k+1));
//   }
//   
//   q=-dmultinorm(N,K,Nk,lprob);
//   
     q=bnw_llikN(K,n,s,nk,Nprior);
     if(q > bound){
       Rprintf("Warning: Rejection sampling bound log(C)=%f exceeded\n",bound);
       Rprintf("         by drawn value of %f.\n",q);
       Rprintf("         Resetting the value to 110 percent of the draw.\n");
       i=0;
       fact=1.1;
       bound=q+log(fact);
       for(k=0;k<(*K);k++){
         Nmle[k]=Nprior[k];
       }
     }
// Rprintf("i=%d Ntotprior=%d lbound=%f q=%f\n",i,Ntotprior, lbound, q);
     if(bound+log(unif_rand()) < q){
      prob[Ntotpriori]=prob[Ntotpriori]+1;
      i++;
      if( ((10*i) % Mi)==0 && Mi > 500){
       Rprintf("Sampled %d from %d\n", i, Mi);}
     }
     if(q > NtotMLE[Ntotpriori]){NtotMLE[Ntotpriori]=q;}
    }
    for(k=0;k<Nlen;k++){
        prob[k]=prob[k]/Md;
    }
//  *bound=exp(lbound)/fact;
    *lbound=bound-log(fact);
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
