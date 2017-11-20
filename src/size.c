/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "size.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void bnw_llik(int *K, int *n, int *s, int *nk, double *Nk, double *llik){
    int i, k;
    double totsize, ll, Nkk;
    totsize=0.;
    ll=0.;
    for(k=0;k<(*K);k++){
      Nkk=Nk[k];
      if(Nkk>0.){
        totsize+=(k+1)*Nkk;
        ll+=lgammafn(Nkk+1.)-lgammafn(Nkk-nk[k]+1.);
      }
    }
    for(i=0;i<(*n);i++){
      ll+=log(s[i]/totsize);
      totsize-=s[i];
    }
    *llik=ll;
}
double bnw_llikN(int *K, int *n, int *s, int *nk, int *Nk){
// K: number of classes
// n: sample size 
// s: sample in sequential order (vector) 
// nk: sample counts (vector) 
// Nk: population counts (vector) 
    int i, k, Nkk, totsize;
    double ll;
    ll=0.;
    totsize=0;
    for(k=0;k<(*K);k++){
      if(Nk[k]<nk[k]) return(-1000000.);
    }
    for(k=0;k<(*K);k++){
      Nkk=Nk[k];
      if(Nkk>0){
        totsize+=(k+1)*Nkk;
        ll+=lgammafn(Nkk+1.)-lgammafn(Nkk-nk[k]+1.);
      }
    }
    for(i=0;i<(*n);i++){
//	  if(totsize<=0) return(-1000000.);
	  ll+=log(s[i]/((double)totsize));
	  totsize-=s[i];
    }
    return(ll);
}
double bnw_llikNf(int *K, int *n, int *s, int *nk, int *Nk){
// No check
// K: number of classes
// n: sample size 
// s: sample in sequential order (vector) 
// nk: sample counts (vector) 
// Nk: population counts (vector) 
    int i, k, Nkk, totsize;
    double ll,dtotsize;
    ll=0.;
    totsize=0;
    for(k=0;k<(*K);k++){
      Nkk=Nk[k];
      if(Nkk>0){
        totsize+=(k+1)*Nkk;
        ll+=lgammafn(Nkk+1.)-lgammafn(Nkk-nk[k]+1.);
      }
    }
    dtotsize=(double)totsize;
    for(i=0;i<(*n);i++){
      ll+=log(s[i]/dtotsize);
      dtotsize-=s[i];
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
	    double *qprob,
	    int *M, double *unpos){
    int i, k, Ki, Ni, Mi, Nni;
    int nvalidhtn, nvalidhtd, nvalidthistime, minvalid;
    double mhtni, vhtni, mhtdi, vhtdi, htn, htd, q, p;
    double mhtn, vhtn, mhtd, vhtd;
    double Nd, lcardN, Kd, Md, lM;
    double *lprob = (double *) malloc(sizeof(double) * (*K));
    double *lqprob = (double *) malloc(sizeof(double) * (*K));
    double *htnv = (double *) malloc(sizeof(double) * (*M));
    double *htdv = (double *) malloc(sizeof(double) * (*M));

    Ni=(int)(*N);
    Ki=(int)(*K);
    Mi=(int)(*M);
    Nd=(double)(*N);
    Md=(double)(*M);
    Kd=(double)(*K);
    lM=log(*M);
//  This is the number of populations of size N with K classes
//  lcardN=lgammafn(Nd+Kd)-lgammafn(Kd)-lgammafn(Nd+1.);
    Nni=Ni;
    for(k=0;k<Ki;k++){
      Nni-=nk[k];
    }
    lcardN=lgammafn(Nni+Kd)-lgammafn(Kd)-lgammafn(Nni+1.);

    for(k=0;k<Ki;k++){
      lprob[k]=log(prob[k]);
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
    minvalid=(int)trunc(0.05*Mi);
    if(minvalid < 100) minvalid=100;
    while(nvalidhtn < minvalid){
     for(i=0;i<Mi;i++){
      rmultinom(Nni, qprob, Ki, Nk);
// Rprintf("N=%d\n",*N);
//for(k=0;k<Ki;k++){
// Rprintf(" %d",Nk[k]);
//}
// Rprintf("\n");
//   Nkk=0;
//   for(k=0;k<(Ki);k++){
//    Nkk+=(Nk[k]*(k+1));
//   }
      q=dmultinorm(&Nni,K,Nk,lqprob);
//    htnv[i]=bnw_unposN(N,K,n,s,nk,Nk,lprob)-q;
      for(k=0;k<Ki;k++){
       Nk[k]=Nk[k]+nk[k];
      }
      p=dmultinorm(N,K,Nk,lprob);
      htnv[i]=bnw_llikNf(K,n,s,nk,Nk)+p-q;
      htdv[i]=-q;
//  Rprintf("i=%d unposN=%f q=%f p=%f\n",i,htnv[i],q,p);
     }
     mhtni=0.;
     mhtdi=0.;
     nvalidthistime=0;
     for(i=0;i<Mi;i++){
       if(htnv[i] > -90000.){
         mhtni+=htnv[i]/Md;
         nvalidthistime++;
       }
       mhtdi+=htdv[i]/Md;
     }
     mhtni=Md*mhtni/nvalidthistime;
     vhtni=0.;
     vhtdi=0.;
     for(i=0;i<Mi;i++){
      if(htnv[i] > -90000.){
        vhtni+=(htnv[i]-mhtni)*(htnv[i]-mhtni)/nvalidthistime;
      }
      vhtdi+=(htdv[i]-mhtdi)*(htdv[i]-mhtdi)/Md;
     }
//   mhtni=log(mhtni)+htnq1-q1-log(mhtdi)+q1+lcardN;
     mhtn=(nvalidhtn*mhtn+nvalidthistime*mhtni)/(nvalidhtn+nvalidthistime);
     vhtn=(nvalidhtn*vhtn+nvalidthistime*vhtni)/(nvalidhtn+nvalidthistime);
     mhtd=(nvalidhtd*mhtd+Mi*mhtdi)/(nvalidhtd+Mi);
     vhtd=(nvalidhtd*vhtd+Mi*vhtdi)/(nvalidhtd+Mi);
     nvalidhtn+=nvalidthistime;
     nvalidhtd+=Mi;
//   Rprintf("N=%d minvalid=%d nvalidhtn=%d nvalidhtd=%d\n",
//    *N,minvalid,nvalidhtn, nvalidhtd);
    }
//  htn=mhtn+0.5*vhtn+log(1.*nvalidhtn);
    htn=mhtn+0.5*vhtn+log(1.*nvalidhtd);
    htd=mhtd+0.5*vhtd+log(1.*nvalidhtd);
//  Rprintf("mhtn=%f vhtn=%f nvalidhtn=%d\n",mhtn,vhtn,nvalidhtn);
//  Rprintf("mhtd=%f vhtd=%f nvalidhtd=%d\n",mhtd,vhtd,nvalidhtd);
//  htd=mhtd+0.5*vhtd+lM;
//  htn=mhtn+0.5*vhtn;
//  htd=mhtd+0.5*vhtd;
    *unpos=htn-htd+lcardN;
//  Rprintf("N=%d Pct Valid=%f\n",Ni,(nvalidhtn*1.)/(1.*nvalidhtd));
//  Rprintf("htn=%f htd=%f lcardN=%f\n",htn,htd,lcardN);
//  htn=htn-htd;
//  htn=log(htn)-log(htd)+lcardN;
//    cpos=log(cpos/(1.*(*M)));
    PutRNGstate();  /* Disable RNG before returning */
    free(lprob);
    free(lqprob);
    free(htnv);
    free(htdv);
//    Rprintf("K=%d lcardN=%f cpos=%f\n",Ki, lcardN,cpos);
// Rprintf("N=%d #pop ties=%d #mean ties %f\n",(*N),Nkk,(Nkk*1.)/(*N));
}
void bnw_NCbound(int *N, int *K, int *n, int *s, int *nk, int *Nk,
	    double *prob,
	    double *qprob,
	    int *M, double *unpos){
    int i, k, Ki, Ni, Mi, Nni;
    double mhtn;
    double mhtni, vhtni, htn;
    double Nd, Kd, Md, lM;
    double *lprob = (double *) malloc(sizeof(double) * (*K));
    double *htnv = (double *) malloc(sizeof(double) * (*M));

    Ni=(int)(*N);
    Ki=(int)(*K);
    Mi=(int)(*M);
    Nd=(double)(*N);
    Md=(double)(*M);
    Kd=(double)(*K);
    lM=log(*M);

    Nni=Ni;
    for(k=0;k<Ki;k++){
      Nni-=nk[k];
      lprob[k]=log(prob[k]);
    }

    GetRNGstate();  /* R function enabling uniform RNG */

    mhtn=0.;
    for(i=0;i<Mi;i++){
/*    Generate a random draw from the population of sizes */
      rmultinom(Nni, qprob, Ki, Nk);
//    htn=bnw_unposN(N,K,n,s,nk,Nk,lprob);
      for(k=0;k<Ki;k++){
       Nk[k]=Nk[k]+nk[k];
      }
      htn=dmultinorm(N,K,Nk,lprob);
//  Rprintf("i=%d Nni=%d htn=%f\n",i,Nni,htn);
      htn=bnw_llikNf(K,n,s,nk,Nk)+htn;
//    htn=bnw_llikN(K,n,s,nk,Nk);
    if(htn > -90000.){
        htnv[i]=htn;
//      mhtn+=exp(0.0001*htn)/Md;
    }else{
        htnv[i]=-10000.;
    }
    }
     mhtni=0.;
     for(i=0;i<Mi;i++){
       if(htnv[i] > -90000.){
         mhtni+=htnv[i]/Md;
       }
     }
     vhtni=0.;
     for(i=0;i<Mi;i++){
      if(htnv[i] > -90000.){
        vhtni+=(htnv[i]-mhtni)*(htnv[i]-mhtni)/Md;
      }
     }
    htn=mhtni+0.5*vhtni;
    *unpos=htn;
//  Rprintf("N=%d Pct Valid=%f\n",Ni,mhtn);
    PutRNGstate();  /* Disable RNG before returning */
    free(lprob);
    free(htnv);
}
void bnw_stocdiscrete(int *N, int *K, int *n, int *s, int *nk, int *Nk,
	    double *qprob,
	    int *M, double *mllik){
    int i, k, Ki, Ni, Nni, Mi;
    double llik, lbound;
    Ki=(int)(*K);
    Ni=(int)(*N);
    Mi=(int)(*M);
    int *Nmle = (int *) malloc(sizeof(int) * Ki);

    Nni=Ni;
    for(k=0;k<Ki;k++){
     Nni-=nk[k];
    }
    GetRNGstate();  /* R function enabling uniform RNG */

    lbound=-1000000.;
    for(i=0;i<Mi;i++){
/*    Generate a random draw from the population of sizes */
      rmultinom(Nni, qprob, Ki, Nk);
      for(k=0;k<Ki;k++){
       Nk[k]=Nk[k]+nk[k];
      }
      llik=bnw_llikNf(K,n,s,nk,Nk);
//     for(k=0;k<Ki;k++){
//if(Nk[k]<0)  Rprintf("Error: i=%d llik=%f lbound=%f Nk[k]=%d\n",i,llik,lbound,Nk[k]);
//     }
//  Rprintf("i=%d llik=%f lbound=%f\n",i,llik,lbound);
      if(llik > lbound){
       for(k=0;k<Ki;k++){
         Nmle[k]=Nk[k];
       }
       lbound=llik;
//  Rprintf("i=%d lbound=%f\n",i,lbound);
      }
    }
    for(k=0;k<Ki;k++){
      Nk[k]=Nmle[k];
    }
    *mllik=lbound;
    PutRNGstate();  /* Disable RNG before returning */
    free(Nmle);
}
double bnw_unposN(int *N, int *K, int *n, int *s, int *nk, int *Nk,
		  double *lprob){
    int i, k, Nkk, totsize;
    double ll;
    for(k=0;k<(*K);k++){
      if(Nk[k]<nk[k]) return(-1000000.);
    }
    totsize=0;
    ll=lgammafn((*N)+1.);
    for(k=0;k<*K;k++){
	  Nkk=Nk[k];
	  if(Nkk>0){
	    totsize+=(k+1)*Nkk;
	    ll+=Nkk*lprob[k]-lgammafn(Nkk-nk[k]+1.);
//   ll+=lgammafn(Nkk+1.)-lgammafn(Nkk-nk[k]+1.);
//   ll+=Nkk*lprob[k]-lgammafn(Nkk+1.);
	  }
    }
    for(i=0;i<(*n);i++){
	  ll+=log(s[i]/(1.*totsize));
	  totsize-=s[i];
    }
    return(ll);
}
void bnw_stocdiscreteimpute(int *N, int *K, int *n, int *s, int *nk, int *Nk,
	    double *qprob,
	    int *nsim, int *M, double *mllik){
    int i, j, k, l, Ki, Ni, Nni, Mi, nsims, ns;
    double llik, lbound;
    ns=(int)(*n);
    nsims=(int)(*nsim);
    Ki=(int)(*K);
    Ni=(int)(*N);
    Mi=(int)(*M);
    int *si = (int *) malloc(sizeof(int) * ns);
    int *nki = (int *) malloc(sizeof(int) * Ki);
    int *Nmle = (int *) malloc(sizeof(int) * Ki);
    int *Nnis = (int *) malloc(sizeof(int) * nsims);
    double *qprobi = (double *) malloc(sizeof(double) * Ki);

    for(l=0;l<nsims;l++){
     Nni=Ni;
     for(k=0;k<Ki;k++){
      Nni-=nk[k];
     }
     Nnis[l]=Nni;
    }
    GetRNGstate();  /* R function enabling uniform RNG */

    lbound=-1000000.;
    for(i=0;i<Mi;i++){
      llik=0.;
      for(l=0;l<nsims;l++){
       for(j=0;j<ns;j++){
        si[j]=s[l*ns+j];
       }
       for(k=0;k<Ki;k++){
        nki[k]=nk[l*Ki+k];
        qprobi[k]=qprob[l*Ki+k];
       }
/*     Generate a random draw from the population of sizes */
       Nni=Nnis[l];
       rmultinom(Nni, qprobi, Ki, Nk);
       for(k=0;k<Ki;k++){
        Nk[k]=Nk[k]+nki[k];
       }
       llik+=bnw_llikNf(K,n,si,nki,Nk);
      }
      llik/=nsims;
//  Rprintf("i=%d llik=%f lbound=%f\n",i,llik,lbound);
      if(llik > lbound){
       for(k=0;k<Ki;k++){
         Nmle[k]=Nk[k];
       }
       lbound=llik;
      }
    }
    for(k=0;k<Ki;k++){
      Nk[k]=Nmle[k];
    }
    *mllik=lbound;
    PutRNGstate();  /* Disable RNG before returning */
//  Rprintf("i=%d llik=%f lbound=%f\n",i,llik,lbound);
    free(si);
    free(nki);
    free(qprobi);
    free(Nmle);
    free(Nnis);
}
