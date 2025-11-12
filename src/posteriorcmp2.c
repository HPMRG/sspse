/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "posteriorcmp.h"
#include "posteriorcmp2.h"
#include "cmp.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gcmp2 (int *pop12,
            int *pop21,
            int *nk,
            int *K,
            int *n1i,
            int *n2i,
            int *n12i,
            int *samplesize, int *warmup, int *interval,
            double *mu, double *dfmu,
            double *sigma, double *dfsigma,
            double *lnlam, double *nu,
            int *Npi,
            double *lnlamproposal,
            double *nuproposal,
            int *N, int *maxN,
            double *sample,
            double *posu,
            double *lpriorm,
            int *warmuptheta,
            int *verbose
            ) {
  int dimsample, Np;
  int step, staken, getone=1, intervalone=1, verboseMHcmp = 0;
  int i, n, Ni, Ki, isamp, iinterval, isamplesize, iwarmup;
  int n1, n2, n12, unrecap;
  double mui, sigmai, lnlami, nui, dsamp, sigma2i;
  double ddfmu, ddfsigma, dnuproposal;
  int tU1, tU2, sizei, imaxN, imaxm, give_log0=0, give_log1=1;
  double r1, r2, gammart, pis, Nd;
  double temp;
  double errval=0.0000000001, lzcmp;

  GetRNGstate();  /* R function enabling uniform RNG */

  n1=(*n1i);
  n2=(*n2i);
  n12=(*n12i);
  Ni=(*N);
  Ki=(*K);
  imaxN=(*maxN);
  n = n1 + n2 - n12; /* The number unique people seen */
  imaxm=imaxN-n;
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iwarmup=(*warmup);
  Np=(*Npi);
  ddfmu=(*dfmu);
  ddfsigma=(*dfsigma);
  dnuproposal=(*nuproposal);

  dimsample=5+Np;

  double *pi = (double *) malloc(sizeof(double) * Ki);
  int *u1 = (int *) malloc(sizeof(int) * imaxN);
  int *u2 = (int *) malloc(sizeof(int) * imaxN);
  int *b1 = (int *) malloc(sizeof(int) * n1);
  int *b2 = (int *) malloc(sizeof(int) * (n2+1));
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  int *Nkpos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdegi = (double *) malloc(sizeof(double) * (Np+1));
  double *psample = (double *) malloc(sizeof(double) * (Np+1));
  double *lnlamsample = (double *) malloc(sizeof(double) * isamplesize);
  double *nusample = (double *) malloc(sizeof(double) * isamplesize);

  for (i=0; i<Ki; i++){
    nk[i]=0;
  }
  unrecap=0;
  for (i=0; i<n; i++){
    if((pop12[i] >0) && (pop12[i] <= Ki)){ u1[i]=pop12[i];}
    if( pop12[i]==0){ u1[i]=1;}
    if( pop12[i]>Ki){ u1[i]=Ki;}
    nk[u1[i]-1]=nk[u1[i]-1]+1;
    unrecap+=u1[i];
  }
  for (i=0; i<n2; i++){
    if((pop21[i] >0) && (pop21[i] <= Ki)){ u2[i]=pop21[i];}
    if( pop21[i]==0){ u2[i]=1;}
    if( pop21[i]>Ki){ u2[i]=Ki;}
    unrecap-=u2[i];
  }

  // b is the cumulative version of d, which is uobs
  // so b1 is cumulative unit sizes for first list
  // b2 is cumulative unit sizes for second list
  b1[n1-1]=u1[n1-1];
  for (i=(n1-2); i>=0; i--){
    b1[i]=b1[i+1]+u1[i];
  }
  b2[n2]=0;
  if(n2 > 0){
    b2[n2-1]=u2[n2-1];
    for (i=(n2-2); i>=0; i--){
      b2[i]=b2[i+1]+u2[i];
    }
  }
  for (i=0; i<Ki; i++){
     Nk[i]=nk[i];
     Nkpos[i]=0;
     posu[i]=0.;
  }
  for (i=n1; i<imaxN; i++){
    u1[i]=u1[(int)trunc(10*unif_rand()+n1-10)];
  }
  // tU1 is the total unit sizes from the first list
  tU1=0;
  for (i=n1; i<Ni; i++){
    tU1+=u1[i];
  }
  for (i=n2; i<imaxN; i++){
    u2[i]=u2[(int)trunc(10*unif_rand()+n2-10)];
  }
  tU2=unrecap;
  for (i=n2; i<Ni; i++){
    tU2+=u2[i];
  }
  /* Draw initial phis */
  r1=0.;
  for (i=0; i<n1; i++){
    r1+=(exp_rand()/(tU1+b1[i]));
  }
  r2=0.;
  for (i=0; i<n2; i++){
    r2+=(exp_rand()/(tU2+b2[i]));
  }

  for (i=0; i<Np; i++){
     psample[i] = 0.01;
  }
  lnlamsample[0] = (*lnlam);
  nusample[0] = (*nu);

  isamp = 0;
  step = -iwarmup;
  while (isamp < isamplesize) {
    /* Draw new theta */
    /* but less often than the other full conditionals */
    if (step == -iwarmup || step % 10 == 0) {
     MHcmptheta(Nk,K,mu,dfmu,sigma,dfsigma,lnlamproposal,nuproposal,
           &Ni, &Np, psample,
           lnlamsample, nusample, &getone, &staken, warmuptheta, &intervalone,
           &verboseMHcmp);
    }

    for (i=0; i<Np; i++){
      pdegi[i] = psample[i];
    }
    lnlami=lnlamsample[0];
    nui=nusample[0];

    /* Draw phis */
    // tU1 is the total unit sizes from first list
    tU1=0;
    for (i=n1; i<Ni; i++){
      tU1+=u1[i];
    }
    tU2=unrecap;
    for (i=n2; i<Ni; i++){
      tU2+=u2[i];
    }
    r1=0.;
    for (i=0; i<n1; i++){
      r1+=(exp_rand()/(tU1+b1[i]));
    }
    r2=0.;
    for (i=0; i<n2; i++){
      r2+=(exp_rand()/(tU2+b2[i]));
    }

    /* Draw new N */

    /* First find the degree distribution */
    pis=0.;
    lzcmp = zcmp(exp(lnlami), nui, errval, Ki, give_log1);
    if(lzcmp < -100000.0){continue;}
    pi[Np]=cmp(Np+1,lnlami,nui,lzcmp,give_log0);
    pis+=pi[Np];
    for (i=Np+1; i<Ki; i++){
      pi[i]=pi[i-1]*exp(lnlami-nui*log((double)(i+1)));
      pis+=pi[i];
    }
//  pis=1.-exp(-lzcmp);
    for (i=0; i<Ki; i++){
      pi[i]/=pis;
    }
    pis=1.;
    for (i=0; i<Np; i++){
      pis-=pdegi[i];
    }
    for (i=0; i<Ki; i++){
      pi[i]*=pis;
    }
    // !!!!! Why this? For non-parametric piece
    for (i=0; i<Np; i++){
      pi[i]=pdegi[i];
    }

    // Now computes mean and s.d. from log-lambda and nu
    mui=0.0;
    sigma2i=0.0;
    for (i=0; i<Ki; i++){
      mui+=pi[i]*(i+1);
      sigma2i+=pi[i]*(i+1)*(i+1);
    }
    sigma2i=sigma2i-mui*mui;
    sigmai = sqrt(sigma2i);


    gammart=0.;
    for (i=0; i<Ki; i++){
      gammart+=(exp(-(r1 + r2)*(i+1))*pi[i]);
    }
    gammart=log(gammart);
    temp = -100000000.0;
    // N = m + n
    // Compute (log) P(m | \theta and data and \Psi)
    for (i=0; i<imaxm; i++){
      lpm[i]=i*gammart+lgamma(n+i+1.)-lgamma(i+1.);
//    Add in the (log) prior on m: P(m)
      lpm[i]+=lpriorm[i];
      if(lpm[i] > temp) temp = lpm[i];
    }
    for (i=0; i<imaxm; i++){
      lpm[i]=exp(lpm[i]-temp);
    }
    for (i=1; i<imaxm; i++){
      lpm[i]=lpm[i-1]+lpm[i];
    }
    temp = lpm[imaxm-1] * unif_rand();
    for (Ni=0; Ni<imaxm; Ni++){
      if(temp <= lpm[Ni]) break;
    }
    // Add back the sample size
    Ni += n;
    if(Ni > imaxN) Ni = imaxN;


    /* Draw unseen sizes */
    for (i=0; i<Ki; i++){
      Nk[i]=nk[i];
    }
    // Set up pi to be cumulative for random draws
    for (i=1; i<Ki; i++){
      pi[i]=pi[i-1]+pi[i];
    }
    for (i=n; i<Ni; i++){
      /* Propose unseen size for unit i */
      /* Use rejection sampling */
      sizei=1000000;
      while(sizei >= Ki){
       sizei=1000000;
       while(log(1.0-unif_rand()) > -(r1+r2)*sizei){
        /* Now propose unseen size for unit i */
        /* In the next two lines a sizei is chosen */
        /* with parameters lnlami and nui */
        temp = unif_rand();
        for (sizei=1; sizei<=Ki; sizei++){
          if(temp <= pi[sizei-1]) break;
        }
       }
      }
      u1[i]=sizei;
      u2[i]=sizei;
      Nk[sizei-1]=Nk[sizei-1]+1;
    }
    if (step > 0 && step % iinterval == 0) {
      /* record statistics for posterity */
      Nd=(double)Ni;
      sample[isamp*dimsample  ]=Nd;
      sample[isamp*dimsample+1]=mui;
      sample[isamp*dimsample+2]=sigmai;
      sample[isamp*dimsample+3]=(double)(Nk[0]);
      temp=0.0;
      for (i=0; i<Ki; i++){
        temp+=(i+1.0)*Nk[i];
      }
      sample[isamp*dimsample+4]=temp;
      for (i=0; i<Np; i++){
        sample[isamp*dimsample+5+i]=pdegi[i];
      }
      for (i=0; i<Ki; i++){
        Nkpos[i]=Nkpos[i]+Nk[i];
        posu[i]+=((Nk[i]*1.)/Nd);
      }
      isamp++;
      if (*verbose && isamplesize % isamp == 0) Rprintf("Taken %d samples...\n", isamp);
    }
    step++;
  }
  dsamp=((double)isamp);
  for (i=0; i<Ki; i++){
    nk[i]=Nkpos[i];
    posu[i]/=dsamp;
  }
  for (i=0; i<n1; i++){
     pop12[i]=u1[i];
  }
  for (i=0; i<n2; i++){
     pop21[i]=u2[i];
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(pi);
  free(u1);
  free(u2);
  free(b1);
  free(b2);
  free(Nk);
  free(Nkpos);
  free(lpm);
  free(pdegi);
  free(psample);
  free(lnlamsample);
  free(nusample);
}
