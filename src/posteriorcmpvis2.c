/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "posteriorcmpvis.h"
#include "posteriorcmpvis2.h"
#include "cmp.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gcmpvis2 (int *pop12, int *pop21,
            int *nk, 
            int *K, 
            int *n1, 
            int *n2, 
            int *n0, 
            int *samplesize, int *burnin, int *interval,
            double *mu, double *dfmu, 
            double *sigma, double *dfsigma,
            double *lnlam, double *nu,
            double *beta0muprior, double *beta0sigmaprior, 
            double *beta1muprior, double *beta1sigmaprior, 
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            int *Npi,
            int *srd, 
            int *numrec, 
            double *rectime,
            int *srd2, 
            int *numrec2, 
            double *rectime2,
            int *rc, 
            int *maxcoupons,
            double *muproposal, 
            double *nuproposal, 
            double *beta0proposal, double *beta1proposal, 
            double *lmemmuproposal, double *memnuproposal,
            int *N, int *maxN, 
            double *sample, 
            int *vsample, 
            int *vsample2, 
            double *posu, 
            double *lpriorm, 
            int *burnintheta,
            int *burninbeta,
            int *verbose
                         ) {
  int dimsample, Np;
  int ni0, ni1, ni2, itemp;
  int step, staken, getone=1, intervalone=1, verboseMHcmp = 0;
  int i, j, k, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  int umax, unrecap;
  double mui, sigmai, dsamp, nui, lnlami, sigma2i;
  double dbeta0, dbeta1;
  double dlmemmu, dmemnu;
  double beta0i, beta1i, lmemmui, memnui;
  int tU1, tU2, sizei, imaxN, imaxm, give_log0=0, give_log1=1;
  double r1, r2, gammart, pis, Nd;
  double temp, temp2, uprob1, uprob2;
  double rtprob, lliki;
  int maxc;
  double errval=0.0000000001, lzcmp;

  GetRNGstate();  /* R function enabling uniform RNG */

  ni1=(*n1);
  ni2=(*n2);
  ni0=(*n0);
  ni = ni1 + ni2 - ni0; /* The number unique people seen */
  Ni=(*N);
  Ki=(*K);
  imaxN=(*maxN);
  imaxm=imaxN-ni;
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburnin=(*burnin);
  Np=(*Npi);
  dbeta0=(*beta0muprior);
  dbeta1=(*beta1muprior);
  dlmemmu=(*lmemmu);
  dmemnu=(*memnu);
  maxc=(*maxcoupons);

  dimsample=5+Np+4;

  double *pi = (double *) malloc(sizeof(double) * Ki);
  double *pd = (double *) malloc(sizeof(double) * Ki);
  double *pdm = (double *) malloc(sizeof(double) * Ki);
  double *pd2 = (double *) malloc(sizeof(double) * Ki);
  int *d1 = (int *) malloc(sizeof(int) * imaxN);
  int *d2 = (int *) malloc(sizeof(int) * imaxN);
  int *b1 = (int *) malloc(sizeof(int) * ni1);
  int *b2 = (int *) malloc(sizeof(int) * (ni2+1));
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  int *Nkpos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdegi = (double *) malloc(sizeof(double) * (Np+1));
  double *psample = (double *) malloc(sizeof(double) * (Np+1));
  double *lnlamsample = (double *) malloc(sizeof(double));
  double *nusample = (double *) malloc(sizeof(double));
  double *beta0sample = (double *) malloc(sizeof(double));
  double *beta1sample = (double *) malloc(sizeof(double));
  double *lmemmusample = (double *) malloc(sizeof(double));
  double *memnusample = (double *) malloc(sizeof(double));

  for (i=0; i<Ki; i++){
    nk[i]=0;
  }
  unrecap=0;
  uprob1=ni;
  for (i=0; i<ni; i++){
    if((pop12[i] >0) && (pop12[i] <= Ki)){ d1[i]=pop12[i];}
    if( pop12[i]==0){ d1[i]=1;}
    if( pop12[i]>Ki){ d1[i]=Ki; uprob1--;}
    nk[d1[i]-1]=nk[d1[i]-1]+1;
    unrecap+=d1[i];
  }
  uprob1/=ni;
  uprob1 = 0.5 + uprob1/2.0;
  uprob2=ni2+1;
  for (i=0; i<ni2; i++){
    if((pop21[i] >0) && (pop21[i] <= Ki)){ d2[i]=pop21[i];}
    if( pop21[i]==0){ d2[i]=1;}
    if( pop21[i]>Ki){ d2[i]=Ki; uprob2--;}
    unrecap-=d2[i];
  }

  uprob2/=ni2+1;
  uprob2 = 0.5 + uprob2/2.0;
  b1[ni1-1]=d1[ni1-1];
  for (i=(ni1-2); i>=0; i--){
    b1[i]=b1[i+1]+d1[i];
  }
  b2[ni2]=0;
  if(ni2 > 0){
    b2[ni2-1]=d2[ni2-1];
    for (i=(ni2-2); i>=0; i--){
      b2[i]=b2[i+1]+d2[i];
    }
  }
  for (i=0; i<Ki; i++){
     Nk[i]=nk[i];
     Nkpos[i]=0;
     posu[i]=0.;
  }
  for (i=ni1; i<imaxN; i++){
    d1[i]=d1[(int)trunc(10*unif_rand()+ni1-10)];
  }
  tU1=0;
  for (i=ni1; i<Ni; i++){
    tU1+=d1[i];
  }
  for (i=ni2; i<imaxN; i++){
    d2[i]=d2[(int)trunc(10*unif_rand()+ni2-10)];
  }
  tU2=unrecap;
  for (i=ni2; i<Ni; i++){
    tU2+=d2[i];
  }
// Rprintf("d1[ni1] %d d2[ni2] %d tU1 %d tU2 %d\n", d1[ni1],d2[ni2],tU1,tU2);
  /* Draw initial phis */
  r1=0.;
  for (i=0; i<ni1; i++){
    r1+=(exp_rand()/(tU1+b1[i]));
  }
  r2=0.;
  for (i=0; i<ni2; i++){
    r2+=(exp_rand()/(tU2+b2[i]));
  }

  for (i=0; i<Np; i++){
     psample[i] = 0.01;
  }
  beta0sample[0] = dbeta0;
  beta1sample[0] = dbeta1;
  lmemmusample[0] = dlmemmu;
  memnusample[0] = dmemnu;
  lnlamsample[0] = (*lnlam);
  nusample[0] = (*nu);

  isamp = 0;
  step = -iburnin;
  while (isamp < isamplesize) {

    /* Draw new theta */
    /* but less often than the other full conditionals */
    if (step == -iburnin || step==(10*(step/10))) { 
     MHcmptheta(Nk,K,mu,dfmu,sigma,dfsigma,muproposal,nuproposal,
       &Ni, &Np, psample,
       lnlamsample, nusample, &getone, &staken, burnintheta, &intervalone, 
       &verboseMHcmp);

     for (i=0; i<Np; i++){
      pdegi[i] = psample[i];
     }
     lnlami=lnlamsample[0];
     nui=nusample[0];
//if(nui > 4.0 || lnlami > 4.5) Rprintf("lnlami %f nui %f dfmu %f\n", lnlami, nui, (*dfmu));
    }
//Rprintf("Finished MHcmpthetavis : isamp %d\n", isamp);

    /* Compute the unit distribution (given the new theta = (mu, sigma)) */
    pis=0.;
    lzcmp = zcmp(exp(lnlami), nui, errval, Ki, give_log1);
    if(lzcmp < -100000.0){continue;}
    pi[Np]=cmp(Np+1,lnlami,nui,lzcmp,give_log0);
//Rprintf("lnlami %f nui %f lzcmp %f pi %f\n", lnlami, nui, lzcmp, pi[Np]);
    for (i=Np+1; i<Ki; i++){
      pi[i]=pi[i-1]*exp(lnlami-nui*log((double)(i+1)));
    }
//  Rprintf("isamp %d pis %f\n", isamp, pis);
    pis=1.-exp(-lzcmp);
    for (i=0; i<Ki; i++){
      pi[i]/=pis;
//Rprintf("i %d pi %f pi0 %f\n", i, pi[i], pi0[i], pis, pis0);
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
    temp=0.;
    for (umax=1; umax<=Ki; umax++){
      temp+=pi[umax-1];
      if(temp > uprob1) break;
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

//  if (step == -iburnin || step==(10*(step/10))) { 
//    Rprintf("mean unit size = %f s.d. = %f\n", mui, sigmai);
//  }
//  Rprintf("burninbeta: %d\n", (*burninbeta));
    /* Draw new beta using MCMC */
    if (step == -iburnin || step==(10*(step/10))) { 
     MHcmpbeta2(d1,d2,n1,n2,K,beta0muprior,beta0sigmaprior,beta1muprior,beta1sigmaprior,
       lmemmu,memdfmu,memnu,memdfnu,srd,numrec,rectime,srd2,numrec2,rectime2,maxcoupons,
       beta0proposal,beta1proposal,
       lmemmuproposal,memnuproposal,
       beta0sample, beta1sample,lmemmusample,memnusample,
       &getone, &staken, burninbeta, &intervalone, 
       &verboseMHcmp);
     beta0i=beta0sample[0];
     beta1i=beta1sample[0];
     lmemmui=lmemmusample[0];
     memnui=memnusample[0];
//Rprintf("Finished MHcmpbeta2 : isamp %d\n", isamp);
//    Rprintf("mean lmemmui = %f\n", temp2);
    }

    /* Draw true degrees (sizes) of the first RDS sample based on the reported degrees*/
    // First reset counts
    for (i=0; i<Ki; i++){
     nk[i]=0;
    }
    for (j=0; j<ni1; j++){
      if(srd[j] <= Ki){
      temp = beta0i + beta1i*rectime[j];
      rtprob = exp(temp)/(1.0+exp(temp));
//    Multiply by the PoissonLogNormal PMF for observation
      for (i=0; i<Ki; i++){
//     if((numrec[j] <= (i+1)) && ((maxc-1) <= (i+1))){
       if(numrec[j] <= (i+1)){
        lliki=0.0;
        if(((i+1) <= maxc)|(numrec[j]<maxc)){
          lliki += dbinom(numrec[j],(i+1),rtprob,give_log1);
        }else{
          lliki += log(1.0-pbinom(maxc-1.0,(i+1),rtprob,give_log0,give_log0));
        }
//      if(srd[j] < 0) Rprintf("srd: %d\n", srd[j]);
        if(srd[j]>=0){
//        Use CMP localized
//        temp= lmemmui + lmemrci*rc[j] + log((double)(i+1));
          temp= lmemmui + log((double)(i+1));
          pis=0.;
          lzcmp = zcmp(exp(temp), memnui, errval, Ki, give_log1);
          if(lzcmp < -100000.0){Rprintf("badlzcmp ");continue;}
          pd2[0]=cmp(1,temp,memnui,lzcmp,give_log0);
          for (k=1; k<Ki; k++){
            pd2[k]=pd2[k-1]*exp(temp-memnui*log((double)(k+1)));
          }
          pis=1.-exp(-lzcmp);
          for (k=0; k<Ki; k++){
            pd2[k]/=pis;
          }
          if(srd[i] <= Ki){
            lliki += log(pd2[srd[j]-1]);
          }
        }
        pd[i]=pi[i]*exp(lliki);
        pdm[i]=exp(lliki);
       }else{
        pd[i]=0.0;
        pdm[i]=0.0;
       }
      }
      // Set up pd to be cumulative for the random draws
      temp2 = pdm[0];
      temp = pdm[0];
      for (i=1; i<Ki; i++){
       pd[i]=pd[i-1]+pd[i];
       temp2+=pdm[i]*(i+1);
       temp+=pdm[i];
      }
      temp2/=temp;
      if(pd[Ki-1]<0.00000000001){
       for (i=0; i<Ki; i++){
        pd[i]=pi[i];
       }
       for (i=1; i<Ki; i++){
        pd[i]=pd[i-1]+pd[i];
       }
      }
      /* Draw unit size for the observed degree */
      /* Now propose the true size for unit i based on reported size and disease status */
      /* In the next three lines a sizei is chosen */
      temp = pd[Ki-1] * unif_rand();
      for (sizei=1; sizei<=Ki; sizei++){
        if(temp <= pd[sizei-1]) break;
      }
      }else{
      /* Deal with the outliers */
        sizei=umax;
      }

      nk[sizei-1]=nk[sizei-1]+1;
      d1[j]=sizei;
     } 
     // Rebuild b1
     b1[ni1-1]=d1[ni1-1];
     for (i=(ni1-2); i>=0; i--){
      b1[i]=b1[i+1]+d1[i];
     }
     /* End of imputed unit sizes for observed in the first RDS*/

//Rprintf("Finished d1 : isamp %d\n", isamp);
    /* Draw true degrees (sizes) of the second RDS sample based on the reported degrees*/
    // First reset counts
//  for (i=0; i<Ki; i++){
//   nk2[i]=0;
//  }
    itemp = 0;
    for (j=0; j<ni2; j++){
      if(srd2[j] <= Ki){
      temp = beta0i + beta1i*rectime2[j];
      rtprob = exp(temp)/(1.0+exp(temp));
//    Multiply by the PoissonLogNormal PMF for observation
      for (i=0; i<Ki; i++){
//     if((numrec2[j] <= (i+1)) && ((maxc-1) <= (i+1))){
       if(numrec2[j] <= (i+1)){
        lliki=0.0;
        if(((i+1) <= maxc)|(numrec2[j]<maxc)){
          lliki += dbinom(numrec2[j],(i+1),rtprob,give_log1);
        }else{
          lliki += log(1.0-pbinom(maxc-1.0,(i+1),rtprob,give_log0,give_log0));
        }
        if(srd2[j]>=0){
//        Use CMP localized
          temp= lmemmui + log((double)(i+1));
          pis=0.;
          lzcmp = zcmp(exp(temp), memnui, errval, Ki, give_log1);
          if(lzcmp < -100000.0){Rprintf("badlzcmp ");continue;}
          pd2[0]=cmp(1,temp,memnui,lzcmp,give_log0);
          for (k=1; k<Ki; k++){
            pd2[k]=pd2[k-1]*exp(temp-memnui*log((double)(k+1)));
          }
          pis=1.-exp(-lzcmp);
          for (k=0; k<Ki; k++){
            pd2[k]/=pis;
          }
          if(srd2[i] <= Ki){
            lliki += log(pd2[srd2[j]-1]);
          }
        }
        pd[i]=pi[i]*exp(lliki);
       }else{
        pd[i]=0.0;
       }
//Rprintf("srd[j] %d i %d pd[i] %f\n", srd[j],i, pd[i]);
//Rprintf("srd[j] %d numrec2 %f i %d lliki %f\n", srd[j], numrec2[j], i+1, lliki);
      }
      // Set up pd to be cumulative for the random draws
//if(j==6) Rprintf("pd[i]: ");
      for (i=1; i<Ki; i++){
//if(j==6) Rprintf("%f ", pd[i]);
       pd[i]=pd[i-1]+pd[i];
      }
//if(j==6) Rprintf("\n");
//if(j==6)Rprintf("beta0i %f beta1i %f lmemmui %f memnui %f rtprob %f pd[Ki-1] %f\n", beta0i,beta1i,lmemmui,memnui,rtprob,pd[Ki-1]);
      if(pd[Ki-1]<0.00000000001){
// Rprintf("fixed bad pd[Ki-1] %f\n", pd[Ki-1]);
       for (i=0; i<Ki; i++){
        pd[i]=pi[i];
       }
       for (i=1; i<Ki; i++){
        pd[i]=pd[i-1]+pd[i];
       }
      }
      /* Draw unit size for the observed degree */
      /* Now propose the true size for unit i based on reported size and disease status */
      /* In the next three lines a sizei is chosen */
      temp = pd[Ki-1] * unif_rand();
      for (sizei=1; sizei<=Ki; sizei++){
        if(temp <= pd[sizei-1]) break;
      }
      }else{
      /* Deal with the outliers */
        sizei=umax;
      }

      d2[j]=sizei;
      if(rc[j]==0){
        nk[sizei-1]=nk[sizei-1]+1;
        d1[ni1 + itemp]=sizei;
        itemp++;
      }
//Rprintf("j %d dis %d sizei %d pd[Ki-1] %f\n", j, ddis, sizei, pd[Ki-1]);
     } 
     itemp=0;
     for (i=0; i<ni1; i++){
       itemp+=d1[i];
     }
     for (i=0; i<ni2; i++){
       if(rc[i]!=0) itemp-=d2[i];
     }
     unrecap=itemp;
     // Rebuild b2
     b2[ni2]=0;
     if(ni2 > 0){
       b2[ni2-1]=d2[ni2-1];
       for (i=(ni2-2); i>=0; i--){
         b2[i]=b2[i+1]+d2[i];
       }
     }
     /* End of imputed unit sizes for observed in second RDS sample */

     /* End of overall loop started at "Draw new beta using MCMC" */

    /* Draw phis */
    tU1=0;
    for (i=ni1; i<Ni; i++){
      tU1+=d1[i];
    }
    tU2=unrecap;
    for (i=ni2; i<Ni; i++){
      tU2+=d2[i];
    }
    r1=0.;
    for (i=0; i<ni1; i++){
      r1+=exp_rand()/(tU1+b1[i]);
    }
    r2=0.;
    for (i=0; i<ni2; i++){
      r2+=exp_rand()/(tU2+b2[i]);
    }

    /* Draw new N */

    gammart=0.;
    for (i=0; i<Ki; i++){
      gammart+=(exp(-(r1+r2)*(i+1))*pi[i]);
    }
    gammart=log(gammart);
    temp = -100000000.0;
    // N = m + n
    // Compute (log) P(m | \theta and data and \Psi)
    for (i=0; i<imaxm; i++){
      lpm[i]=i*gammart+lgamma(ni+i+1.)-lgamma(i+1.);
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
    Ni += ni;
    if(Ni > imaxN) Ni = imaxN;
//Rprintf("Finished Ni : isamp %d\n", isamp);


    /* Draw unseen sizes */
    for (i=0; i<Ki; i++){
//Rprintf("i %d nk[i] %f\n", i, nk[i]/(1.0*ni));
      Nk[i]=nk[i];
    }
    // Set up pi to be cumulative for random draws
    for (i=1; i<Ki; i++){
//    Rprintf("i %d pi[i] %f\n", i, pi[i]);
      pi[i]=pi[i-1]+pi[i];
    }
    for (i=ni; i<Ni; i++){
      /* Propose unseen size for unit i */
      /* Use rejection sampling */
      sizei=1000000;
      while(sizei >= Ki){
       sizei=1000000;
       while(log(1.0-unif_rand()) > -(r1+r2)*sizei){
        /* Now propose unseen size for unit i */
        /* In the next two lines a sizei is chosen */
        /* with parameters mui and nui */
        temp = unif_rand();
//      gammart = pi[Ki-1] * unif_rand();
        for (sizei=1; sizei<=Ki; sizei++){
          if(temp <= pi[sizei-1]) break;
        }
//      Rprintf("sizei %d pi[Ki-1] %f gammart %f\n", sizei, pi[Ki-1],gammart);
       }
      }
//    if(sizei >= Ki){sizei=Ki-1;}
      d1[i]=sizei;
      d2[i]=sizei;
//    if((sizei <= 0) | (sizei > Ki-1)) Rprintf("sizei %d r %f\n", sizei,r);
      Nk[sizei-1]=Nk[sizei-1]+1;
    }
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
      Nd=(double)Ni;
      sample[isamp*dimsample  ]=Nd;
//if(sigmai > 4.0 || mui > 4.5) Rprintf("sample: %f %f\n", mui,sigmai);
// Rprintf("sample: %f %f\n", mui,sigmai);
      sample[isamp*dimsample+1]=mui;
      sample[isamp*dimsample+2]=sigmai;
      sample[isamp*dimsample+3]=(double)(Nk[0]);
      temp=0.0;
      for (i=0; i<Ki; i++){
        temp+=(i+1.0)*Nk[i];
      }
      sample[isamp*dimsample+4]=temp;
      sample[isamp*dimsample+5]=beta0i;
      sample[isamp*dimsample+6]=beta1i;
      sample[isamp*dimsample+7]=lmemmui;
      sample[isamp*dimsample+8]=memnui;
      for (i=0; i<Np; i++){
        sample[isamp*dimsample+9+i]=pdegi[i];
      }
      for (i=0; i<Ki; i++){
        Nkpos[i]=Nkpos[i]+Nk[i];
        posu[i]+=((Nk[i]*1.)/Nd);
      }
      for (i=0; i<ni1; i++){
         vsample[isamp*ni1+i]=d1[i];
      }
      for (i=0; i<ni2; i++){
        vsample2[isamp*ni2+i]=d2[i];
      }
      isamp++;
      if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
//    if (*verbose) Rprintf("Taken %d samples...\n", isamp);
 Rprintf("Taken %d samples...\n", isamp);
// Rprintf("d1[ni1-1] %d d2[ni2-1] %d d1[ni1+10] %d tU1 %d tU2 %d\n", d1[ni1-1],d2[ni2-1],d1[ni1+10],tU1,tU2);
    }
    step++;
  }
  dsamp=((double)isamp);
  for (i=0; i<Ki; i++){
    nk[i]=Nkpos[i];
    posu[i]/=dsamp;
  }
  for (i=0; i<ni1; i++){
     pop12[i]=d1[i];
  }
  for (i=0; i<ni2; i++){
     pop21[i]=d2[i];
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(pi);
  free(pd);
  free(pd2);
  free(d1);
  free(d2);
  free(psample);
  free(pdegi);
  free(b1);
  free(b2);
  free(Nk);
  free(Nkpos);
  free(lpm);
  free(lnlamsample);
  free(nusample);
  free(beta0sample);
  free(beta1sample);
  free(lmemmusample);
  free(memnusample);
}

void MHcmpbeta2 (int *d1, int *d2, int *n1, int *n2, int *K,
            double *beta0, double *beta0s, double *beta1, double *beta1s, 
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            int *srd, 
            int *numrec, 
            double *rectime,
            int *srd2, 
            int *numrec2, 
            double *rectime2,
            int *maxcoupons,
            double *beta0proposal, double *beta1proposal, 
            double *lmemmuproposal, double *memnuproposal, 
            double *beta0sample, double *beta1sample,
            double *lmemmusample, double *memnusample,
            int *samplesize, int *staken, int *burninbeta, int *interval,
            int *verbose
         ) {
  int Ki, maxc, ni1, ni2;
  int step, taken, give_log1=1, give_log0=0;
  int i, k, isamp, iinterval, isamplesize, iburninbeta;
  double ip, cutoff;
  double temp, rtprob;
  double beta0star, beta1star, beta0i, beta1i;
  double qi, qstar, lliki, llikstar;
  double lmemmustar, memnustar, lmemmui, memnui;
  double rmemnui, rmemnustar;
  double pibeta0star, pibeta0i;
  double pibeta1star, pibeta1i;
  double pimemstar, pimemi;
  double dbeta0, dbeta0s, dbeta1, dbeta1s;
  double dlmemmu, dmemdfmu, dmemdfnu, rmemdfmu;
  double dmemnu, dmemnur;
  double dbeta0proposal, dbeta1proposal;
  double dlmemmuproposal, dmemnuproposal;
  double pis, errval=0.0000000001, lzcmp;

//Rprintf("burninbeta: %d\n", (*burninbeta));

  Ki=(*K);
  double *pd = (double *) malloc(sizeof(double) * Ki);

//GetRNGstate();  /* R function enabling uniform RNG */

  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburninbeta=(*burninbeta);
  dbeta0=(*beta0);
  dbeta0s=(*beta0s);
  dbeta1=(*beta1);
  dbeta1s=(*beta1s);
  dlmemmu=(*lmemmu);
  dmemnu=(*memnu);
  dmemnur=sqrt(dmemnu);
  dmemdfmu=(*memdfmu);
  rmemdfmu=sqrt(dmemdfmu);
  dmemdfnu=(*memdfnu);
  dbeta0proposal=(*beta0proposal);
  dbeta1proposal=(*beta1proposal);
  dlmemmuproposal=(*lmemmuproposal);
  dmemnuproposal=(*memnuproposal);

  // First set starting values
  isamp = taken = 0;
  step = -iburninbeta;
  ni1 =(*n1);
  ni2 =(*n2);
  maxc=(*maxcoupons);
  beta0i = beta0sample[0];
  beta1i = beta1sample[0];
  lmemmui = lmemmusample[0];
  memnui = memnusample[0];

  // Compute initial current lik
  lliki = 0.0;
  for (i=0; i<ni1; i++){
    temp = beta0i + beta1i*rectime[i];
    rtprob = exp(temp)/(1.0+exp(temp));
//  if((numrec[i] <= d1[i]) && ((maxc-1) <= d1[i])){
    if(numrec[i] <= d1[i]){
     if((d1[i] <= maxc)|(numrec[i]<maxc)){
      lliki += dbinom(numrec[i],d1[i],rtprob,give_log1);
     }else{
      lliki += log(1.0-pbinom(maxc-1.0,d1[i],rtprob,give_log0,give_log0));
     }
     if(srd[i]>=0){
//    Use CMP localized
      temp= lmemmui + log((double)(d1[i]));
      pis=0.;
      lzcmp = zcmp(exp(temp), memnui, errval, Ki, give_log1);
      if(lzcmp < -100000.0){Rprintf("badlzcmp ");continue;}
      pd[0]=cmp(1,temp,memnui,lzcmp,give_log0);
      for (k=1; k<Ki; k++){
        pd[k]=pd[k-1]*exp(temp-memnui*log((double)(k+1)));
      }
      pis=1.-exp(-lzcmp);
      for (k=0; k<Ki; k++){
        pd[k]/=pis;
//  Rprintf(" %f",pd[k]);
      }
//  Rprintf("\n %d %f\n",i,lliki);
      if(srd[i] <= Ki){
        lliki += log(pd[srd[i]-1]);
      }
     }
    }else{
     lliki = -100000.0; 
    }
  }
  for (i=0; i<ni2; i++){
    temp = beta0i + beta1i*rectime2[i];
    rtprob = exp(temp)/(1.0+exp(temp));
//  if((numrec2[i] <= d2[i]) && ((maxc-1) <= d2[i])){
    if(numrec2[i] <= d2[i]){
     if((d2[i] <= maxc)|(numrec2[i]<maxc)){
      lliki += dbinom(numrec2[i],d2[i],rtprob,give_log1);
     }else{
      lliki += log(1.0-pbinom(maxc-1.0,d2[i],rtprob,give_log0,give_log0));
     }
     if(srd2[i]>=0){
//    Use CMP localized
      temp= lmemmui + log((double)(d2[i]));
      pis=0.;
      lzcmp = zcmp(exp(temp), memnui, errval, Ki, give_log1);
      if(lzcmp < -100000.0){Rprintf("badlzcmp ");continue;}
      pd[0]=cmp(1,temp,memnui,lzcmp,give_log0);
      for (k=1; k<Ki; k++){
        pd[k]=pd[k-1]*exp(temp-memnui*log((double)(k+1)));
      }
      pis=1.-exp(-lzcmp);
      for (k=0; k<Ki; k++){
        pd[k]/=pis;
      }
      if(srd2[i] <= Ki){
        lliki += log(pd[srd2[i]-1]);
      }
     }
    }else{
     lliki = -100000.0; 
    }
  }
  if(!isfinite(lliki)) lliki = -100000.0; 

//  Rprintf("%d of %d MH samples taken in %i steps %i; cutoff=%f\n",
//   isamp, isamplesize, step, iburninbeta, cutoff);
//Rprintf("numrec =%d rectime =%f\n",numrec[3],rectime[3]);
//Rprintf("numrec2=%d rectime2=%f\n",numrec2[3],rectime2[3]);
//Rprintf("New call: lliki=%f lmemmui=%f memnui=%f rtprob=%f\n",lliki,lmemmui,memnui,rtprob);

  // Compute initial prior
  pibeta0i = dnorm(beta0i, dbeta0, dbeta0s, give_log1);
  pibeta1i = dnorm(beta1i, dbeta1, dbeta1s, give_log1);
  rmemnui  = sqrt(memnui);
  pimemi = dnorm(lmemmui, dlmemmu, rmemnui/rmemdfmu, give_log1);
  pimemi = pimemi+dsclinvchisq(memnui, dmemdfnu, dmemnu);

  qi = dnorm(log(memnui/memnui)/dmemnuproposal,0.,1.,give_log1)
       -log(dmemnuproposal*memnui);

  // Now do the MCMC updates (starting with the burnin updates)
  while (isamp < isamplesize && step < 1000) {
    /* Propose new beta */
    beta0star = rnorm(beta0i, dbeta0proposal);
    beta1star = rnorm(beta1i, dbeta1proposal);
    /* Propose new memnu and lmemmu */
    lmemmustar = rnorm(lmemmui, dlmemmuproposal);
//  memnustar = rnorm(memnui, dmemnuproposal);
    memnustar = memnui*exp(rnorm(0., dmemnuproposal));
    rmemnustar = sqrt(memnustar);

// for (i=0; i<ni1; i++){
//    Rprintf("%f %i %f\n",numrec[i],d1[i],rectime[i]);
// }
  
    llikstar = 0.0;
    for (i=0; i<ni1; i++){
      temp = beta0star + beta1star*rectime[i];
      rtprob = exp(temp)/(1.0+exp(temp));
//    if((numrec[i] <= d1[i]) && ((maxc-1) <= d1[i])){
      if(numrec[i] <= d1[i]){
       if((d1[i] <= maxc)|(numrec[i]<maxc)){
        llikstar += dbinom(numrec[i],d1[i],rtprob,give_log1);
       }else{
        llikstar += log(1.0-pbinom(maxc-1.0,d1[i],rtprob,give_log0,give_log0));
       }
       if(srd[i]>=0){
//      Use CMP localized
        temp= lmemmustar + log((double)(d1[i]));
        pis=0.;
        lzcmp = zcmp(exp(temp), memnustar, errval, Ki, give_log1);
        if(lzcmp < -100000.0){Rprintf("badlzcmp ");continue;}
        pd[0]=cmp(1,temp,memnustar,lzcmp,give_log0);
        for (k=1; k<Ki; k++){
          pd[k]=pd[k-1]*exp(temp-memnustar*log((double)(k+1)));
        }
        pis=1.-exp(-lzcmp);
        for (k=0; k<Ki; k++){
          pd[k]/=pis;
//  Rprintf(" %f",pd[k]);
        }
//  Rprintf("\n %d srd %d pd %f llikstar %f temp %f\n",i,srd[i],pd[srd[i]-1],llikstar, temp);
        if(srd[i] <= Ki){
          llikstar += log(pd[srd[i]-1]);
        }
//  Rprintf("\n %d srd %d pd %f %f\n",i,srd[i],pd[srd[i]-1],llikstar);
//  Rprintf("\n %d %f\n",i,llikstar);
       }
      }else{
       llikstar = -100000.0; 
      }
    }
    for (i=0; i<ni2; i++){
      temp = beta0star + beta1star*rectime2[i];
      rtprob = exp(temp)/(1.0+exp(temp));
//    if((numrec2[i] <= d2[i]) && ((maxc-1) <= d2[i])){
      if(numrec2[i] <= d2[i]){
       if((d2[i] <= maxc)|(numrec2[i]<maxc)){
        llikstar += dbinom(numrec2[i],d2[i],rtprob,give_log1);
       }else{
        llikstar += log(1.0-pbinom(maxc-1.0,d2[i],rtprob,give_log0,give_log0));
       }
       if(srd2[i]>=0){
//      Use CMP localized
        temp= lmemmustar + log((double)(d2[i]));
        pis=0.;
        lzcmp = zcmp(exp(temp), memnustar, errval, Ki, give_log1);
        if(lzcmp < -100000.0){Rprintf("badlzcmp ");continue;}
        pd[0]=cmp(1,temp,memnustar,lzcmp,give_log0);
        for (k=1; k<Ki; k++){
          pd[k]=pd[k-1]*exp(temp-memnustar*log((double)(k+1)));
        }
        pis=1.-exp(-lzcmp);
        for (k=0; k<Ki; k++){
          pd[k]/=pis;
//  Rprintf(" %f",pd[k]);
        }
        if(srd2[i] <= Ki){
          llikstar += log(pd[srd2[i]-1]);
        }
//  Rprintf("\n %d srd %d pd %f %f\n",i,srd2[i],pd[srd2[i]-1],llikstar);
       }
      }else{
       llikstar = -100000.0; 
      }
    }
    if(!isfinite(llikstar)) llikstar = -100000.0; 

    /* Calculate pieces of the prior. */
    pibeta0star = dnorm(beta0star, dbeta0, dbeta0s, give_log1);
    pibeta1star = dnorm(beta1star, dbeta1, dbeta1s, give_log1);
    rmemnustar  = sqrt(memnustar);
    pimemstar = dnorm(lmemmustar, dlmemmu, rmemnustar/rmemdfmu, give_log1);
    pimemstar = pimemstar+dsclinvchisq(memnustar, dmemdfnu, dmemnu);

    qstar = dnorm(log(memnustar/memnui)/dmemnuproposal,0.,1.,give_log1)
                  -log(dmemnuproposal*memnustar);

    /* Calculate ratio */
    ip =      pibeta0star-pibeta0i;
    ip = ip + pibeta1star-pibeta1i;
    ip = ip + pimemstar - pimemi;

// Rprintf("pibeta0star=%f pibeta0i=%f pibeta1star=%f pibeta1i=%f pimemstar=%f pimemi=%f\n",pibeta0star,pibeta0i,pibeta1star,pibeta1i,pimemstar,pimemi);

    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
    then let the MH probability equal min{exp(cutoff), 1.0}.
    But we'll do it in log space instead.  */
//  if (*verbose)
//  Rprintf("Now proposing %d MH steps %f ip1...\n", step, ip);
//  cutoff = ip + lliki-llikstar;
    cutoff = ip + llikstar - lliki + qi - qstar;
      
//  Rprintf("Proposed: cutoff=%f ip=%f llikstar=%f lliki= %f qi=%f qstar=%f\n", cutoff, ip, llikstar, lliki, qi, qstar);
//  Rprintf("Proposed: beta0i=%f beta0star=%f beta1s=%f beta1star=%f\n", beta0i,beta0star,beta1i,beta1star);

//  if (*verbose)
//      Rprintf("Now proposing %d MH steps %f cutoff...\n", step, cutoff);

    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
// Rprintf("Accepted: cutoff=%f ip=%f llikstar=%f lliki=%f qi=%f qstar=%f\n", cutoff, ip, llikstar, lliki, qi, qstar);
      /* Make proposed changes */
      beta0i = beta0star;
      beta1i = beta1star;
      lmemmui = lmemmustar;
      memnui = memnustar;
      lliki = llikstar;
      qi = qstar;
      pibeta0i = pibeta0star;
      pibeta1i = pibeta1star;
      pimemi = pimemstar;
      taken++;
      if (step > 0 && step==(iinterval*(step/iinterval))) { 
        /* record statistics for posterity */
        beta0sample[isamp]=beta0i;
        beta1sample[isamp]=beta1i;
        lmemmusample[isamp]=lmemmui;
        memnusample[isamp]=memnui;
        isamp++;
        if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d MH samples...\n", isamp);
      }
    }
//  Rprintf("%d of %d MH samples taken in %d steps %d; cutoff=%f\n",
//   isamp, isamplesize, step, iburninbeta, cutoff);
    step++;
  }
//PutRNGstate();  /* Disable RNG before returning */
  free(pd);
  /*Check for interrupts (if recursion is taking way too long...)*/
  R_CheckUserInterrupt();
//Rprintf("Done: %d MH samples taken with %d steps\n", taken, step);
  *staken = taken;
}
