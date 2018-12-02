/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "posteriorcmpvis.h"
#include "cmp.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gcmpvis (int *pop,
            int *nk,
            int *K,
            int *n,
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
            int *maxcoupons,
            double *lnlamproposal,
            double *nuproposal,
            double *beta0proposal, double *beta1proposal,
            double *lmemmuproposal, double *memnuproposal,
            int *N, int *maxN,
            double *sample,
            int *vsample,
            double *ppos,
            double *lpriorm,
            int *burnintheta,
            int *burninbeta,
            int *verbose
                         ) {
  int dimsample, Np;
  int step, staken, getone=1, intervalone=1, verboseMHcmp = 0;
  int i, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  int j, k;
  double mui, sigmai, lnlami, nui, dsamp;
  double sigma2i;
  double dbeta0, dbeta1;
  double dlmemmu, dmemnu;
  double beta0i, beta1i, lmemmui, memnui;
  int tU, sizei, imaxN, imaxm, give_log0=0, give_log1=1;
  double r, gammart, pis, Nd;
  double temp;
  double rtprob, lliki;
  int maxc;
  double errval=0.0000000001, lzcmp;

  GetRNGstate();  /* R function enabling uniform RNG */

  ni=(*n);
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
  double *pd2 = (double *) malloc(sizeof(double) * Ki);
  int *d = (int *) malloc(sizeof(int) * imaxN);
  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  int *Nkpos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdegi = (double *) malloc(sizeof(double) * (Np+1));
  double *psample = (double *) malloc(sizeof(double) * (Np+1));
  double *lnlamsample = (double *) malloc(sizeof(double));
  double *beta0sample = (double *) malloc(sizeof(double));
  double *beta1sample = (double *) malloc(sizeof(double));
  double *nusample = (double *) malloc(sizeof(double));
  double *lmemmusample = (double *) malloc(sizeof(double));
  double *memnusample = (double *) malloc(sizeof(double));

  for (i=0; i<ni; i++){
    if((pop[i]>0) && (pop[i] <= Ki)){ d[i]=pop[i];}
    if( pop[i]==0){ d[i]=1;}
    if( pop[i]>Ki){ d[i]=Ki;}
  }
  b[ni-1]=d[ni-1];
  for (i=(ni-2); i>=0; i--){
    b[i]=b[i+1]+d[i];
  }
  for (i=0; i<Ki; i++){
     Nk[i]=nk[i];
     Nkpos[i]=0;
     ppos[i]=0.;
  }
  tU=0;
  for (i=ni; i<Ni; i++){
    d[i]=d[(int)trunc(10.0*unif_rand()+ni-10)];
    tU+=d[i];
  }
  /* Draw initial phis */
  r=0.;
  for (i=0; i<ni; i++){
    r+=(exp_rand()/(tU+b[i]));
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
     MHcmptheta(Nk,K,mu,dfmu,sigma,dfsigma,lnlamproposal,nuproposal,
       &Ni, &Np, psample,
       lnlamsample, nusample, &getone, &staken, burnintheta, &intervalone, 
       &verboseMHcmp);

     lnlami=lnlamsample[0];
     nui=nusample[0];

     for (i=0; i<Np; i++){
      pdegi[i] = psample[i];
     }
    }

    /* Compute the unit distribution (given the new theta = (lnlam, nu)) */
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

    // Now computes mean and s.d. from log-lambda and nu
    mui=0.0;
    sigma2i=0.0;
    for (i=0; i<Ki; i++){
      mui+=pi[i]*(i+1);
      sigma2i+=pi[i]*(i+1)*(i+1);
    }
    sigma2i=sigma2i-mui*mui;
    sigmai = sqrt(sigma2i);

    /* Draw new beta using MCMC */
    if (step == -iburnin || step==(20*(step/20))) { 
//  if (step == -iburnin) { 
     MHcmpmem(d,n,K,beta0muprior,beta0sigmaprior,beta1muprior,beta1sigmaprior,
       lmemmu,memdfmu,memnu,memdfnu,srd,numrec,rectime,maxcoupons,
       beta0proposal,beta1proposal,
       lmemmuproposal,memnuproposal,
       beta0sample, beta1sample,lmemmusample,memnusample,
       &getone, &staken, burninbeta, &intervalone, 
       &verboseMHcmp);
     beta0i=beta0sample[0];
     beta1i=beta1sample[0];
     lmemmui=lmemmusample[0];
     memnui=memnusample[0];
//   beta0i= -3.063687833  ;
//   beta1i=  0.001548059 ;
//   lmemmui= 1.380479339 ;
//   memnui= 0.913915196 ;
//  }
//  Rprintf("step %d lmemmui: %f memnui %f beta0i %f beta1i %f\n",step,lmemmui,memnui,beta0i,beta1i);

    /* Draw true degrees (sizes) based on the reported degrees*/
    // First reset counts
    for (i=0; i<Ki; i++){
     nk[i]=0;
    }
    for (j=0; j<ni; j++){
      temp = beta0i + beta1i*rectime[j];
      rtprob = exp(temp)/(1.0+exp(temp));
//    Multiply by the PoissonLogNormal PMF for observation
      for (i=0; i<Ki; i++){
       // Next to exclude unit sizes inconsistent with the number of recruits
//     if((numrec[j] <= (i+1)) & ((maxc-1) <= (i+1))){
       if((numrec[j] <= (i+1))){
        lliki=0.0;
        if(((i+1) <= maxc)|(numrec[j]<maxc)){
          lliki += dbinom(numrec[j],(i+1),rtprob,give_log1);
        }else{
          lliki += log(1.0-pbinom(maxc-1.0,(i+1),rtprob,give_log0,give_log0));
        }
        if(srd[j]>=0){
//        lliki += log(poilog(srd[j],log((double)(i+1))-lmemmui,memnui));
//        Use CMP localized
          temp=lmemmui + log((double)(i+1));
          pis=0.;
          lzcmp = zcmp(exp(temp), memnui, errval, Ki, give_log1);
          if(lzcmp < -100000.0){Rprintf("badlzcmp ");continue;}
          pd2[0]=cmp(1,temp,memnui,lzcmp,give_log0);
          for (k=1; k<Ki; k++){
            pd2[k]=pd2[k-1]*exp(temp-memnui*log((double)(k+1)));
          }
          pis=1.-exp(-lzcmp);
//if(j==6) Rprintf("pd2[i]: %f ",pis);
          for (k=0; k<Ki; k++){
            pd2[k]/=pis;
//if(j==6) Rprintf("%f ", pd2[i]);
          }
//if(j==6) Rprintf("\n");
//if(j==6) Rprintf("pd2[srd]: %f\n",pd2[srd[j]-1]);
          lliki += log(pd2[srd[j]-1]);
        }
        pd[i]=pi[i]*exp(lliki);
       }else{
        pd[i]=0.0;
       }
//     pd[i]=pi[i];
//Rprintf("srd[j] %d i %d pd[i] %f\n", srd[j],i, pd[i]);
//Rprintf("srd[j] %d numrec %f i %d lliki %f\n", srd[j], numrec[j], i+1, lliki);
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
//     Rprintf("fixed bad pd[Ki-1] %f\n", pd[Ki-1]);
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

      nk[sizei-1]=nk[sizei-1]+1;
      d[j]=sizei;

//    if(d[j] < numrec[j]) Rprintf("Warning: j %d d[j] %d numrec[j] %d maxc %d: %f %f %f %f\n",j,d[j],numrec[j],maxc,pd[0],pd[1],pd[2],pd[3]);
//Rprintf("j %d dis %d sizei %d pd[Ki-1] %f\n", j, ddis, sizei, pd[Ki-1]);
     } 
     // Rebuild b
     b[ni-1]=d[ni-1];
     for (i=(ni-2); i>=0; i--){
      b[i]=b[i+1]+d[i];
     }
// Rprintf("j %d d[j] %d pd[Ki-1] %f\n", j, d[j], pd[Ki-1]);
//
     /* End of imputed unit sizes for observed */
    }

//  }
// Rprintf("step %d b[ni-2] %d; End of imputed unit sizes for observed\n", step, b[ni-2]);

    /* Draw phis */
    tU=0;
    for (i=ni; i<Ni; i++){
      tU+=d[i];
    }
    r=0.;
    for (i=0; i<ni; i++){
//    phi[i]=(tU+b[i])*exp_rand();
      r+=exp_rand()/(tU+b[i]);
    }

    /* Draw new N */

    gammart=0.;
    for (i=0; i<Ki; i++){
      gammart+=(exp(-r*(i+1))*pi[i]);
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
       while(log(1.0-unif_rand()) > -r*sizei){
        /* Now propose unseen size for unit i */
        /* In the next two lines a sizei is chosen */
        /* with parameters lnlami and nui */
        temp = unif_rand();
//      gammart = pi[Ki-1] * unif_rand();
        for (sizei=1; sizei<=Ki; sizei++){
          if(temp <= pi[sizei-1]) break;
        }
//      Rprintf("sizei %d pi[Ki-1] %f gammart %f\n", sizei, pi[Ki-1],gammart);
       }
      }
//    if(sizei >= Ki){sizei=Ki-1;}
      d[i]=sizei;
//    if((sizei <= 0) | (sizei > Ki-1)) Rprintf("sizei %d r %f\n", sizei,r);
      Nk[sizei-1]=Nk[sizei-1]+1;
    }
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
      Nd=(double)Ni;
      sample[isamp*dimsample  ]=Nd;
//if(nui > 4.0 || lnlami > 4.5) Rprintf("sample: %f %f\n", lnlami,nui);
// Rprintf("sample: step %d %f %f %d %d %d\n", step, lnlami,nui,Np,isamp*dimsample+8,isamp*ni+ni-1);
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
        ppos[i]+=((Nk[i]*1.)/Nd);
      }
      for (i=0; i<ni; i++){
        vsample[isamp*ni+i]=d[i];
      }
      isamp++;
      if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
//    if (*verbose) Rprintf("Taken %d samples...\n", isamp);
    }
    step++;
  }
  dsamp=((double)isamp);
  for (i=0; i<Ki; i++){
    nk[i]=Nkpos[i];
    ppos[i]/=dsamp;
  }
  for (i=0; i<ni; i++){
     pop[i]=d[i];
  }
 // Rprintf("ni %d Ki %d\n", ni, Ki);
  PutRNGstate();  /* Disable RNG before returning */
  free(pi);
  free(pd);
  free(pd2);
  free(d);
  free(b);
  free(Nk);
  free(Nkpos);
  free(lpm);
  free(pdegi);
  free(psample);
  free(lnlamsample);
  free(nusample);
  free(beta0sample);
  free(beta1sample);
  free(lmemmusample);
  free(memnusample);
}

void MHcmpmem (int *d, int *n, int *K,
            double *beta0, double *beta0sd, double *beta1, double *beta1sd, 
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            int *srd, 
            int *numrec, 
            double *rectime,
            int *maxcoupons,
            double *beta0proposal, double *beta1proposal, 
            double *lmemmuproposal, double *memnuproposal, 
            double *beta0sample, double *beta1sample,
            double *lmemmusample, double *memnusample,
            int *samplesize, int *staken, int *burnin, int *interval,
            int *verbose
         ) {
  int Ki, maxc, ni;
  int step, taken, give_log1=1, give_log0=0;
  int i, k, isamp, iinterval, isamplesize, iburnin;
  double ip, cutoff;
  double temp, rtprob;
  double beta0star, beta1star, beta0i, beta1i;
  double qi, qstar, lliki, llikstar;
  double lmemmustar, memnustar, lmemmui, memnui;
  double memnu2i, memnu2star;
  double pibeta0star, pibeta0i;
  double pibeta1star=0.0, pibeta1i=0.0;
  double pimemstar, pimemi;
  double dbeta0, dbeta0sd, dbeta1, dbeta1sd;
  double dlmemmu, dmemdfmu, dmemdfnu, rmemdfmu;
  double dmemnu, dmemnu2;
  double dbeta0proposal, dbeta1proposal;
  double dlmemmuproposal, dmemnuproposal;
  double pis, errval=0.0000000001, lzcmp;

  Ki=(*K);
  double *pd = (double *) malloc(sizeof(double) * Ki);

//GetRNGstate();  /* R function enabling uniform RNG */

  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburnin=(*burnin);
  dbeta0=(*beta0);
  dbeta0sd=(*beta0sd);
  dbeta1=(*beta1);
  dbeta1sd=(*beta1sd);
  dlmemmu=(*lmemmu);
  dmemnu=(*memnu);
  dmemnu2=dmemnu*dmemnu;
  dmemdfmu=(*memdfmu);
  rmemdfmu=sqrt(dmemdfmu);
  dmemdfnu=(*memdfnu);
  dbeta0proposal=(*beta0proposal);
  dbeta1proposal=(*beta1proposal);
  dlmemmuproposal=(*lmemmuproposal);
  dmemnuproposal=(*memnuproposal);

  // First set starting values
  isamp = taken = 0;
  step = -iburnin;
  ni =(*n);
  maxc=(*maxcoupons);
  beta0i = beta0sample[0];
  beta1i = beta1sample[0];
  lmemmui = lmemmusample[0];
  memnui = memnusample[0];
  memnu2i = memnui*memnui;

  // Compute initial current lik
  lliki = 0.0;
  for (i=0; i<ni; i++){
    temp = beta0i + beta1i*rectime[i];
    rtprob = exp(temp)/(1.0+exp(temp));
//  if((numrec[i] <= d[i]) && ((maxc-1) <= d[i])){
    if(numrec[i] <= d[i]){
     if((d[i] <= maxc)|(numrec[i]<maxc)){
       lliki += dbinom(numrec[i],d[i],rtprob,give_log1);
     }else{
       lliki += log(1.0-pbinom(maxc-1.0,d[i],rtprob,give_log0,give_log0));
     }
     if(srd[i]>=0){
//    lliki += log(poilog(srd[i],log((double)(d[i]))-lmemmui,memnui));
//    Use CMP localized
      temp=lmemmui+log(d[i]);
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
      lliki += log(pd[srd[i]-1]);
     }
    }else{
     lliki = -100000.0; 
    }
  }
  if(!isfinite(lliki)) lliki = -100000.0; 

//Rprintf("New call: lliki=%f lmemmui=%f memnui=%f rtprob=%f\n",lliki,lmemmui,memnui,rtprob);

  // Compute initial prior
  pibeta0i = dnorm(beta0i, dbeta0, dbeta0sd, give_log1);
  if(dbeta1sd > 0.0) pibeta1i = dnorm(beta1i, dbeta1, dbeta1sd, give_log1);
  pimemi = dnorm(lmemmui, dlmemmu, memnui/rmemdfmu, give_log1);
  pimemi = pimemi+dsclinvchisq(memnu2i, dmemdfnu, dmemnu2);

//Rprintf("dmemdfmu %f dmemdfnu %f ddfmu %f ddfsigma %f\n", memdfmu, memdfnu, dfmu, dfsigma);
//Rprintf("dmemdfmu %f dmemdfnu %f rmemdfmu %f\n", dmemdfmu, dmemdfnu, rmemdfmu);

  // Now do the MCMC updates (starting with the burnin updates)
  while (isamp < isamplesize && step < iburnin) {
    /* Propose new beta */
    beta0star = rnorm(beta0i, dbeta0proposal);
    if(dbeta1sd > 0.0){
      beta1star = rnorm(beta1i, dbeta1proposal);
    }else{
      beta1star = beta1i;
    }
    /* Propose new memnu and lmemmu */
    lmemmustar = rnorm(lmemmui, dlmemmuproposal);
//  memnustar = rnorm(memnui, dmemnuproposal);
    memnu2star = memnu2i*exp(rnorm(0., dmemnuproposal));
    memnustar = sqrt(memnu2star);

// for (i=0; i<ni; i++){
//    Rprintf("%f %i %f\n",numrec[i],d[i],rectime[i]);
// }
  
    llikstar = 0.0;
    for (i=0; i<ni; i++){
      temp = beta0star + beta1star*rectime[i];
      rtprob = exp(temp)/(1.0+exp(temp));
//    Rprintf("%f %f %f %f\n",llikstar,lmemmustar,memnustar,rtprob);
//    if((numrec[i] <= d[i]) && ((maxc-1) <= d[i])){
      if(numrec[i] <= d[i]){
       if((d[i] <= maxc)|(numrec[i]<maxc)){
        llikstar += dbinom(numrec[i],d[i],rtprob,give_log1);
//      Rprintf("< %f %f %f %f\n",numrec[i],rtprob[i],u,dbinom(numrec[i],u,rtprob[i],give_log0));
       }else{
        llikstar += log(1.0-pbinom(maxc-1.0,d[i],rtprob,give_log0,give_log0));
//      Rprintf("< %f\n",1.0-pbinom(maxc-1.0,u,rtprob[i],give_log0,give_log0));
       }
       if(srd[i]>=0){
//      llikstar += log(poilog(srd[i],log(d[i])-lmemmustar,exp(memnustar)));
//      Use CMP localized
        temp=lmemmustar+log(d[i]);
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
        }
        llikstar += log(pd[srd[i]-1]);
//  Rprintf("llikstar: %i %f %f %f %f %f %f %f\n",i, llikstar,lmemmustar,memnustar,beta0star,beta1star,temp,rtprob);
       }
      }else{
       llikstar = -100000.0; 
      }
    }
    if(!isfinite(llikstar)) llikstar = -100000.0; 

//  Rprintf("llikstar: %f %f %f %f\n",llikstar,lmemmustar,memnustar,rtprob);

//    Rprintf("dn %f\n", poilog(srd[0],log(d[0])-lmemmustar,exp(memnustar)));

//      Rprintf("dn %f\n", log(temp));
//  if(!isnan(temp)) llikstar += temp; 

    /* Calculate pieces of the prior. */
    pibeta0star = dnorm(beta0star, dbeta0, dbeta0sd, give_log1);
    if(dbeta1sd > 0.0) pibeta1star = dnorm(beta1star, dbeta1, dbeta1sd, give_log1);
    pimemstar = dnorm(lmemmustar, dlmemmu, memnustar/rmemdfmu, give_log1);
    pimemstar = pimemstar+dsclinvchisq(memnu2star, dmemdfnu, dmemnu2);

    qi = dnorm(log(memnui/memnustar)/dmemnuproposal,0.,1.,give_log1)
         -log(dmemnuproposal*memnui);

    qstar = dnorm(log(memnustar/memnui)/dmemnuproposal,0.,1.,give_log1)
                  -log(dmemnuproposal*memnustar);

    /* Calculate ratio */
    ip =      pibeta0star-pibeta0i;
    if(dbeta1sd > 0.0) ip = ip + pibeta1star-pibeta1i;
    ip = ip + pimemstar - pimemi;

//  Rprintf("pibeta0star=%f pibeta0i=%f %f %f %f %f\n",pibeta0star,pibeta0i,pibeta1star,pibeta1i,pimemstar,pimemi);

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
//    Rprintf("Now proposing %d MH steps %f cutoff...\n", step, cutoff);

    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
//  Rprintf("Accepted: cutoff=%f ip=%f llikstar=%f lliki=%f qi=%f qstar=%f\n",cutoff, ip,llikstar, lliki, qi, qstar);
      /* Make proposed changes */
      beta0i = beta0star;
      beta1i = beta1star;
      lmemmui = lmemmustar;
      memnui = memnustar;
      memnu2i = memnu2star;
      lliki = llikstar;
      qi = qstar;
      pibeta0i = pibeta0star;
      if(dbeta1sd > 0.0) pibeta1i = pibeta1star;
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
//  Rprintf("%d of %d MH samples taken in %d steps; cutoff=%f\n", isamp, samplesize, step, cutoff);
    step++;
  }
//PutRNGstate();  /* Disable RNG before returning */
  free(pd);
  /*Check for interrupts (if recursion is taking way too long...)*/
  R_CheckUserInterrupt();
//Rprintf("Done: %d MH samples taken with %d steps\n", taken, step);
  *staken = taken;
}
