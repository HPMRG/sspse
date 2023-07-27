/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "posteriorcmpwpvis_binomial.h"
#include "cmp.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gcmpwpvis_binomial (int *pop,
            int *K,
            int *n,
            int *samplesize, int *burnin, int *interval,
            double *mu, double *dfmu,
            double *sigma, double *dfsigma,
            double *lnlam, double *nu,
            double *beta0muprior, double *beta0sigmaprior,
            double *betatmuprior, double *betatsigmaprior,
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            double *memod,
            int *Npi,
            int *srd,
            int *numrec,
            double *rectime,
            int *maxcoupons,
            double *lnlamproposal,
            double *nuproposal,
            double *beta0proposal, double *betatproposal,
            double *lmemmuproposal, double *memnuproposal,
            int *N, int *maxN,
            double *sample,
            int *vsample,
            double *posu,
            double *posd,
            double *lpriorm,
            int *burnintheta,
            int *burninbeta,
            int *verbose
                         ) {
  int dimsample, Np;
  int step, staken, getone=1, intervalone=1, verboseMHcmp = 0;
  int i, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  int j, k;
  int umax;
  double alpha, pnb, rnb;
  double mui, sigmai, lnlami, nui, dsamp;
  double sigma2i;
  double dbeta0, dbetat;
  double dlmemmu, dmemnu;
  double beta0i, betati, lmemmui, memnui;
  double memmui;
  int tU, sizei, imaxN, imaxm, give_log0=0, give_log1=1;
  double r, gammart, pis, pis2, Nd;
  double temp, uprob;
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
  dbetat=(*betatmuprior);
  dlmemmu=(*lmemmu);
  dmemnu=(*memnu);
  alpha=(*memod);
  maxc=(*maxcoupons);
  
  dimsample=5+Np+4;
  pnb=(alpha-1.)/alpha;

  double *pi = (double *) malloc(sizeof(double) * Ki);
  double *pd = (double *) malloc(sizeof(double) * Ki);
  double *pd2 = (double *) malloc(sizeof(double) * (10*Ki));
  int *u = (int *) malloc(sizeof(int) * imaxN);
  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  int *nk = (int *) malloc(sizeof(int) * Ki);
  int *Nkpos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdegi = (double *) malloc(sizeof(double) * (Np+1));
  double *psample = (double *) malloc(sizeof(double) * (Np+1));
  double *lnlamsample = (double *) malloc(sizeof(double));
  double *nusample = (double *) malloc(sizeof(double));
  double *beta0sample = (double *) malloc(sizeof(double));
  double *betatsample = (double *) malloc(sizeof(double));
  double *lmemmusample = (double *) malloc(sizeof(double));
  double *memnusample = (double *) malloc(sizeof(double));

  for (i=0; i<Ki; i++){
    nk[i]=0;
  }
  uprob=ni;
  for (i=0; i<ni; i++){
    if((pop[i]>0) && (pop[i] <= Ki)){ u[i]=pop[i];}
    if( pop[i]==0){ u[i]=1;}
    if( pop[i]>Ki){ u[i]=Ki; uprob--;}
    nk[u[i]-1]=nk[u[i]-1]+1;
  }
  uprob/=ni;
  uprob = 0.5 + uprob/2.0;
// Rprintf("uprob %f\n", uprob);
  b[ni-1]=u[ni-1];
  for (i=(ni-2); i>=0; i--){
    b[i]=b[i+1]+u[i];
  }
  for (i=0; i<Ki; i++){
     Nk[i]=nk[i];
     Nkpos[i]=0;
     posu[i]=0.;
     posd[i]=0.;
  }
  for (i=ni; i<imaxN; i++){
    u[i]=u[(int)trunc(10*unif_rand()+ni-10)];
  }
  tU=0;
  for (i=ni; i<Ni; i++){
    tU+=u[i];
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
  betatsample[0] = dbetat;
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
    pis+=pi[Np];
    for (i=Np+1; i<Ki; i++){
      pi[i]=pi[i-1]*exp(lnlami-nui*log((double)(i+1)));
      pis+=pi[i];
    }
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

    /* Draw new beta using MCMC */
    if (step == -iburnin || step==(20*(step/20))) { 
     MHwpmem_binomial(u,n,K,beta0muprior,beta0sigmaprior,betatmuprior,betatsigmaprior,
       lmemmu,memdfmu,memnu,memdfnu,srd,numrec,rectime,maxcoupons,
       beta0proposal,betatproposal,
       lmemmuproposal,memnuproposal,
       beta0sample, betatsample,lmemmusample,memnusample,
       &getone, &staken, burninbeta, &intervalone, 
       &verboseMHcmp);
     beta0i=beta0sample[0];
     betati=betatsample[0];
     lmemmui=lmemmusample[0];
     memnui=memnusample[0];

    /* Draw true unit sizes based on the reported degrees*/
    // First reset counts
    for (i=0; i<Ki; i++){
     nk[i]=0;
    }
    umax = 0;
    for (j=0; j<ni; j++){
     if(srd[j] <= (10*Ki)){
      temp = beta0i + betati*rectime[j];
      rtprob = exp(temp)/(1.0+exp(temp));
//    Multiply by the Conway-Maxwell-Poisson PMF for observation
      for (i=0; i<Ki; i++){
       // Next to exclude unit sizes inconsistent with the number of recruits
       if((numrec[j] <= (i+1))){
        lliki=0.0;
        if(((i+1) <= maxc)|(numrec[j]<maxc)){
          lliki += dbinom(numrec[j],(i+1),rtprob,give_log1);
        }else{
          lliki += log(1.0-pbinom(maxc-1.0,(i+1),rtprob,give_log0,give_log0));
        }
        if(srd[j]>=0){
//        Use WP for localized
          memmui = exp(lmemmui)*(i+1.);
          rnb=memmui/(alpha-1.);
          pd2[0]= exp(-fabs(memmui-1.)/sqrt(memnui));
          pis=pd2[0];
          for (k=1; k<(10*Ki); k++){
            pd2[k]= pd2[k-1]*(k+rnb)*pnb*exp((fabs(k-memmui)-fabs(k+1-memmui))/sqrt(memnui))/((double)(k+1));
            pis+=pd2[k];
          }
          for (k=0; k<(10*Ki); k++){
            pd2[k]/=pis;
          }
          pis2=0.;
          for (k=Ki; k<(10*Ki); k++){
            pis2+=pd2[k];
          }
          if(srd[j] <= 10*Ki){
            lliki += log(pd2[srd[j]-1]);
          }else{
            lliki += log(pis2);
          }
        }
        pd[i]=pi[i]*exp(lliki)*(i+1.);
       }else{
        pd[i]=0.0;
       }
      }
    pis=0.;
    for (i=0; i<Ki; i++){
      pis+=pd[i];
    }
    for (i=0; i<Ki; i++){
      pd[i]/=pis;
    }
      for (i=1; i<Ki; i++){
       pd[i]=pd[i-1]+pd[i];
      }
      /* Draw unit size for the observed degree */
      /* Now propose the true size for unit i based on reported size and disease status */
      /* In the next three lines a sizei is chosen */
      temp = pd[Ki-1] * unif_rand();
      for (sizei=1; sizei<=Ki; sizei++){
        if(temp <= pd[sizei-1]) break;
      }
      nk[sizei-1]=nk[sizei-1]+1;
      u[j]=sizei;

    }
     } 

     /* Deal with the outliers */
     umax=nk[Ki-1];
     sizei=Ki;
     for (j=0; j<ni; j++){
      if(srd[j] > 10*Ki){
       while((umax==0) & (sizei > 1)){
         sizei--;
         umax=nk[sizei-1];
       }
       u[j]=sizei;
       umax--;
      }
     }
     for (j=0; j<ni; j++){
      if(srd[j] > 10*Ki){
       nk[u[j]-1]=nk[u[j]-1]+1;
      }
     }

     // Rebuild b
     b[ni-1]=u[ni-1];
     for (i=(ni-2); i>=0; i--){
      b[i]=b[i+1]+u[i];
     }
     /* End of imputed unit sizes for observed */
    }

    /* Draw phis */
    tU=0;
    for (i=ni; i<Ni; i++){
      tU+=u[i];
    }
    r=0.;
    for (i=0; i<ni; i++){
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
      Nk[i]=nk[i];
    }
    // Set up pi to be cumulative for random draws
    for (i=1; i<Ki; i++){
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
        for (sizei=1; sizei<=Ki; sizei++){
          if(temp <= pi[sizei-1]) break;
        }
       }
      }
      u[i]=sizei;
      Nk[sizei-1]=Nk[sizei-1]+1;
    }
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
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
      sample[isamp*dimsample+5]=beta0i;
      sample[isamp*dimsample+6]=betati;
      sample[isamp*dimsample+7]=lmemmui;
      sample[isamp*dimsample+8]=memnui;
      for (i=0; i<Np; i++){
        sample[isamp*dimsample+9+i]=pdegi[i];
      }
      for (i=0; i<Ki; i++){
        Nkpos[i]=Nkpos[i]+Nk[i];
        posu[i]+=((Nk[i]*1.)/Nd);
      }
      for (i=0; i<ni; i++){
        vsample[isamp*ni+i]=u[i];
        if((srd[i]>0) && (srd[i] <= 10*Ki)){ posd[srd[i]-1]+=(1./Nd); }
      }
//    Record the predicted degrees (not unit sizes)
      for (i=ni; i<Ni; i++){
        memmui = exp(lmemmui)*u[i];
        rnb=memmui/(alpha-1.);
        pd[0]= exp(-fabs(memmui-1.)/sqrt(memnui));
        pis=pd[0];
        for (k=1; k<Ki; k++){
          pd[k]= pd[k-1]*(k+rnb)*pnb*exp((fabs(k-memmui)-fabs(k+1-memmui))/sqrt(memnui))/((double)(k+1));
          pis+=pd[k];
        }
        for (k=0; k<Ki; k++){
          pd[k]/=pis;
          posd[k]+=(pd[k]/Nd);
        }
      }
      isamp++;
      if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
    }
    step++;
  }
  dsamp=((double)isamp);
  for (i=0; i<Ki; i++){
    nk[i]=Nkpos[i];
    posu[i]/=dsamp;
    posd[i]/=dsamp;
  }
  for (i=0; i<ni; i++){
     pop[i]=u[i];
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(pi);
  free(pd);
  free(pd2);
  free(u);
  free(b);
  free(nk);
  free(Nk);
  free(Nkpos);
  free(lpm);
  free(pdegi);
  free(psample);
  free(lnlamsample);
  free(nusample);
  free(beta0sample);
  free(betatsample);
  free(lmemmusample);
  free(memnusample);
}

void MHwpmem_binomial (int *u, int *n, int *K,
            double *beta0, double *beta0sd, double *betat, double *betatsd, 
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            int *srd, 
            int *numrec, 
            double *rectime,
            int *maxcoupons,
            double *beta0proposal, double *betatproposal, 
            double *lmemmuproposal, double *memnuproposal, 
            double *beta0sample, double *betatsample,
            double *lmemmusample, double *memnusample,
            int *samplesize, int *staken, int *burnin, int *interval,
            int *verbose
         ) {
  int Ki, maxc, ni;
  int step, taken, give_log1=1, give_log0=0;
  int i, k, isamp, iinterval, isamplesize, iburnin;
  double ip, cutoff;
  double temp, rtprob;
  double beta0star, betatstar, beta0i, betati;
  double qi, qstar, lliki, llikstar;
  double lmemmustar, memnustar, lmemmui, memnui;
  double memmui, memmustar;
  double rmemnui, rmemnustar;
  double pibeta0star, pibeta0i;
  double pibetatstar=0.0, pibetati=0.0;
  double pimemstar, pimemi;
  double dbeta0, dbeta0sd, dbetat, dbetatsd;
  double dlmemmu, dmemdfmu, dmemdfnu, rmemdfmu;
  double dmemnu, dmemnur;
  double dbeta0proposal, dbetatproposal;
  double dlmemmuproposal, dmemnuproposal;
  double pis, pis2;
  double alpha=25., pnb, rnb;

  Ki=(*K);
  double *pd = (double *) malloc(sizeof(double) * (10*Ki));
  pnb=(alpha-1.)/alpha;

//GetRNGstate();  /* R function enabling uniform RNG */

  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburnin=(*burnin);
  dbeta0=(*beta0);
  dbeta0sd=(*beta0sd);
  dbetat=(*betat);
  dbetatsd=(*betatsd);
  dlmemmu=(*lmemmu);
  dmemnu=(*memnu);
  dmemnur=sqrt(dmemnu);
  dmemdfmu=(*memdfmu);
  rmemdfmu=sqrt(dmemdfmu);
  dmemdfnu=(*memdfnu);
  dbeta0proposal=(*beta0proposal);
  dbetatproposal=(*betatproposal);
  dlmemmuproposal=(*lmemmuproposal);
  dmemnuproposal=(*memnuproposal);

  // First set starting values
  isamp = taken = 0;
  step = -iburnin;
  ni =(*n);
  maxc=(*maxcoupons);
  beta0i = beta0sample[0];
  betati = betatsample[0];
  lmemmui = lmemmusample[0];
  memnui = memnusample[0];
  rmemnui = sqrt(memnui);

  // Compute initial current lik
  lliki = 0.0;
  for (i=0; i<ni; i++){
    temp = beta0i + betati*rectime[i];
    rtprob = exp(temp)/(1.0+exp(temp));
    if(numrec[i] <= u[i]){
     if((u[i] <= maxc)|(numrec[i]<maxc)){
       lliki += dbinom(numrec[i],u[i],rtprob,give_log1);
     }else{
       lliki += log(1.0-pbinom(maxc-1.0,u[i],rtprob,give_log0,give_log0));
     }
     if(srd[i]>=0){
//    Use WP for localized
      memmui = exp(lmemmui)*u[i];
      rnb=memmui/(alpha-1.);
      pd[0]= memmui*exp(-fabs(memmui-1.)/sqrt(memnui));
      pis=pd[0];
      for (k=1; k<(10*Ki); k++){
        pd[k]= pd[k-1]*(k+rnb)*pnb*exp((fabs(k-memmui)-fabs(k+1-memmui))/sqrt(memnui))/((double)(k+1));
        pis+=pd[k];
      }
      for (k=0; k<(10*Ki); k++){
        pd[k]/=pis;
      }
      pis2=0.;
      for (k=Ki; k<(10*Ki); k++){
        pis2+=pd[k];
      }
      if(srd[i] <= 10*Ki){
        lliki += log(pd[srd[i]-1]);
      }else{
        lliki += log(pis2);
      }
      lliki += log((double)(u[i]));
     }
    }
  }
  if(!isfinite(lliki)) lliki = -100000.0; 

  // Compute initial prior
  pibeta0i = dnorm(beta0i, dbeta0, dbeta0sd, give_log1);
  if(dbetatsd > 0.0) pibetati = dnorm(betati, dbetat, dbetatsd, give_log1);
  pimemi = dnorm(lmemmui, dlmemmu, rmemnui/rmemdfmu, give_log1);
  pimemi = pimemi+dsclinvchisq(memnui, dmemdfnu, dmemnu);

  // Now do the MCMC updates (starting with the burnin updates)
  while (isamp < isamplesize && step < iburnin) {
    /* Propose new beta */
    beta0star = rnorm(beta0i, dbeta0proposal);
    if(dbetatsd > 0.0){
      betatstar = rnorm(betati, dbetatproposal);
    }else{
      betatstar = betati;
    }
    /* Propose new memnu and lmemmu */
    lmemmustar = rnorm(lmemmui, dlmemmuproposal);
// VIP Remember this next line hold the optimism fixed at 1!!! VIP
    lmemmustar = lmemmui;

//  lmemmustar = 0.0;
//  rmemnustar = rnorm(memnui, dmemnuproposal);
    memnustar = memnui*exp(rnorm(0., dmemnuproposal));
    rmemnustar = sqrt(memnustar);

    llikstar = 0.0;
    for (i=0; i<ni; i++){
      temp = beta0star + betatstar*rectime[i];
      rtprob = exp(temp)/(1.0+exp(temp));
      if(numrec[i] <= u[i]){
       if((u[i] <= maxc)|(numrec[i]<maxc)){
        llikstar += dbinom(numrec[i],u[i],rtprob,give_log1);
       }else{
        llikstar += log(1.0-pbinom(maxc-1.0,u[i],rtprob,give_log0,give_log0));
       }
       if(srd[i]>=0){
        memmustar = exp(lmemmustar)*u[i];
        rnb=memmustar/(alpha-1.);
        pd[0]= exp(-fabs(memmustar-1.)/sqrt(memnustar));
        pis=pd[0];
        for (k=1; k<(10*Ki); k++){
          pd[k]= pd[k-1]*(k+rnb)*pnb*exp((fabs(k-memmustar)-fabs(k+1-memmustar))/sqrt(memnustar))/((double)(k+1));
          pis+=pd[k];
        }
        for (k=0; k<(10*Ki); k++){
          pd[k]/=pis;
        }
        pis2=0.;
        for (k=Ki; k<(10*Ki); k++){
          pis2+=pd[k];
        }
        if(srd[i] <= 10*Ki){
          llikstar += log(pd[srd[i]-1]);
        }else{
          llikstar += log(pis2);
        }
        llikstar += log((double)(u[i]));
       }
      }
    }
    if(!isfinite(llikstar)) llikstar = -100000.0; 

    /* Calculate pieces of the prior. */
    pibeta0star = dnorm(beta0star, dbeta0, dbeta0sd, give_log1);
    if(dbetatsd > 0.0) pibetatstar = dnorm(betatstar, dbetat, dbetatsd, give_log1);
    pimemstar = dnorm(lmemmustar, dlmemmu, rmemnustar/rmemdfmu, give_log1);
    pimemstar = pimemstar+dsclinvchisq(memnustar, dmemdfnu, dmemnu);
//  pimemstar = dsclinvchisq(memnustar, dmemdfnu, dmemnu);

    qi = dnorm(log(memnui/memnustar)/dmemnuproposal,0.,1.,give_log1)
         -log(dmemnuproposal*memnui);

    qstar = dnorm(log(memnustar/memnui)/dmemnuproposal,0.,1.,give_log1)
                  -log(dmemnuproposal*memnustar);

    /* Calculate ratio */
    ip =      pibeta0star-pibeta0i;
    if(dbetatsd > 0.0) ip = ip + pibetatstar-pibetati;
    ip = ip + pimemstar - pimemi;

    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
    then let the MH probability equal min{exp(cutoff), 1.0}.
    But we'll do it in log space instead.  */
    cutoff = ip + llikstar - lliki + qi - qstar;
      
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed changes */
      beta0i = beta0star;
      betati = betatstar;
      lmemmui = lmemmustar;
      memnui = memnustar;
      rmemnui = rmemnustar;
      lliki = llikstar;
      qi = qstar;
      pibeta0i = pibeta0star;
      if(dbetatsd > 0.0) pibetati = pibetatstar;
      pimemi = pimemstar;
      taken++;
      if (step > 0 && step==(iinterval*(step/iinterval))) { 
        /* record statistics for posterity */
        beta0sample[isamp]=beta0i;
        betatsample[isamp]=betati;
        lmemmusample[isamp]=lmemmui;
        memnusample[isamp]=memnui;
        isamp++;
        if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d MH samples...\n", isamp);
      }
    }
    step++;
  }
//PutRNGstate();  /* Disable RNG before returning */
  free(pd);
  /*Check for interrupts (if recursion is taking way too long...)*/
  R_CheckUserInterrupt();
  *staken = taken;
}
