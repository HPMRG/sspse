/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "mcmc.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gspps (int *pop,
            int *nk, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *maxN, 
            double *sample, 
            int *burnintheta,
	    int *fVerbose
			 ) {
  int step, staken, getone=1, intervalone=1, fVerboseMHpln = 0;
  int i, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  double mu, mu0i, sigma;
  double dkappa0, ddf0, dmu0, dsigma0, dmuproposal, dsigmaproposal;
  int tU, popi, imaxN;
  double r, gammart, pis;

  GetRNGstate();  /* R function enabling uniform RNG */

  ni=(*n);
  Ni=(*N);
  Ki=(*K);
  imaxN=(*maxN);
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburnin=(*burnin);
  dkappa0=(*kappa0);
  ddf0=(*df0);
  dsigma0=(*sigma0);
  dmu0=(*mu0);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  double *pi = (double *) malloc(sizeof(double) * Ki);
  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  int *Nkpos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxN);
  double *musample = (double *) malloc(sizeof(double));
  double *sigmasample = (double *) malloc(sizeof(double));

  b[ni-1]=pop[ni-1];
  for (i=(ni-2); i>=0; i--){
    b[i]=b[i+1]+pop[i];
  }
  for (i=0; i<Ki; i++){
     Nk[i]=nk[i];
     Nkpos[i]=isamplesize*nk[i];
  }
  tU=0;
  for (i=ni; i<Ni; i++){
    tU+=pop[i];
  }
  /* Draw phis */
  r=0.;
  for (i=0; i<ni; i++){
    r+=(exp_rand()/(tU+b[i]));
  }

  musample[0] = dmu0;
  sigmasample[0] = dsigma0;

  isamp = 0;
  step = -iburnin;
  while (isamp < isamplesize) {
    /* Draw new theta */
    MHpln(Nk,K,mu0,kappa0,sigma0,df0,muproposal,sigmaproposal,&Ni,
	  musample, sigmasample, &getone, &staken, burnintheta, &intervalone, 
	  &fVerboseMHpln);

    mu0i=musample[0];
    sigma=sigmasample[0];

    /* Draw new N */
    pis=0.;
    for (i=0; i<Ki; i++){
      pi[i]=poilog(i+1,mu0i,sigma);
      pis+=pi[i];
    }
    gammart=0.;
    for (i=0; i<Ki; i++){
      gammart+=(exp(-r*(i+1))*pi[i]);
    }
    gammart=log(gammart);
    tU = -1000000000;
    for (i=0; i<imaxN; i++){
      lpm[i]=i*gammart+lgamma(ni+i+1.)-lgamma(i+1.);
      if(lpm[i] > tU) tU = lpm[i];
    }
    for (i=0; i<imaxN; i++){
      lpm[i]=exp(lpm[i]-tU);
    }
    for (i=1; i<imaxN; i++){
      lpm[i]=lpm[i-1]+lpm[i];
    }
    gammart = lpm[imaxN-1] * unif_rand();
    Ni = 0;
    while(gammart > lpm[Ni]){Ni++;}
    Ni+= ni;
    if(Ni >= imaxN) Ni = imaxN-1;
		    
    /* Draw phis */
    tU=0;
    for (i=ni; i<Ni; i++){
      tU+=pop[i];
    }
    r=0.;
    for (i=0; i<ni; i++){
//    phi[i]=(tU+b[i])*exp_rand();
      r+=exp_rand()/(tU+b[i]);
    }

    /* Draw unseen sizes */
    for (i=0; i<Ki; i++){
      Nk[i]=nk[i];
    }
    for (i=ni; i<Ni; i++){
      /* Propose unseen size for unit i */
      mu = exp(rnorm(mu0i, sigma));
      /* Use rejection sampling */
      popi=rpois(mu);
      while((popi == 0) || (log(1.0-unif_rand()) > -r*popi)){
       mu = exp(rnorm(mu0i, sigma));
       popi=rpois(mu);
      }
      if(popi > Ki){popi=Ki;}
      pop[i]=popi;
      Nk[popi-1]=Nk[popi-1]+1;
    }
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
      sample[isamp*4  ]=(double)(Ni);
      sample[isamp*4+1]=mu0i;
      sample[isamp*4+2]=sigma;
      sample[isamp*4+3]=(double)(Nk[0]);
      for (i=0; i<Ki; i++){
        Nkpos[i]=Nkpos[i]+Nk[i];
      }
      isamp++;
      if (*fVerbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
    }
    step++;
  }
  for (i=0; i<Ki; i++){
    nk[i]=Nkpos[i];
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(pi);
  free(b);
  free(Nk);
  free(Nkpos);
  free(lpm);
  free(musample);
  free(sigmasample);
}

void gsppsN (int *pop,
            int *nk, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, 
            int *sample, 
            int *burnintheta,
	    int *fVerbose
			 ) {
  int step, staken, getone=1, intervalone=1, fVerboseMHpln = 0;
  int i, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  double mu;
  double dkappa0, ddf0, dmu0, dsigma0, dmuproposal, dsigmaproposal;
  int tU, popi;
  double r;

  GetRNGstate();  /* R function enabling uniform RNG */

  ni=(*n);
  Ni=(*N);
  Ki=(*K);
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburnin=(*burnin);
  dkappa0=(*kappa0);
  ddf0=(*df0);
  dsigma0=(*sigma0);
  dmu0=(*mu0);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  double *phi = (double *) malloc(sizeof(double) * ni);
  double *musample = (double *) malloc(sizeof(double));
  double *sigmasample = (double *) malloc(sizeof(double));
  b[ni-1]=pop[ni-1];
  for (i=(ni-2); i>=0; i--){
    b[i]=b[i+1]+pop[i];
  }
  for (i=0; i<Ki; i++){
     Nk[i]=nk[i];
  }

  isamp = 0;
  step = -iburnin;
  while (isamp < isamplesize) {
    /* Draw new theta */
    MHpln(Nk,K,mu0,kappa0,sigma0,df0,muproposal,sigmaproposal,N,
	  musample, sigmasample, &getone, &staken, burnintheta, &intervalone, 
	  &fVerboseMHpln);
    tU=0;
    for (i=ni; i<Ni; i++){
      tU+=pop[i];
    }
    /* Draw phis */
    r=0.;
    for (i=0; i<ni; i++){
//    phi[i]=(tU+b[i])*exp_rand();
      phi[i]=exp_rand()/(tU+b[i]);
      r+=phi[i];
    }
    for (i=0; i<Ki; i++){
      Nk[i]=nk[i];
    }
    /* Draw unseen sizes */
    for (i=ni; i<Ni; i++){
      /* Propose unseen size for unit i */
      mu = exp(rnorm(musample[0], sigmasample[0]));
      /* Use rejection sampling */
      popi=rpois(mu);
      while(log(1.0-unif_rand()) > -r*popi){
       popi=rpois(mu);
      }
      if(popi > Ki){popi=Ki;}
      pop[i]=popi;
      Nk[popi]=Nk[popi]+1;
//    if (*fVerbose) Rprintf("Ni %d ni %d mu %f pop[i] %d r %f\n", Ni, ni, mu, pop[i], r);
    }
		    
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
      if (*fVerbose) Rprintf("isamp %d pop[501] %d\n", isamp, pop[501]);
      for (i=0; i<Ni; i++){
        sample[isamp*Ni+i]=pop[i];
      }
      isamp++;
      if (*fVerbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
    }
    step++;
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(b);
  free(Nk);
  free(phi);
  free(musample);
  free(sigmasample);
}

void MHpln (int *nk, int *K,
	    double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, 
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *fVerbose
			 ) {
  int step, taken, give_log=1;
  int i, Ki, Ni, isamp, iinterval, isamplesize, iburnin;
  double ip, cutoff;
  double mustar, mui, lp;
  double pis, pstars;
  double sigmastar, sigmai, sigma2star, sigma2i, qsigma2star, qsigma2i;
  double pithetastar, pithetai;
  double dkappa0, rkappa0, ddf0, dmu0;
  double dsigma0, dsigma20, dmuproposal, dsigmaproposal;

  GetRNGstate();  /* R function enabling uniform RNG */

  Ki=(*K);
  double *pstar = (double *) malloc(sizeof(double) * Ki);
  double *pi = (double *) malloc(sizeof(double) * Ki);

  Ni=(*N);
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburnin=(*burnin);
  dkappa0=(*kappa0);
  rkappa0=sqrt(dkappa0);
  ddf0=(*df0);
  dsigma0=(*sigma0);
  dsigma20=(dsigma0*dsigma0);
  dmu0=(*mu0);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);
  isamp = taken = 0;
  step = -iburnin;
  mui = musample[0];
  sigmai = sigmasample[0];
  sigma2i  = sigmai*sigmai;
  pithetai = dnorm(mui, dmu0, sigmai/rkappa0, give_log);
  pithetai = pithetai+dsclinvchisq(sigma2i, ddf0, dsigma20);
  pis=0.;
  for (i=0; i<Ki; i++){
    pi[i] = poilog(i+1,mui,sigmai);
    pis+=pi[i];
  }
  while (isamp < isamplesize) {
    /* Propose new theta */
    mustar = rnorm(mui, dmuproposal);
    sigma2star = sigma2i*exp(rnorm(0., dsigmaproposal));
    sigmastar = sqrt(sigma2star);

//  Rprintf("%f %f %f %f %f\n", mustar, dmu0, sigma2star, dkappa0, sigma2i);
    /* Calculate pieces of the posterior. */
    qsigma2star = dnorm(log(sigma2star/sigma2i)/dsigmaproposal,0.,1.,give_log)
                  -log(dsigmaproposal*sigma2star);
    pithetastar = dnorm(mustar, dmu0, sigmastar/rkappa0, give_log);
    pithetastar = pithetastar+dsclinvchisq(sigma2star, ddf0, dsigma20);
    qsigma2i = dnorm(log(sigma2i/sigma2star)/dsigmaproposal,0.,1.,give_log)
               -log(dsigmaproposal*sigma2i);

    /* Calculate ratio */
    ip = pithetastar-pithetai;
//  if (*fVerbose){
//    Rprintf("pithetastar %f pithetai %f qsigma2star %f qsigmai %f\n",
//           pithetastar,pithetai,qsigma2star,qsigma2i);
//  }
    pstars=0.;
    for (i=0; i<Ki; i++){
      pstar[i] = poilog(i+1,mustar,sigmastar);
      pstars+=pstar[i];
    }
    for (i=0; i<Ki; i++){
//  Rprintf("%d %f %f\n", i,poilog(s[i],mustar,sigmastar),
//    poilog(s[i],mui,sigmai));
      lp = log(pis*pstar[i]/(pstars*pi[i]));
      if((lp > -100.) && (lp<100.)){ip += (nk[i]*lp);}
//    Rprintf("%d %f\n", i, log(poilog(s[i],mustar,sigmastar)/poilog(s[i],mui,sigmai)));
    }
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
    then let the MH probability equal min{exp(cutoff), 1.0}.
    But we'll do it in log space instead.  */
//  if (*fVerbose)
//    Rprintf("Now proposing %d MH steps %f ip1...\n", step, ip);
    cutoff = ip + qsigma2i-qsigma2star;
      
//  if (*fVerbose)
//    Rprintf("Now proposing %d MH steps %f cutoff...\n", step, cutoff);

    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed changes */
      sigmai = sigmastar;
      mui    = mustar;
      sigma2i = sigma2star;
      qsigma2i = qsigma2star;
      pithetai = pithetastar;
      pis=pstars;
      for (i=0; i<Ki; i++){
        pi[i] = pstar[i];
      }
      taken++;
      if (step > 0 && step==(iinterval*(step/iinterval))) { 
        /* record statistics for posterity */
        musample[isamp]=mui;
        sigmasample[isamp]=sigmai;
        isamp++;
        if (*fVerbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d MH samples...\n", isamp);
      }
    }
    step++;
  }
  free(pi);
  free(pstar);
  PutRNGstate();  /* Disable RNG before returning */
  *staken = taken;
}

void MHpriorpln (double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *fVerbose
			 ) {
  int step, taken, give_log=1;
  int isamp, iinterval, isamplesize, iburnin;
  double ip, cutoff;
  double mustar, mui;
  double sigmastar, sigmai, sigma2star, sigma2i, qsigma2star, qsigma2i;
  double pithetastar, pithetai;
  double dkappa0, rkappa0, ddf0, dmu0;
  double dsigma0, dsigma20, dmuproposal, dsigmaproposal;

  GetRNGstate();  /* R function enabling uniform RNG */

  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburnin=(*burnin);
  dkappa0=(*kappa0);
  rkappa0=sqrt(dkappa0);
  ddf0=(*df0);
  dsigma0=(*sigma0);
  dsigma20=(dsigma0*dsigma0);
  dmu0=(*mu0);
  dmuproposal=(*muproposal);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);
  isamp = taken = 0;
  step = -iburnin;
  mui = dmu0;
  sigma2i = dsigma20;
  sigmai  = sqrt(sigma2i);
  pithetai = dnorm(mui, dmu0, sigmai/rkappa0, give_log);
  pithetai = pithetai+dsclinvchisq(sigma2i, ddf0, dsigma20);
  while (isamp < isamplesize) {
    /* Propose new theta */
    mustar = rnorm(mui, dmuproposal);
    sigma2star = sigma2i*exp(rnorm(0., dsigmaproposal));
    sigmastar = sqrt(sigma2star);

    /* Calculate pieces of the posterior. */
    qsigma2star = dnorm(log(sigma2star/sigma2i)/dsigmaproposal,0.,1.,give_log)
                  -log(dsigmaproposal*sigma2star);
    pithetastar = dnorm(mustar, dmu0, sigmastar/rkappa0, give_log);
    pithetastar = pithetastar+dsclinvchisq(sigma2star, ddf0, dsigma20);
    qsigma2i = dnorm(log(sigma2i/sigma2star)/dsigmaproposal,0.,1.,give_log)
               -log(dsigmaproposal*sigma2i);

    /* Calculate ratio */
    ip = pithetastar-pithetai;
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
    then let the MH probability equal min{exp(cutoff), 1.0}.
    But we'll do it in log space instead.  */
    cutoff = ip + qsigma2i-qsigma2star;
      
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
//    if (*fVerbose) {
//     sigmai = exp(mustar+0.5*sigmastar);
//     sigmai = exp(mustar+0.5*sigmastar);
//     if(sigmai > 12.){
//     Rprintf("step %d mu %f pithetastar %f pithetai %f qsigma2star %f qsigma2i %f mustar %f sigmastar %f cutoff %f\n", 
//		           step, sigmai, pithetastar, pithetai, qsigma2star, qsigma2i, mustar, sigmastar, cutoff);
//			   }
//			   }
      /* Make proposed changes */
      sigmai = sigmastar;
      mui    = mustar;
      sigma2i = sigma2star;
      pithetai = pithetastar;
//    if (*fVerbose) Rprintf("step %d mui %f sigmai %f\n", step, mui, sigmai);
      taken++;
      if (step > 0 && step==(iinterval*(step/iinterval))) { 
        /* record statistics for posterity */
        musample[isamp]=mui;
        sigmasample[isamp]=sigmai;
        isamp++;
        if (*fVerbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d MH samples...\n", isamp);
      }
    }
    step++;
  }
  PutRNGstate();  /* Disable RNG before returning */
  *staken = taken;
}
