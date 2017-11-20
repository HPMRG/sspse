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
	    int *Npi,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *maxN, 
            double *sample, 
            double *ppos, 
            int *burnintheta,
	    int *verbose
			 ) {
  int step, staken, getone=1, intervalone=1, verboseMHpln = 0;
  int dimsample, Np;
  int i, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  double mu, mui, sigmai;
  double dkappa0, ddf0, dmu0, dsigma0, dmuproposal, dsigmaproposal;
  int tU, popi, imaxN, imaxm;
  double r, gammart, pis;

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
  dkappa0=(*kappa0);
  ddf0=(*df0);
  dsigma0=(*sigma0);
  dmu0=(*mu0);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  dimsample=4+Np;

  double *pi = (double *) malloc(sizeof(double) * Ki);
  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  int *Nkpos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdegi = (double *) malloc(sizeof(double) * (Np+1));
  double *psample = (double *) malloc(sizeof(double) * (Np+1));
  double *musample = (double *) malloc(sizeof(double));
  double *sigmasample = (double *) malloc(sizeof(double));

  b[ni-1]=pop[ni-1];
  for (i=(ni-2); i>=0; i--){
    b[i]=b[i+1]+pop[i];
  }
  for (i=0; i<Ki; i++){
     Nk[i]=nk[i];
     Nkpos[i]=0;
     ppos[i]=0.;
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

  for (i=0; i<Np; i++){
     psample[i] = 0.01;
  }
  musample[0] = dmu0;
  sigmasample[0] = dsigma0;

  isamp = 0;
  step = -iburnin;
  while (isamp < isamplesize) {
    /* Draw new theta */
    MHplnorig(Nk,K,mu0,kappa0,sigma0,df0,muproposal,sigmaproposal,
          &Ni, &Np, psample,
	  musample, sigmasample, &getone, &staken, burnintheta, &intervalone, 
	  &verboseMHpln);

    for (i=0; i<Np; i++){
      pdegi[i] = psample[i];
    }
    mui=musample[0];
    sigmai=sigmasample[0];

    /* First find the degree distribution */
    pis=0.;
    for (i=Np; i<Ki; i++){
      pi[i]=poilog(i+1,mui,sigmai);
//    pis+=pi[i];
    }
    pis=1.-poilog(0,mui,sigmai);
    for (i=Np; i<Ki; i++){
      pi[i]/=pis;
    }
    pis=1.;
    for (i=0; i<Np; i++){
      pis-=pdegi[i];
    }
    for (i=0; i<Ki; i++){
      pi[i]*=pis;
    }
    for (i=0; i<Np; i++){
      pi[i]=pdegi[i];
    }
    gammart=0.;
    for (i=0; i<Ki; i++){
      gammart+=(exp(-r*(i+1))*pi[i]);
    }
    gammart=log(gammart);
    tU = -1000000000;
    for (i=0; i<imaxm; i++){
      lpm[i]=i*gammart+lgamma(ni+i+1.)-lgamma(i+1.);
      if(lpm[i] > tU) tU = (int)lpm[i];
    }
    for (i=0; i<imaxm; i++){
      lpm[i]=exp(lpm[i]-tU);
    }
    for (i=1; i<imaxm; i++){
      lpm[i]=lpm[i-1]+lpm[i];
    }
    gammart = lpm[imaxm-1] * unif_rand();
    for (Ni=0; Ni<imaxm; Ni++){
      if(gammart <= lpm[Ni]) break;
    }
    Ni+= ni;
    if(Ni > imaxN) Ni = imaxN;
		    
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
      /* Use rejection sampling */
      popi=0;
      while((popi == 0) || (log(1.0-unif_rand()) > -r*popi)){
       mu = exp(rnorm(mui, sigmai));
       if(mu < 5.*Ki){
         popi=(int)rpois(mu);
         if(popi < 0){popi=0;}
       }else{
    if (*verbose) Rprintf("mu > 5.*Ki; mu %f K %d\n", mu, Ki);
         popi=0;
       }
      }
      if(popi > Ki){popi=Ki;}
      pop[i]=popi;
      Nk[popi-1]=Nk[popi-1]+1;
    }
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
      sample[isamp*dimsample  ]=(double)(Ni);
      sample[isamp*dimsample+1]=mui;
      sample[isamp*dimsample+2]=sigmai;
      sample[isamp*dimsample+3]=(double)(Nk[0]);
      for (i=0; i<Np; i++){
        sample[isamp*dimsample+4+i]=pdegi[i];
      }
      for (i=0; i<Ki; i++){
        Nkpos[i]=Nkpos[i]+Nk[i];
        ppos[i]+=((Nk[i]*1.)/Ni);
      }
      isamp++;
      if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
      if (*verbose) Rprintf("Taken %d samples...\n", isamp);
    }
    step++;
  }
  for (i=0; i<Ki; i++){
    nk[i]=Nkpos[i];
    ppos[i]/=isamp;
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(pi);
  free(b);
  free(Nk);
  free(Nkpos);
  free(lpm);
  free(pdegi);
  free(psample);
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
	    int *Npi,
            double *muproposal, 
            double *sigmaproposal, 
            int *N,
            int *sample, 
            int *burnintheta,
	    int *verbose
			 ) {
  int Np;
  int step, staken, getone=1, intervalone=1, verboseMHpln = 0;
  int i, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  double mu;
  double dkappa0, ddf0, dmu0, dsigma0, dmuproposal, dsigmaproposal;
  int tU, popi;
  double r;

  GetRNGstate();  /* R function enabling uniform RNG */

  ni=(*n);
  Ni=(*N);
  Ki=(*K);
  Np=(*Npi);
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
  double *psample = (double *) malloc(sizeof(double) * (Np+1));
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
    MHplnorig(Nk,K,mu0,kappa0,sigma0,df0,muproposal,sigmaproposal,N,
          &Np, psample,
	  musample, sigmasample, &getone, &staken, burnintheta, &intervalone, 
	  &verboseMHpln);
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
      popi=0;
      while((popi == 0) || (log(1.0-unif_rand()) > -r*popi)){
       if(mu < 5.*Ki){
        popi=(int)rpois(mu);
        if(popi < 0){popi=0;}
       }else{
        popi=0;
       }
      }
      if(popi > Ki){popi=Ki;}
      pop[i]=popi;
      Nk[popi-1]=Nk[popi-1]+1;
//    if (*verbose) Rprintf("Ni %d ni %d mu %f pop[i] %d r %f\n", Ni, ni, mu, pop[i], r);
    }
		    
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
      if (*verbose) Rprintf("isamp %d pop[501] %d\n", isamp, pop[501]);
      for (i=0; i<Ni; i++){
        sample[isamp*Ni+i]=pop[i];
      }
      isamp++;
      if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
    }
    step++;
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(b);
  free(Nk);
  free(phi);
  free(psample);
  free(musample);
  free(sigmasample);
}

void MHplnorig (int *Nk, int *K,
	    double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Npi, double *psample,
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *verbose
			 ) {
  int Np;
  int step, taken, give_log=1;
  int i, Ki, Ni, isamp, iinterval, isamplesize, iburnin;
  double ip, cutoff;
  double mustar, mui, lp;
  double pis, pstars;
  double sigmastar, sigmai, sigma2star, sigma2i, qsigma2star, qsigma2i;
  double pithetastar, pithetai;
  double dkappa0, rkappa0, ddf0, dmu0, logK;
  double dsigma0, dsigma20, dmuproposal, dsigmaproposal;

  GetRNGstate();  /* R function enabling uniform RNG */

  Ki=(*K);
  Np=(*Npi);
  double *pstar = (double *) malloc(sizeof(double) * Ki);
  double *pi = (double *) malloc(sizeof(double) * Ki);
  double *odegstar = (double *) malloc(sizeof(double) * (Np+1));
  double *odegi = (double *) malloc(sizeof(double) * (Np+1));
  double *pdegstar = (double *) malloc(sizeof(double) * (Np+1));
  double *pdegi = (double *) malloc(sizeof(double) * (Np+1));

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
  logK=log((double)Ki);
  isamp = taken = 0;
  step = -iburnin;
  pis=1.;
  for (i=0; i<Np; i++){
    pdegi[i] = psample[i];
    pis-=pdegi[i];
  }
  for (i=0; i<Np; i++){
    odegi[i] = log(pdegi[i]/pis);
  }
  mui = musample[0];
  sigmai = sigmasample[0];
  sigma2i  = sigmai*sigmai;
  pithetai = dnorm(mui, dmu0, sigmai/rkappa0, give_log);
  pithetai = pithetai+dsclinvchisq(sigma2i, ddf0, dsigma20);
  pis=0.;
  for (i=Np; i<Ki; i++){
   pi[i]=poilog(i+1,mui,sigmai);
// pis+=pi[i];
  }
  pis=1.-poilog(0,mui,sigmai);
  for (i=Np; i<Ki; i++){
   pi[i]/=pis;
  }
  pis=1.;
  for (i=0; i<Np; i++){
    pis-=pdegi[i];
  }
  for (i=0; i<Np; i++){
    pi[i]=pdegi[i];
  }
  for (i=Np; i<Ki; i++){
    pi[i]*=pis;
  }
  while (isamp < isamplesize) {
    /* Propose new theta */
    /* Now the degree distribution model parameters */
    for (i=0; i<Np; i++){
      odegstar[i] = rnorm(odegi[i], dmuproposal);
    }
    /* Convert from odds to probabilities */
    pis=1.;
    for (i=0; i<Np; i++){
      pdegstar[i] = exp(odegstar[i]);
      pis+=pdegstar[i];
    }
    for (i=0; i<Np; i++){
      pdegstar[i]/=pis;
    }
    /* Now the degree distribution (log) mean and s.d. parameters */
    mustar = rnorm(mui, dmuproposal);
    sigma2star = sigma2i*exp(rnorm(0., dsigmaproposal));
    sigmastar = sqrt(sigma2star);
    /* Check for magnitude */

  if(sigma2star > 1000) Rprintf("%f %f %f %f %f\n", mustar, dmu0, sigma2star, dkappa0, sigma2i);
    /* Calculate pieces of the posterior. */
    qsigma2star = dnorm(log(sigma2star/sigma2i)/dsigmaproposal,0.,1.,give_log)
                  -log(dsigmaproposal*sigma2star);
    pithetastar = dnorm(mustar, dmu0, sigmastar/rkappa0, give_log);
    pithetastar = pithetastar+dsclinvchisq(sigma2star, ddf0, dsigma20);
    qsigma2i = dnorm(log(sigma2i/sigma2star)/dsigmaproposal,0.,1.,give_log)
               -log(dsigmaproposal*sigma2i);

    /* Calculate ratio */
    ip = pithetastar-pithetai;
    pstars=0.;
    for (i=Np; i<Ki; i++){
      pstar[i]=poilog(i+1,mustar,sigmastar);
//    pstars+=pstar[i];
    }
    pstars=1.-poilog(0,mustar,sigmastar);
    for (i=Np; i<Ki; i++){
     pstar[i]/=pstars;
    }
    pstars=1.;
    for (i=0; i<Np; i++){
      pstars-=pdegstar[i];
    }
    for (i=0; i<Np; i++){
      pstar[i]=pdegstar[i];
    }
    for (i=Np; i<Ki; i++){
      pstar[i]*=pstars;
    }
    for (i=0; i<Ki; i++){
     if(Nk[i]>0){
      lp = log(pstar[i]/pi[i]);
      if(fabs(lp)<100.){ip += (Nk[i]*lp);}
     }
//    Rprintf("%d %f\n", i, log(poilog(s[i],mustar,sigmastar)/poilog(s[i],mui,sigmai)));
    }
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
    then let the MH probability equal min{exp(cutoff), 1.0}.
    But we'll do it in log space instead.  */
//  if (*verbose)
//    Rprintf("Now proposing %d MH steps %f ip1...\n", step, ip);
    cutoff = ip + qsigma2i-qsigma2star;
      
//  if (*verbose)
//    Rprintf("Now proposing %d MH steps %f cutoff...\n", step, cutoff);

    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed changes */
      for (i=0; i<Np; i++){
        odegi[i] = odegstar[i];
        pdegi[i] = pdegstar[i];
      }
      mui    = mustar;
      sigmai = sigmastar;
      sigma2i = sigma2star;
      qsigma2i = qsigma2star;
      pithetai = pithetastar;
      for (i=0; i<Ki; i++){
        pi[i] = pstar[i];
      }
      taken++;
      if (step > 0 && step==(iinterval*(step/iinterval))) { 
        /* record statistics for posterity */
        musample[isamp]=mui;
        sigmasample[isamp]=sigmai;
        for (i=0; i<Np; i++){
          psample[i]=pdegi[i];
        }
        isamp++;
        if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d MH samples...\n", isamp);
      }
    }
    step++;
  }
  free(pi);
  free(pstar);
  free(odegstar);
  free(odegi);
  free(pdegstar);
  free(pdegi);
  PutRNGstate();  /* Disable RNG before returning */
  *staken = taken;
}

void MHpriorpln (double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *verbose
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
//    if (*verbose) {
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
//    if (*verbose) Rprintf("step %d mui %f sigmai %f\n", step, mui, sigmai);
      taken++;
      if (step > 0 && step==(iinterval*(step/iinterval))) { 
        /* record statistics for posterity */
        musample[isamp]=mui;
        sigmasample[isamp]=sigmai;
        isamp++;
        if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d MH samples...\n", isamp);
      }
    }
    step++;
  }
  PutRNGstate();  /* Disable RNG before returning */
  *staken = taken;
}
