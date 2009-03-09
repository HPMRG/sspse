/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "dis.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gsppsdis (int *pop, int *dis,
            int *nk0, int *nk1, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu0, double *mu1, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *maxN, 
            double *sample, 
            int *burnintheta,
            int *fVerbose
              ) {
  int step, staken, getone=1, intervalone=1, fVerboseMHdis = 0;
  int i, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  double mu, mu0i, mu1i, pbeta, beta, sigma;
  double dkappa0, ddf0, dmu0, dmu1, dsigma0, dmuproposal, dsigmaproposal;
  int tU, popi, imaxN, itotdis0, itotdis;
  double r, gamma0rt, gamma1rt, p0is, p1is;

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
  dmu1=(*mu1);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  double *p0i = (double *) malloc(sizeof(double) * Ki);
  double *p1i = (double *) malloc(sizeof(double) * Ki);
  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk0 = (int *) malloc(sizeof(int) * Ki);
  int *Nk0pos = (int *) malloc(sizeof(int) * Ki);
  int *Nk1 = (int *) malloc(sizeof(int) * Ki);
  int *Nk1pos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxN);
  double *musample = (double *) malloc(sizeof(double) * 2);
  double *betasample = (double *) malloc(sizeof(double));
  double *sigmasample = (double *) malloc(sizeof(double));

  b[ni-1]=pop[ni-1];
  for (i=(ni-2); i>=0; i--){
    b[i]=b[i+1]+pop[i];
  }
  for (i=0; i<Ki; i++){
     Nk0[i]=nk0[i];
     Nk0pos[i]=isamplesize*nk0[i];
     Nk1[i]=nk1[i];
     Nk1pos[i]=isamplesize*nk1[i];
  }
  tU=0;
  for (i=ni; i<Ni; i++){
    tU+=pop[i];
  }
  /* Draw phis */
  r=0.;
  itotdis0=0;
  for (i=0; i<ni; i++){
    r+=(exp_rand()/(tU+b[i]));
    itotdis0+=dis[i];
  }
  itotdis=itotdis0;
  for (i=ni; i<Ni; i++){
    itotdis+=dis[i];
  }

  betasample[0] = -1.386294;
  musample[0] = dmu0;
  musample[1] = dmu1;
  sigmasample[0] = dsigma0;

  isamp = 0;
  step = -iburnin;
  while (isamp < isamplesize) {
    /* Draw new theta */
    MHdis(Nk0,Nk1,&itotdis,K,mu0,mu1,kappa0,sigma0,df0,muproposal,sigmaproposal,
	  &Ni, musample, betasample, sigmasample, &getone, &staken, 
	  burnintheta, &intervalone, 
	  &fVerboseMHdis);

    beta=betasample[0];
    pbeta=exp(beta)/(1.+exp(beta));
    mu0i=musample[0];
    mu1i=musample[1];
    sigma=sigmasample[0];

    /* Draw new N */
    p0is=0.;
    p1is=0.;
    for (i=0; i<Ki; i++){
      p0i[i]=poilog(i+1,mu0i,sigma);
      p1i[i]=poilog(i+1,mu1i,sigma);
      p0is+=p0i[i];
      p1is+=p1i[i];
    }
    gamma0rt=0.;
    gamma1rt=0.;
    for (i=0; i<Ki; i++){
      gamma0rt+=(exp(-r*(i+1))*p0i[i]);
      gamma1rt+=(exp(-r*(i+1))*p1i[i]);
    }
//    gamma0rt=log(gamma0rt);
//    gamma1rt=log(gamma1rt);
    gamma0rt=log((1.-pbeta)*gamma0rt+pbeta*gamma1rt);
    tU = -1000000000;
    for (i=0; i<imaxN; i++){
     lpm[i]=lgamma(ni+i+1.)-lgamma(i+1)+i*gamma0rt;
     if(lpm[i] > tU) tU = lpm[i];
    }
    for (i=0; i<imaxN; i++){
      lpm[i]=exp(lpm[i]-tU);
    }
    for (i=1; i<imaxN; i++){
      lpm[i]=lpm[i-1]+lpm[i];
    }
    gamma0rt = lpm[imaxN-1] * unif_rand();
    Ni = 0;
    while(gamma0rt > lpm[Ni]){Ni++;}
//  if (*fVerbose) Rprintf("Ni %d lpm[imaxN-1] %f lpm[Ni] %f\n", Ni, lpm[imaxN-1],
//  lpm[Ni]);
//  }
    Ni+= ni;
    if(Ni >= imaxN) Ni = imaxN-1;

//  if (*fVerbose) Rprintf("step %d Ni %d itotdis %d beta %f mu0 %f mu1 %f s0 %f r %f\n",
//  step, Ni, itotdis, betasample[0], musample[0], musample[1], sigmasample[0], r);

    /* Draw phis */
    tU=0;
    for (i=ni; i<Ni; i++){
      tU+=pop[i];
    }
    r=0.;
    for (i=0; i<ni; i++){
      r+=exp_rand()/(tU+b[i]);
    }

    /* Draw unseen sizes */
    for (i=0; i<Ki; i++){
      Nk0[i]=nk0[i];
      Nk1[i]=nk1[i];
    }
    for (i=ni; i<Ni; i++){
      /* Use rejection sampling */
      popi=0;
      while((popi == 0) || (log(1.0-unif_rand()) > -r*popi)){
        /* First propose unseen disease status for unit i */
        if(unif_rand() < pbeta){
          dis[i]=1;
        }else{
          dis[i]=0;
        }
        /* Now propose unseen size for unit i based on disease status */
        mu = exp(rnorm(musample[dis[i]], sigma));
        popi=rpois(mu);
      }
      if(popi > Ki){popi=Ki;}
      pop[i]=popi;
      if(dis[i]==1){
        Nk1[popi-1]=Nk1[popi-1]+1;
      }else{
        Nk0[popi-1]=Nk0[popi-1]+1;
      }
    }
    itotdis=itotdis0;
    for (i=ni; i<Ni; i++){
      itotdis+=dis[i];
    }
		    
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
//    if (*fVerbose) Rprintf("isamp %d pop[501] %d\n", isamp, pop[501]);
      sample[isamp*7  ]=(double)(Ni);
      sample[isamp*7+1]=mu0i;
      sample[isamp*7+2]=mu1i;
      sample[isamp*7+3]=sigma;
      sample[isamp*7+4]=(double)(Nk0[0]);
      sample[isamp*7+5]=beta;
      sample[isamp*7+6]=(double)(itotdis);
      for (i=0; i<Ki; i++){
        Nk0pos[i]=Nk0pos[i]+Nk0[i];
        Nk1pos[i]=Nk1pos[i]+Nk1[i];
      }
      isamp++;
      if (*fVerbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
//    if (*fVerbose) Rprintf("r %f gammart %f\n", r, gammart);
//    if (*fVerbose) Rprintf("Ni %d lpm[0] %f imaxN %d\n", Ni, lpm[0], imaxN);
    }
    step++;
  }
  for (i=0; i<Ki; i++){
    nk0[i]=Nk0pos[i];
    nk1[i]=Nk1pos[i];
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(p0i);
  free(p1i);
  free(b);
  free(Nk0);
  free(Nk0pos);
  free(Nk1);
  free(Nk1pos);
  free(lpm);
  free(musample);
  free(betasample);
  free(sigmasample);
}

void MHdis (int *Nk0, int *Nk1, int *totdis, int *K,
	    double *mu0, double *mu1, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, 
            double *musample, double *betasample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *fVerbose
			 ) {
  int step, taken, give_log=1;
  int i, Ki, Ni, isamp, iinterval, isamplesize, iburnin, itotdis;
  double ip, cutoff;
  double mu0star, mu1star, mu0i, mu1i, lp;
  double pbeta, betastar, betai;
  double p0is, p1is, p0stars, p1stars;
  double sigmastar, sigmai, sigma2star, sigma2i, qsigma2star, qsigma2i;
  double pithetastar, pithetai;
  double dkappa0, rkappa0, ddf0, dmu0, dmu1;
  double dsigma0, dsigma20, dmuproposal, dsigmaproposal;

  GetRNGstate();  /* R function enabling uniform RNG */

  Ki=(*K);
  double *p0star = (double *) malloc(sizeof(double) * Ki);
  double *p0i = (double *) malloc(sizeof(double) * Ki);
  double *p1star = (double *) malloc(sizeof(double) * Ki);
  double *p1i = (double *) malloc(sizeof(double) * Ki);

  Ni=(*N);
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburnin=(*burnin);
  dkappa0=(*kappa0);
  rkappa0=sqrt(dkappa0);
  ddf0=(*df0);
  dsigma0=(*sigma0);
  dsigma20=(dsigma0*dsigma0);
  itotdis=(*totdis);
  dmu0=(*mu0);
  dmu1=(*mu1);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  isamp = taken = 0;
  step = -iburnin;
  betai = betasample[0];
  mu0i = musample[0];
  mu1i = musample[1];
  sigmai = sigmasample[0];
  sigma2i  = sigmai*sigmai;

  pithetai = dnorm(mu0i, dmu0, sigmai/rkappa0, give_log);
  pithetai = pithetai+dnorm(mu1i, dmu1, sigmai/rkappa0, give_log);
  pithetai = pithetai+dsclinvchisq(sigma2i, ddf0, dsigma20);
  p0is=0.;
  p1is=0.;
  for (i=0; i<Ki; i++){
    p0i[i] = poilog(i+1,mu0i,sigmai);
    p1i[i] = poilog(i+1,mu1i,sigmai);
    p0is+=p0i[i];
    p1is+=p1i[i];
  }
  while (isamp < isamplesize) {
//  Rprintf("step %d Ni %d itotdis %d isamp %d\n", step, Ni, itotdis, isamp);
    /* Propose new theta */
    betastar = rnorm(betai, dmuproposal);
    mu0star = rnorm(mu0i, dmuproposal);
    mu1star = rnorm(mu1i, dmuproposal);
    sigma2star = sigma2i*exp(rnorm(0., dsigmaproposal));
    sigmastar = sqrt(sigma2star);

//  Rprintf("%f %f %f %f %f\n", mu0star, dmu0, dmu1, sigma2star, dkappa0, sigma2i);
    /* Calculate pieces of the posterior. */
    qsigma2star = dnorm(log(sigma2star/sigma2i)/dsigmaproposal,0.,1.,give_log)
                  -log(dsigmaproposal*sigma2star);
    pithetastar = dnorm(mu0star, dmu0, sigmastar/rkappa0, give_log);
    pithetastar = pithetastar+dnorm(mu1star, dmu1, sigmastar/rkappa0, give_log);
    pithetastar = pithetastar+dsclinvchisq(sigma2star, ddf0, dsigma20);
    qsigma2i = dnorm(log(sigma2i/sigma2star)/dsigmaproposal,0.,1.,give_log)
               -log(dsigmaproposal*sigma2i);

    /* Calculate ratio */
    ip = pithetastar-pithetai;
    /* Add the disease status */
    pbeta = exp(betastar)/(1.+exp(betastar));
//    Rprintf("pbeta %f betastar %f\n", pbeta, betastar);
    ip+=(itotdis*log(pbeta)+(Ni-itotdis)*log(1.-pbeta));
    pbeta = exp(betai)/(1.+exp(betai));
    ip-=(itotdis*log(pbeta)+(Ni-itotdis)*log(1.-pbeta));
//    Rprintf("pbeta %f betai %f betastar %f\n", pbeta, betai, betastar);
//    Rprintf("pithetastar %f pithetai %f qsigma2star %f qsigmai %f\n",
//           pithetastar,pithetai,qsigma2star,qsigma2i);
    p0stars=0.;
    p1stars=0.;
    for (i=0; i<Ki; i++){
      p0star[i] = poilog(i+1,mu0star,sigmastar);
      p0stars+=p0star[i];
      p1star[i] = poilog(i+1,mu1star,sigmastar);
      p1stars+=p1star[i];
    }
    for (i=0; i<Ki; i++){
//  Rprintf("%d %f %f\n", i,poilog(s[i],mu0star,sigmastar),
//    poilog(s[i],mu0i,sigmai));
      lp = log(p0is*p0star[i]/(p0stars*p0i[i]));
      if((lp > -100.) && (lp<100.)){ip += (Nk0[i]*lp);}
//    Rprintf("%f %f ", lp, ip);
      lp = log(p1is*p1star[i]/(p1stars*p1i[i]));
      if((lp > -100.) && (lp<100.)){ip += (Nk1[i]*lp);}
    }
//    Rprintf("%f %f\n", lp, ip);
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
    then let the MH probability equal min{exp(cutoff), 1.0}.
    But we'll do it in log space instead.  */
    cutoff = ip + qsigma2i-qsigma2star;
      
//    Rprintf("Now proposing %d MH steps %f cutoff...\n", step, cutoff);

    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed changes */
      betai = betastar;
      mu0i    = mu0star;
      mu1i    = mu1star;
      sigmai = sigmastar;
      sigma2i = sigma2star;
      qsigma2i = qsigma2star;
      pithetai = pithetastar;
      p0is=p0stars;
      p1is=p1stars;
      for (i=0; i<Ki; i++){
        p0i[i] = p0star[i];
        p1i[i] = p1star[i];
      }
      taken++;
//  if (*fVerbose)
//    Rprintf("Taken %d MH steps...\n", taken);
//    }
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
      sigmasample[isamp]=sigmai;
      musample[2*isamp]=mu0i;
      musample[2*isamp+1]=mu1i;
      betasample[isamp]=betai;
      isamp++;
      if (*fVerbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
    }
    }
    step++;
  }
  free(p0i);
  free(p0star);
  free(p1i);
  free(p1star);
  PutRNGstate();  /* Disable RNG before returning */
  *staken = taken;
}

void MHpriordis (double *mu0, double *mu1, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            double *musample, double *betasample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *fVerbose
			 ) {
  int step, taken, give_log=1;
  int isamp, iinterval, isamplesize, iburnin;
  double ip, cutoff;
  double mu0star, mu1star, mu0i, mu1i;
  double betastar, betai;
  double sigmastar, sigmai, sigma2star, sigma2i, qsigma2star, qsigma2i;
  double pithetastar, pithetai;
  double dkappa0, rkappa0, ddf0, dmu0, dmu1;
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
  dmu1=(*mu1);
  dmuproposal=(*muproposal);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);
  isamp = taken = 0;
  step = -iburnin;
  mu0i = dmu0;
  mu1i = dmu1;
  betai = -1.386294;
  sigma2i = dsigma20;
  sigmai  = sqrt(sigma2i);
  pithetai = dnorm(mu0i, dmu0, sigmai/rkappa0, give_log);
  pithetai = pithetai+dnorm(mu1i, dmu1, sigmai/rkappa0, give_log);
  pithetai = pithetai+dsclinvchisq(sigma2i, ddf0, dsigma20);
  while (isamp < isamplesize) {
    /* Propose new theta */
    betastar = rnorm(betai, dmuproposal);
    mu0star = rnorm(mu0i, dmuproposal);
    mu1star = rnorm(mu1i, dmuproposal);
    sigma2star = sigma2i*exp(rnorm(0., dsigmaproposal));
    sigmastar = sqrt(sigma2star);

    /* Calculate pieces of the posterior. */
    qsigma2star = dnorm(log(sigma2star/sigma2i)/dsigmaproposal,0.,1.,give_log)
                  -log(dsigmaproposal*sigma2star);
    pithetastar = dnorm(mu0star, dmu0, sigmastar/rkappa0, give_log);
    pithetastar = pithetastar+dnorm(mu1star, dmu1, sigmastar/rkappa0, give_log);
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
      /* Make proposed changes */
      betai = betastar;
      mu0i    = mu0star;
      mu1i    = mu1star;
      sigmai = sigmastar;
      sigma2i = sigma2star;
      pithetai = pithetastar;
//    if (*fVerbose) Rprintf("step %d mu0i %f sigmai %f\n", step, mu0i, sigmai);
      taken++;
      if (step > 0 && step==(iinterval*(step/iinterval))) { 
        /* record statistics for posterity */
        musample[2*isamp]=mu0i;
        musample[2*isamp+1]=mu1i;
        betasample[isamp]=betai;
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
