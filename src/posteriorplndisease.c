/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "posteriorplndisease.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gplndisease (int *pop, int *dis,
            int *nk0, int *nk1, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu0, double *mu1, double *kappa0, 
            double *sigma0,  double *sigma1, double *df0,
	    int *Np0, int *Np1,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *maxN, 
            double *sample, 
            double *p0pos, double *p1pos, 
            double *ppos, 
	    double *lpriorm,
            int *burnintheta,
            int *verbose
              ) {
  int dimsample, Nnp0, Nnp1;
  int step, staken, getone=1, intervalone=1, verboseMHdisease = 0;
  int i, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  double mu0i, mu1i, pbeta, beta, sigma0i, sigma1i, dsamp;
  double dkappa0, ddf0, dmu0, dmu1, dsigma0, dsigma1, dmuproposal, dsigmaproposal;
  int tU, popi, imaxN, imaxm, itotdis0, itotdis;
  double r, gamma0rt, gamma1rt, p0is, p1is, Nd;

  GetRNGstate();  /* R function enabling uniform RNG */

  ni=(*n);
  Ni=(*N);
  Ki=(*K);
  imaxN=(*maxN);
  imaxm=imaxN-ni;
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburnin=(*burnin);
  Nnp0=(*Np0);
  Nnp1=(*Np1);
  dkappa0=(*kappa0);
  ddf0=(*df0);
  dsigma0=(*sigma0);
  dsigma1=(*sigma1);
  dmu0=(*mu0);
  dmu1=(*mu1);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  dimsample=8+Nnp0+Nnp1;

  double *p0i = (double *) malloc(sizeof(double) * Ki);
  double *p1i = (double *) malloc(sizeof(double) * Ki);
  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk0 = (int *) malloc(sizeof(int) * Ki);
  int *Nk0pos = (int *) malloc(sizeof(int) * Ki);
  int *Nk1 = (int *) malloc(sizeof(int) * Ki);
  int *Nk1pos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdeg0i = (double *) malloc(sizeof(double) * Nnp0);
  double *pdeg1i = (double *) malloc(sizeof(double) * Nnp1);
  double *psample = (double *) malloc(sizeof(double) * (Nnp0+Nnp1));
  double *musample = (double *) malloc(sizeof(double) * 2);
  double *betasample = (double *) malloc(sizeof(double));
  double *sigmasample = (double *) malloc(sizeof(double) * 2);

  b[ni-1]=pop[ni-1];
  for (i=(ni-2); i>=0; i--){
    b[i]=b[i+1]+pop[i];
  }
  for (i=0; i<Ki; i++){
     Nk0[i]=nk0[i];
     Nk0pos[i]=0;
     Nk1[i]=nk1[i];
     Nk1pos[i]=0;
     p0pos[i]=0.;
     p1pos[i]=0.;
     ppos[i]=0.;
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
  for (i=0; i<Nnp0; i++){
     psample[i    ] = 0.01;
  }
  for (i=0; i<Nnp1; i++){
     psample[i+Nnp0]=0.01;
  }
  musample[0] = dmu0;
  musample[1] = dmu1;
  sigmasample[0] = dsigma0;
  sigmasample[1] = dsigma1;

  isamp = 0;
  step = -iburnin;
  while (isamp < isamplesize) {
    /* Draw new theta */
    MHplndisease(Nk0,Nk1,&itotdis,K,mu0,mu1,kappa0,sigma0,sigma1,df0,
          muproposal, sigmaproposal,
	  &Ni, &Nnp0, &Nnp1, psample, 
	  musample, betasample, sigmasample, &getone, &staken, 
	  burnintheta, &intervalone, 
	  &verboseMHdisease);

    beta=betasample[0];
    pbeta=exp(beta)/(1.+exp(beta));
    for (i=0; i<Nnp0; i++){
      pdeg0i[i] = psample[i];
    }
    for (i=0; i<Nnp1; i++){
      pdeg1i[i] = psample[i+Nnp0];
    }
    mu0i=musample[0];
    mu1i=musample[1];
    sigma0i=sigmasample[0];
    sigma1i=sigmasample[1];

    /* Draw new N */

    /* First find the degree distribution */
    p0is=0.;
    p1is=0.;
    for (i=Nnp0; i<Ki; i++){
      p0i[i]=poilog(i+1,mu0i,sigma0i);
    }
    p0is=1.-poilog(0,mu0i,sigma0i);
    for (i=Nnp1; i<Ki; i++){
      p1i[i]=poilog(i+1,mu1i,sigma1i);
//    p1is+=p1i[i];
    }
    p1is=1.-poilog(0,mu1i,sigma1i);
    for (i=0; i<Ki; i++){
      p0i[i]=p0i[i]/p0is;
      p1i[i]=p1i[i]/p1is;
    }
    p0is=1.;
    p1is=1.;
    for (i=0; i<Nnp0; i++){
      p0is-=pdeg0i[i];
    }
    for (i=0; i<Nnp1; i++){
      p1is-=pdeg1i[i];
    }
    for (i=0; i<Ki; i++){
      p0i[i]=p0i[i]*p0is;
      p1i[i]=p1i[i]*p1is;
    }
    for (i=0; i<Nnp0; i++){
      p0i[i]=pdeg0i[i];
    }
    for (i=0; i<Nnp1; i++){
      p1i[i]=pdeg1i[i];
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
    // N = m + n
    // Compute (log) P(m | \theta and data and \Psi)
    for (i=0; i<imaxm; i++){
     lpm[i]=lgamma(ni+i+1.)-lgamma(i+1.)+i*gamma0rt;
     //  Add in the (log) prior on m: P(m)
     lpm[i]=lpm[i]+lpriorm[i];
     if(lpm[i] > tU) tU = (int)lpm[i];
    }
    for (i=0; i<imaxm; i++){
      lpm[i]=exp(lpm[i]-tU);
    }
    for (i=1; i<imaxm; i++){
      lpm[i]=lpm[i-1]+lpm[i];
    }
    gamma0rt = lpm[imaxm-1] * unif_rand();
    for (Ni=0; Ni<imaxm; Ni++){
      if(gamma0rt <= lpm[Ni]) break;
    }
//  if (*verbose) Rprintf("Ni %d lpm[imaxm-1] %f lpm[Ni] %f\n", Ni, lpm[imaxm-1],
//  lpm[Ni]);
//  }
    Ni += ni;
    if(Ni > imaxN) Ni = imaxN;

//  if (*verbose) Rprintf("step %d Ni %d itotdis %d beta %f mu0 %f mu1 %f s0 %f r %f\n",
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
    // Set up pi for random draws
    for (i=1; i<Ki; i++){
      p0i[i]=p0i[i-1]+p0i[i];
      p1i[i]=p1i[i-1]+p1i[i];
    }
    for (i=ni; i<Ni; i++){
      /* Propose unseen size for unit i */
      /* Use rejection sampling */
      popi=1000000;
      while(popi >= Ki){
       popi=1000000;
       while((log(1.0-unif_rand()) > -r*popi)){
        /* First propose unseen disease status for unit i */
        if(unif_rand() < pbeta){
          dis[i]=1;
          /* Now propose unseen size for unit i based on disease status */
          /* In the next two lines a popi is chosen */
          /* with parameters mu1i and sigma1i */
          gamma1rt = p1i[Ki-1] * unif_rand();
          for (popi=1; popi<=Ki; popi++){
            if(gamma1rt <= p1i[popi-1]) break;
          }
        }else{
          dis[i]=0;
          /* Now propose unseen size for unit i based on non-disease status */
          /* In the next two lines a popi is chosen */
          /* with parameters mu0i and sigma0i */
          gamma0rt = p0i[Ki-1] * unif_rand();
          for (popi=1; popi<=Ki; popi++){
            if(gamma0rt <= p0i[popi-1]) break;
          }
        }
       }
      }
//    if(popi >= Ki){popi=Ki-1;}
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
//    if (*verbose) Rprintf("isamp %d pop[501] %d\n", isamp, pop[501]);
      Nd=(double)Ni;
      sample[isamp*dimsample  ]=Nd;
      sample[isamp*dimsample+1]=mu0i;
      sample[isamp*dimsample+2]=mu1i;
      sample[isamp*dimsample+3]=sigma0i;
      sample[isamp*dimsample+4]=sigma1i;
      sample[isamp*dimsample+5]=(double)(Nk0[0]+Nk1[0]);
      sample[isamp*dimsample+6]=beta;
      sample[isamp*dimsample+7]=(double)(itotdis);
      for (i=0; i<Nnp0; i++){
        sample[isamp*dimsample+8+i]=pdeg0i[i];
      }
      for (i=0; i<Nnp1; i++){
        sample[isamp*dimsample+8+Nnp0+i]=pdeg1i[i];
      }
//    N0d=0.;
      for (i=0; i<Ki; i++){
        Nk0pos[i]=Nk0pos[i]+Nk0[i];
        Nk1pos[i]=Nk1pos[i]+Nk1[i];
        p0pos[i]+=(Nk0[i]/Nd);
        p1pos[i]+=(Nk1[i]/Nd);
        ppos[i]+=(Nk0[i]+Nk1[i])/Nd;
//
//      N0d+=Nk0[i];
//      N1d=Ni-N0d;
      }
      isamp++;
      if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
//    if (*verbose) Rprintf("r %f gammart %f\n", r, gammart);
//    if (*verbose) Rprintf("Ni %d lpm[0] %f imaxm %d\n", Ni, lpm[0], imaxm);
    }
    step++;
  }
  dsamp=((double)isamp);
  for (i=0; i<Ki; i++){
    nk0[i]=Nk0pos[i];
    nk1[i]=Nk1pos[i];
    p0pos[i]=p0pos[i]/dsamp;
    p1pos[i]=p1pos[i]/dsamp;
    ppos[i]=ppos[i]/dsamp;
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(psample);
  free(pdeg0i);
  free(pdeg1i);
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

void MHplndisease (int *Nk0, int *Nk1, int *totdis, int *K,
	    double *mu0, double *mu1, double *kappa0, 
            double *sigma0,  double *sigma1, double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Np0, int *Np1, double *psample,
            double *musample, double *betasample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *verbose
			 ) {
  int Nnp0, Nnp1;
  int step, taken, give_log1=1;
  int i, Ki, Ni, isamp, iinterval, isamplesize, iburnin, itotdis;
  double ip, cutoff;
  double mu0star, mu1star, mu0i, mu1i, lp;
  double pbeta, betastar, betai;
  double p0is, p1is, p0stars, p1stars;
  double sigma0star, sigma1star, sigma0i, sigma1i;
  double sigma02star, sigma12star, sigma02i, sigma12i;
  double qsigma02star, qsigma12star, qsigma02i, qsigma12i;
  double pithetastar, pithetai;
  double dkappa0, rkappa0, ddf0, dmu0, dmu1;
  double dsigma0, dsigma1, dsigma02, dsigma12, dmuproposal, dsigmaproposal;

  GetRNGstate();  /* R function enabling uniform RNG */

  Ki=(*K);
  Nnp0=(*Np0);
  Nnp1=(*Np1);
  double *p0star = (double *) malloc(sizeof(double) * Ki);
  double *p0i = (double *) malloc(sizeof(double) * Ki);
  double *p1star = (double *) malloc(sizeof(double) * Ki);
  double *p1i = (double *) malloc(sizeof(double) * Ki);
  double *odeg0star = (double *) malloc(sizeof(double) * Nnp0);
  double *odeg0i = (double *) malloc(sizeof(double) * Nnp0);
  double *odeg1star = (double *) malloc(sizeof(double) * Nnp1);
  double *odeg1i = (double *) malloc(sizeof(double) * Nnp1);
  double *pdeg0star = (double *) malloc(sizeof(double) * Nnp0);
  double *pdeg0i = (double *) malloc(sizeof(double) * Nnp0);
  double *pdeg1star = (double *) malloc(sizeof(double) * Nnp1);
  double *pdeg1i = (double *) malloc(sizeof(double) * Nnp1);

  Ni=(*N);
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iburnin=(*burnin);
  dkappa0=(*kappa0);
  rkappa0=sqrt(dkappa0);
  ddf0=(*df0);
  dsigma0=(*sigma0);
  dsigma1=(*sigma1);
  dsigma02=(dsigma0*dsigma0);
  dsigma12=(dsigma1*dsigma1);
  itotdis=(*totdis);
  dmu0=(*mu0);
  dmu1=(*mu1);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  // First set starting values
  isamp = taken = 0;
  step = -iburnin;
  betai = betasample[0];
  p0is=1.;
  p1is=1.;
  for (i=0; i<Nnp0; i++){
    pdeg0i[i] = psample[i];
    p0is-= pdeg0i[i];
  }
  for (i=0; i<Nnp1; i++){
    pdeg1i[i] = psample[i+Nnp0];
    p1is-= pdeg1i[i];
  }
  for (i=0; i<Nnp0; i++){
    odeg0i[i] = log(pdeg0i[i]/p0is);
  }
  for (i=0; i<Nnp1; i++){
    odeg1i[i] = log(pdeg1i[i]/p1is);
  }
  mu0i = musample[0];
  mu1i = musample[1];
  sigma0i = sigmasample[0];
  sigma1i = sigmasample[1];
  sigma02i = sigma0i*sigma0i;
  sigma12i = sigma1i*sigma1i;

  pithetai = dnorm(mu0i, dmu0, sigma0i/rkappa0, give_log1);
  pithetai = pithetai+dnorm(mu1i, dmu1, sigma1i/rkappa0, give_log1);
  pithetai = pithetai+dsclinvchisq(sigma02i, ddf0, dsigma02);
  pithetai = pithetai+dsclinvchisq(sigma12i, ddf0, dsigma12);
  p0is=0.;
  p1is=0.;
  for (i=Nnp0; i<Ki; i++){
    p0i[i]=poilog(i+1,mu0i,sigma0i);
//  p0is+=p0i[i];
  }
  p0is=1.-poilog(0,mu0i,sigma0i);
  for (i=Nnp1; i<Ki; i++){
    p1i[i]=poilog(i+1,mu1i,sigma1i);
//  p1is+=p1i[i];
  }
  p1is=1.-poilog(0,mu1i,sigma1i);
  for (i=0; i<Ki; i++){
    p0i[i]=p0i[i]/p0is;
    p1i[i]=p1i[i]/p1is;
  }
  p0is=1.;
  p1is=1.;
  for (i=0; i<Nnp0; i++){
    p0is-=pdeg0i[i];
  }
  for (i=0; i<Nnp1; i++){
    p1is-=pdeg1i[i];
  }
  for (i=0; i<Ki; i++){
    p0i[i]=p0i[i]*p0is;
    p1i[i]=p1i[i]*p1is;
  }
  for (i=0; i<Nnp0; i++){
    p0i[i]=pdeg0i[i];
  }
  for (i=0; i<Nnp1; i++){
    p1i[i]=pdeg1i[i];
  }

  // Now do the MCMC updates (starting with the burnin updates)
  while (isamp < isamplesize) {
//  Rprintf("step %d Ni %d itotdis %d isamp %d\n", step, Ni, itotdis, isamp);
    /* Propose new theta */
    /* Start with the disease status parameters */
    betastar = rnorm(betai, dmuproposal);
    pbeta = exp(betastar)/(1.+exp(betastar));
    /* Now the degree distribution model parameters */
    for (i=0; i<Nnp0; i++){
      odeg0star[i] = rnorm(odeg0i[i], dmuproposal);
    }
    for (i=0; i<Nnp1; i++){
      odeg1star[i] = rnorm(odeg1i[i], dmuproposal);
    }
    /* Convert from odds to probabilities */
    p0is=1.;
    p1is=1.;
    for (i=0; i<Nnp0; i++){
      pdeg0star[i] = exp(odeg0star[i]);
      p0is+=pdeg0star[i];
    }
    for (i=0; i<Nnp1; i++){
      pdeg1star[i] = exp(odeg1star[i]);
      p1is+=pdeg1star[i];
    }
    for (i=0; i<Nnp0; i++){
      pdeg0star[i]/=p0is;
    }
    for (i=0; i<Nnp1; i++){
      pdeg1star[i]/=p1is;
    }
    /* Now the degree distribution (log) mean adn s.d. parameters */
    mu0star = rnorm(mu0i, dmuproposal);
    mu1star = rnorm(mu1i, dmuproposal);
    sigma02star = sigma02i*exp(rnorm(0., dsigmaproposal));
    sigma0star = sqrt(sigma02star);
    sigma12star = sigma12i*exp(rnorm(0., dsigmaproposal));
    sigma1star = sqrt(sigma12star);

//  Rprintf("%f %f %f %f %f\n", mu0star, dmu0, dmu1, sigma2star, dkappa0, sigma2i);
    /* Calculate pieces of the posterior. */
    qsigma02star = dnorm(log(sigma02star/sigma02i)/dsigmaproposal,0.,1.,give_log1)
                  -log(dsigmaproposal*sigma02star);
    qsigma12star = dnorm(log(sigma12star/sigma12i)/dsigmaproposal,0.,1.,give_log1)
                  -log(dsigmaproposal*sigma12star);
    pithetastar = dnorm(mu0star, dmu0, sigma0star/rkappa0, give_log1);
    pithetastar = pithetastar+dnorm(mu1star, dmu1, sigma1star/rkappa0, give_log1);
    pithetastar = pithetastar+dsclinvchisq(sigma02star, ddf0, dsigma02);
    pithetastar = pithetastar+dsclinvchisq(sigma12star, ddf0, dsigma12);
    qsigma02i = dnorm(log(sigma02i/sigma02star)/dsigmaproposal,0.,1.,give_log1)
               -log(dsigmaproposal*sigma02i);
    qsigma12i = dnorm(log(sigma12i/sigma12star)/dsigmaproposal,0.,1.,give_log1)
               -log(dsigmaproposal*sigma12i);

    /* Calculate ratio */
    ip = pithetastar-pithetai;
    /* Add the disease status */
//    Rprintf("pbeta %f betastar %f\n", pbeta, betastar);
    ip+=(itotdis*log(pbeta)+(Ni-itotdis)*log(1.-pbeta));
    pbeta = exp(betai)/(1.+exp(betai));
    ip-=(itotdis*log(pbeta)+(Ni-itotdis)*log(1.-pbeta));
//    Rprintf("pbeta %f betai %f betastar %f\n", pbeta, betai, betastar);
//    Rprintf("pithetastar %f pithetai %f qsigma2star %f qsigmai %f\n",
//           pithetastar,pithetai,qsigma2star,qsigma2i);
    p0stars=0.;
    p1stars=0.;
    for (i=Nnp0; i<Ki; i++){
      p0star[i]=poilog(i+1,mu0star,sigma0star);
//    p0stars+=p0star[i];
    }
    p0stars=1.-poilog(0,mu0star,sigma0star);
    for (i=Nnp1; i<Ki; i++){
      p1star[i]=poilog(i+1,mu1star,sigma1star);
//    p1stars+=p1star[i];
    }
    p1stars=1.-poilog(0,mu1star,sigma1star);
    for (i=Nnp0; i<Ki; i++){
      p0star[i]/=p0stars;
    }
    for (i=Nnp1; i<Ki; i++){
      p1star[i]/=p1stars;
    }
    p0stars=1.;
    p1stars=1.;
    for (i=0; i<Nnp0; i++){
      p0stars-=pdeg0star[i];
    }
    for (i=0; i<Nnp1; i++){
      p1stars-=pdeg1star[i];
    }
    for (i=Nnp0; i<Ki; i++){
      p0star[i]*=p0stars;
    }
    for (i=Nnp1; i<Ki; i++){
      p1star[i]*=p1stars;
    }
    for (i=0; i<Nnp0; i++){
      p0star[i]=pdeg0star[i];
    }
    for (i=0; i<Nnp1; i++){
      p1star[i]=pdeg1star[i];
    }
    for (i=0; i<Ki; i++){
     if(Nk0[i]>0){
      lp = log(p0star[i]/p0i[i]);
      if(fabs(lp) < 100.){ip += (Nk0[i]*lp);}
     }
     if(Nk1[i]>0){
      lp = log(p1star[i]/p1i[i]);
      if(fabs(lp) < 100.){ip += (Nk1[i]*lp);}
     }
    }
//    Rprintf("%f %f\n", lp, ip);
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
    then let the MH probability equal min{exp(cutoff), 1.0}.
    But we'll do it in log space instead.  */
    cutoff = ip + qsigma02i-qsigma02star + qsigma12i-qsigma12star;
      
//    Rprintf("Now proposing %d MH steps %f cutoff...\n", step, cutoff);

    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed changes */
      betai = betastar;
      for (i=0; i<Nnp0; i++){
        odeg0i[i] = odeg0star[i];
        pdeg0i[i] = pdeg0star[i];
      }
      for (i=0; i<Nnp1; i++){
        odeg1i[i] = odeg1star[i];
        pdeg1i[i] = pdeg1star[i];
      }
      mu0i    = mu0star;
      mu1i    = mu1star;
      sigma0i = sigma0star;
      sigma1i = sigma1star;
      sigma02i = sigma02star;
      sigma12i = sigma12star;
      qsigma02i = qsigma02star;
      qsigma12i = qsigma12star;
      pithetai = pithetastar;
      for (i=0; i<Ki; i++){
        p0i[i] = p0star[i];
        p1i[i] = p1star[i];
      }
      taken++;
//  if (*verbose)
//    Rprintf("Taken %d MH steps...\n", taken);
//    }
      if (step > 0 && step==(iinterval*(step/iinterval))) { 
        /* record statistics for posterity */
        sigmasample[2*isamp]=sigma0i;
        sigmasample[2*isamp+1]=sigma1i;
        musample[2*isamp]=mu0i;
        musample[2*isamp+1]=mu1i;
        betasample[isamp]=betai;
        for (i=0; i<Nnp0; i++){
          psample[i    ]=pdeg0i[i];
        }
        for (i=0; i<Nnp1; i++){
          psample[i+Nnp0]=pdeg1i[i];
        }
        isamp++;
        if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d samples...\n", isamp);
      }
    }
    step++;
  }
  free(p0i);
  free(p0star);
  free(p1i);
  free(p1star);
  free(odeg0i);
  free(odeg0star);
  free(odeg1i);
  free(odeg1star);
  free(pdeg0i);
  free(pdeg0star);
  free(pdeg1i);
  free(pdeg1star);
  PutRNGstate();  /* Disable RNG before returning */
  *staken = taken;
}
