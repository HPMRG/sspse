/*******************************************************************/
/* Computation of the log-likelihood and marginal likelihood of size*/
/*******************************************************************/

#include "likelihoodcmp.h"
#include "cmp.h"
//#include "cmp179.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void lcmp (int *pop,
            int *nk, 
            int *K, 
            int *n, 
            int *samplesize, int *warmup, int *interval,
            double *mu, double *kappa, 
            double *sigma,  double *df,
	    int *Npi,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *maxN, 
            double *sample, 
            double *ppos, 
            double *lpriorm, 
            int *warmuptheta,
	    int *verbose
			 ) {
  int dimsample, Np;
  int step, staken, getone=1, intervalone=1, verboseMHcmp = 0;
  int i, ni, Ni, Ki, isamp, iinterval, isamplesize, iwarmup;
  double mui, sigmai, dsamp;
  double dkappa, ddf, dmu, dsigma, dmuproposal, dsigmaproposal;
  int tU, sizei, imaxN, imaxm, give_log0=0, give_log1=1;
  double r, gammart, pis, Nd;
  double temp;
  double errval=0.000000001, lzcmp;

  GetRNGstate();  /* R function enabling uniform RNG */

  ni=(*n);
  Ni=(*N);
  Ki=(*K);
  imaxN=(*maxN);
  imaxm=imaxN-ni;
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iwarmup=(*warmup);
  Np=(*Npi);
  dkappa=(*kappa);
  ddf=(*df);
  dsigma=(*sigma);
  dmu=(*mu);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  dimsample=4+Np;

  double *pi = (double *) malloc(sizeof(double) * Ki);
  double *pd = (double *) malloc(sizeof(double) * Ki);
  int *d = (int *) malloc(sizeof(int) * imaxN);
  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  int *Nkpos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdegi = (double *) malloc(sizeof(double) * (Np+1));
  double *psample = (double *) malloc(sizeof(double) * (Np+1));
  double *musample = (double *) malloc(sizeof(double) * isamplesize);
  double *sigmasample = (double *) malloc(sizeof(double) * isamplesize);

  for (i=0; i<ni; i++){
    if((pop[i]>0) && (pop[i] <= Ki)){ d[i]=pop[i];}
    if(pop[i]==0){ d[i]=1;}
    if(pop[i]>Ki){ d[i]=Ki;}
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
    tU+=pop[i];
  }
  /* Draw initial phis */
  r=0.;
  for (i=0; i<ni; i++){
    r+=(exp_rand()/(tU+b[i]));
  }

  for (i=0; i<Np; i++){
     psample[i] = 0.01;
  }
  musample[0] = dmu;
  sigmasample[0] = dsigma;

  isamp = 0;
  step = -iwarmup;
  while (isamp < isamplesize) {
    /* Draw new theta */
    /* but less often than the other full conditionals */
    if (step == -iwarmup || step==(10*(step/10))) { 
     MHlcmp(Nk,K,mu,kappa,sigma,df,muproposal,sigmaproposal,
           &Ni, &Np, psample,
	   musample, sigmasample, &getone, &staken, warmuptheta, &intervalone, 
	   &verboseMHcmp);
    }

    for (i=0; i<Np; i++){
      pdegi[i] = psample[i];
    }
    mui=musample[0];
    sigmai=sigmasample[0];

    /* Draw new N */

    /* First find the degree distribution */
    pis=0.;
    lzcmp = zcmp(exp(mui), sigmai, errval, Ki, give_log1);
    if(lzcmp < -100000.0){continue;}
    pi[Np]=cmp(Np+1,mui,sigmai,lzcmp,give_log0);
    for (i=Np+1; i<Ki; i++){
      pi[i]=pi[i-1]*exp(mui-sigmai*log((double)(i+1)));
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
        /* with parameters mui and sigmai */
	temp = unif_rand();
//      gammart = pi[Ki-1] * unif_rand();
        for (sizei=1; sizei<=Ki; sizei++){
          if(temp <= pi[sizei-1]) break;
        }
//      Rprintf("sizei %d pi[Ki-1] %f gammart %f\n", sizei, pi[Ki-1],gammart);
       }
      }
//    if(sizei >= Ki){sizei=Ki-1;}
      pop[i]=sizei;
//    if((sizei <= 0) | (sizei > Ki-1)) Rprintf("sizei %d r %f\n", sizei,r);
      Nk[sizei-1]=Nk[sizei-1]+1;
    }
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
      Nd=(double)Ni;
      sample[isamp*dimsample  ]=Nd;
      sample[isamp*dimsample+1]=mui;
      sample[isamp*dimsample+2]=sigmai;
      sample[isamp*dimsample+3]=(double)(Nk[0]);
      for (i=0; i<Np; i++){
        sample[isamp*dimsample+4+i]=pdegi[i];
      }
      for (i=0; i<Ki; i++){
        Nkpos[i]=Nkpos[i]+Nk[i];
        ppos[i]+=((Nk[i]*1.)/Nd);
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
  PutRNGstate();  /* Disable RNG before returning */
  free(pi);
  free(pd);
  free(psample);
  free(pdegi);
  free(b);
  free(d);
  free(Nk);
  free(Nkpos);
  free(lpm);
  free(musample);
  free(sigmasample);
}

void MHlcmp (int *Nk, int *K,
	    double *mu, double *kappa, 
            double *sigma,  double *df,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Npi, double *psample,
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *warmup, int *interval,
	    int *verbose
			 ) {
  int Np;
  int step, taken, give_log1=1, give_log0=0;
  int i, Ki, Ni, isamp, iinterval, isamplesize, iwarmup;
  double ip, cutoff;
  double mustar, mui, lp;
  double pis, pstars;
  double sigmastar, sigmai, sigma2star, sigma2i, qsigma2star, qsigma2i;
//double pithetastar, pithetai;
  double dkappa, rkappa, ddf, dmu;
  double dsigma, dsigma2, dmuproposal, dsigmaproposal;
  double errval=0.000000001, lzcmp;

  GetRNGstate();  /* R function enabling uniform RNG */

  Ki=(*K);
  Np=(*Npi);
  double *pstar = (double *) malloc(sizeof(double) * Ki);
  double *pi = (double *) malloc(sizeof(double) * Ki);
  double *odegstar = (double *) malloc(sizeof(double) * Np);
  double *odegi = (double *) malloc(sizeof(double) * Np);
  double *pdegstar = (double *) malloc(sizeof(double) * Np);
  double *pdegi = (double *) malloc(sizeof(double) * Np);

  Ni=(*N);
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iwarmup=(*warmup);
  dkappa=(*kappa);
  rkappa=sqrt(dkappa);
  ddf=(*df);
  dsigma=(*sigma);
  dsigma2=(dsigma*dsigma);
  dmu=(*mu);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  // First set starting values
  isamp = taken = 0;
  step = -iwarmup;
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
//pithetai = dnorm(mui, dmu, sigmai/rkappa, give_log1);
//pithetai = pithetai+dsclinvchisq(sigma2i, ddf, dsigma2);
//    Rprintf("mui %f sigmai %f lzcmp %f\n", mui, sigmai, lzcmp);
  pis=0.;
  lzcmp = zcmp(exp(mui), sigmai, errval, Ki, give_log1);
  pi[Np]=cmp(Np+1,mui,sigmai,lzcmp,give_log0);
  for (i=Np+1; i<Ki; i++){
    pi[i]=pi[i-1]*exp(mui-sigmai*log((double)(i+1)));
  }
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
    pi[i]=pi[i]*pis;
  }
  for (i=0; i<Np; i++){
    pi[i]=pdegi[i];
  }
  // Bin last group
  pis=1.;
  for (i=0; i<(Ki-1); i++){
    pis-=pi[i];
  }
  pi[Ki-1]=pis;

  // Now do the MCMC updates (starting with the warmup updates)
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

//  if(sigma2star > 1000) Rprintf("%f %f %f %f %f\n", mustar, dmu, sigma2star, dkappa, sigma2i);
    /* Calculate pieces of the likelihood. */
    qsigma2star = dnorm(log(sigma2star/sigma2i)/dsigmaproposal,0.,1.,give_log1)
                  -log(dsigmaproposal*sigma2star);
//  pithetastar = dnorm(mustar, dmu, sigmastar/rkappa, give_log1);
//  pithetastar = pithetastar+dsclinvchisq(sigma2star, ddf, dsigma2);
    qsigma2i = dnorm(log(sigma2i/sigma2star)/dsigmaproposal,0.,1.,give_log1)
               -log(dsigmaproposal*sigma2i);

    /* Calculate ratio */
//  ip = pithetastar-pithetai;
    ip = 0.0;
//    Rprintf("mustar %f sigmastar %f lzcmp %f\n", mustar, sigmastar, lzcmp);
    pstars=0.;
    lzcmp = zcmp(exp(mustar), sigmastar, errval, Ki, give_log1);
    pstar[Np]=cmp(Np+1,mustar,sigmastar,lzcmp,give_log0);
    for (i=Np+1; i<Ki; i++){
      pstar[i]=pstar[i-1]*exp(mustar-sigmastar*log((double)(i+1)));
    }
    pstars=1.-exp(-lzcmp);
    for (i=0; i<Ki; i++){
      pstar[i]/=pstars;
//Rprintf("i %d pstar %f pi0 %f\n", i, pstar[i], pi0[i], pis, pis0);
    }
    pstars=1.;
    for (i=0; i<Np; i++){
      pstars-=pdegstar[i];
    }
    for (i=Np; i<Ki; i++){
      pstar[i]*=pstars;
    }
    for (i=0; i<Np; i++){
      pstar[i]=pdegstar[i];
    }
    // Bin last group
    pstars=1.;
    for (i=0; i<(Ki-1); i++){
      pstars-=pstar[i];
    }
    pstar[Ki-1]=pstars;

    for (i=0; i<Ki; i++){
     if(Nk[i]>0){
      lp = log(pstar[i]/pi[i]);
      if(fabs(lp) < 100.){ip += (Nk[i]*lp);}
     }
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
//    pithetai = pithetastar;
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
  free(odegi);
  free(odegstar);
  free(pdegi);
  free(pdegstar);
  PutRNGstate();  /* Disable RNG before returning */
  *staken = taken;
}
