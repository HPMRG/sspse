/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "mcmc.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void MetropolisHastings (double *s,
		         double *mu0, double *kappa0, 
                         double *sigma20,  double *df0,
                         double *sigmaproposal, 
                         int *N, 
			 double *musample, double *sigmasample,
			 int *nsteps, int *staken, int *burnin, int *interval,
			 int fVerbose
			 ) {
  int step, taken;
  int i, Ni, isamp, iinterval, isteps, iburnin;
  double ip, cutoff;
  double mustar, mui, lp;
  double sigmastar, sigmai, sigma2star, sigma2i, qsigma2star, qsigma2i;
  double pithetastar, pithetai;
  double dkappa0, ddf0, dmu0, dsigma20, dsigmaproposal;

  GetRNGstate();  /* R function enabling uniform RNG */

  Ni=(*N);
  isteps=(*nsteps);
  iinterval=(*interval);
  iburnin=(*burnin);
  dkappa0=(*kappa0);
  ddf0=(*df0);
  dsigma20=(*sigma20);
  dmu0=(*mu0);
  dsigmaproposal=(*sigmaproposal);
  isamp = taken = 0;
  step = -iburnin;
/*  if (fVerbose)
    Rprintf("Now proposing %d MH steps... ", *nsteps); */
  mui = dmu0;
  sigma2i = dsigma20;
  while (isamp < isteps) {
    /* Propose new theta */
    mustar = rnorm(mui, dsigmaproposal);
    sigma2star = rsclinvchisq(ddf0,sigma2i);

//  Rprintf("%f %f %f %f %f\n", mustar, dmu0, sigma2star, dkappa0, sigma2i);
    /* Calculate pieces of the posterior. */
    qsigma2star = dsclinvchisq(sigma2star,ddf0, dsigma20,1);
    qsigma2i    = dsclinvchisq(sigma2i   ,ddf0, dsigma20,1);
    pithetastar = dnorm(mustar, dmu0, sigma2star/dkappa0,1);
    pithetai    = dnorm(mui   , dmu0, sigma2i/dkappa0,   1);

    sigmai    = sqrt(sigma2i   );
    sigmastar = sqrt(sigma2star);
    /* Calculate ratio */
    ip = pithetastar-pithetai;
  if (fVerbose)
//    Rprintf("Now proposing %d MH steps %f ip0...\n", step, ip);
    for (i=0; i<Ni; i++){
//  Rprintf("%d %f %f\n", i,poilog(s[i],mustar,sigmastar),
//    poilog(s[i],mui,sigmai));
      lp = log(poilog(s[i],mustar,sigmastar)/poilog(s[i],mui,sigmai));
      if(lp > -100000){ip += lp;}
//    Rprintf("%d %f\n", i, log(poilog(s[i],mustar,sigmastar)/poilog(s[i],mui,sigmai)));
    }
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
    then let the MH probability equal min{exp(cutoff), 1.0}.
    But we'll do it in log space instead.  */
//  if (fVerbose)
//    Rprintf("Now proposing %d MH steps %f ip1...\n", step, ip);
    cutoff = ip + qsigma2i-qsigma2star;
      
//  if (fVerbose)
//    Rprintf("Now proposing %d MH steps %f cutoff...\n", step, cutoff);

    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed chnages */
      sigmai = sigmastar;
      mui    = mustar;
      taken++;
//  if (fVerbose)
//    Rprintf("Taken %d MH steps...\n", taken);
//    }
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
      musample[isamp]=mui;
      sigmasample[isamp]=sigmai;
      isamp++;
      if (fVerbose) Rprintf("Taken %d MH steps...\n", isamp);
    }
    }
    step++;
  }
  PutRNGstate();  /* Disable RNG before returning */
  *staken = taken;
}
