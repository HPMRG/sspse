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
            int *samplesize, int *burnin, int *interval,
            double *mu, double *kappa, 
            double *sigma,  double *df,
	    int *Npi,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *maxN, 
            double *sample, 
            double *ppos, 
            double *lpriorm, 
            int *burnintheta,
	    double *lambdad,
	    double *nud,
	    int *verbose
			 ) {
  int dimsample, Np;
  int step, staken, getone=1, intervalone=1, verboseMHcmp = 0;
  int i, j, compute, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  double mui, sigmai, dsamp;
  double dkappa, ddf, dmu, dsigma, dmuproposal, dsigmaproposal;
  int tU, sizei, imaxN, imaxm, give_log0=0, give_log1=1;
  int maxpop;
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
  iburnin=(*burnin);
  Np=(*Npi);
  dkappa=(*kappa);
  ddf=(*df);
  dsigma=(*sigma);
  dmu=(*mu);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  dimsample=4+Np;

  double *pi = (double *) malloc(sizeof(double) * Ki);
  double *pd = (double *) malloc(sizeof(double) * ni);
  int *d = (int *) malloc(sizeof(int) * ni);
  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  int *Nkpos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdegi = (double *) malloc(sizeof(double) * (Np+1));
  double *psample = (double *) malloc(sizeof(double) * (Np+1));
  double *musample = (double *) malloc(sizeof(double));
  double *sigmasample = (double *) malloc(sizeof(double));

  maxpop=0;
  for (i=0; i<ni; i++){
    if((pop[i]>0) && (pop[i] <= Ki)){ d[i]=pop[i];}
    if(pop[i]==0){ d[i]=1;}
    if(pop[i]>Ki){ d[i]=Ki;}
    if(pop[i]>maxpop){maxpop=pop[i];}
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
  step = -iburnin;
  while (isamp < isamplesize) {
    /* Draw new theta */
    /* but less often than the other full conditionals */
    if (step == -iburnin || step==(10*(step/10))) { 
     MHlcmp(Nk,K,mu,kappa,sigma,df,muproposal,sigmaproposal,
           &Ni, &Np, psample,
	   musample, sigmasample, &getone, &staken, burnintheta, &intervalone, 
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
    Ni = 0;
    while(temp > lpm[Ni]){Ni++;}
    // Add back the sample size
    Ni += ni;
    if(Ni > imaxN) Ni = imaxN;
		    
    if((fabs(lambdad[0])>0.0000001) | (fabs(nud[0])>0.0000001)){
Rprintf("No! lambdad[0] %f nud[0] %f\n", lambdad[0], nud[0]);
    for (i=0; i<Ki; i++){
      nk[i]=0;
    }

    /* Draw true degrees (sizes) based on the reported degrees*/
    /* First find the reported degree distribution */
    for (j=0; j<=maxpop; j++){
//Rprintf("j %d pop[j] %d\n", j, pop[j]);
     compute=0;
     for (i=0; i<ni; i++){if(pop[i]==(j)){compute=1;}}
     if(compute==1){
//    Next four lines for cmp reporting distribution
//    ?? Should it be cmp(j+1,...) or cmp(j,...)??
      for (i=0; i<Ki; i++){
       lzcmp = zcmp(exp(lambdad[i]),nud[i], errval, Ki, give_log1);
//     pd[i]=pi[i]*cmp(j,lambdad[i],nud[i],lzcmp,give_log0);
       pd[i]=pi[i]*cmp(j+1,lambdad[i],nud[i],lzcmp,give_log0);
      }
//     Next seven lines for proportional reporting distribution
//       for (i=0; i<Ki; i++){
//        pd[i]   = pgamma(2.0*lambdad[i]/((j)+0.5),1.0,1.0,0,0);
//        if((j)>0){
//         pd[i] -= pgamma(2.0*lambdad[i]/((j)-0.5),1.0,1.0,0,0);
//        }
if((pd[i]<0.0 ) | (pd[i]>1.0)){ Rprintf("j %d pop[j] %d i %d pd[i] %f\n", j, pop[j],i, pd[i]);
 Rprintf("i %d pi[i] %f, gammart %f\n", i, pi[i],  gammart);
 }
//       if(j==75 & isamp == 4){
////      for (i=0; i<100; i++){
//Rprintf("j %d dis %d i %d l[i] %f pd[i] %f\n", j, ddis, i, lambdad[i], pd[i]);
//}// }
//      pd[i]=p1i[i]*pd[i];
//       }
      // Set up pd to be cumulative for the random draws
      for (i=1; i<Ki; i++){
       pd[i]=pd[i-1]+pd[i];
if((pd[i]<0.0 ) | (pd[i]>1.0)){ Rprintf("j %d pop[j] %d i %d pd[i] %f\n", j, pop[j],i, pd[i]);}
      }
      /* Draw unobserved degrees sizes */
      for (i=0; i<ni; i++){
       if(pop[i]==(j)){
        /* Now propose the true size for unit i based on reported size and disease status */
        /* In the next three lines a sizei is chosen */
        sizei=1;
        temp = pd[Ki-1] * unif_rand();
        while(temp > pd[sizei-1]){sizei++;}
        nk[sizei-1]=nk[sizei-1]+1;
        d[i]=sizei;
//Rprintf("j %d dis %d sizei %d pd[Ki-1] %f\n", j, ddis, sizei, pd[Ki-1]);
       }
      }
     } //compute
    } //for j
    b[ni-1]=d[ni-1];
    for (i=(ni-2); i>=0; i--){
      b[i]=b[i+1]+d[i];
    }
// Rprintf("j %d d[j] %d pd[Ki-1] %f\n", j, d[j], pd[Ki-1]);
    }

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
        sizei = 1;
	temp = unif_rand();
//      gammart = pi[Ki-1] * unif_rand();
        while(temp > pi[sizei-1]){sizei++;}
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
  free(psample);
  free(pdegi);
  free(b);
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
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *verbose
			 ) {
  int Np;
  int step, taken, give_log1=1, give_log0=0;
  int i, Ki, Ni, isamp, iinterval, isamplesize, iburnin;
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
  iburnin=(*burnin);
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

  // Now do the MCMC updates (starting with the burnin updates)
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
