/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "posteriorcmp.h"
#include "cmp.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gcmp (int *pop,
            int *nk, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu, double *dfmu, 
            double *sigma, double *dfsigma,
            int *Npi,
            double *lnlamproposal, 
            double *nuproposal, 
            int *N, int *maxN, 
            double *sample, 
            double *ppos, 
            double *lpriorm, 
            int *burnintheta,
            int *verbose
                         ) {
  int dimsample, Np;
  int step, staken, getone=1, intervalone=1, verboseMHcmp = 0;
  int i, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  double lnlami, nui, dsamp;
  double dmu, dsigma;
  int tU, sizei, imaxN, imaxm, give_log0=0, give_log1=1;
  int maxpop;
  double r, gammart, pis, Nd;
  double temp;
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
  dsigma=(*sigma);
  dmu=(*mu);

  dimsample=5+Np;

  double *pi = (double *) malloc(sizeof(double) * Ki);
  int *d = (int *) malloc(sizeof(int) * ni);
  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  int *Nkpos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdegi = (double *) malloc(sizeof(double) * (Np+1));
  double *psample = (double *) malloc(sizeof(double) * (Np+1));
  double *lnlamsample = (double *) malloc(sizeof(double));
  double *nusample = (double *) malloc(sizeof(double));

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
  lnlamsample[0] = dmu;
  nusample[0] = dsigma;

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

     for (i=0; i<Np; i++){
      pdegi[i] = psample[i];
     }
     lnlami=lnlamsample[0];
     nui=nusample[0];
//if(nui > 4.0 || lnlami > 4.5) Rprintf("lnlami %f nui %f dfmu %f\n", lnlami, nui, (*dfmu));
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
      pop[i]=sizei;
//    if((sizei <= 0) | (sizei > Ki-1)) Rprintf("sizei %d r %f\n", sizei,r);
      Nk[sizei-1]=Nk[sizei-1]+1;
    }
    if (step > 0 && step==(iinterval*(step/iinterval))) { 
      /* record statistics for posterity */
      Nd=(double)Ni;
      sample[isamp*dimsample  ]=Nd;
//if(nui > 4.0 || lnlami > 4.5) Rprintf("sample: %f %f\n", lnlami,nui);
// Rprintf("sample: %f %f\n", lnlami,nui);
      sample[isamp*dimsample+1]=lnlami;
      sample[isamp*dimsample+2]=nui;
      sample[isamp*dimsample+3]=(double)(Nk[0]);
      temp=0.0;
      for (i=0; i<Ki; i++){
        temp+=(i+1.0)*Nk[i];
      }
      sample[isamp*dimsample+4]=temp;
      for (i=0; i<Np; i++){
        sample[isamp*dimsample+5+i]=pdegi[i];
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
  free(d);
  free(psample);
  free(pdegi);
  free(b);
  free(Nk);
  free(Nkpos);
  free(lpm);
  free(lnlamsample);
  free(nusample);
}

void MHcmptheta (int *Nk, int *K,
            double *mu, double *dfmu, 
            double *sigma,  double *dfsigma,
            double *lnlamproposal, 
            double *nuproposal, 
            int *N, int *Npi, double *psample,
            double *lnlamsample, double *nusample,
            int *samplesize, int *staken, int *burnintheta, int *interval,
            int *verbose
         ) {
  int Np;
  int step, taken, give_log1=1, give_log0=0;
  int i, Ki, Ni, isamp, iinterval, isamplesize, iburnintheta;
  double ip, cutoff;
  double mui, mustar, lnlamstar, lnlami, lp;
  double pis, pstars;
  double sigmastar, sigmai, sigma2star, sigma2i, qnustar, qnui;
  double nustar, nui;
  double pithetastar, pithetai;
  double ddfmu, rdfmu, ddfsigma, dmu;
  double dsigma, dsigma2, dlnlamproposal, dnuproposal;
  double errval=0.0000000001, lzcmp;

//GetRNGstate();  /* R function enabling uniform RNG */

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
  iburnintheta=(*burnintheta);
  ddfmu=(*dfmu);
  rdfmu=sqrt(ddfmu);
  ddfsigma=(*dfsigma);
  dsigma=(*sigma);
  dsigma2=(dsigma*dsigma);
  dmu=(*mu);
  dnuproposal=(*nuproposal);
  dlnlamproposal=(*lnlamproposal);

  // First set starting values
  isamp = taken = 0;
  step = -iburnintheta;
  pis=1.;
  for (i=0; i<Np; i++){
    pdegi[i] = psample[i];
    pis-=pdegi[i];
  }
  for (i=0; i<Np; i++){
    odegi[i] = log(pdegi[i]/pis);
  }
  lnlami = lnlamsample[0];
  nui = nusample[0];
// if(nui > 4.0 || lnlami > 4.5) Rprintf("%f %f\n", lnlami,nui);
// Rprintf("%f %f\n", lnlami,nui);
//    Rprintf("mui %f nui %f lzcmp %f\n", mui, nui, lzcmp);
  pis=0.;
  lzcmp = zcmp(exp(lnlami), nui, errval, 2*Ki, give_log1);
//Rprintf("lnlami %f nui %f lzcmp %f\n", lnlami, nui, lzcmp);
  pi[Np]=cmp(Np+1,lnlami,nui,lzcmp,give_log0);
  for (i=Np+1; i<Ki; i++){
    pi[i]=pi[i-1]*exp(lnlami-nui*log((double)(i+1)));
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

  // Now computes mean and s.d. from log-lambda and nu
  mui=0.0;
  sigma2i=0.0;
  for (i=0; i<Ki; i++){
    mui+=pi[i]*(i+1);
    sigma2i+=pi[i]*(i+1)*(i+1);
  }
  sigma2i=sigma2i-mui*mui;

  sigmai  = sqrt(sigma2i);
  pithetai = dnorm(mui, dmu, sigmai/rdfmu, give_log1);
  pithetai = pithetai+dsclinvchisq(sigma2i, ddfsigma, dsigma2);

  // Now do the MCMC updates (starting with the burnin updates)
  while (isamp < isamplesize) {
    /* Propose new theta */
    /* Now the degree distribution model parameters */
    for (i=0; i<Np; i++){
      odegstar[i] = rnorm(odegi[i], dlnlamproposal);
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
    lnlamstar = rnorm(lnlami, dlnlamproposal);
    nustar = nui*exp(rnorm(0., dnuproposal));
    /* Check for magnitude */

    pstars=0.;
    lzcmp = zcmp(exp(lnlamstar), nustar, errval, 2*Ki, give_log1);
//if(nustar > 4.0 || lnlamstar > 4.5)  Rprintf("lnlamstar %f nustar %f lzcmp %f\n", lnlamstar, nustar, lzcmp);
//    if(nustar > 4.0 || lnlamstar > 4.5){step++;continue;}
    pstar[Np]=cmp(Np+1,lnlamstar,nustar,lzcmp,give_log0);
//  Rprintf("lnlamstar %f nustar %f lzcmp %f pstar %f\n", lnlamstar, nustar, lzcmp,pstar[Np]);
    for (i=Np+1; i<Ki; i++){
      pstar[i]=pstar[i-1]*exp(lnlamstar-nustar*log((double)(i+1)));
    }
    pstars=1.-exp(-lzcmp);
    for (i=0; i<Ki; i++){
      pstar[i]/=pstars;
//if(pstar[Np] < 0.00001){Rprintf("i %d pstar %f pi %f\n", i, pstar[i], pi[i]);}
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

    // Now compute mean and s.d. from log-lambda and nu
//  pstars=0.;
    mustar=0.0;
    sigma2star=0.0;
    for (i=0; i<Ki; i++){
      mustar+=pstar[i]*(i+1);
      sigma2star+=pstar[i]*(i+1)*(i+1);
//    pstars+=pstar[i];
    }
    sigma2star=sigma2star-mustar*mustar;
//  Rprintf("Check: %f\n", pstars);

    sigmastar  = sqrt(sigma2star);

//  if(nustar > 4.0 || lnlamstar > 4.5) Rprintf("%f %f %f %f %f\n", lnlamstar, dlnlam, nustar, ddflnlam, nui);
    /* Calculate pieces of the posterior. */
    qnustar = dnorm(log(nustar/nui)/dnuproposal,0.,1.,give_log1)
                  -log(dnuproposal*nustar);
    pithetastar = dnorm(mustar, dmu, sigmastar/rdfmu, give_log1);
    pithetastar = pithetastar+dsclinvchisq(sigma2star, ddfsigma, dsigma2);
    qnui = dnorm(log(nui/nustar)/dnuproposal,0.,1.,give_log1)
               -log(dnuproposal*nui);

    /* Calculate ratio */
    ip = pithetastar-pithetai;

    for (i=0; i<Ki; i++){
//if(pstar[Np] < 0.00001){Rprintf("i %d pstar %f pi %f\n", i, pstar[i], pi[i]);}
     if(Nk[i]>0){
      lp = log(pstar[i]/pi[i]);
      if(fabs(lp) < 100.){ip += (Nk[i]*lp);}
     }
    }
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
    then let the MH probability equal min{exp(cutoff), 1.0}.
    But we'll do it in log space instead.  */
//  if (*verbose)
//  Rprintf("Now proposing %d MH steps %f ip1...\n", step, ip);
    cutoff = ip + qnui-qnustar;
      
// Rprintf("mui %f mustar %f dmu %f ip %f pithetai, %f pithetastar %f cutoff %f\n", mui, mustar, dmu, ip, pithetai, pithetastar, cutoff);
//  if (*verbose)
//    Rprintf("Now proposing %d MH steps %f cutoff...\n", step, cutoff);

    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed changes */
      for (i=0; i<Np; i++){
        odegi[i] = odegstar[i];
        pdegi[i] = pdegstar[i];
      }
      lnlami    = lnlamstar;
      nui = nustar;
      qnui = qnustar;
      pithetai = pithetastar;
      for (i=0; i<Ki; i++){
        pi[i] = pstar[i];
      }
      taken++;
      if (step > 0 && step==(iinterval*(step/iinterval))) { 
        /* record statistics for posterity */
        lnlamsample[isamp]=lnlami;
        nusample[isamp]=nui;
        for (i=0; i<Np; i++){
          psample[i]=pdegi[i];
        }
        isamp++;
        if (*verbose && isamplesize==(isamp*(isamplesize/isamp))) Rprintf("Taken %d MH samples...\n", isamp);
//      Rprintf("Taken %d MH samples...\n", isamp);
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
//PutRNGstate();  /* Disable RNG before returning */
  /*Check for interrupts (if recursion is taking way too long...)*/
  R_CheckUserInterrupt();
  *staken = taken;
}
