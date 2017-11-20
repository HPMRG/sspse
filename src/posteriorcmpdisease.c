/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "posteriorcmpdisease.h"
#include "cmp.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gcmpdisease (int *pop, int *dis,
            int *nk0, int *nk1, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu0, double *mu1, double *kappa0, 
            double *sigma0,  double *sigma1, double *df0,
	    int *Np0i, int *Np1i,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *maxN, 
            double *sample, 
            double *p0pos, double *p1pos, 
            double *ppos, 
	    double *lpriorm,
            int *burnintheta,
	    double *lambdad,
	    double *nud,
            int *verbose
              ) {
  int dimsample, Np0, Np1;
  int step, staken, getone=1, intervalone=1, verboseMHdisease = 0;
  int i, j, compute, ni, Ni, Ki, isamp, iinterval, isamplesize, iburnin;
  double mu0i, mu1i, pbeta, beta, sigma0i, sigma1i, dsamp;
  double dkappa0, ddf0, dmu0, dmu1, dsigma0, dsigma1, dmuproposal, dsigmaproposal;
  int tU, sizei, imaxN, imaxm, itotdis0, itotdis, give_log0=0, give_log1=1;
  int maxpop, ddis;
  double r, gamma0rt, gamma1rt, p0is, p1is, Nd;
  double gammart, temp;
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
  Np0=(*Np0i);
  Np1=(*Np1i);
  dkappa0=(*kappa0);
  ddf0=(*df0);
  dsigma0=(*sigma0);
  dsigma1=(*sigma1);
  dmu0=(*mu0);
  dmu1=(*mu1);
  dsigmaproposal=(*sigmaproposal);
  dmuproposal=(*muproposal);

  dimsample=8+Np0+Np1;

  double *p0i = (double *) malloc(sizeof(double) * Ki);
  double *p1i = (double *) malloc(sizeof(double) * Ki);
  double *pd = (double *) malloc(sizeof(double) * Ki);
  int *d = (int *) malloc(sizeof(int) * ni);
  int *b = (int *) malloc(sizeof(int) * ni);
  int *Nk0 = (int *) malloc(sizeof(int) * Ki);
  int *Nk0pos = (int *) malloc(sizeof(int) * Ki);
  int *Nk1 = (int *) malloc(sizeof(int) * Ki);
  int *Nk1pos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdeg0i = (double *) malloc(sizeof(double) * Np0);
  double *pdeg1i = (double *) malloc(sizeof(double) * Np1);
  double *psample = (double *) malloc(sizeof(double) * (Np0+Np1));
  double *musample = (double *) malloc(sizeof(double) * 2);
  double *betasample = (double *) malloc(sizeof(double));
  double *sigmasample = (double *) malloc(sizeof(double) * 2);

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
  /* Draw initial phis */
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
  for (i=0; i<Np0; i++){
     psample[i    ] = 0.01;
  }
  for (i=0; i<Np1; i++){
     psample[i+Np0]=0.01;
  }
  musample[0] = dmu0;
  musample[1] = dmu1;
  sigmasample[0] = dsigma0;
  sigmasample[1] = dsigma1;

  isamp = 0;
  step = -iburnin;
  while (isamp < isamplesize) {

    /* Draw new theta */
    /* but less often than the other full conditionals */
    if (step == -iburnin || step==(10*(step/10))) { 
// Rprintf("Ni %d itotdis %d K %d Nk0[0] %d Nk1[0] %d mu0 %f mu1 %f sigma0 %f sigma1\n", Ni, itotdis, *K, Nk0[0], Nk1[0],*mu0,*mu1,*sigma0,*sigma1);
    MHcmpdisease(Nk0,Nk1,&itotdis,K,mu0,mu1,kappa0,sigma0,sigma1,df0,
          muproposal, sigmaproposal,
	  &Ni, &Np0, &Np1, psample, 
	  musample, betasample, sigmasample, &getone, &staken, 
	  burnintheta, &intervalone, 
	  &verboseMHdisease);
    }

    beta=betasample[0];
    pbeta=exp(beta)/(1.+exp(beta));
    for (i=0; i<Np0; i++){
      pdeg0i[i] = psample[i];
    }
    for (i=0; i<Np1; i++){
      pdeg1i[i] = psample[i+Np0];
    }
    mu0i=musample[0];
    mu1i=musample[1];
    sigma0i=sigmasample[0];
    sigma1i=sigmasample[1];

    /* Draw new N */

    /* First find the degree distribution */
    p0is=0.;
    p1is=0.;
    lzcmp = zcmp(exp(mu0i), sigma0i, errval, 2*Ki, give_log1);
    if(lzcmp < -100000.0){continue;}
    p0i[Np0]=cmp(Np0+1,mu0i,sigma0i,lzcmp,give_log0);
    for (i=Np0+1; i<Ki; i++){
//    p0i[i]=cmp(i+1,mu0i,sigma0i,lzcmp,give_log0);
      p0i[i]=p0i[i-1]*exp(mu0i-sigma0i*log((double)(i+1)));
// Rprintf("i %d ss %e c %e lzcmp %e\n", i, p0i[i], cmp(i+1,mu0i,sigma0i,lzcmp,give_log0),lzcmp);
    }
//  p0is=1.-cmp(0,mu0i,sigma0i,lzcmp,give_log0);
    p0is=1.-exp(-lzcmp);
    lzcmp = zcmp(exp(mu1i), sigma1i, errval, 2*Ki, give_log1);
    if(lzcmp < -100000.0){continue;}
    p1i[Np1]=cmp(Np1+1,mu1i,sigma1i,lzcmp,give_log0);
    for (i=Np1+1; i<Ki; i++){
//    p1i[i]=cmp(i+1,mu1i,sigma1i,lzcmp,give_log0);
      p1i[i]=p1i[i-1]*exp(mu1i-sigma1i*log((double)(i+1)));
//    p1is+=p1i[i];
    }
//  p1is=1.-cmp(0,mu1i,sigma1i,lzcmp,give_log0);
    p1is=1.-exp(-lzcmp);
    for (i=0; i<Ki; i++){
      p0i[i]/=p0is;
      p1i[i]/=p1is;
    }
    p0is=1.;
    p1is=1.;
    for (i=0; i<Np0; i++){
      p0is-=pdeg0i[i];
    }
    for (i=0; i<Np1; i++){
      p1is-=pdeg1i[i];
    }
    for (i=0; i<Ki; i++){
      p0i[i]*=p0is;
      p1i[i]*=p1is;
    }
    for (i=0; i<Np0; i++){
      p0i[i]=pdeg0i[i];
    }
    for (i=0; i<Np1; i++){
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
    gammart=log((1.-pbeta)*gamma0rt+pbeta*gamma1rt);
    temp = -100000000.0;
    // N = m + n
    // Compute (log) P(m | \theta and data and \Psi)
    for (i=0; i<imaxm; i++){
     lpm[i]=lgamma(ni+i+1.)-lgamma(i+1.)+i*gammart;
     //  Add in the (log) prior on m: P(m)
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
//  if (*verbose) Rprintf("Ni %d lpm[imaxm-1] %f lpm[Ni] %f\n", Ni, lpm[imaxm-1],
//  lpm[Ni]);
//  }
    // Add back the sample size
    Ni += ni;
    if(Ni > imaxN) Ni = imaxN;

//  if (*verbose) Rprintf("step %d Ni %d itotdis %d beta %f mu0 %f mu1 %f s0 %f r %f\n",
//  step, Ni, itotdis, betasample[0], musample[0], musample[1], sigmasample[0], r);

// Rprintf("lambdad[0] %f\n", lambdad[0]);
// Rprintf("nud[0] %f\n", nud[0]);
    if((fabs(lambdad[0])>0.0000001) | (fabs(nud[0])>0.0000001)){
// Rprintf("lambdad[0] %f\n", lambdad[0]);
    for (i=0; i<Ki; i++){
      nk0[i]=0;
      nk1[i]=0;
    }

    /* Draw true degrees (sizes) based on the reported degrees*/
    /* First find the reported degree distribution */
    for (j=0; j<=maxpop; j++){
    for (ddis=0; ddis<2; ddis++){
//Rprintf("j %d pop[j] %d\n", j, pop[j]);
     compute=0;
     for (i=0; i<ni; i++){if((pop[i]==(j)) && (dis[i]==ddis)){compute=1;}}
     if(compute==1){
      if(ddis==1){
//    Next four lines for cmp reporting distribution
//    ?? Should it be cmp(j+1,...) or cmp(j,...)??
//    I think j=1 as pd[i], p1i[i], etc, refer to degree i+1 
//    as all degrees must be positive.
      for (i=0; i<Ki; i++){
       lzcmp = zcmp(exp(lambdad[i]),nud[i], errval, 2*Ki, give_log1);
       pd[i]=p1i[i]*cmp(j+1,lambdad[i],nud[i],lzcmp,give_log0);
      }
//     Next seven lines for proportional reporting distribution
//       for (i=0; i<Ki; i++){
//        pd[i]   = pgamma(2.0*lambdad[i]/((j)+0.5),1.0,1.0,0,0);
//        if((j)>0){
//         pd[i] -= pgamma(2.0*lambdad[i]/((j)-0.5),1.0,1.0,0,0);
//        }
// if((pd[i]<0.0 ) | (pd[i]>1.0)){ Rprintf("j %d pop[j] %d i %d pd[i] %f\n", j, pop[j],i, pd[i]);
//  Rprintf("i %d p1i[i] %f, gamma0rt %f gamma1rt %f \n", i, p1i[i],  gamma0rt,  gamma1rt);
//  }
//       if(j==75 & isamp == 4){
////      for (i=0; i<100; i++){
//Rprintf("j %d dis %d i %d l[i] %f pd[i] %f\n", j, ddis, i, lambdad[i], pd[i]);
//}// }
//      pd[i]=p1i[i]*pd[i];
//       }
      }else{
//      Next four lines for cmp reporting distribution
       for (i=0; i<Ki; i++){
        lzcmp = zcmp(exp(lambdad[i]),nud[i], errval, 2*Ki, give_log1);
        pd[i]=p0i[i]*cmp(j,lambdad[i],nud[i],lzcmp,give_log0);
       }
//     Next seven lines for proportional reporting distribution
//       for (i=0; i<Ki; i++){
//        pd[i]   = pgamma(2.0*lambdad[i]/((j)+0.5),1.0,1.0,0,0);
//        if((j)>0){
//         pd[i] -= pgamma(2.0*lambdad[i]/((j)-0.5),1.0,1.0,0,0);
//        }
//      pd[i]=p0i[i]*pd[i];
// Rprintf("i %d pd[i] %f, p0i[i] %f, pop[i] %d lambdad[i] %f nud[i] %f \n", i, log(pd[i]), p0i[i], pop[i], lambdad[i],nud[i]);
//       }
      }
      // Set up pd to be cumulative for the random draws
      for (i=1; i<Ki; i++){
       pd[i]=pd[i-1]+pd[i];
// if((pd[i]<0.0 ) | (pd[i]>1.0)){ Rprintf("j %d pop[j] %d i %d pd[i] %f\n", j, pop[j],i, pd[i]);}
      }
      /* Draw unobserved degrees sizes */
      for (i=0; i<ni; i++){
       if((pop[i]==(j)) && (dis[i]==ddis)){
        /* Now propose the true size for unit i based on reported size and disease status */
        /* In the next three lines a sizei is chosen */
        temp = pd[Ki-1] * unif_rand();
        for (sizei=1; sizei<=Ki; sizei++){
          if(temp <= pd[sizei-1]) break;
        }
        if(dis[i]==1){
          nk1[sizei-1]=nk1[sizei-1]+1;
        }else{
          nk0[sizei-1]=nk0[sizei-1]+1;
        }
        d[i]=sizei;
//Rprintf("j %d dis %d sizei %d pd[Ki-1] %f\n", j, ddis, sizei, pd[Ki-1]);
       }
      }
     } //compute
    }} //for j and ddis
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
      r+=exp_rand()/(tU+b[i]);
    }

    /* Draw unseen sizes */
    for (i=0; i<Ki; i++){
      Nk0[i]=nk0[i];
      Nk1[i]=nk1[i];
    }
    // Set up p0i and p1i to be cumulative for the random draws
    for (i=1; i<Ki; i++){
      p0i[i]=p0i[i-1]+p0i[i];
      p1i[i]=p1i[i-1]+p1i[i];
    }
    for (i=ni; i<Ni; i++){
      /* Propose unseen size for unit i */
      /* Use rejection sampling */
      sizei=1000000;
      while(sizei >= Ki){
       sizei=1000000;
       while(log(1.0-unif_rand()) > -r*sizei){
        /* First propose unseen disease status for unit i */
        if(unif_rand() < pbeta){
          dis[i]=1;
          /* Now propose unseen size for unit i based on disease status */
          /* In the next two lines a sizei is chosen */
          /* with parameters mu1i and sigma1i */
//        temp = p1i[Ki-1] * unif_rand();
          temp = unif_rand();
          for (sizei=1; sizei<=Ki; sizei++){
            if(temp <= p1i[sizei-1]) break;
          }
        }else{
          dis[i]=0;
          /* Now propose unseen size for unit i based on non-disease status */
          /* In the next two lines a sizei is chosen */
          /* with parameters mu0i and sigma0i */
//        temp = p0i[Ki-1] * unif_rand();
          temp = unif_rand();
          for (sizei=1; sizei<=Ki; sizei++){
            if(temp <= p0i[sizei-1]) break;
          }
        }
       }
      }
//    if(sizei >= Ki){sizei=Ki-1;}
      pop[i]=sizei;
      if(dis[i]==1){
        Nk1[sizei-1]=Nk1[sizei-1]+1;
      }else{
        Nk0[sizei-1]=Nk0[sizei-1]+1;
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
      for (i=0; i<Np0; i++){
        sample[isamp*dimsample+8+i]=pdeg0i[i];
      }
      for (i=0; i<Np1; i++){
        sample[isamp*dimsample+8+Np0+i]=pdeg1i[i];
      }
//    N0d=0.;
      for (i=0; i<Ki; i++){
        Nk0pos[i]=Nk0pos[i]+Nk0[i];
        Nk1pos[i]=Nk1pos[i]+Nk1[i];
        p0pos[i]+=((double)(Nk0[i])/Nd);
        p1pos[i]+=((double)(Nk1[i])/Nd);
        ppos[i]+=((double)(Nk0[i]+Nk1[i]))/Nd;
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
// Organize and return
  dsamp=((double)isamp);
  for (i=0; i<Ki; i++){
    nk0[i]=Nk0pos[i];
    nk1[i]=Nk1pos[i];
    p0pos[i]=p0pos[i]/dsamp;
    p1pos[i]=p1pos[i]/dsamp;
    ppos[i]=ppos[i]/dsamp;
  }
  for (i=0; i<ni; i++){
     pop[i]=d[i];
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(psample);
  free(pdeg0i);
  free(pdeg1i);
  free(pd);
  free(p0i);
  free(p1i);
  free(b);
  free(d);
  free(Nk0);
  free(Nk0pos);
  free(Nk1);
  free(Nk1pos);
  free(lpm);
  free(musample);
  free(betasample);
  free(sigmasample);
}

void MHcmpdisease (int *Nk0, int *Nk1, int *totdis, int *K,
	    double *mu0, double *mu1, double *kappa0, 
            double *sigma0,  double *sigma1, double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Np0i, int *Np1i, double *psample,
            double *musample, double *betasample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *verbose
			 ) {
  int Np0, Np1;
  int step, taken, give_log0=0, give_log1=1;
  int i, Ki, Ni, isamp, iinterval, isamplesize, iburnin, itotdis;
  double ip, cutoff;
  double mu0star, mu1star, mu0i, mu1i, lp;
  double pbetastar, pbeta, betastar, betai;
  double p0is, p1is, p0stars, p1stars;
  double sigma0star, sigma1star, sigma0i, sigma1i;
  double sigma02star, sigma12star, sigma02i, sigma12i;
  double qsigma02star, qsigma12star, qsigma02i, qsigma12i;
  double pithetastar, pithetai;
  double dkappa0, rkappa0, ddf0, dmu0, dmu1;
  double dsigma0, dsigma1, dsigma02, dsigma12, dmuproposal, dsigmaproposal;
  double errval=0.000000001, lzcmp;

  GetRNGstate();  /* R function enabling uniform RNG */

  Ki=(*K);
  Np0=(*Np0i);
  Np1=(*Np1i);
  double *p0star = (double *) malloc(sizeof(double) * Ki);
  double *p0i = (double *) malloc(sizeof(double) * Ki);
  double *p1star = (double *) malloc(sizeof(double) * Ki);
  double *p1i = (double *) malloc(sizeof(double) * Ki);
  double *odeg0star = (double *) malloc(sizeof(double) * Np0);
  double *odeg0i = (double *) malloc(sizeof(double) * Np0);
  double *odeg1star = (double *) malloc(sizeof(double) * Np1);
  double *odeg1i = (double *) malloc(sizeof(double) * Np1);
  double *pdeg0star = (double *) malloc(sizeof(double) * Np0);
  double *pdeg0i = (double *) malloc(sizeof(double) * Np0);
  double *pdeg1star = (double *) malloc(sizeof(double) * Np1);
  double *pdeg1i = (double *) malloc(sizeof(double) * Np1);

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
  for (i=0; i<Np0; i++){
    pdeg0i[i] = psample[i];
    p0is-= pdeg0i[i];
  }
  for (i=0; i<Np1; i++){
    pdeg1i[i] = psample[i+Np0];
    p1is-= pdeg1i[i];
  }
  for (i=0; i<Np0; i++){
    odeg0i[i] = log(pdeg0i[i]/p0is);
  }
  for (i=0; i<Np1; i++){
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
  lzcmp = zcmp(exp(mu0i), sigma0i, errval, 2*Ki, give_log1);
  p0i[Np0]=cmp(Np0+1,mu0i,sigma0i,lzcmp,give_log0);
  for (i=Np0+1; i<Ki; i++){
//  p0i[i]=cmp(i+1,mu0i,sigma0i,lzcmp,give_log0);
    p0i[i]=p0i[i-1]*exp(mu0i-sigma0i*log((double)(i+1)));
//  p0is+=p0i[i];
  }
//p0is=1.-cmp(0,mu0i,sigma0i,lzcmp,give_log0);
  p0is=1.-exp(-lzcmp);
  lzcmp = zcmp(exp(mu1i), sigma1i, errval, 2*Ki, give_log1);
  p1i[Np1]=cmp(Np1+1,mu1i,sigma1i,lzcmp,give_log0);
  for (i=Np1+1; i<Ki; i++){
//  p1i[i]=cmp(i+1,mu1i,sigma1i,lzcmp,give_log0);
    p1i[i]=p1i[i-1]*exp(mu1i-sigma1i*log((double)(i+1)));
//  p1is+=p1i[i];
  }
//p1is=1.-cmp(0,mu1i,sigma1i,lzcmp,give_log0);
  p1is=1.-exp(-lzcmp);
  for (i=0; i<Ki; i++){
    p0i[i]=p0i[i]/p0is;
    p1i[i]=p1i[i]/p1is;
  }
  p0is=1.;
  p1is=1.;
  for (i=0; i<Np0; i++){
    p0is-=pdeg0i[i];
  }
  for (i=0; i<Np1; i++){
    p1is-=pdeg1i[i];
  }
  for (i=0; i<Ki; i++){
    p0i[i]=p0i[i]*p0is;
    p1i[i]=p1i[i]*p1is;
  }
  for (i=0; i<Np0; i++){
    p0i[i]=pdeg0i[i];
  }
  for (i=0; i<Np1; i++){
    p1i[i]=pdeg1i[i];
  }
  // Bin last group
  p0is=1.;
  p1is=1.;
  for (i=0; i<(Ki-1); i++){
    p0is-=p0i[i];
    p1is-=p1i[i];
  }
  p0i[Ki-1]=p0is;
  p1i[Ki-1]=p1is;

  // Now do the MCMC updates (starting with the burnin updates)
  while (isamp < isamplesize) {
//  Rprintf("step %d Ni %d itotdis %d isamp %d\n", step, Ni, itotdis, isamp);
    /* Propose new theta */
    /* Start with the disease status parameters */
    betastar = rnorm(betai, dmuproposal);
    pbetastar = exp(betastar)/(1.+exp(betastar));
    /* Now the degree distribution model parameters */
    for (i=0; i<Np0; i++){
      odeg0star[i] = rnorm(odeg0i[i], dmuproposal);
    }
    for (i=0; i<Np1; i++){
      odeg1star[i] = rnorm(odeg1i[i], dmuproposal);
    }
    /* Convert from odds to probabilities */
    p0is=1.;
    p1is=1.;
    for (i=0; i<Np0; i++){
      pdeg0star[i] = exp(odeg0star[i]);
      p0is+=pdeg0star[i];
    }
    for (i=0; i<Np1; i++){
      pdeg1star[i] = exp(odeg1star[i]);
      p1is+=pdeg1star[i];
    }
    for (i=0; i<Np0; i++){
      pdeg0star[i]/=p0is;
    }
    for (i=0; i<Np1; i++){
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
//    Rprintf("pbetastar %f betastar %f\n", pbetastar, betastar);
    ip+=(itotdis*log(pbetastar)+(Ni-itotdis)*log(1.-pbetastar));
    pbeta = exp(betai)/(1.+exp(betai));
    ip-=(itotdis*log(pbeta)+(Ni-itotdis)*log(1.-pbeta));
//    Rprintf("pbeta %f betai %f betastar %f\n", pbeta, betai, betastar);
//    Rprintf("pithetastar %f pithetai %f qsigma2star %f qsigmai %f\n",
//           pithetastar,pithetai,qsigma2star,qsigma2i);
    p0stars=0.;
    p1stars=0.;
    lzcmp = zcmp(exp(mu0star), sigma0star, errval, 2*Ki, give_log1);
    if(lzcmp < -100000.0){continue;}
    p0star[Np0]=cmp(Np0+1,mu0star,sigma0star,lzcmp,give_log0);
    for (i=Np0+1; i<Ki; i++){
//    p0star[i]=cmp(i+1,mu0star,sigma0star,lzcmp,give_log0);
      p0star[i]=p0star[i-1]*exp(mu0star-sigma0star*log((double)(i+1)));
//    p0stars+=p0star[i];
    }
//  p0stars=1.-cmp(0,mu0star,sigma0star,lzcmp,give_log0);
    p0stars=1.-exp(-lzcmp);
    lzcmp = zcmp(exp(mu1star), sigma1star, errval, 2*Ki, give_log1);
    if(lzcmp < -100000.0){continue;}
    p1star[Np1]=cmp(Np1+1,mu1star,sigma1star,lzcmp,give_log0);
    for (i=Np1+1; i<Ki; i++){
//    p1star[i]=cmp(i+1,mu1star,sigma1star,lzcmp,give_log0);
      p1star[i]=p1star[i-1]*exp(mu1star-sigma1star*log((double)(i+1)));
//    p1stars+=p1star[i];
    }
//  p1stars=1.-cmp(0,mu1star,sigma1star,lzcmp,give_log0);
    p1stars=1.-exp(-lzcmp);
    for (i=Np0; i<Ki; i++){
      p0star[i]/=p0stars;
    }
    for (i=Np1; i<Ki; i++){
      p1star[i]/=p1stars;
    }
    p0stars=1.;
    p1stars=1.;
    for (i=0; i<Np0; i++){
      p0stars-=pdeg0star[i];
    }
    for (i=0; i<Np1; i++){
      p1stars-=pdeg1star[i];
    }
    for (i=Np0; i<Ki; i++){
      p0star[i]*=p0stars;
    }
    for (i=Np1; i<Ki; i++){
      p1star[i]*=p1stars;
    }
    for (i=0; i<Np0; i++){
      p0star[i]=pdeg0star[i];
    }
    for (i=0; i<Np1; i++){
      p1star[i]=pdeg1star[i];
    }
    // Bin last group
    p0stars=1.;
    p1stars=1.;
    for (i=0; i<(Ki-1); i++){
      p0stars-=p0star[i];
      p1stars-=p1star[i];
    }
    p0star[Ki-1]=p0stars;
    p1star[Ki-1]=p1stars;

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
      
//  Rprintf("Now proposing %d MH steps %f cutoff...\n", step, cutoff);

    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed changes */
      betai = betastar;
      for (i=0; i<Np0; i++){
        odeg0i[i] = odeg0star[i];
        pdeg0i[i] = pdeg0star[i];
      }
      for (i=0; i<Np1; i++){
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
        musample[2*isamp]=mu0i;
        musample[2*isamp+1]=mu1i;
        sigmasample[2*isamp]=sigma0i;
        sigmasample[2*isamp+1]=sigma1i;
        betasample[isamp]=betai;
        for (i=0; i<Np0; i++){
          psample[i    ]=pdeg0i[i];
        }
        for (i=0; i<Np1; i++){
          psample[i+Np0]=pdeg1i[i];
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
  /*Check for interrupts (if recursion is taking way too long...)*/
  R_CheckUserInterrupt();
  *staken = taken;
}
