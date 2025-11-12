/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "posteriorcmpvis.h"
#include "posteriorcmpwpvis2.h"
#include "cmp.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gcmpwpvis2 (int *pop12, int *pop21,
            int *K,
            int *n1i,
            int *n2i,
            int *n12i,
            int *samplesize, int *warmup, int *interval,
            double *mu, double *dfmu,
            double *sigma, double *dfsigma,
            double *lnlam, double *nu,
            double *beta0muprior, double *beta0sigmaprior,
            double *betatmuprior, double *betatsigmaprior,
            double *betaumuprior, double *betausigmaprior,
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            double *memod,
            int *Npi,
            int *srd,
            int *numrec,
            double *rectime,
            int *srd2,
            int *numrec2,
            double *rectime2,
            int *rc,
            int *maxcoupons,
            double *lnlamproposal,
            double *nuproposal,
            double *beta0proposal, double *betatproposal, double *betauproposal,
            double *lmemmuproposal, double *memnuproposal,
            int *N, int *maxN,
            double *sample,
            int *vsample,
            int *vsample2,
            double *posu,
            double *posd,
            double *lpriorm,
            int *warmuptheta,
            int *warmupbeta,
            int *verbose
                         ) {
  int dimsample, Np;
  int n12, n1, n2, itemp;
  int step, staken, getone=1, intervalone=1, verboseMHcmp = 0;
  int i, j, k, n, Ni, Ki, isamp, iinterval, isamplesize, iwarmup;
  int umax, unrecap;
  double alpha, pnb, rnb;
  double memmui;
  double mui, sigmai, dsamp, nui, lnlami, sigma2i;
  double dbeta0, dbetat, dbetau;
  double dlmemmu, dmemnu, dmemmu;
  double beta0i, betati, betaui, lmemmui, memnui;
  int tU1, tU2, sizei, imaxN, imaxm, give_log0=0, give_log1=1;
  double r1, r2, gammart, pis, pis2, Nd;
  double temp, temp2, uprob1, uprob2;
  double lliki;
  int maxc;
  double errval=0.0000000001, lzcmp;

  GetRNGstate();  /* R function enabling uniform RNG */

  n1=(*n1i);
  n2=(*n2i);
  n12=(*n12i);
  n = n1 + n2 - n12; /* The number unique people seen */
  Ni=(*N);
  Ki=(*K);
  imaxN=(*maxN);
  imaxm=imaxN-n;
  isamplesize=(*samplesize);
  iinterval=(*interval);
  iwarmup=(*warmup);
  Np=(*Npi);
  dbeta0=(*beta0muprior);
  dbetat=(*betatmuprior);
  dbetau=(*betaumuprior);
  dlmemmu=(*lmemmu);
  dmemmu=exp(*lmemmu);
  dmemnu=(*memnu);
  alpha=(*memod);
  maxc=(*maxcoupons);

  dimsample=5+Np+4+1;
  pnb=(alpha-1.)/alpha;

  double *pi = (double *) malloc(sizeof(double) * Ki);
  double *pd = (double *) malloc(sizeof(double) * Ki);
  double *pdm = (double *) malloc(sizeof(double) * Ki);
  double *pd2 = (double *) malloc(sizeof(double) * (10*Ki));
  int *u1 = (int *) malloc(sizeof(int) * imaxN);
  int *u2 = (int *) malloc(sizeof(int) * imaxN);
  int *b1 = (int *) malloc(sizeof(int) * n1);
  int *b2 = (int *) malloc(sizeof(int) * n2);
  int *nk = (int *) malloc(sizeof(int) * Ki);
  int *Nk = (int *) malloc(sizeof(int) * Ki);
  int *Nkpos = (int *) malloc(sizeof(int) * Ki);
  double *lpm = (double *) malloc(sizeof(double) * imaxm);
  double *pdegi = (double *) malloc(sizeof(double) * (Np+1));
  double *psample = (double *) malloc(sizeof(double) * (Np+1));
  double *lnlamsample = (double *) malloc(sizeof(double) * isamplesize);
  double *nusample = (double *) malloc(sizeof(double) * isamplesize);
  double *beta0sample = (double *) malloc(sizeof(double) * isamplesize);
  double *betatsample = (double *) malloc(sizeof(double) * isamplesize);
  double *betausample = (double *) malloc(sizeof(double) * isamplesize);
  double *lmemmusample = (double *) malloc(sizeof(double) * isamplesize);
  double *memnusample = (double *) malloc(sizeof(double) * isamplesize);

  for (i=0; i<Ki; i++){
    nk[i]=0;
  }
  unrecap=0;
  uprob1=n;
// Rprintf("n12 %d n1 %d n2 %d n %d\n", n12, n1, n2,n);
  for (i=0; i<n1; i++){
    srd[i] = (int)round( srd[i] / dmemmu);
    if(srd[i]<1){ srd[i]=1;}
  }
  for (i=0; i<n; i++){
 Rprintf("i %d pop12[i] %d dmemmu %f round(pop12[i] / dmemmu) %d\n", i, pop12[i], dmemmu, (int)round(pop12[i] / dmemmu));
    if((pop12[i] >0) && (pop12[i] <= Ki*dmemmu)){ u1[i]=(int)round(pop12[i] / dmemmu);}
    if( pop12[i]==0){ u1[i]=1;}
    if( pop12[i]>Ki*dmemmu){ u1[i]=Ki; uprob1--;}
    nk[u1[i]-1]=nk[u1[i]-1]+1;
    unrecap+=u1[i];
  }
  // VIP setting the mem.optimism.prior = mem.optimism to 1
  // as the network.sizes have been divided by it
  dlmemmu=0.0;
  *lmemmu=0.0;
  //
  uprob1/=n;
  uprob1 = 0.5 + uprob1/2.0;
  uprob2=n2+1;
  for (i=0; i<n2; i++){
    srd2[i] = (int)round( srd2[i] / dmemmu);
    if(srd2[i]<1){ srd2[i]=1;}
  }
  for (i=0; i<n2; i++){
 Rprintf("i %d pop21[i] %d dmemmu %f (int)round(pop21[i] / dmemmu) %d\n", i, pop21[i], dmemmu, (int)round(pop21[i] / dmemmu));
    if((pop21[i] >0) && (pop21[i] <= Ki*dmemmu)){ u2[i]=(int)round(pop21[i] / dmemmu);}
    if( pop21[i]==0){ u2[i]=1;}
    if( pop21[i]>Ki*dmemmu){ u2[i]=Ki; uprob2--;}
    unrecap-=u2[i];
    if(rc[i]!=0) unrecap-=u2[i];
  }

  uprob2/=n2+1;
  uprob2 = 0.5 + uprob2/2.0;
  b1[n1-1]=u1[n1-1];
  for (i=(n1-2); i>=0; i--){
    b1[i]=b1[i+1]+u1[i];
  }
  if(n2 > 0){
    b2[n2-1]=u2[n2-1];
    for (i=(n2-2); i>=0; i--){
      b2[i]=b2[i+1]+u2[i];
    }
  }
  for (i=0; i<Ki; i++){
     Nk[i]=nk[i];
     Nkpos[i]=0;
     posu[i]=0.;
  }
  for (i=0; i<10*Ki; i++){
     posd[i]=0.;
  }
  for (i=n1; i<imaxN; i++){
    u1[i]=u1[(int)trunc(10*unif_rand()+n1-10)];
  }
  tU1=0;
  for (i=n1; i<Ni; i++){
    tU1+=u1[i];
  }
  for (i=n2; i<imaxN; i++){
    u2[i]=u2[(int)trunc(10*unif_rand()+n2-10)];
  }
  tU2=unrecap;
  for (i=n2; i<Ni; i++){
    tU2+=u2[i];
  }
  /* Draw initial phis */
  r1=0.;
  for (i=0; i<n1; i++){
    r1+=(exp_rand()/(tU1+b1[i]));
  }
  r2=0.;
  for (i=0; i<n2; i++){
    r2+=(exp_rand()/(tU2+b2[i]));
  }

  for (i=0; i<Np; i++){
     psample[i] = 0.01;
  }
  beta0sample[0] = dbeta0;
  betatsample[0] = dbetat;
  betausample[0] = dbetau;
  lmemmusample[0] = dlmemmu;
  memnusample[0] = dmemnu;
  lnlamsample[0] = (*lnlam);
  nusample[0] = (*nu);

  isamp = 0;
  step = -iwarmup;
  while (isamp < isamplesize) {
// Rprintf("isamp %d step %d lnlamproposal %f nuproposa %f \n", isamp, step, lnlamproposal,nuproposal);

    /* Draw new theta */
    /* but less often than the other full conditionals */
    if (step == -iwarmup || step % 10 == 0) {
     MHcmptheta(Nk,K,mu,dfmu,sigma,dfsigma,lnlamproposal,nuproposal,
       &Ni, &Np, psample,
       lnlamsample, nusample, &getone, &staken, warmuptheta, &intervalone,
       &verboseMHcmp);
     }

// Rprintf("lnlamsample[0] %f\n", lnlamsample[0]);
     lnlami=lnlamsample[0];
     nui=nusample[0];

     for (i=0; i<Np; i++){
      pdegi[i] = psample[i];
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
 // pis=1.-exp(-lzcmp);
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
    if (step == -iwarmup || step % 10 == 0) {
     MHwpmem2(u1,u2,&n1,&n2,K,beta0muprior,beta0sigmaprior,betatmuprior,betatsigmaprior,betaumuprior,betausigmaprior,
       lmemmu,memdfmu,memnu,memdfnu,memod,srd,numrec,rectime,srd2,numrec2,rectime2,maxcoupons,
       beta0proposal, betatproposal, betauproposal,
       lmemmuproposal, memnuproposal,
       beta0sample, betatsample, betausample, lmemmusample,memnusample,
       &getone, &staken, warmupbeta, &intervalone,
       &verboseMHcmp);
     beta0i=beta0sample[0];
     betati=betatsample[0];
     betaui=betausample[0];
     lmemmui=lmemmusample[0];
     memnui=memnusample[0];
    }

    /* Draw true unit (sizes) of the first RDS sample based on the reported degrees*/
    // First reset counts
    for (i=0; i<Ki; i++){
     nk[i]=0;
    }
    umax=0;
    for (j=0; j<n1; j++){
      if(srd[j] <= (10*Ki)){
//     Multiply by the Conway-Maxwell-Poisson PMF for observation
       for (i=0; i<Ki; i++){
        // Next to exclude unit sizes inconsistent with the number of recruits
        temp = beta0i + betati*log(rectime[j]) + betaui*log(i+1.0);
         lliki=0.0;
         if(numrec[j]<maxc){
           lliki += dpois(numrec[j],exp(temp),give_log1);
         }else{
 	  lliki += log(1.0-ppois(maxc-1.0,exp(temp),give_log0,give_log0));
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
          if(srd[i] <= 10*Ki){
            lliki += log(pd2[srd[j]-1]);
          }else{
            lliki += log(pis2);
          }
        }
        pd[i]=pi[i]*exp(lliki)*(i+1.);
        pdm[i]=exp(lliki);
      }
      // Set up pd to be cumulative for the random draws
      temp2 = pdm[0];
      temp = pdm[0];
      for (i=1; i<Ki; i++){
       pd[i]=pd[i-1]+pd[i];
       temp2+=pdm[i]*(i+1);
       temp+=pdm[i];
      }
      temp2/=temp;
      /* Draw unit size for the observed degree */
      /* Now propose the true size for unit i based on reported size and disease status */
      /* In the next three lines a sizei is chosen */
      temp = pd[Ki-1] * unif_rand();
      for (sizei=1; sizei<=Ki; sizei++){
        if(temp <= pd[sizei-1]) break;
      }
      nk[sizei-1]=nk[sizei-1]+1;
      u1[j]=sizei;

     }
     }
     /* Deal with the outliers */
     umax=nk[Ki-1];
     sizei=Ki;
     for (j=0; j<n1; j++){
      if(srd[j] > 10*Ki){
 Rprintf("j %d srd[j] %d nk[j] %d\n", j, srd[j], nk[j]);
       while((umax==0) & (sizei > 1)){
         sizei--;
         umax=nk[sizei-1];
       }
       u1[j]=sizei;
       umax--;
      }
     }
     for (j=0; j<n1; j++){
      if(srd[j] > 10*Ki){
       nk[u1[j]-1]=nk[u1[j]-1]+1;
      }
     }

     // Rebuild b1
     b1[n1-1]=u1[n1-1];
     for (i=(n1-2); i>=0; i--){
      b1[i]=b1[i+1]+u1[i];
     }
     /* End of imputed unit sizes for observed in the first RDS*/

    /* Draw true degrees (sizes) of the second RDS sample based on the reported degrees*/
    for (j=0; j<n2; j++){
     if(srd2[j] <= (10*Ki)){
//    Multiply by the Conway-Maxwell-Poisson PMF for observation
      for (i=0; i<Ki; i++){
       // Next to exclude unit sizes inconsistent with the number of recruits
       temp = beta0i + betati*log(rectime2[j]) + betaui*log(i+1.0);
        lliki=0.0;
        if(numrec2[j]<maxc){
          lliki += dpois(numrec2[j],exp(temp),give_log1);
        }else{
	  lliki += log(1.0-ppois(maxc-1.0,exp(temp),give_log0,give_log0));
        }
        if(srd2[j]>=0){
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
          if(srd2[i] <= 10*Ki){
            lliki += log(pd2[srd2[j]-1]);
          }else{
            lliki += log(pis2);
          }
        }
        pd[i]=pi[i]*exp(lliki)*(i+1.);
      }
      // Set up pd to be cumulative for the random draws
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
      u2[j]=sizei;

     }
    }

     /* Deal with the outliers */
     umax=nk[Ki-1];
     sizei=Ki;
     for (j=0; j<n2; j++){
      if(srd2[j] > 10*Ki){
       while((umax==0) & (sizei > 1)){
         sizei--;
         umax=nk[sizei-1];
       }
       u2[j]=sizei;
       umax--;
      }
     }

     itemp=0;
     for (j=0; j<n2; j++){
      if(rc[j]==0){
        nk[u2[j]-1]=nk[u2[j]-1]+1;
        u1[n1 + itemp]=u2[j];
        itemp++;
      }
     }

     itemp=0;
     for (i=0; i<n1; i++){
       itemp+=u1[i];
     }
     for (i=0; i<n2; i++){
       if(rc[i]!=0) itemp-=u2[i];
     }
     unrecap=itemp;
     // Rebuild b2
     if(n2 > 0){
       b2[n2-1]=u2[n2-1];
       for (i=(n2-2); i>=0; i--){
         b2[i]=b2[i+1]+u2[i];
       }
     }
     /* End of imputed unit sizes for observed in second RDS sample */

     /* End of overall loop started at "Draw new beta using MCMC" */

    /* Draw phis */
    tU1=0;
    for (i=n1; i<Ni; i++){
      tU1+=u1[i];
    }
    tU2=unrecap;
    for (i=n2; i<Ni; i++){
      tU2+=u2[i];
    }
    r1=0.;
    for (i=0; i<n1; i++){
      r1+=exp_rand()/(tU1+b1[i]);
    }
    r2=0.;
    for (i=0; i<n2; i++){
      r2+=exp_rand()/(tU2+b2[i]);
    }

    /* Draw new N */

    gammart=0.;
    for (i=0; i<Ki; i++){
      gammart+=(exp(-(r1+r2)*(i+1))*pi[i]);
    }
    gammart=log(gammart);
    temp = -100000000.0;
    // N = m + n
    // Compute (log) P(m | \theta and data and \Psi)
    for (i=0; i<imaxm; i++){
      lpm[i]=i*gammart+lgamma(n+i+1.)-lgamma(i+1.);
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
    Ni += n;
    if(Ni > imaxN) Ni = imaxN;


    /* Draw unseen sizes */
    for (i=0; i<Ki; i++){
      Nk[i]=nk[i];
    }
    // Set up pi to be cumulative for random draws
    for (i=1; i<Ki; i++){
      pi[i]=pi[i-1]+pi[i];
    }
    for (i=n; i<Ni; i++){
      /* Propose unseen size for unit i */
      /* Use rejection sampling */
      sizei=1000000;
      while(sizei >= Ki){
       sizei=1000000;
       while(log(1.0-unif_rand()) > -(r1+r2)*sizei){
        /* Now propose unseen size for unit i */
        /* In the next two lines a sizei is chosen */
        /* with parameters lnlami and nui */
        temp = unif_rand();
        for (sizei=1; sizei<=Ki; sizei++){
          if(temp <= pi[sizei-1]) break;
        }
       }
      }
      u1[i]=sizei;
      Nk[sizei-1]=Nk[sizei-1]+1;
    }
    for (i=n; i<Ni; i++){
      /* Propose unseen size for unit i */
      /* Use rejection sampling */
      sizei=1000000;
      while(sizei >= Ki){
       sizei=1000000;
       while(log(1.0-unif_rand()) > -(r1+r2)*sizei){
        /* Now propose unseen size for unit i */
        /* In the next two lines a sizei is chosen */
        /* with parameters lnlami and nui */
        temp = unif_rand();
        for (sizei=1; sizei<=Ki; sizei++){
          if(temp <= pi[sizei-1]) break;
        }
       }
      }
      u2[i]=sizei;
    }
    if (step > 0 && step % iinterval == 0) {
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
      sample[isamp*dimsample+7]=betaui;
      sample[isamp*dimsample+8]=lmemmui;
      sample[isamp*dimsample+9]=memnui;
      for (i=0; i<Np; i++){
        sample[isamp*dimsample+10+i]=pdegi[i];
      }
      for (i=0; i<Ki; i++){
        Nkpos[i]=Nkpos[i]+Nk[i];
        posu[i]+=(((double)(Nk[i]))/Nd);
      }
      for (i=0; i<n1; i++){
        vsample[isamp*n1+i]=u1[i];
        if((srd[i]>0) && (srd[i] <= 10*Ki)){ posd[srd[i]-1]+=(1./Nd); }
      }
      for (i=0; i<n2; i++){
        vsample2[isamp*n2+i]=u2[i];
      }
//    Record the predicted degrees (not unit sizes)
      for (i=n1; i<Ni; i++){
        memmui = exp(lmemmui)*u1[i];
        rnb=memmui/(alpha-1.);
        pd2[0]= exp(-fabs(memmui-1.)/sqrt(memnui));
        pis=pd2[0];
        for (k=1; k<10*Ki; k++){
          pd2[k]= pd2[k-1]*(k+rnb)*pnb*exp((fabs(k-memmui)-fabs(k+1-memmui))/sqrt(memnui))/((double)(k+1));
          pis+=pd2[k];
        }
        for (k=0; k<10*Ki; k++){
          pd2[k]/=pis;
          posd[k]+=(pd2[k]/Nd);
        }
      }
      isamp++;
      if (*verbose && isamplesize % isamp == 0) Rprintf("Taken %d samples...\n", isamp);
    }
    step++;
  }
  dsamp=((double)isamp);
  for (i=0; i<Ki; i++){
    nk[i]=Nkpos[i];
    posu[i]/=dsamp;
  }
  for (i=0; i<Ki; i++){
 Rprintf("i %d posu[i] %f nk[i] %d Nk[i] %d\n", i, posu[i], nk[i], Nk[i]);
  }
  for (i=0; i<n1; i++){
     pop12[i]=u1[i];
  }
  for (i=0; i<n2; i++){
     pop21[i]=u2[i];
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(pi);
  free(pd);
  free(pdm);
  free(pd2);
  free(u1);
  free(u2);
  free(b1);
  free(b2);
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
  free(betausample);
  free(lmemmusample);
  free(memnusample);
}

void MHwpmem2 (int *u1, int *u2, int *n1i, int *n2i, int *K,
            double *beta0, double *beta0s, double *betat, double *betats, double *betau, double *betaus,
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            double *memod,
            int *srd,
            int *numrec,
            double *rectime,
            int *srd2,
            int *numrec2,
            double *rectime2,
            int *maxcoupons,
            double *beta0proposal, double *betatproposal, double *betauproposal,
            double *lmemmuproposal, double *memnuproposal,
            double *beta0sample, double *betatsample, double *betausample,
            double *lmemmusample, double *memnusample,
            int *samplesize, int *staken, int *warmupbeta, int *interval,
            int *verbose
         ) {
  int Ki, maxc, n1, n2;
  int step, taken, give_log1=1, give_log0=0;
  int i, k, isamp, iinterval, isamplesize, iwarmupbeta;
  double ip, cutoff;
  double temp;
  double beta0star, betatstar, beta0i, betati;
  double betaustar, betaui;
  double qi, qstar, lliki, llikstar;
  double lmemmustar, lmemmui, memmui, memnui;
  double memmustar, memnustar;
  double rmemnui, rmemnustar;
  double pibeta0star, pibeta0i;
  double pibetatstar, pibetati;
  double pibetaustar, pibetaui;
  double pimemstar, pimemi;
  double dbeta0, dbeta0s, dbetat, dbetats;
  double dbetau, dbetaus;
  double dlmemmu, dmemdfmu, dmemdfnu, rmemdfmu;
  double dmemnu, dmemnur;
  double dbeta0proposal, dbetatproposal, dbetauproposal;
  double dlmemmuproposal, dmemnuproposal;
  double pis, pis2;
  double alpha, pnb, rnb;

  Ki=(*K);
  double *pd = (double *) malloc(sizeof(double) * (10*Ki));

//GetRNGstate();  /* R function enabling uniform RNG */

  isamplesize=(*samplesize);
  iinterval=(*interval);
  iwarmupbeta=(*warmupbeta);
  dbeta0=(*beta0);
  dbeta0s=(*beta0s);
  dbetat=(*betat);
  dbetats=(*betats);
  dbetau=(*betau);
  dbetaus=(*betaus);
  dlmemmu=(*lmemmu);
  dmemnu=(*memnu);
  dmemnur=sqrt(dmemnu);
  dmemdfmu=(*memdfmu);
  rmemdfmu=sqrt(dmemdfmu);
  dmemdfnu=(*memdfnu);
  dbeta0proposal=(*beta0proposal);
  dbetatproposal=(*betatproposal);
  dbetauproposal=(*betauproposal);
  dlmemmuproposal=(*lmemmuproposal);
  dmemnuproposal=(*memnuproposal);
  alpha=(*memod);

  pnb=(alpha-1.)/alpha;

  // First set starting values
  isamp = taken = 0;
  step = -iwarmupbeta;
  n1 =(*n1i);
  n2 =(*n2i);
  maxc=(*maxcoupons);
  beta0i = beta0sample[0];
  betati = betatsample[0];
  betaui = betausample[0];
  lmemmui = lmemmusample[0];
  memnui = memnusample[0];
  rmemnui  = sqrt(memnui);

  // Compute initial current lik
  lliki = 0.0;
  for (i=0; i<n1; i++){
     temp = beta0i + betati*log(rectime[i]) + betaui*log(u1[i]);
     if(numrec[i]<maxc){
       lliki += dpois(numrec[i],exp(temp),give_log1);
     }else{
       lliki += log(1.0-ppois(maxc-1.0,exp(temp),give_log0,give_log0));
     }
     if(srd[i]>=0){
//    Use WP localized
      memmui = exp(lmemmui)*(u1[i]);
      rnb=memmui/(alpha-1.);
      pd[0]= exp(-fabs(memmui-1.)/sqrt(memnui));
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
      lliki += log((double)(u1[i]));
    }
  }
  for (i=0; i<n2; i++){
     temp = beta0i + betati*log(rectime2[i]) + betaui*log(u2[i]);
     if(numrec2[i]<maxc){
       lliki += dpois(numrec2[i],exp(temp),give_log1);
     }else{
       lliki += log(1.0-ppois(maxc-1.0,exp(temp),give_log0,give_log0));
     }
     if(srd2[i]>=0){
//    Use WP localized
      memmui = exp(lmemmui)*(u2[i]);
      rnb=memmui/(alpha-1.);
      pd[0]= exp(-fabs(memmui-1.)/sqrt(memnui));
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
      if(srd2[i] <= 10*Ki){
        lliki += log(pd[srd2[i]-1]);
      }else{
        lliki += log(pis2);
      }
      lliki += log((double)(u2[i]));
     }
  }
  if(!isfinite(lliki)) lliki = -100000.0;

  // Compute initial prior
  pibeta0i = dnorm(beta0i, dbeta0, dbeta0s, give_log1);
  if(dbetats > 0.0) pibetati = dnorm(betati, dbetat, dbetats, give_log1);
  if(dbetaus > 0.0) pibetaui = dnorm(betaui, dbetau, dbetaus, give_log1);
  pimemi = dnorm(lmemmui, dlmemmu, rmemnui/rmemdfmu, give_log1);
  pimemi = pimemi+dsclinvchisq(memnui, dmemdfnu, dmemnu);

//qi = dnorm(log(memnui/memnui)/dmemnuproposal,0.,1.,give_log1)
//     -log(dmemnuproposal*memnui);

  // Now do the MCMC updates (starting with the warmup updates)
  while (isamp < isamplesize && step < iwarmupbeta) {
    /* Propose new beta */
    beta0star = rnorm(beta0i, dbeta0proposal);
    if(dbetats > 0.0){
      betatstar = rnorm(betati, dbetatproposal);
    }else{
      betatstar = betati;
    }
    if(dbetaus > 0.0){
      betaustar = rnorm(betaui, dbetauproposal);
    }else{
      betaustar = betaui;
    }
    /* Propose new memnu and lmemmu */
    lmemmustar = rnorm(lmemmui, dlmemmuproposal);
    // VIP Remember this next line hold the optimism fixed at
    // its initial value !!! VIP
    lmemmustar = lmemmui;

//  memnustar = rnorm(memnui, dmemnuproposal);
    memnustar = memnui*exp(rnorm(0., dmemnuproposal));
    rmemnustar = sqrt(memnustar);

    llikstar = 0.0;
    for (i=0; i<n1; i++){
       temp = beta0star + betatstar*log(rectime[i]) + betaustar*log(u1[i]);
       if(numrec[i]<maxc){
         llikstar += dpois(numrec[i],exp(temp),give_log1);
       }else{
         llikstar += log(1.0-ppois(maxc-1.0,exp(temp),give_log0,give_log0));
       }
       if(srd[i]>=0){
//      Use WP localized
        memmustar = exp(lmemmustar)*(u1[i]);
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
        llikstar += log((double)(u1[i]));

       }
    }
    for (i=0; i<n2; i++){
       temp = beta0star + betatstar*log(rectime2[i]) + betaustar*log(u2[i]);
       if(numrec2[i]<maxc){
         llikstar += dpois(numrec2[i],exp(temp),give_log1);
       }else{
         llikstar += log(1.0-ppois(maxc-1.0,exp(temp),give_log0,give_log0));
       }
       if(srd2[i]>=0){
//      Use WP localized
        memmustar = exp(lmemmustar)*(u2[i]);
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
        if(srd2[i] <= 10*Ki){
          llikstar += log(pd[srd2[i]-1]);
        }else{
          llikstar += log(pis2);
        }
        llikstar += log((double)(u2[i]));
       }
    }
    if(!isfinite(llikstar)) llikstar = -100000.0;

    /* Calculate pieces of the prior. */
    pibeta0star = dnorm(beta0star, dbeta0, dbeta0s, give_log1);
    if(dbetats > 0.0) pibetatstar = dnorm(betatstar, dbetat, dbetats, give_log1);
    if(dbetaus > 0.0) pibetaustar = dnorm(betaustar, dbetau, dbetaus, give_log1);
    rmemnustar  = sqrt(memnustar);
    pimemstar = dnorm(lmemmustar, dlmemmu, rmemnustar/rmemdfmu, give_log1);
    pimemstar = pimemstar+dsclinvchisq(memnustar, dmemdfnu, dmemnu);

    qi = dnorm(log(memnui/memnustar)/dmemnuproposal,0.,1.,give_log1)
                  -log(dmemnuproposal*memnui);

    qstar = dnorm(log(memnustar/memnui)/dmemnuproposal,0.,1.,give_log1)
                  -log(dmemnuproposal*memnustar);

    /* Calculate ratio */
    ip =      pibeta0star-pibeta0i;
    if(dbetats > 0.0) ip = ip + pibetatstar-pibetati;
    if(dbetaus > 0.0) ip = ip + pibetaustar-pibetaui;
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
      betaui = betaustar;
      lmemmui = lmemmustar;
      memnui = memnustar;
      rmemnui = rmemnustar;
      lliki = llikstar;
      qi = qstar;
      pibeta0i = pibeta0star;
      if(dbetats > 0.0) pibetati = pibetatstar;
      if(dbetaus > 0.0) pibetaui = pibetaustar;
      pimemi = pimemstar;
      taken++;
      if (step > 0 && step==(iinterval*(step/iinterval))) {
        /* record statistics for posterity */
        beta0sample[isamp]=beta0i;
        betatsample[isamp]=betati;
        betausample[isamp]=betaui;
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
