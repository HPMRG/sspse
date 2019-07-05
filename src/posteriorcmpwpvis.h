#ifndef POSTERIORCMPWPVIS_H
#define POSTERIORCMPWPVIS_H

void gcmpwpvis (int *pop,
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu, double *dfmu, 
            double *sigma, double *dfsigma,
	    double *lnlam, double *nu,
            double *beta0muprior, double *beta0sigmaprior, 
            double *beta1muprior, double *beta1sigmaprior, 
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            double *memod,
            int *Npi,
            int *srd, 
            int *numrec, 
            double *rectime,
            int *maxcoupons,
            double *lnlamproposal, 
            double *nuproposal, 
            double *beta0proposal, double *beta1proposal, 
            double *lmemmuproposal, double *memnuproposal,
            int *N, int *maxN, 
            double *sample, 
            int *vsample, 
            double *posu, 
            double *posd, 
            double *lpriorm, 
            int *burnintheta,
            int *burninbeta,
            int *verbose
         );

void MHwpmem (int *d, int *n, int *K,
            double *beta0, double *beta0sd, double *beta1, double *beta1sd, 
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            int *srd, 
            int *numrec, 
            double *rectime,
            int *maxcoupons,
            double *beta0proposal, double *beta1proposal, 
            double *lmemmuproposal, double *memnuproposal, 
            double *beta0sample, double *beta1sample,
            double *lmemmusample, double *memnusample,
            int *samplesize, int *staken, int *burnin, int *interval,
            int *verbose
         );

void MHcmptheta (int *Nk, int *K,
            double *mu, double *dfmu, 
            double *sigma,  double *dfsigma,
            double *lnlamproposal, 
            double *nuproposal, 
            int *N, int *Npi, double *psample,
            double *lnlamsample, double *nusample,
            int *samplesize, int *staken, int *burnintheta, int *interval,
            int *verbose
         );

#endif /* POSTERIORCMPWPVIS_H */
