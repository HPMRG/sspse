#ifndef POSTERIORCMPVIS_H
#define POSTERIORCMPVIS_H

void gcmpvis (int *pop,
            int *K,
            int *n,
            int *samplesize, int *burnin, int *interval,
            double *mu, double *dfmu,
            double *sigma, double *dfsigma,
            double *lnlam, double *nu,
            double *beta0muprior, double *beta0sigmaprior,
            double *betatmuprior, double *betatsigmaprior,
            double *betaumuprior, double *betausigmaprior,
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            int *Npi,
            int *srd,
            int *numrec,
            double *rectime,
            int *maxcoupons,
            double *lnlamproposal,
            double *nuproposal,
            double *beta0proposal, double *betatproposal, double *betauproposal,
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

void MHcmpmem (int *u, int *n, int *K,
            double *beta0, double *beta0sd, double *betat, double *betatsd, double *betau, double *betausd,
            double *lmemmu, double *memdfmu,
            double *memnu, double *memdfnu,
            int *srd, 
            int *numrec, 
            double *rectime,
            int *maxcoupons,
            double *beta0proposal, double *betatproposal, double *betauproposal, 
            double *lmemmuproposal, double *memnuproposal, 
            double *beta0sample, double *betatsample, double *betausample,
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

#endif /* POSTERIORCMPVIS_H */
