#ifndef POSTERIORCMPVIS_H
#define POSTERIORCMPVIS_H

void gcmpvis (int *pop,
            int *nk, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu, double *dfmu, 
            double *sigma, double *dfsigma,
            double *beta0muprior, double *beta0sigmaprior, 
            double *beta1muprior, double *beta1sigmaprior, 
            double *memlnopt, double *memdflnopt,
            double *memnu, double *memdfnu,
            int *Npi,
            int *srd, 
            double *numrec, 
            double *rectime,
            int *maxcoupons,
            double *lnlamproposal, 
            double *nuproposal, 
            double *beta0proposal, double *beta1proposal, 
            double *memlnoptproposal, double *memnuproposal,
            int *N, int *maxN, 
            double *sample, 
            int *vsample, 
            double *ppos, 
            double *lpriorm, 
            int *burnintheta,
            int *burninbeta,
            int *verbose
         );

void MHcmpmem (int *d, int *n, int *K,
            double *beta0, double *beta0sd, double *beta1, double *beta1sd, 
            double *memlnopt, double *memdflnopt,
            double *memnu, double *memdfnu,
            int *srd, 
            double *numrec, 
            double *rectime,
            int *maxcoupons,
            double *beta0proposal, double *beta1proposal, 
            double *memlnoptproposal, double *memnuproposal, 
            double *beta0sample, double *beta1sample,
            double *memlnoptsample, double *memnusample,
            int *samplesize, int *staken, int *burnin, int *interval,
            int *verbose
         );

void MHcmptheta (int *Nk, int *K,
            double *mu, double *dfmu, 
            double *sigma,  double *dfsigma,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Npi, double *psample,
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
            int *verbose
         );

#endif /* POSTERIORCMPVIS_H */
