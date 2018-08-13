#ifndef POSTERIORCMPVIS_H
#define POSTERIORCMPVIS_H

void gcmpvis (int *pop,
            int *nk, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu, double *dfmu, 
            double *sigma,  double *dfsigma,
            double *beta0muprior, double *beta0sigmaprior, 
            double *beta1muprior, double *beta1sigmaprior, 
            double *memmu,  double *memdfmu,
            double *memsd,  double *memdfsd,
            int *Npi,
            int *srd, 
            double *numrec, 
            double *rectime,
            int *maxcoupons,
            double *muproposal, 
            double *sigmaproposal, 
            double *beta0proposal, double *beta1proposal, 
            double *memmuproposal, double *memsdproposal,
            int *N, int *maxN, 
            double *sample, 
            int *vsample, 
            double *ppos, 
            double *lpriorm, 
            int *burnintheta,
            int *burninbeta,
            int *verbose
         );

void MHcmpbeta (int *pop, int *n, int *K,
            double *beta0, double *beta0s, double *beta1, double *beta1s, 
            double *memmu, double *memdfmu,
            double *memsd, double *memdfsd,
            int *srd, 
            double *numrec, 
            double *rectime,
            int *maxcoupons,
            double *beta0proposal, double *beta1proposal, 
            double *memmuproposal, double *memsdproposal, 
            double *beta0sample, double *beta1sample,
            double *memmusample, double *memsdsample,
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
