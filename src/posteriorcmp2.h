#ifndef POSTERIORCMP2_H
#define POSTERIORCMP2_H

void gcmp2 (int *pop12,
            int *pop21, 
            int *nk, 
            int *K, 
            int *n1, 
            int *n2, 
            int *n12, 
            int *samplesize, int *warmup, int *interval,
            double *mu, double *dfmu, 
            double *sigma, double *dfsigma,
	    double *lnlam, double *nu,
            int *Npi,
            double *lnlamproposal, 
            double *nuproposal, 
            int *N, int *maxN, 
            double *sample, 
            double *posu, 
            double *lpriorm, 
            int *warmuptheta,
            int *verbose
           );
#endif /* POSTERIORCMP2_H */
