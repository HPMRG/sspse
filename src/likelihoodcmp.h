#ifndef LIKELIHOODCMP_H
#define LIKELIHOODCMP_H

void lcmp (int *pop,
           int *nk,
           int *K,
           int *n,
           int *samplesize, int *warmup, int *interval,
           double *mu, double *kappa,
           double *sigma,  double *df,
           int *Npi,
           double *muproposal,
           double *sigmaproposal,
           int *N, int *maxN,
           double *sample,
           double *ppos,
           double *lpriorm,
           int *warmuptheta,
           int *verbose
			 );
void MHlcmp (int *Nk, int *K,
            double *mu, double *kappa, 
            double *sigma,  double *df,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Npi, double *psample,
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *warmup, int *interval,
	    int *verbose
			 );
#endif /* LIKELIHOODCMP_H */
