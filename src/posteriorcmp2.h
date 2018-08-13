#ifndef POSTERIORCMP2_H
#define POSTERIORCMP2_H

void gcmp2 (int *pop12,
           int *pop21,
           int *nk,
           int *K,
           int *n1,
           int *n2,
           int *n0,
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
           int *verbose
			 );
#endif /* POSTERIORCMP2_H */
