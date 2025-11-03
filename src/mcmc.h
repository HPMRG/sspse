#ifndef MCMC_H
#define MCMC_H

/* I define Scale-Inverse-Chi-Square distribution Y with a scale parameter "scale"
   and degrees of freedom "df" as

   Y=scale*df/X,
   where X~Chi-sq(df).
   
   Here, E[Y]=scale*df/(df-2)

   These functions evaluate the density and generate deviates from this distribution.
*/

#define dsclinvchisqboth(x,df,scale,give_log) (give_log ? (dchisq((df)*((double)(scale))/((double)(x)),df,1)+log(((double)(scale))*(df)/(double)((x)*(x)))) : (dchisq((df)*((double)(scale))/((double)(x)),df,0)*(df)*((double)(scale))/(double)((x)*(x))))

#define dsclinvchisq(x,df,scale) (dchisq((df)*((double)(scale))/((double)(x)),df,1)+log(((double)(scale))*(df)/(double)((x)*(x))))

#define rsclinvchisq(df,scale) ((scale)*(df)/(rchisq(df)))

double poilog(int x, double mu, double sig);
void gspps (int *pop,
            int *nk, 
            int *K, 
            int *n, 
            int *samplesize, int *warmup, int *interval,
            double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            int *Npi,
	    double *muproposal, 
            double *sigmaproposal, 
            int *N, int *maxN, 
            double *sample, 
	    double *ppos,
            int *warmuptheta,
	    int *verbose
			 );
void gsppsN (int *pop,
            int *nk, 
            int *K, 
            int *n, 
            int *samplesize, int *warmup, int *interval,
            double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            int *Npi,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, 
            int *sample, 
            int *warmuptheta,
	    int *verbose
			 );
void MHplnorig (int *Nk, int *K,
	    double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Npi, double *psample,
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *warmup, int *interval,
	    int *verbose
			 );
void MHpriorpln (double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *warmup, int *interval,
	    int *verbose
			 );
#endif /* MCMC_H */
