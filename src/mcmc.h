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

double poilog(int x, double my, double sig);
void gspps (int *pop,
            int *nk, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            int *Npi,
	    double *muproposal, 
            double *sigmaproposal, 
            int *N, int *maxN, 
            double *sample, 
            int *burnintheta,
	    int *fVerbose
			 );
void gsppsN (int *pop,
            int *nk, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            int *Npi,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, 
            int *sample, 
            int *burnintheta,
	    int *fVerbose
			 );
void MHpln (int *Nk, int *K,
	    double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Npi, double *psample,
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *fVerbose
			 );
void MHpriorpln (double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *fVerbose
			 );
void getinclC (int *N,
            int *pop,
            double *size, 
            int *K, 
            int *n, 
            int *samplesize,
            int *Nk,
	    int *fVerbose
			 );
static void ProbSampleNoReplace(int n, double *p, int *perm, int nans, int *ans);
#endif /* MCMC_H */
