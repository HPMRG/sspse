#ifndef POSTERIORPLN_H
#define POSTERIORPLN_H

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
void gpln (int *pop,
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
	    double *lpriorm,
            int *warmuptheta,
	    int *verbose
			 );
void MHpln (int *Nk, int *K,
	    double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Npi, double *psample,
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *warmup, int *interval,
	    int *verbose
			 );
#endif /* POSTERIORPLN_H */
