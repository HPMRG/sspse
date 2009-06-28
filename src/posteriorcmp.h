#ifndef POSTERIORCMP_H
#define POSTERIORCMP_H

void dcmp_mu (int *x, double *mu, double *sig, int *n, double *err, int *give_log, double *val);

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

double cmp(int x, double llambda, double nu, double lzcmp, int give_log);
void dcmp (int *x, double *lambda, double *nu, int *n, double *err, int *give_log, double *val);
void gcmp (int *pop,
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
	    double *ppos,
	    double *lpriorm,
            int *burnintheta,
	    int *fVerbose
			 );
void MHcmp (int *Nk, int *K,
	    double *mu0, double *kappa0, 
            double *sigma0,  double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Npi, double *psample,
            double *musample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *fVerbose
			 );
double zcmp(double lambda, double nu, double err, int give_log);
void vzcmp(double *lambda, double *nu, double *err, int *give_log, double *out);
void rcmp (int *x, double *lambda, double *nu, int *n, int *K, double *err);
#endif /* POSTERIORCMP_H */
