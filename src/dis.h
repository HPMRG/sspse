#ifndef DIS_H
#define DIS_H

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
void gsppsdis (int *pop, int *dis,
            int *nk0, int *nk1, 
            int *K, 
            int *n, 
            int *samplesize, int *burnin, int *interval,
            double *mu0, double *mu1, double *kappa0, 
            double *sigma0,  double *sigma1, double *df0,
	    int *Np0i, int *Np1i,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *maxN, 
            double *sample, 
            double *p0pos, double *p1pos, 
            int *burnintheta,
	    int *verbose
			 );
void MHdis (int *Nk0, int *Nk1, int *totdis, int *K,
	    double *mu0, double *mu1, double *kappa0, 
            double *sigma0,  double *sigma1, double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            int *N, int *Np0i, int *Np1i, double *psample,
            double *musample, double *betasample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *verbose
			 );
void MHpriordis (double *mu0, double *mu1, double *kappa0, 
            double *sigma0,  double *sigma1, double *df0,
            double *muproposal, 
            double *sigmaproposal, 
            double *musample, double *betasample, double *sigmasample,
            int *samplesize, int *staken, int *burnin, int *interval,
	    int *verbose
			 );
#endif /* DIS_H */
