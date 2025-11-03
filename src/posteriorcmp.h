#ifndef POSTERIORCMP_H
#define POSTERIORCMP_H

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

void gcmp (int *pop,
            int *K, 
            int *n, 
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

void MHcmptheta (int *Nk, int *K,
            double *mu, double *dfmu, 
            double *sigma,  double *dfsigma,
            double *lnlamproposal, 
            double *nuproposal, 
            int *N, int *Npi, double *psample,
            double *lnlamsample, double *nusample,
            int *samplesize, int *staken, int *warmuptheta, int *interval,
            int *verbose
         );

#endif /* POSTERIORCMP_H */
