#ifndef MCMC_H
#define MCMC_H

/* I define Scale-Inverse-Chi-Square distribution Y with a scale parameter "scale"
   and degrees of freedom "df" as

   Y=scale*df/X,
   where X~Chi-sq(df).
   
   Here, E[Y]=scale*df/(df-2)

   These functions evaluate the density and generate deviates from this distribution.
*/

#define dsclinvchisq(x,df,scale,give_log) (give_log ? (dchisq((df)*((double)(scale))/((double)(x)),df,1)+log(((double)(scale))*(df)/(double)((x)*(x)))) : (dchisq((df)*((double)(scale))/((double)(x)),df,0)*(df)*((double)(scale))/(double)((x)*(x))))

#define rsclinvchisq(df,scale) ((scale)*(df)/(rchisq(df)))

double poilog(int x, double my, double sig);
void MetropolisHastings (double *s,
		         double *mu0, double *kappa0, 
                         double *sigma20,  double *df0,
                         double *sigmaproposal, 
			 int *N,
			 double *musample, double *sigmasample,
			 int *nsteps, int *staken, int *burnin, int *interval,
			 int fVerbose
			 );
#endif /* MCMC_H */
