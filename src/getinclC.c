/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "getinclC.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void getinclC (int *N, 
            int *pop,
            double *size, 
            int *K, 
            int *n, 
            int *samplesize,
            int *Nk
		 ) {
  int i, ni, Ni, Ki, isamp, isamplesize;
  // N = populaton size
  // pop = class order of ith member of the pop i=1,...,N
  // size = class of ith member of the pop i=1,...,N i.e. degree
  // K = number of classes (of degrees)
  // n = RDS sample size
  // samplesize = number of w.o.replacement samples to take
  // Nk = class values k=1:K 

  GetRNGstate();  /* R function enabling uniform RNG */

  ni=(*n);
  Ni=(*N);
  Ki=(*K);
  isamplesize=(*samplesize);
  // ni = RDS sample size
  // Ni = populaton size
  // Ki = number of classes (of degrees)
  // isamplesize = number of w.o.replacement samples to take

  int *perm = (int *) malloc(sizeof(int) * Ni);
  int *tperm = (int *) malloc(sizeof(int) * Ni);
  double *tsize = (double *) malloc(sizeof(double) * Ni);
  int *samp = (int *) malloc(sizeof(int) * ni);

  for (i=0; i<Ki; i++){
     Nk[i]=0;
  }
  /* Record element identities */
  for (i = 0; i < Ni; i++)
    perm[i] = i + 1;

  /* Sort probabilities into descending order */
  /* Order element identities in parallel */
  /* perm is the permutation order of the ith element of the pop */
  revsort(size, perm, Ni);

  for(isamp = 0 ; isamp < isamplesize ; isamp++){
    /* Draw new sample */
    for(i = 0 ; i < Ni ; i++){
      tsize[i]=size[i];
      tperm[i]=perm[i];
    }
  /* Sample ni from population with Ni elements with the prob */
  /* of the ith pop in tsize[i] (ordered in descending order given */
  /* by the permutation in perm */
    ProbSampleNoReplace(Ni, tsize, tperm, ni, samp);
//Rprintf("isamp %d %d %d %d\n",isamp,Nk[0],Nk[1],Nk[2]);
//Rprintf("%d %d %d %d %d %d\n",samp[0],samp[1],samp[2],samp[3],samp[4],samp[5]);
    /* Tabulate */
    for(i = 0 ; i < ni ; i++){
      Nk[pop[samp[i]-1]]++;
    }
//  if (isverbose && isamplesize==(10*isamp*(isamplesize/(10*isamp)))){
//    Rprintf("Taken %d samples...\n", isamp);
//  }
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(samp);
  free(tsize);
  free(perm);
  free(tperm);
}
static void ProbSampleNoReplace(int n, double *p, int *perm, int nans, int *ans)
{
    double rT, mass, totalmass;
    int i, j, k, n1;

    /* Record element identities */
//  for (i = 0; i < n; i++)
//    perm[i] = i + 1;

    /* Sort probabilities into descending order */
    /* Order element identities in parallel */
//  revsort(p, perm, n);

    /* Compute the sample */
    totalmass = 1;
    for (i = 0, n1 = n-1; i < nans; i++, n1--) {
	rT = totalmass * unif_rand();
	mass = 0;
	for (j = 0; j < n1; j++) {
	    mass += p[j];
	    if (rT <= mass)
		break;
	}
	ans[i] = perm[j];
	totalmass -= p[j];
	for(k = j; k < n1; k++) {
	    p[k] = p[k + 1];
	    perm[k] = perm[k + 1];
	}
//Rprintf("i %d ans %d\n",i,ans[i]);
    }
}
