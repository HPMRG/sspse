/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "getinclCstacked.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void getinclCstacked (
            int *nbyclass,
            double *size, 
            int *K, 
            int *n, 
            int *samplesize,
            int *Nk
                 ) {
  int i, ni, Ki, isamp, isamplesize;
  // nbyclass = the number of members of class k=1,...,K
  // size = the size (i.e., degree) of the kth class k=1,...,K
  // K = number of classes (of degrees)
  // n = RDS sample size
  // samplesize = number of w.o.replacement samples to take
  // Nk = output: the total number of times a member of class k=1:K is sampled.

  GetRNGstate();  /* R function enabling uniform RNG */

  ni=(*n);
  Ki=(*K);
  isamplesize=(*samplesize);
  // ni = RDS sample size
  // Ki = number of classes (of degrees)
  // isamplesize = number of w.o.replacement samples to take

  int *perm = (int *) malloc(sizeof(int) * Ki);
  int *tperm = (int *) malloc(sizeof(int) * Ki);
  double *tsize = (double *) malloc(sizeof(double) * Ki);
  int *tnbyclass = (int *) malloc(sizeof(int) * Ki);
  int *samp = (int *) malloc(sizeof(int) * ni);

  for (i=0; i<Ki; i++){
     Nk[i]=0;
  }
  /* Record element identities */
  for (i = 0; i < Ki; i++)
    perm[i] = i + 1;

  /* Sort probabilities into descending order */
  /* Order element identities in parallel */
  /* perm is the permutation order of the ith element of the pop */
  revsort(size, perm, Ki);
  /* Order element nbyclass also */
  for(i = 0 ; i < Ki ; i++){
   tnbyclass[i]=nbyclass[i];
  }
  for(i = 0 ; i < Ki ; i++){
    nbyclass[i]=tnbyclass[perm[i]-1];
  }
  for(isamp = 0 ; isamp < isamplesize ; isamp++){
    /* Draw new sample */
    for(i = 0 ; i < Ki ; i++){
      tnbyclass[i]=nbyclass[i];
      tsize[i]=size[i];
      tperm[i]=perm[i];
    }
  /* Sample ni from population with Ni=sum(nbyclass) elements with the prob */
  /* of the ith pop in tsize[i] (ordered in descending order given */
  /* by the permutation in perm */
    ProbSampleNoReplaceStacked(Ki, tnbyclass, tsize, tperm, ni, samp);
//Rprintf("isamp %d %d %d %d\n",isamp,Nk[0],Nk[1],Nk[2]);
//Rprintf("%d %d %d %d %d %d\n",samp[0],samp[1],samp[2],samp[3],samp[4],samp[5]);
    /* Tabulate */
    for(i = 0 ; i < ni ; i++){
//    Nk[pop[samp[i]-1]]++; // edit here as original pop used here
      Nk[samp[i]-1]++;
    }
//  if (isverbose && isamplesize==(10*isamp*(isamplesize/(10*isamp)))){
//    Rprintf("Taken %d samples...\n", isamp);
//  }
  }
  PutRNGstate();  /* Disable RNG before returning */
  free(samp);
  free(tsize);
  free(tnbyclass);
  free(tperm);
  free(perm);
}
static void ProbSampleNoReplaceStacked(int n, int *nbyclass, double *p, int *perm, int nans, int *ans)
{
  // n = number of classes (of degrees)
  // nbyclass = the number of members of class k=1,...,K
  // p = class of ith member of the pop i=1,...,N i.e. degree
  // perm = permutation of class in descending order k=1:n
  // nans = RDS sample size
  // ans = sample values drawn in sequential order i=1:nans
    double rT, mass, totalmass;
    int i, j, nby;

    /* Record element identities */
//  for (i = 0; i < n; i++)
//    perm[i] = i + 1;

    /* Sort probabilities into descending order */
    /* Order element identities in parallel */
//  revsort(p, perm, n);

    /* Compute the sample */
    totalmass = 1.0;
    for (i = 0; i < nans; i++) {
        rT = totalmass * unif_rand();
        mass = 0.0;
        for (j = 0; j < n; j++) {
            mass += p[j];
            if (rT <= mass)
                break;
        }
        ans[i] = perm[j];
//Rprintf("j %d perm[j] %d p[j] %f nbyclass[j] %d totalmass %f rT %f n %d\n",j,perm[j],
//p[j], nbyclass[j], totalmass,rT, n);
//        /* update the reduced probabilities */
        totalmass -= (p[j] / nbyclass[j]);
        p[j] *= (1.0-1.0/nbyclass[j]);
        nbyclass[j]--;
        if(j < n - 1 && p[j] < p[j+1]){
          perm[j] = perm[j+1];
          perm[j+1] = ans[i];
          mass = p[j];
          p[j] = p[j+1];
          p[j+1] = mass;
          nby = nbyclass[j];
          nbyclass[j] = nbyclass[j+1];
          nbyclass[j+1] = nby;
//Rprintf("j %d perm[j] %d perm[j+1] %d p[j] %f p[j+1] %f nbyclass[j] %d totalmass %f rT %f\n",j,perm[j],
//perm[j+1], p[j], p[j+1], nbyclass[j], totalmass,rT);
        }
//Rprintf("i %d ans %d\n",i,ans[i]);
    }
}
