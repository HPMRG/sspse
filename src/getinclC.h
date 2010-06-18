#ifndef GETINCLC_H
#define GETINCLC_H

void getinclC (int *N,
            int *pop,
            double *size, 
            int *K, 
            int *n, 
            int *samplesize,
            int *Nk
		 );
static void ProbSampleNoReplace(int n, double *p, int *perm, int nans, int *ans);
#endif /* GETINCLC_H */
