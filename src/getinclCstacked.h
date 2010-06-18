#ifndef GETINCLCSTACKED_H
#define GETINCLCSTACKED_H

void getinclCstacked (
            int *nbyclass,
            double *size, 
            int *K, 
            int *n, 
            int *samplesize,
            int *Nk
		 );
static void ProbSampleNoReplaceStacked(int n, int *nbyclass, double *p, int *perm, int nans, int *ans);
#endif /* GETINCLCSTACKED_H */
