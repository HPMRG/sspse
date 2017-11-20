#include <stdio.h>
#include <float.h>
#include <R.h>

void ppspolya(double *y, double *size, double *w, int *nin, int *Nin)
{
	int i, j;
	int n = nin[0];
	int N = Nin[0];
        double a;

	GetRNGstate();

	for (i=n; i<N; i++) {
		a = w[i-1] * unif_rand();
                for (j=0; j<N; j++){
                  if(a <= w[j]) break;
                }
		y[i] = y[j];
		size[i] = size[j];
	}

	PutRNGstate();
}
