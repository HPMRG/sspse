#ifndef SIZE_H
#define SIZE_H

void bnw_llik(int *K, int *n, int *s, int *nk, double *Nk, double *llik);
void bnw_s_llik(int *K, int *n, int *s, double *Nk, double *llik);
void dwarC(int *N, double *mu, double *rho, double *pmf);
double ldwar(double *N, double *mu, double *rho);
void bnw_unpos(int *K, int *n, int *s, int *nk, double *Nk, double *mu,
               double *rho, double *unpos);
double bnw_llikN(int *K, int *n, int *s, int *nk, int *Nk);
double dmultinorm(int *N, int *K, int *Nk, double *lprob);
double bnw_unposN(int *N, int *K, int *n, int *s, int *nk, int *Nk, double *lprob);
void bnw_NC(int *N, int *K, int *n, int *s, int *nk, int *Nk,
	    double *prob,
	    double *mu,
  	    double *rho, int *M, double *unpos);
void bnw_mp(int *N, int *lenN, int *K, int *n, int *s, int *nk,
	    double *Cval,
	    double *prob,
	    int *Nprior,
	    int *Nmle,
	    double *mu,
	    double *rho, int *M);

double ldwarint(int *N, double *mu, double *rho);

#endif /* SIZE_H */
