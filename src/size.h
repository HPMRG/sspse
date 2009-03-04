#ifndef SIZE_H
#define SIZE_H

void bnw_llik(int *K, int *n, int *s, int *nk, double *Nk, double *llik);
void bnw_s_llik(int *K, int *n, int *s, double *Nk, double *llik);
void dwarC(int *N, double *mu, double *rho, double *pmf);
double ldwar(double *N, double *mu, double *rho);
void bnw_unpos(int *K, int *n, int *s, int *nk, double *Nk, double *mu,
               double *rho, double *unpos);
double bnw_llikN(int *K, int *n, int *s, int *nk, int *Nk);
double bnw_llikNf(int *K, int *n, int *s, int *nk, int *Nk);
double dmultinorm(int *N, int *K, int *Nk, double *lprob);
double bnw_unposN(int *N, int *K, int *n, int *s, int *nk, int *Nk, double *lprob);
void bnw_NCbound(int *N, int *K, int *n, int *s, int *nk, int *Nk,
	    double *prob,
	    double *qprob,
  	    int *M, double *unpos);
void bnw_stocdiscrete(int *N, int *K, int *n, int *s, int *nk, int *Nk,
	    double *qprob,
  	    int *M, double *mllik);
void bnw_stocdiscreteimpute(int *N, int *K, int *n, int *s, int *nk, int *Nk,
	    double *qprob,
  	    int *nsim, int *M, double *mllik);
void bnw_NC(int *N, int *K, int *n, int *s, int *nk, int *Nk,
	    double *prob,
	    double *qprob,
  	    int *M, double *unpos);
void bnw_NCwar(int *N, int *K, int *n, int *s, int *nk, int *Nk,
            double *prob,
            double *mu,
            double *rho, int *M, double *unpos);
void bnw_mpwar(int *N, int *lenN, int *K, int *n, int *s, int *nk,
	    double *lbound,
	    double *prob,
	    double *NtotMLE,
	    int *Nprior,
	    int *Nmle,
	    double *mu,
	    double *rho, int *M);
void bnw_mp(int *N, int *lenN, int *K, int *n, int *s, int *nk,
	    double *lbound,
	    double *dprob,
	    double *prob,
	    double *NtotMLE,
	    int *Nprior,
	    int *Nmle,
	    int *M);

double ldwarint(int *N, double *mu, double *rho);
#endif /* SIZE_H */
