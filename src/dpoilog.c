
   #include <math.h>
   #include <R.h>
   #include <R_ext/Utils.h>
   #include <R_ext/Applic.h>
   #include <Rinternals.h>


   /*function prototypes */
   double poilog(int x, double my, double sig);
   double bipoilog(int x, int y, double my1, double my2, double sig1, double sig2, double ro);
   double maxf(int x, double my, double sig);
   double upper(int x, double m, double my, double sig);
   double lower(int x, double m, double my, double sig);
   double my_f(double z, int x, double my, double sig, double fac);
   double my_f2(double z,int y,int x,double my1,double my2,double sig1,double sig2,double ro,double fac);
   void my_f_vec(double *z, int n, void *p);
   void my_f2_vec(double *z, int n, void *p);


/* ---------------------------------------------------------------------------*/

  double maxf(int x, double my, double sig)
  {
     double d,z;
     z=0;
     d=100;
     while (d>0.00001) {
       if (x-1-exp(z)-1/sig*(z-my)>0) z=z+d; else z=z-d;
       d=d/2;
     }
     return(z);
   }

/* ---------------------------------------------------------------------------*/


  double upper(int x, double m, double my, double sig)
  {
     double d,z,mf;
     mf = (x-1)*m-exp(m)-0.5/sig*((m-my)*(m-my));
     z = m+20;
     d = 10;
     while (d>0.000001) {
        if ((x-1)*z-exp(z)-0.5/sig*((z-my)*(z-my))-mf+log(1000000)>0) z=z+d; else z=z-d;
        d=d/2;
     }
     return(z);
   }


/* ---------------------------------------------------------------------------*/


  double lower(int x, double m, double my, double sig)
  {
     double d,z,mf;
     mf = (x-1)*m-exp(m)-0.5/sig*((m-my)*(m-my));
     z = m-20;
     d = 10;
     while (d>0.000001) {
        if ((x-1)*z-exp(z)-0.5/sig*((z-my)*(z-my))-mf+log(1000000)>0) z=z-d; else z=z+d;
        d=d/2;
     }
     return(z);
   }


/* ---------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/

  struct my_f_params { int x; double sig; double my; double fac;};
  struct my_f2_params { int x; int y; double sig1; double sig2; double my1; double my2;double ro; double fac;};
/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/

  double my_f(double z, int x, double my, double sig, double fac)
  {
     return exp(z*x-exp(z)-0.5/sig*((z-my)*(z-my))-fac);
  }


  double my_f2(double z,int y,int x,double my1,double my2,double sig1,double sig2,double ro,double fac)
  {
     return( poilog(y,my2+ro*sqrt(sig2/sig1)*(z-my1),sig2*(1-ro*ro)) *
             exp(x*z-exp(z)-fac-0.5/sig1*(z-my1)*(z-my1)) );
  }
/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/

  void my_f_vec(double *z, int n, void *p)
    {
       int i;
       struct my_f_params * params = (struct my_f_params *)p;
       int x      = (params->x);
       double sig = (params->sig);
       double my  = (params->my);
       double fac = (params->fac);
       for (i=0;i<n;i++) z[i]=my_f(z[i],x,my,sig,fac);
       return;
  }


  void my_f2_vec(double *z, int n, void * p)
  {
	 int i;
     struct my_f2_params * params = (struct my_f2_params *)p;
     int x       = (params->x);
     int y       = (params->y);
     double sig1 = (params->sig1);
     double sig2 = (params->sig2);
     double my1  = (params->my1);
     double my2  = (params->my2);
     double ro   = (params->ro);
     double fac  = (params->fac);
     for (i=0;i<n;i++) z[i]=my_f2(z[i],y,x,my1,my2,sig1,sig2,ro,fac);
     return;
   }

/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/


  double poilog(int x, double my, double sig)
  {
     double a,b,m,fac,val;
     double result, abserr;
     int last, neval, ier;
     int lenw;
     int *iwork;
     double *work;
     int limit=100;
     double reltol=0.00001;
     double abstol=0.00001;
     lenw = 4 * limit;
     iwork =   (int *) Calloc(limit, int);
     work = (double *) Calloc(lenw,  double);

     m=maxf(x,my,sig);
     a=lower(x,m,my,sig);
     b=upper(x,m,my,sig);
     fac = lgamma(x+1);

     struct my_f_params p = { x, sig, my, fac };

     Rdqags(my_f_vec, (void *) &p, &a, &b,
            &abstol,&reltol,&result,&abserr,&neval,&ier,
            &limit,&lenw, &last,iwork, work);

     if (ier!=0) error("error in integration\n");

     val = result*(1/sqrt(2*M_PI*sig));
     Free(iwork);
     Free(work);
     return(val);
  }


/* ---------------------------------------------------------------------------*/

  double bipoilog(int x, int y, double my1, double my2, double sig1, double sig2, double ro)
  {

     double a,b,m,fac,val;
     double result, abserr;
	 int last, neval, ier;
	 int lenw;
	 int *iwork;
	 double *work;
	 int limit=100;
	 double reltol=0.00001;
	 double abstol=0.00001;
	 lenw = 4 * limit;
	 iwork =   (int *) Calloc(limit, int);
     work = (double *) Calloc(lenw,  double);

     m=maxf(x,my1,sig1);
     a=lower(x,m,my1,sig1);
     b=upper(x,m,my1,sig1);
     fac = lgamma(x+1);

     struct my_f2_params p = { x, y, sig1, sig2, my1, my2, ro, fac};

     Rdqags(my_f2_vec, (void *) &p, &a, &b,
            &abstol,&reltol,&result,&abserr,&neval,&ier,
            &limit,&lenw, &last,iwork, work);

     if (ier!=0) error("error in integration\n");

     val = result*(1/sqrt(2*M_PI*sig1));
     Free(iwork);
     Free(work);

     return(val);

   }


/* ---------------------------------------------------------------------------*/
    void poilog2 (int *x, int *y, double *my1, double *my2, double *sig1, double *sig2, double *ro, int *nrN, double *val)
     {
       	int i;
        for (i = 0; i < *nrN; i++){
           val[i] = bipoilog(x[i],y[i],*my1,*my2,*sig1,*sig2,*ro);
           R_CheckUserInterrupt();
         }
      }

    void poilog1 (int *x, double *my, double *sig, int *nrN, double *val)
     {
       	int i;
        for (i = 0; i < *nrN; i++){
           val[i] = poilog(x[i],*my,*sig);
           R_CheckUserInterrupt();
         }
      }
/* ---------------------------------------------------------------------------*/


