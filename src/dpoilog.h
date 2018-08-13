   /*function prototypes */
   double poilog(int x, double my, double sig);
   double poilogmaxf(int x, double my, double sig);
   double poilogupper(int x, double m, double my, double sig);
   double poiloglower(int x, double m, double my, double sig);
   double poilogmy_f(double z, int x, double my, double sig, double fac);
   double poilogmy_f2(double z,int y,int x,double my1,double my2,double sig1,double sig2,double ro,double fac);
   void poilogmy_f_vec(double *z, int n, void *p);
   void poilogmy_f2_vec(double *z, int n, void *p);
