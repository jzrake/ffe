#include <math.h>
#include "ffe.h"


#ifndef M_PI
#define M_PI 3.1415926535897926
#endif


static int order = 0;
static int num_bins = 128;
static int array_size = 512;
static double Uintegral(double *Adata, double *Jdata, double *Sdata, double A, int N);


int ffe_make_nonlinear_eqilibrium(struct ffe_sim *sim)
{

  int N = array_size;
  
  double *Adata = (double*) malloc(N * N * sizeof(double));
  double *Jdata = (double*) malloc(N * N * sizeof(double));
  double *Sdata = (double*) malloc(N * N * sizeof(double));
  double Amax = 0.0;
  double Amin = 0.0;
  
  for (int i=0; i<N; ++i) {
    for (int j=0; j<N; ++j) {

#define n (2*order + 1)
#define c cos
#define s sin
#define p pow

      double x = 2 * M_PI * i / N;
      double y = 2 * M_PI * j / N;

      double A0 = c(x) - s(y);
      double A1 = fabs(A0) < 1e-14 ? 1.0 : A0;

      Adata[i*N+j] = pow(A0, n);
      Jdata[i*N+j] = -0.5*n*p(A1,-2+n)*(-4+2*n-n*c(2*x)+n*c(2*y)+4*c(x)*s(y));
      
      if (Adata[i*N+j] > Amax) Amax = Adata[i*N+j];
      if (Adata[i*N+j] < Amin) Amin = Adata[i*N+j];

#undef n
#undef c
#undef s
#undef p
      
    }
  }

  cow_histogram *dSdA_hist = cow_histogram_new();

  cow_histogram_setlower(dSdA_hist, 0, Amin);
  cow_histogram_setupper(dSdA_hist, 0, Amax);
  cow_histogram_setnbins(dSdA_hist, 0, num_bins);  
  cow_histogram_setspacing(dSdA_hist, COW_HIST_SPACING_LINEAR);
  cow_histogram_setbinmode(dSdA_hist, COW_HIST_BINMODE_DENSITY);
  cow_histogram_commit(dSdA_hist);
  
  for (int i=0; i<N; ++i) {
    for (int j=0; j<N; ++j) {
      cow_histogram_addsample1(dSdA_hist, Adata[i*N + j], 1.0);
    }
  }
  cow_histogram_seal(dSdA_hist);
  
  for (int i=0; i<N; ++i) {
    for (int j=0; j<N; ++j) {
      Sdata[i*N + j] = cow_histogram_sample1(dSdA_hist, Adata[i*N + j]);
    }
  }

  double *Atable = (double*) malloc(num_bins * sizeof(double));
  double *Btable = (double*) malloc(num_bins * sizeof(double));

  double U0 = Uintegral(Adata, Jdata, Sdata, 0.0, N);
  
  for (int n=0; n<num_bins; ++n) {

    double A = n * Amax / (num_bins-1);
    double U = Uintegral(Adata, Jdata, Sdata, A, N) - U0;
    
    Atable[n] = A;
    Btable[n] = sqrt(2 * U);
  }

  FILE *outf = fopen("test.dat", "w");
  for (int n=0; n<num_bins; ++n) {
    double A = Atable[n];
    double B = Btable[n];
    fprintf(outf, "%f %f\n", A, B);
  }
  fclose(outf);
  
  cow_histogram_del(dSdA_hist);

  free(Atable);
  free(Btable);
  free(Adata);
  free(Jdata);
  free(Sdata);
  
  return 0;
}



double Uintegral(double *Adata, double *Jdata, double *Sdata, double A, int N)
{

  double U = 0.0; /* \int_{0 < Az(S) < Az}{dS J / S'(A)} */
    
  for (int i=0; i<N; ++i) {
    for (int j=0; j<N; ++j) {

      if (Adata[i*N + j] < A) {
	U += Jdata[i*N +j] / Sdata[i*N + j];
      }
	
    }
  }

  return U;
}
