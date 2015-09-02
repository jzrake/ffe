#include <math.h>
#include "ffe.h"


#ifndef M_PI
#define M_PI 3.1415926535897926
#endif


static double Uintegral(double *Adata, double *Jdata, double *Sdata, double A, int N);


void ffe_nle_null(struct ffe_nle *nle)
{
  nle->order = 0;
  nle->num_bins = 0;
  nle->array_size = 0;
  nle->Atable = NULL;
  nle->Btable = NULL;
}


void ffe_nle_init(struct ffe_nle *nle, int order, int num_bins, int array_size)
{
  int N = array_size;

  nle->order = order;
  nle->num_bins = num_bins;
  nle->array_size = array_size;
  
  double *Adata = (double*) malloc(N * N * sizeof(double));
  double *Jdata = (double*) malloc(N * N * sizeof(double));
  double *Sdata = (double*) malloc(N * N * sizeof(double));
  double Amax = 0.0;
  double Amin = 0.0;
  
  for (int i=0; i<N; ++i) {
    for (int j=0; j<N; ++j) {

#define n nle->order
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




  /* -------------------------------------------------------------------
   * Tabulate the U(A) integral
   */
  double *Atable = (double*) malloc(num_bins * sizeof(double));
  double *Btable = (double*) malloc(num_bins * sizeof(double));
  double *Utable = (double*) malloc(num_bins * sizeof(double));
  double U0 = Uintegral(Adata, Jdata, Sdata, 0.0, N); 
  for (int n=0; n<num_bins; ++n) {
    double A = Amin + n * (Amax - Amin) / (num_bins-1);
    double U = Uintegral(Adata, Jdata, Sdata, A, N) - U0;
    
    Atable[n] = A;
    Utable[n] = U;
  }


  /* -------------------------------------------------------------------
   * Determine the minimum value of the U(A) integral
   */
  double Umin = Utable[0];
  for (int n=0; n<num_bins; ++n) {
    if (Utable[n] < Umin) {
      Umin = Utable[n];
    }
  }


  
  for (int n=0; n<num_bins; ++n) {
    double U = Utable[n];
    double f = 1.0; /* 1.0 gives the smallest mean field value */
    Btable[n] = sqrt(2*(U - f*Umin));
  }


  /* -------------------------------------------------------------------
   * Store the result in the nle struct and clean up
   */
  nle->Atable = Atable;
  nle->Btable = Btable;
  nle->Utable = Utable;
  cow_histogram_del(dSdA_hist);
  free(Adata);
  free(Jdata);
  free(Sdata);
}


void ffe_nle_free(struct ffe_nle *nle)
{
  free(nle->Atable);
  free(nle->Btable);
  free(nle->Utable);
}


void ffe_nle_sample(struct ffe_nle *nle, double x, double y, double B[4])
{

  x *= 2 * M_PI;
  y *= 2 * M_PI;

  int m = nle->order;
  double C = cos(x) - sin(y);
  double A = pow(C, m);
  
  for (int n=0; n<nle->num_bins; ++n) {
    if (nle->Atable[n] < A && A < nle->Atable[n+1]) {

      double A0 = nle->Atable[n];
      double dA = (nle->Atable[n+1] - nle->Atable[n]);
      double dB = (nle->Btable[n+1] - nle->Btable[n]);
      
      B[3] = nle->Btable[n] + (dB / dA) * (A - A0);
      break;
    }
  }
  
  B[1] = -m * cos(y) * pow(C, m-1);
  B[2] = +m * sin(x) * pow(C, m-1);
}


void ffe_nle_write_table(struct ffe_nle *nle, const char *fname)
{
  FILE *outf = fopen(fname, "w");
  if (outf == NULL) return;
  for (int n=0; n<nle->num_bins; ++n) {
    double A = nle->Atable[n];
    double B = nle->Btable[n];
    double U = nle->Utable[n];
    fprintf(outf, "%+12.10e %+12.10e %+12.10e\n", A, B, U);
  }
  fclose(outf);
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
