#include <math.h>
#include <complex.h>
#include "ffe.h"
#include "jsw_rand.h"


#ifndef M_PI
#define M_PI 3.1415926535897926
#endif




static void random_beltrami_field(double x[4], double B[4], int model, int k2, double h);





void initial_data_emwave(struct ffe_sim *sim, double x[4], double E[4], double B[4])
{
  E[1] = sin(2 * M_PI * x[3]);
  E[2] = 0.0;
  E[3] = 0.0;

  B[1] = 0.0;
  B[2] = sin(2 * M_PI * x[3]);
  B[3] = 0.0;
}

void initial_data_alfvenwave(struct ffe_sim *sim, double x[4], double E[4], double B[4])
{
  E[1] = 0.0;
  E[2] = 0.0;
  E[3] = 0.0;

  B[1] = 0.0;
  B[2] = 0.1 * sin(2 * M_PI * x[3]);
  B[3] = 1.0;
}

void initial_data_abc(struct ffe_sim *sim, double x[4], double E[4], double B[4])
{
  double a = sim->abc_coefficients[0];
  double b = sim->abc_coefficients[1];
  double c = sim->abc_coefficients[2];
  double h = sim->fractional_helicity;
  double alpha = sqrt(sim->alpha_squared) * 2 * M_PI;

  E[1] = 0.0;
  E[2] = 0.0;
  E[3] = 0.0;

  B[1] = 0.0;
  B[2] = 0.0;
  B[3] = 0.0;

  B[2] += a * cos(alpha * x[1]);
  B[3] -= a * sin(alpha * x[1]) * h;
  B[1] += 0.0;

  B[3] += b * cos(alpha * x[2]);
  B[1] -= b * sin(alpha * x[2]) * h;
  B[2] += 0.0;

  B[1] += c * cos(alpha * x[3]);
  B[2] -= c * sin(alpha * x[3]) * h;
  B[3] += 0.0;
}

void initial_data_beltrami(struct ffe_sim *sim, double x[4], double E[4], double B[4])
{
  E[1] = 0.0;
  E[2] = 0.0;
  E[3] = 0.0;

  random_beltrami_field(x, B, 0, sim->alpha_squared, sim->fractional_helicity);
}








typedef complex double Complex;
typedef struct fourier_mode
{
  double k[4];
  Complex A[4];
} fourier_mode;




void random_beltrami_field(double x[4], double B[4], int model, int k2, double h)
{
#define RAND jsw_random_double(&rand, -1, 1)
  int m,i,j,k;
  jsw_rand_t rand;
  jsw_seed(&rand, model);


  int k2_sphere = k2;
  int k_cube = floor(sqrt(k2_sphere)) + 1;
  fourier_mode *modes = NULL;
  int num_modes = 0;
  Complex A[4] = {0, 0, 0, 0};
  double amp;

  for (i=-k_cube; i<=k_cube; ++i) {
    for (j=-k_cube; j<=k_cube; ++j) {
      for (k=-k_cube; k<=k_cube; ++k) {
  	fourier_mode M;

  	if (i*i + j*j + k*k != k2_sphere) {
  	  continue;
  	}
  	else {
  	  M.k[0] = 0.0;
  	  M.k[1] = i;
  	  M.k[2] = j;
  	  M.k[3] = k;
  	  M.A[0] = 0.0;
  	  M.A[1] = RAND + RAND*I*h;
  	  M.A[2] = RAND + RAND*I*h;
  	  M.A[3] = RAND + RAND*I*h;
  	  amp = sqrt(M.A[1]*conj(M.A[1]) + M.A[2]*conj(M.A[2]) + M.A[3]*conj(M.A[3]));
  	  M.A[1] /= amp;
  	  M.A[2] /= amp;
  	  M.A[3] /= amp;
  	  num_modes += 1;
  	  modes = (fourier_mode *) realloc(modes, num_modes * sizeof(fourier_mode));
  	  modes[num_modes-1] = M;

  	  /* printf("k[%d] = [%d %d %d] is on shell\n", num_modes, i, j, k); */
  	}
      }
    }
  }


  for (m=0; m<num_modes; ++m) {
    fourier_mode M = modes[m];
    double a = sqrt(M.k[1]*M.k[1] + M.k[2]*M.k[2] + M.k[3]*M.k[3]);
    Complex K[4] = {0, I*M.k[1], I*M.k[2], I*M.k[3]};
    Complex Ikx  = (K[1]*x[1] + K[2]*x[2] + K[3]*x[3]) * 2 * M_PI;
    Complex P[4] = {0, M.A[1], M.A[2], M.A[3]}; /* a times psi */

    Complex T[4] = {0, /* T = K cross (a psi) */
		    K[2]*P[3] - K[3]*P[2],
		    K[3]*P[1] - K[1]*P[3],
		    K[1]*P[2] - K[2]*P[1]};

    Complex S[4] = {0, /* S = K cross T / alpha */
		    (K[2]*T[3] - K[3]*T[2])/a,
		    (K[3]*T[1] - K[1]*T[3])/a,
		    (K[1]*T[2] - K[2]*T[1])/a};

    A[1] += (h*S[1] + T[1])/a * cexp(Ikx);
    A[2] += (h*S[2] + T[2])/a * cexp(Ikx);
    A[3] += (h*S[3] + T[3])/a * cexp(Ikx);
  }


  free(modes);

  B[1] = A[1];
  B[2] = A[2];
  B[3] = A[3];

#undef RAND
}
