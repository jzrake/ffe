#include <math.h>
#include <complex.h>
#include "ffe.h"
#include "jsw_rand.h"


#ifndef M_PI
#define M_PI 3.1415926535897926
#endif




static void random_beltrami_field(double x[4], double B[4], int model, int k2, double h, int rank);
static void random_beltrami_field_kk65641(double x[4], double B[4], int model, int k2, double h, int rank);




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

  int rank = 0;
  rank += sim->Ni > 1;
  rank += sim->Nj > 1;
  rank += sim->Nk > 1;

  if (sim->alpha_squared == 65641) {
    random_beltrami_field_kk65641(x, B, 0, sim->alpha_squared, sim->fractional_helicity, rank);
  }
  else {
    random_beltrami_field(x, B, 0, sim->alpha_squared, sim->fractional_helicity, rank);
  }
}

void initial_data_clayer(struct ffe_sim *sim, double x[4], double E[4], double B[4])
{
  double w = 0.01;
  double dx = 0.0;//5 * w * sin(2 * M_PI * x[2]);
  double Bz = tanh((x[1] - 0.25 + dx)/w) - tanh((x[1] - 0.75 + dx)/w) - 1.0;

  E[1] = 0.0;
  E[2] = 0.0;
  E[3] = 0.0;

  B[1] = 0.0;
  B[2] = sqrt(1.0 - Bz*Bz) * (x[1] > 0.5 ? 1 : -1);
  B[3] = Bz;
}






typedef complex double Complex;
typedef struct fourier_mode
{
  double k[4];
  Complex A[4];
} fourier_mode;




void random_beltrami_field(double x[4], double B[4], int model, int k2, double h, int rank)
{
#define RAND jsw_random_double(&rand, -1, 1)
  int m,d,i,j,k;
  jsw_rand_t rand;
  jsw_seed(&rand, model);


  int k2_sphere = k2;
  int k_cube = floor(sqrt(k2_sphere)) + 1;
  fourier_mode *modes = NULL;
  int num_modes = 0;
  Complex A[4] = {0, 0, 0, 0};


  for (i=1; i<=k_cube; ++i) {
    for (j=-k_cube*(rank >= 2); j<=k_cube*(rank >= 2); ++j) {
      for (k=-k_cube*(rank >= 3); k<=k_cube*(rank >= 3); ++k) {

  	fourier_mode M;
	double phase = RAND * M_PI;


  	if (i*i + j*j + k*k != k2_sphere) {
  	  continue;
  	}
  	else {
  	  //printf("k[%d] = [%d %d %d] is on shell\n", num_modes, i, j, k);

  	  M.k[0] = 0.0;
  	  M.k[1] = i;
  	  M.k[2] = j;
  	  M.k[3] = k;


	  double e0[4] = {0, RAND, RAND, RAND}; /* any vector not parallel to k */
	  double e1[4] = CROSS(M.k, e0);
	  double A1 = sqrt(DOT(e1, e1));

	  for (d=1; d<=3; ++d) {
	    e1[d] /= A1;
	  }

	  M.A[0] = 0.0;
	  for (d=1; d<=3; ++d) {	  
	    M.A[d] = e1[d] * cexp(I * phase);
	  }

  	  num_modes += 1;
  	  modes = (fourier_mode *) realloc(modes, num_modes * sizeof(fourier_mode));
  	  modes[num_modes-1] = M;


	  /* reality condition */
	  for (d=1; d<=3; ++d) {
	    M.k[d] = -M.k[d];
	    M.A[d] = conj(M.A[d]);
	  }
  	  num_modes += 1;
  	  modes = (fourier_mode *) realloc(modes, num_modes * sizeof(fourier_mode));
  	  modes[num_modes-1] = M;
	}
      }
    }
  }

  for (m=0; m<num_modes; ++m) {
    double a = sqrt(k2);
    fourier_mode M = modes[m];
    Complex K[4] = {0, I*M.k[1], I*M.k[2], I*M.k[3]};
    Complex Ikx  = (K[1]*x[1] + K[2]*x[2] + K[3]*x[3]) * 2 * M_PI;
    Complex P[4] = {0, M.A[1], M.A[2], M.A[3]}; /* a times psi */
    Complex H[4] = CROSS(K, P);

    A[1] += (P[1] + h*H[1]/a) * cexp(Ikx);
    A[2] += (P[2] + h*H[2]/a) * cexp(Ikx);
    A[3] += (P[3] + h*H[3]/a) * cexp(Ikx);
  }


  free(modes);

  /* printf("im(A) = %+8.6e %+8.6e %+8.6e\n", */
  /* 	 cimag(A[1]), cimag(A[2]), cimag(A[3])); */

  B[1] = A[1];
  B[2] = A[2];
  B[3] = A[3];

#undef RAND
}



/*
 * Hard-coded high-k^2 IC, the one above is too slow for high k
*/
void random_beltrami_field_kk65641(double x[4], double B[4], int model, int k2, double h, int rank)
{
#define RAND jsw_random_double(&rand, -1, 1)
  int m,d,index=0;
  jsw_rand_t rand;
  jsw_seed(&rand, model);

  int k2_sphere = k2;
  //int k_cube = floor(sqrt(k2_sphere)) + 1;
  fourier_mode *modes = NULL;
  int num_modes = 0;
  Complex A[4] = {0, 0, 0, 0};


  int k_index[] = {
    155, -204, 0,
    155, +204, 0,
    165, -196, 0,
    165, +196, 0,
    196, -165, 0,
    196, +165, 0,
    204, -155, 0,
    204, +155, 0 };

  for (index=0; index<8; ++index) {

    fourier_mode M;
    double phase = RAND * M_PI;

    int i = k_index[index*3 + 0];
    int j = k_index[index*3 + 1];
    int k = k_index[index*3 + 2];

    if (i*i + j*j + k*k != k2_sphere) {
      continue;
    }
    else {
      /* printf("k[%d] = [%d %d %d] is on shell\n", num_modes, i, j, k); */

      M.k[0] = 0.0;
      M.k[1] = i;
      M.k[2] = j;
      M.k[3] = k;


      double e0[4] = {0, RAND, RAND, RAND}; /* any vector not parallel to k */
      double e1[4] = CROSS(M.k, e0);
      double A1 = sqrt(DOT(e1, e1));

      for (d=1; d<=3; ++d) {
	e1[d] /= A1;
      }

      M.A[0] = 0.0;
      for (d=1; d<=3; ++d) {	  
	M.A[d] = e1[d] * cexp(I * phase);
      }

      num_modes += 1;
      modes = (fourier_mode *) realloc(modes, num_modes * sizeof(fourier_mode));
      modes[num_modes-1] = M;


      /* reality condition */
      for (d=1; d<=3; ++d) {
	M.k[d] = -M.k[d];
	M.A[d] = conj(M.A[d]);
      }
      num_modes += 1;
      modes = (fourier_mode *) realloc(modes, num_modes * sizeof(fourier_mode));
      modes[num_modes-1] = M;
    }
  }

  for (m=0; m<num_modes; ++m) {
    double a = sqrt(k2);
    fourier_mode M = modes[m];
    Complex K[4] = {0, I*M.k[1], I*M.k[2], I*M.k[3]};
    Complex Ikx  = (K[1]*x[1] + K[2]*x[2] + K[3]*x[3]) * 2 * M_PI;
    Complex P[4] = {0, M.A[1], M.A[2], M.A[3]}; /* a times psi */
    Complex H[4] = CROSS(K, P);

    A[1] += (P[1] + h*H[1]/a) * cexp(Ikx);
    A[2] += (P[2] + h*H[2]/a) * cexp(Ikx);
    A[3] += (P[3] + h*H[3]/a) * cexp(Ikx);
  }


  free(modes);

  /* printf("im(A) = %+8.6e %+8.6e %+8.6e\n", */
  /* 	 cimag(A[1]), cimag(A[2]), cimag(A[3])); */

  B[1] = A[1];
  B[2] = A[2];
  B[3] = A[3];

#undef RAND
}
