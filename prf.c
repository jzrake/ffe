#include "ffe.h"

void ffe_perf(struct ffe_sim *sim)
{
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];

  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);

  void *C = cow_start_clock();

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);
    double divB = D1(B,1) + D2(B,2) + D3(B,3);
    double rotB[4] = {0, D2(B,3) - D3(B,2), D3(B,1) - D1(B,3), D1(B,2) - D2(B,1)};
    double lplB[4] = {0, KOV(B,1), KOV(B,2), KOV(B,3)};

    double foo = divB*divB + DOT(rotB, rotB) + DOT(lplB, lplB);

    if (foo < 0.0) {
      printf("QUOI?\n");
    }
  }

  int local_grid_size = cow_domain_getnumlocalzonesinterior(sim->domain,
							    COW_ALL_DIMS);
  double seconds = cow_stop_clock(C);
  double kzps = 1e-3 * local_grid_size / seconds / sim->omp_num_threads;
  printf("[perf] %6.4e seconds (%3.2f kzps)\n", seconds, kzps);
}
