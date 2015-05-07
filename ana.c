#include <stdio.h>
#include <math.h>
#include "ffe.h"



/*
 * Carry out measurement diagnostic
 * =====================================================================
 */
void ffe_sim_measure(struct ffe_sim *sim, struct ffe_measure *meas)
{
#define GLB_AVG(x) x = cow_domain_dblsum(sim->domain, x) / Nt

  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int Nt = cow_domain_getnumglobalzones(sim->domain, COW_ALL_DIMS);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];

  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);

  meas->electric_energy = 0.0;
  meas->magnetic_energy = 0.0;
  meas->magnetic_monopole = 0.0;
  meas->magnetic_helicity = 0.0; /* TODO */

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);

    double divB = D1(B,1) + D2(B,2) + D3(B,3);    

    double EE = DOT(&E[m], &E[m]);
    double BB = DOT(&B[m], &B[m]);

    meas->electric_energy += 0.5 * EE;
    meas->magnetic_energy += 0.5 * BB;
    meas->magnetic_monopole += fabs(divB);

  }

  GLB_AVG(meas->electric_energy);
  GLB_AVG(meas->magnetic_energy);
  GLB_AVG(meas->magnetic_monopole);

#undef GLB_AVG
}



void ffe_sim_analyze(struct ffe_sim *sim, char *filename)
{
  cow_domain *domain = sim->domain;
  cow_dfield *electric = sim->electric[0];
  cow_dfield *magnetic = sim->magnetic[0];
  cow_dfield *vecpoten = cow_dfield_new();
  cow_dfield *jcurrent = cow_dfield_new();


  char gname[1024];
  char nname[1024];

  snprintf(gname, 1024, "spectra-%06d", sim->status.iteration);
  snprintf(nname, 1024, "%12.10e", sim->status.time_simulation);


  /* Data fields setup */
  /* ---------------------------------------------------- */
  cow_dfield_setdomain(vecpoten, domain);
  cow_dfield_setname(vecpoten, "vector_potential");
  cow_dfield_addmember(vecpoten, "A1");
  cow_dfield_addmember(vecpoten, "A2");
  cow_dfield_addmember(vecpoten, "A3");
  cow_dfield_commit(vecpoten);

  cow_dfield_setdomain(jcurrent, domain);
  cow_dfield_setname(jcurrent, "electric_current");
  cow_dfield_addmember(jcurrent, "J1");
  cow_dfield_addmember(jcurrent, "J2");
  cow_dfield_addmember(jcurrent, "J3");
  cow_dfield_commit(jcurrent);



  /* Histograms setup */
  /* ---------------------------------------------------- */
  cow_histogram *Pe = cow_histogram_new();
  cow_histogram_setlower(Pe, 0, 1);
  cow_histogram_setupper(Pe, 0, 8192);
  cow_histogram_setnbins(Pe, 0, 4096);
  cow_histogram_setspacing(Pe, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Pe, "electric");
  cow_histogram_setfullname(Pe, nname);


  cow_histogram *Pb = cow_histogram_new();
  cow_histogram_setlower(Pb, 0, 1);
  cow_histogram_setupper(Pb, 0, 8192);
  cow_histogram_setnbins(Pb, 0, 4096);
  cow_histogram_setspacing(Pb, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Pb, "magnetic");
  cow_histogram_setfullname(Pb, nname);


  cow_histogram *Hr = cow_histogram_new();
  cow_histogram_setlower(Hr, 0, 1);
  cow_histogram_setupper(Hr, 0, 8192);
  cow_histogram_setnbins(Hr, 0, 4096);
  cow_histogram_setspacing(Hr, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Hr, "helicity-real");
  cow_histogram_setfullname(Hr, nname);


  cow_histogram *Hi = cow_histogram_new();
  cow_histogram_setlower(Hi, 0, 1);
  cow_histogram_setupper(Hi, 0, 8192);
  cow_histogram_setnbins(Hi, 0, 4096);
  cow_histogram_setspacing(Hi, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Hi, "helicity-imag");
  cow_histogram_setfullname(Hi, nname);


  cow_fft_inversecurl(magnetic, vecpoten);
  cow_fft_curl       (magnetic, jcurrent);


  cow_fft_pspecvecfield(magnetic, Pb);
  cow_fft_pspecvecfield(electric, Pe);
  cow_fft_helicityspec(magnetic, Hr, Hi);

  int Nt = cow_domain_getnumglobalzones(sim->domain, COW_ALL_DIMS);
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  double *A = (double*) cow_dfield_getdatabuffer(vecpoten);
  double *B = (double*) cow_dfield_getdatabuffer(magnetic);

  double htot = 0.0;

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {
    int m = INDV(i,j,k);
    htot += DOT(&A[m], &B[m]);
  }

  htot = cow_domain_dblsum(domain, htot) / Nt;
  printf("[main] total magnetic helicity: %8.6e\n", htot);



  if (filename) {

    cow_histogram_dumphdf5(Pe, filename, gname);
    cow_histogram_dumphdf5(Pb, filename, gname);
    cow_histogram_dumphdf5(Hr, filename, gname);
    cow_histogram_dumphdf5(Hi, filename, gname);

    if (0) { /* write derived fields */
      cow_dfield_write(magnetic, filename);
      cow_dfield_write(vecpoten, filename);
      cow_dfield_write(jcurrent, filename);
    }

  }

  cow_histogram_del(Pb);
  cow_histogram_del(Pe);
  cow_histogram_del(Hr);
  cow_histogram_del(Hi);

  cow_dfield_del(vecpoten);
  cow_dfield_del(jcurrent);
}
