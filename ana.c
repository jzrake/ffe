#include <stdio.h>
#include <math.h>
#include "ffe.h"

#ifndef M_PI
#define M_PI 3.1415926535897926
#endif



int ffe_measure_fscanf(struct ffe_measure *meas, FILE *F)
{
  return fscanf(F, "%lf %lf %lf %lf %lfe\n",
		&meas->time_simulation,
		&meas->electric_energy,
		&meas->magnetic_energy,
		&meas->magnetic_helicity,
		&meas->magnetic_monopole);
}



int ffe_measure_fprintf(struct ffe_measure *meas, FILE *F)
{
  return fprintf(F, "%+12.10e %+12.10e %+12.10e %+12.10e %+12.10e\n",
		 meas->time_simulation,
		 meas->electric_energy,
		 meas->magnetic_energy,
		 meas->magnetic_helicity,
		 meas->magnetic_monopole);
}



/*
 * Carry out measurement diagnostic
 * =====================================================================
 */
void ffe_sim_measure(struct ffe_sim *sim, struct ffe_measure *meas)
{
#define GLB_AVG(x) x = cow_domain_dblsum(sim->domain, x) / Nt

  long long Nt = cow_domain_getnumglobalzones(sim->domain, COW_ALL_DIMS);
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];

  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);

  meas->time_simulation = sim->status.time_simulation;
  meas->electric_energy = 0.0;
  meas->magnetic_energy = 0.0;
  meas->magnetic_monopole = 0.0;
  /* meas->magnetic_helicity = 0.0; we assume it was already calculated */

  FOR_ALL_INTERIOR_NO_THREAD(Ni, Nj, Nk) {

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



void ffe_sim_analyze(struct ffe_sim *sim, struct ffe_measure *meas, char *filename)
{
  cow_domain *domain = sim->domain;
  cow_dfield *electric = sim->electric[0];
  cow_dfield *magnetic = sim->magnetic[0];
  cow_dfield *vecpoten = cow_dfield_new();
  cow_dfield *jcurrent = cow_dfield_new();


  char gname[1024];
  char nname[1024];
  double sim_alpha_value = sqrt(sim->alpha_squared);


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
  cow_histogram_setupper(Pe, 0, sim->max_pspec_bin);
  cow_histogram_setnbins(Pe, 0, sim->num_pspec_bins);
  cow_histogram_setspacing(Pe, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Pe, "electric");
  cow_histogram_setfullname(Pe, nname);


  cow_histogram *Pb = cow_histogram_new();
  cow_histogram_setlower(Pb, 0, 1);
  cow_histogram_setupper(Pb, 0, sim->max_pspec_bin);
  cow_histogram_setnbins(Pb, 0, sim->num_pspec_bins);
  cow_histogram_setspacing(Pb, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Pb, "magnetic");
  cow_histogram_setfullname(Pb, nname);


  cow_histogram *Hr = cow_histogram_new();
  cow_histogram_setlower(Hr, 0, 1);
  cow_histogram_setupper(Hr, 0, sim->max_pspec_bin);
  cow_histogram_setnbins(Hr, 0, sim->num_pspec_bins);
  cow_histogram_setspacing(Hr, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Hr, "helicity");
  cow_histogram_setfullname(Hr, nname);


  cow_histogram *alpha_hist = cow_histogram_new();
  cow_histogram_setlower(alpha_hist, 0, -8 * sim_alpha_value); /* educated guess of extent */
  cow_histogram_setupper(alpha_hist, 0, +8 * sim_alpha_value);
  cow_histogram_setnbins(alpha_hist, 0, sim->num_pspec_bins);
  cow_histogram_setspacing(alpha_hist, COW_HIST_SPACING_LINEAR);
  cow_histogram_setbinmode(alpha_hist, COW_HIST_BINMODE_COUNTS);
  cow_histogram_setfullname(alpha_hist, "alpha-hist");
  cow_histogram_setnickname(alpha_hist, "alpha-hist");
  cow_histogram_setdomaincomm(alpha_hist, domain);
  cow_histogram_commit(alpha_hist);


  cow_histogram *mup_hist = cow_histogram_new(); /* 1 - cos(theta) */
  cow_histogram_setlower(mup_hist, 0, 1e-8);
  cow_histogram_setupper(mup_hist, 0, 1.0);
  cow_histogram_setnbins(mup_hist, 0, 512);
  cow_histogram_setspacing(mup_hist, COW_HIST_SPACING_LOG);
  cow_histogram_setbinmode(mup_hist, COW_HIST_BINMODE_DENSITY);
  cow_histogram_setfullname(mup_hist, "mup-hist");
  cow_histogram_setnickname(mup_hist, "mup-hist");
  cow_histogram_setdomaincomm(mup_hist, domain);
  cow_histogram_commit(mup_hist);


  cow_histogram *mum_hist = cow_histogram_new(); /* 1 + cos(theta) */
  cow_histogram_setlower(mum_hist, 0, 1e-8);
  cow_histogram_setupper(mum_hist, 0, 1.0);
  cow_histogram_setnbins(mum_hist, 0, 512);
  cow_histogram_setspacing(mum_hist, COW_HIST_SPACING_LOG);
  cow_histogram_setbinmode(mum_hist, COW_HIST_BINMODE_DENSITY);
  cow_histogram_setfullname(mum_hist, "mum-hist");
  cow_histogram_setnickname(mum_hist, "mum-hist");
  cow_histogram_setdomaincomm(mum_hist, domain);
  cow_histogram_commit(mum_hist);



  cow_fft_inversecurl(magnetic, vecpoten);
  cow_fft_curl       (magnetic, jcurrent);


  cow_fft_pspecvecfield(magnetic, Pb);
  cow_fft_pspecvecfield(electric, Pe);
  cow_fft_helicityspec(magnetic, Hr, 'b');

  long long Nt = cow_domain_getnumglobalzones(sim->domain, COW_ALL_DIMS);
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  double *E = (double*) cow_dfield_getdatabuffer(electric);
  double *A = (double*) cow_dfield_getdatabuffer(vecpoten);
  double *B = (double*) cow_dfield_getdatabuffer(magnetic);
  double *J = (double*) cow_dfield_getdatabuffer(jcurrent);

  double htot = 0.0;
  double mtot = 0.0;
  double utot = 0.0;
  double etot = 0.0;
  double ftot = 0.0; /* Lorentz force squared */
  
  FOR_ALL_INTERIOR_NO_THREAD(Ni, Nj, Nk) {
    int m = INDV(i,j,k);

    double EE = DOT(&E[m], &E[m]);
    double BB = DOT(&B[m], &B[m]);
    double AB = DOT(&A[m], &B[m]);
    double JB = DOT(&J[m], &B[m]);
    double JJ = DOT(&J[m], &J[m]);
    double F[4] = CROSS(&J[m], &B[m]);
    double FF = DOT(F, F);
    
    cow_histogram_addsample1(alpha_hist, JB/BB/(2*M_PI), 1.0);
    cow_histogram_addsample1(mup_hist, 1 - JB / sqrt(BB * JJ), 1.0);
    cow_histogram_addsample1(mum_hist, 1 + JB / sqrt(BB * JJ), 1.0);

    htot += AB;
    mtot += JB;
    utot += BB * 0.5;
    etot += EE * 0.5;
    ftot += FF;
  }

  cow_histogram_seal(alpha_hist);
  cow_histogram_seal(mup_hist);
  cow_histogram_seal(mum_hist);

  htot = cow_domain_dblsum(domain, htot) / Nt;
  mtot = cow_domain_dblsum(domain, mtot) / Nt;
  utot = cow_domain_dblsum(domain, utot) / Nt;
  etot = cow_domain_dblsum(domain, etot) / Nt;
  ftot = cow_domain_dblsum(domain, ftot) / Nt;
  
  printf("[meas] Hm=%8.6e U=%8.6e a=%8.6e F=%4.3f\n", htot, utot+etot,
	 2*utot/htot/(2*M_PI),
	 sqrt(ftot)/(2*utot)); /* This quantity is |J cross B| / B^2, has
				  dimensions of inverse length */


  if (meas) {
    meas->magnetic_helicity = htot;
  }


  if (filename) {

    cow_histogram_dumphdf5(Pe, filename, gname);
    cow_histogram_dumphdf5(Pb, filename, gname);
    cow_histogram_dumphdf5(Hr, filename, gname);
    cow_histogram_dumphdf5(alpha_hist, filename, gname);
    cow_histogram_dumphdf5(mup_hist, filename, gname);
    cow_histogram_dumphdf5(mum_hist, filename, gname);

    if (0) { /* write derived fields */
      cow_dfield_write(vecpoten, filename);
      cow_dfield_write(jcurrent, filename);
    }

  }

  cow_histogram_del(Pb);
  cow_histogram_del(Pe);
  cow_histogram_del(Hr);
  cow_histogram_del(alpha_hist);
  cow_histogram_del(mup_hist);
  cow_histogram_del(mum_hist);

  cow_dfield_del(vecpoten);
  cow_dfield_del(jcurrent);
}
