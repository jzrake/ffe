#include <math.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h> /* mkdir */
#include "ffe.h"


/*
 * Evalute the current based on E, B, and gradients
 * =====================================================================
 */
static void ffe_ohms_law(enum FfeSimParameter flag_ohms_law,
			 double damping_timescale,
			 double E[4], double rotE[4], double divE,
			 double B[4], double rotB[4], double divB, double J[4])
{
  switch (flag_ohms_law) {

  case FFE_OHMS_LAW_VACUUM:
    J[1] = 0.0;
    J[2] = 0.0;
    J[3] = 0.0;
    break;
  case FFE_OHMS_LAW_FORCE_FREE: /* eqn 11: Pfeiffer (2013) */
    {
      double S[4] = CROSS(E, B);
      double B2 = DOT(B, B);
      double Jb = DOT(B, rotB) - DOT(E, rotE);
      double Js = divE;

      if (B2 < 1e-12) B2 = 1e-12; /* Don't divide by zero! */

      J[1] = (B[1] * Jb + S[1] * Js) / B2;
      J[2] = (B[2] * Jb + S[2] * Js) / B2;
      J[3] = (B[3] * Jb + S[3] * Js) / B2;

      if (damping_timescale > 0.0) {
	J[1] += E[1] / damping_timescale;
	J[2] += E[2] / damping_timescale;
	J[3] += E[3] / damping_timescale;
      }

    }
    break;
  }
}



static void enforce_EB_constraints(double E[4], double B[4])
{
    /*
     * Subtract out the component of E parallel to B
     */
    double BB = DOT(B, B);
    double EB = DOT(E, B);


    E[1] -= EB/BB * B[1];
    E[2] -= EB/BB * B[2];
    E[3] -= EB/BB * B[3];




    /*
     * Cap the electric field vector to ensure E <= B
     */
    double EE = DOT(E, E);

    if (EE > BB) {

      double f = sqrt(BB/EE);

      E[1] *= f;
      E[2] *= f;
      E[3] *= f;

    }
}


/*
 * Initialze a new simulation instance. Parameters that need to be initialized
 * before this call:
 *
 * -> Ni, Nj, Nk
 * -> initial_data
 * -> flag_ohms_law
 * -> time_between_checkpoints
 * -> time_final
 * -> output_directory
 * -> eps_parameter
 * -> cfl_parameter
 * -> alpha_squared
 * -> fractional_helicity
 *
 * =====================================================================
 */
void ffe_sim_init(struct ffe_sim *sim)
{
  sim->domain = cow_domain_new();


  sim->grid_spacing[0] = 0.0;
  sim->grid_spacing[1] = 1.0 / sim->Ni;
  sim->grid_spacing[2] = 1.0 / sim->Nj;
  sim->grid_spacing[3] = 1.0 / sim->Nk;



  cow_domain_setndim(sim->domain, (sim->Ni>1) + (sim->Nj>1) + (sim->Nk>1));
  cow_domain_setsize(sim->domain, 0, sim->Ni);
  cow_domain_setsize(sim->domain, 1, sim->Nj);
  cow_domain_setsize(sim->domain, 2, sim->Nk);
  cow_domain_setguard(sim->domain, FFE_NG);
  cow_domain_commit(sim->domain);

  for (int n=0; n<6; ++n) {
    sim->electric[n] = cow_dfield_new();
    sim->magnetic[n] = cow_dfield_new();
    sim->psifield[n] = cow_dfield_new();

    cow_dfield_setname(sim->electric[n], "electric");
    cow_dfield_setname(sim->magnetic[n], "magnetic");
    cow_dfield_setname(sim->psifield[n], "psifield");

    cow_dfield_setdomain(sim->electric[n], sim->domain);
    cow_dfield_setdomain(sim->magnetic[n], sim->domain);
    cow_dfield_setdomain(sim->psifield[n], sim->domain);

    cow_dfield_addmember(sim->electric[n], "E1");
    cow_dfield_addmember(sim->electric[n], "E2");
    cow_dfield_addmember(sim->electric[n], "E3");
    cow_dfield_addmember(sim->magnetic[n], "B1");
    cow_dfield_addmember(sim->magnetic[n], "B2");
    cow_dfield_addmember(sim->magnetic[n], "B3");

    cow_dfield_addmember(sim->psifield[n], "P");

    cow_dfield_commit(sim->electric[n]);
    cow_dfield_commit(sim->magnetic[n]);
    cow_dfield_commit(sim->psifield[n]);
  }
}



/*
 * Finalize a simulation instance
 * =====================================================================
 */
void ffe_sim_free(struct ffe_sim *sim)
{

  for (int n=0; n<6; ++n) {
    cow_dfield_del(sim->electric[n]);
    cow_dfield_del(sim->magnetic[n]);
    cow_dfield_del(sim->psifield[n]);
  }

  cow_domain_del(sim->domain);
}



int ffe_sim_problem_setup(struct ffe_sim *sim, const char *problem_name)
{
  if (problem_name == NULL) {
    printf("1. emwave\n");
    printf("2. alfvenwave\n");
    printf("3. abc\n");
    printf("4. beltrami\n");
    printf("5. clayer\n");
    return 0;
  }
  else if (!strcmp(problem_name, "emwave")) {
    sim->initial_data = initial_data_emwave;
    sim->flag_ohms_law = FFE_OHMS_LAW_VACUUM;
    return 0;
  }
  else if (!strcmp(problem_name, "alfvenwave")) {
    sim->initial_data = initial_data_alfvenwave;
    sim->flag_ohms_law = FFE_OHMS_LAW_FORCE_FREE;
    return 0;
  }
  else if (!strcmp(problem_name, "abc")) {
    sim->initial_data = initial_data_abc;
    sim->flag_ohms_law = FFE_OHMS_LAW_FORCE_FREE;
    return 0;
  }
  else if (!strcmp(problem_name, "beltrami")) {
    sim->initial_data = initial_data_beltrami;
    sim->flag_ohms_law = FFE_OHMS_LAW_FORCE_FREE;
    return 0;
  }
  else if (!strcmp(problem_name, "clayer")) {
    sim->initial_data = initial_data_clayer;
    sim->flag_ohms_law = FFE_OHMS_LAW_FORCE_FREE;
    sim->alpha_squared = 16384;
    return 0;
  }
  else {
    return 1;
  }
}



/*
 * Evaluate initial data
 * =====================================================================
 */
void ffe_sim_initial_data(struct ffe_sim *sim)
{
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  int ti = cow_dfield_getstride(sim->psifield[0], 0);
  int tj = cow_dfield_getstride(sim->psifield[0], 1);
  int tk = cow_dfield_getstride(sim->psifield[0], 2);
  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);
  double *P = cow_dfield_getdatabuffer(sim->psifield[0]);

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {    

    int m = INDV(i,j,k);
    int n = INDS(i,j,k);

    double x[4] = {0,
    		   cow_domain_positionatindex(sim->domain, 0, i),
    		   cow_domain_positionatindex(sim->domain, 1, j),
    		   cow_domain_positionatindex(sim->domain, 2, k)};

    sim->initial_data(sim, x, &E[m], &B[m]);
    P[n] = 0.0;
  }

  cow_dfield_syncguard(sim->electric[0]);
  cow_dfield_syncguard(sim->magnetic[0]);
  cow_dfield_syncguard(sim->psifield[0]);
}



/*
 * Advance the simulation by one Runge-Kutta substep
 * =====================================================================
 */
void ffe_sim_advance_rk(struct ffe_sim *sim, int RKstep)
{
  double RKparam_array[5] = {0.0, 0.5, 0.5, 1.0};
  double RKparam = RKparam_array[RKstep];

  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  int ti = cow_dfield_getstride(sim->psifield[0], 0);
  int tj = cow_dfield_getstride(sim->psifield[0], 1);
  int tk = cow_dfield_getstride(sim->psifield[0], 2);

  double *E0 = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B0 = cow_dfield_getdatabuffer(sim->magnetic[0]);
  double *P0 = cow_dfield_getdatabuffer(sim->psifield[0]);
  double *E  = cow_dfield_getdatabuffer(sim->electric[1]);
  double *B  = cow_dfield_getdatabuffer(sim->magnetic[1]);
  double *P  = cow_dfield_getdatabuffer(sim->psifield[1]);

  double *dtE = cow_dfield_getdatabuffer(sim->electric[RKstep+2]);
  double *dtB = cow_dfield_getdatabuffer(sim->magnetic[RKstep+2]);
  double *dtP = cow_dfield_getdatabuffer(sim->psifield[RKstep+2]);

  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];
  double dt = sim->status.time_step;

  double ch2 = 1.00; /* Dedner wave speed (squared) */
  double tau = 0.02; /* Dedner damping time */
  double KO_const = 0.0;

  switch (FFE_DISSIPATION_ORDER) {
  case 2: KO_const = -1./4 ; break;
  case 4: KO_const = -1./16; break;
  case 6: KO_const = -1./64; break;
  }


  /* ===========================================================================
   * Fill in the n-th (n=0,1,2,3) field register
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);
    int n = INDS(i,j,k);

    E[m+1] = E0[m+1] + dt * RKparam * dtE[m+1];
    E[m+2] = E0[m+2] + dt * RKparam * dtE[m+2];
    E[m+3] = E0[m+3] + dt * RKparam * dtE[m+3];

    B[m+1] = B0[m+1] + dt * RKparam * dtB[m+1];
    B[m+2] = B0[m+2] + dt * RKparam * dtB[m+2];
    B[m+3] = B0[m+3] + dt * RKparam * dtB[m+3];

    P[n] = P0[n] + dt * RKparam * dtP[n];

    enforce_EB_constraints(&E[m], &B[m]);
  }


  cow_dfield_syncguard(sim->electric[1]);
  cow_dfield_syncguard(sim->magnetic[1]);
  cow_dfield_syncguard(sim->psifield[1]);




  /* ===========================================================================
   * Fill in the n-th (n=0,1,2,3) time-derivative register, reading from the
   * [1] field register
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);
    int n = INDS(i,j,k);

    double divE = D1(E,1) + D2(E,2) + D3(E,3);
    double divB = D1(B,1) + D2(B,2) + D3(B,3);    
    double rotE[4] = {0, D2(E,3) - D3(E,2), D3(E,1) - D1(E,3), D1(E,2) - D2(E,1)};
    double rotB[4] = {0, D2(B,3) - D3(B,2), D3(B,1) - D1(B,3), D1(B,2) - D2(B,1)};




    /* Electric and magnetic current evaluation */
    double J[4] = { 0, 0, 0, 0 };
    double F[4] = { 0, 0, 0, 0 };

    ffe_ohms_law(sim->flag_ohms_law, sim->damping_timescale,
		 &E[m], rotE, divE, &B[m], rotB, divB, J);

    F[1] = -S1(P); /* (minus) grad-psi = magnetic current */
    F[2] = -S2(P);
    F[3] = -S3(P);


    /* Maxwell's equations */
    for (int d=1; d<=3; ++d) {
      dtE[m+d] = +rotB[d] - J[d];
      dtB[m+d] = -rotE[d] + F[d];
    }



    /* Dedner */
    dtP[n] = -ch2 * divB - P[n] / tau;




    /* Kreiss-Oliger dissipation */
    if (sim->kreiss_oliger_mode == 'f') {
      double lplE[4] = {0, KOV(E,1), KOV(E,2), KOV(E,3)};
      double lplB[4] = {0, KOV(B,1), KOV(B,2), KOV(B,3)};
      double eps = sim->eps_parameter;

      for (int d=1; d<=3; ++d) {
	dtE[m+d] -= eps * KO_const * lplE[d];
	dtB[m+d] -= eps * KO_const * lplB[d];
      }

      dtP[n] -= eps * KO_const * KOS(P);
    }



    /* Hyperbolicity terms, eqn 48-49: Pfeiffer (2013) */
    if (sim->pfeiffer_terms == 't') {

      double B2 = DOT(&B[m], &B[m]);
      double gradEdotB[4] = {0, 0, 0, 0};

      double S[4]   = CROSS(&E[m], &B[m]);
      double ct1[4] = CROSS(&E[m], gradEdotB);
      double ct3[4] = CROSS(&B[m], gradEdotB);
      double ct2[4] = {0, S[1] * divB, S[2] * divB, S[3] * divB};

      double gamma1 = 0.0;
      double gamma2 = 1.0;
      double gamma3 = 1.0;

      if (B2 < 1e-12) B2 = 1e-12;

      for (int d=1; d<=3; ++d) {
	dtE[m+d] -= gamma1 / B2 * ct1[d];
	dtB[m+d] -= gamma2 / B2 * ct2[d] + gamma3 / B2 * ct3[d];
      }
    }
  }
}



/*
 * Apply Kreiss-Oliger operator to subtract high frequencies
 * =====================================================================
 */
void ffe_sim_kreiss_oliger(struct ffe_sim *sim)
{
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  int ti = cow_dfield_getstride(sim->psifield[0], 0);
  int tj = cow_dfield_getstride(sim->psifield[0], 1);
  int tk = cow_dfield_getstride(sim->psifield[0], 2);
  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);
  double *P = cow_dfield_getdatabuffer(sim->psifield[0]);
  double *dE = cow_dfield_getdatabuffer(sim->electric[1]);
  double *dB = cow_dfield_getdatabuffer(sim->magnetic[1]);
  double *dP = cow_dfield_getdatabuffer(sim->psifield[1]);

  double KO_const = 0.0;

  switch (FFE_DISSIPATION_ORDER) {
  case 2: KO_const = -1./4 ; break;
  case 4: KO_const = -1./16; break;
  case 6: KO_const = -1./64; break;
  }

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);
    int n = INDS(i,j,k);

    for (int d=1; d<=3; ++d) {
      dE[m+d] = KOV(E,d);
      dB[m+d] = KOV(B,d);
    }

    dP[n] = KOS(P);
  }

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);
    int n = INDS(i,j,k);

    double eps = sim->eps_parameter;

    for (int d=1; d<=3; ++d) {
      E[m+d] -= eps * KO_const * dE[m+d];
      B[m+d] -= eps * KO_const * dB[m+d];
    }

    P[n] -= eps * KO_const * dP[n];

    enforce_EB_constraints(&E[m], &B[m]);
  }


  cow_dfield_syncguard(sim->electric[0]);
  cow_dfield_syncguard(sim->magnetic[0]);
  cow_dfield_syncguard(sim->psifield[0]);
}



/*
 * Average Runge-Kutta substeps to complete a full time step
 * =====================================================================
 */
void ffe_sim_average_rk(struct ffe_sim *sim)
{
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);

  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  int ti = cow_dfield_getstride(sim->psifield[0], 0);
  int tj = cow_dfield_getstride(sim->psifield[0], 1);
  int tk = cow_dfield_getstride(sim->psifield[0], 2);

  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);
  double *P = cow_dfield_getdatabuffer(sim->psifield[0]);

  double *dtE0 = cow_dfield_getdatabuffer(sim->electric[2]);
  double *dtB0 = cow_dfield_getdatabuffer(sim->magnetic[2]);
  double *dtP0 = cow_dfield_getdatabuffer(sim->psifield[2]);
  double *dtE1 = cow_dfield_getdatabuffer(sim->electric[3]);
  double *dtB1 = cow_dfield_getdatabuffer(sim->magnetic[3]);
  double *dtP1 = cow_dfield_getdatabuffer(sim->psifield[3]);
  double *dtE2 = cow_dfield_getdatabuffer(sim->electric[4]);
  double *dtB2 = cow_dfield_getdatabuffer(sim->magnetic[4]);
  double *dtP2 = cow_dfield_getdatabuffer(sim->psifield[4]);
  double *dtE3 = cow_dfield_getdatabuffer(sim->electric[5]);
  double *dtB3 = cow_dfield_getdatabuffer(sim->magnetic[5]);
  double *dtP3 = cow_dfield_getdatabuffer(sim->psifield[5]);

  double dt = sim->status.time_step;


  /* ===========================================================================
   * Average the RK substeps, write result into register 0
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);
    int n = INDS(i,j,k);

    E[m+1] += dt/6 * (dtE0[m+1] + 2*dtE1[m+1] + 2*dtE2[m+1] + dtE3[m+1]);
    E[m+2] += dt/6 * (dtE0[m+2] + 2*dtE1[m+2] + 2*dtE2[m+2] + dtE3[m+2]);
    E[m+3] += dt/6 * (dtE0[m+3] + 2*dtE1[m+3] + 2*dtE2[m+3] + dtE3[m+3]);

    B[m+1] += dt/6 * (dtB0[m+1] + 2*dtB1[m+1] + 2*dtB2[m+1] + dtB3[m+1]);
    B[m+2] += dt/6 * (dtB0[m+2] + 2*dtB1[m+2] + 2*dtB2[m+2] + dtB3[m+2]);
    B[m+3] += dt/6 * (dtB0[m+3] + 2*dtB1[m+3] + 2*dtB2[m+3] + dtB3[m+3]);

    P[n] += dt/6 * (dtP0[n] + 2*dtP1[n] + 2*dtP2[n] + dtP3[n]);

    enforce_EB_constraints(&E[m], &B[m]);
  }


  cow_dfield_syncguard(sim->electric[0]);
  cow_dfield_syncguard(sim->magnetic[0]);
  cow_dfield_syncguard(sim->psifield[0]);
}



/*
 * Advance the simulation by one full iteration
 * =====================================================================
 */
void ffe_sim_advance(struct ffe_sim *sim)
{
  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];
  double dt = sim->cfl_parameter * MIN3(dx,dy,dz);

  sim->status.time_step = dt;

  ffe_sim_advance_rk(sim, 0);
  ffe_sim_advance_rk(sim, 1);
  ffe_sim_advance_rk(sim, 2);
  ffe_sim_advance_rk(sim, 3);
  ffe_sim_average_rk(sim);

  if (sim->kreiss_oliger_mode == 'c') {
    ffe_sim_kreiss_oliger(sim);
  }

  sim->status.iteration += 1;
  sim->status.time_simulation += dt;
}





void ffe_sim_write_checkpoint(struct ffe_sim *sim, const char *base_name)
{
  char chkpt_name[1024];
  if (base_name == NULL) {
    snprintf(chkpt_name, 1024, "%s/chkpt.%04d.h5",
	     sim->output_directory,
	     sim->status.checkpoint_number);
  }
  else {
    snprintf(chkpt_name, 1024, "%s/chkpt.%s.h5",
	     sim->output_directory,
	     base_name);
  }

  cow_domain_setcollective(sim->domain, sim->io_use_collective);
  cow_domain_setchunk(sim->domain, sim->io_use_chunked);
  cow_domain_setalign(sim->domain,
		      sim->io_align_threshold * 1024,
		      sim->io_disk_block_size * 1024);

  cow_dfield_write(sim->electric[0], chkpt_name);
  cow_dfield_write(sim->magnetic[0], chkpt_name);
  cow_dfield_write(sim->psifield[0], chkpt_name);

  if (strchr(sim->write_derived_fields, 'J') != NULL) {
    cow_dfield *jcurrent = cow_dfield_new();
    cow_dfield_setdomain(jcurrent, sim->domain);
    cow_dfield_setname(jcurrent, "electric_current");
    cow_dfield_addmember(jcurrent, "J1");
    cow_dfield_addmember(jcurrent, "J2");
    cow_dfield_addmember(jcurrent, "J3");
    cow_dfield_commit(jcurrent);
    cow_fft_curl(sim->magnetic[0], jcurrent);
    cow_dfield_write(jcurrent, chkpt_name);
    cow_dfield_del(jcurrent);
  }

  if (strchr(sim->write_derived_fields, 'A') != NULL) {
    cow_dfield *vecpoten = cow_dfield_new();
    cow_dfield_setdomain(vecpoten, sim->domain);
    cow_dfield_setname(vecpoten, "vector_potential");
    cow_dfield_addmember(vecpoten, "A1");
    cow_dfield_addmember(vecpoten, "A2");
    cow_dfield_addmember(vecpoten, "A3");
    cow_dfield_commit(vecpoten);
    cow_fft_inversecurl(sim->magnetic[0], vecpoten);
    cow_dfield_write(vecpoten, chkpt_name);
    cow_dfield_del(vecpoten);
  }

  if (cow_domain_getcartrank(sim->domain) == 0) {
    read_write_status(sim, chkpt_name, 'w');
    read_write_sim(sim, chkpt_name, 'w');
  }
}



int main(int argc, char **argv)
{
  cow_init(0, NULL, 0);


  int invalid_cfg = 0;
  int restarted_run = 0;
  char logfile_name[1024];
  char anlfile_name[1024];
  struct ffe_sim sim;
  struct ffe_measure measure;

  memset(sim.output_directory, '\0', 1024);
  memset(sim.problem_name, '\0', 1024);
  memset(sim.write_derived_fields, '\0', 1024);

  strcpy(sim.output_directory, ".");
  strcpy(sim.write_derived_fields, "J");

  sim.Ni = 128;
  sim.Nj = 128;
  sim.Nk = 1;
  sim.time_final = 1.0;
  sim.time_between_checkpoints = 1.0;
  sim.cfl_parameter = 0.10;
  sim.eps_parameter = 0.50; /* [0-1] */
  sim.kreiss_oliger_mode = 'c'; /* TODO: correct means of doing 'fine' */
  sim.pfeiffer_terms = 'f';
  sim.damping_timescale = -1.0;
  sim.alpha_squared = 1.0;
  sim.fractional_helicity = 1.0; /* [0-1] */
  sim.abc_coefficients[0] = 1.0;
  sim.abc_coefficients[1] = 1.0;
  sim.abc_coefficients[2] = 0.0;
  sim.measure_cadence = 1;
  sim.analyze_cadence = 100;
  sim.num_pspec_bins = 4096;
  sim.max_pspec_bin = 8192;
  sim.io_use_collective = 1;
  sim.io_use_chunked = 1;
  sim.io_align_threshold = 1; /* means no threshold */
  sim.io_disk_block_size = 1; /* means no alignment */




  /*
   * Print a help message
   * ===================================================================
   */
  printf("\nForce-free electrodynamics solver\n");
  printf("Jonathan Zrake, Stanford University (2015)\n");

  if (argc == 1) {
    printf("usage: ffe <problem-name> [tmax=1.0] [N=16,16,16]\n");
    printf("problems are:\n");
    ffe_sim_problem_setup(NULL, NULL);
    cow_finalize();
    return 0;
  }
  else {
    strncpy(sim.problem_name, argv[1], 1024);
  }



  /*
   * Set up the problem defaults
   * ===================================================================
   */
  if (strstr(argv[1], ".h5") != 0) {
   
    invalid_cfg += read_write_sim(&sim, argv[1], 'r');
    invalid_cfg += read_write_status(&sim, argv[1], 'r');

    if (invalid_cfg == 0) {
      restarted_run = 1;
    }
  }
  else {

    /* set up a fresh status struct */
    sim.status.iteration = 0;
    sim.status.checkpoint_number = 0;
    sim.status.time_simulation = 0.0;
    sim.status.time_step = 0.0;
    sim.status.time_last_checkpoint = 0;
    sim.status.kzps = 0.0;

  }



  if (invalid_cfg) {
    cow_finalize();
    return 0;
  }
  else if (ffe_sim_problem_setup(&sim, sim.problem_name)) {
    printf("[ffe] error: unkown problem name: '%s', choose one of\n",
	   sim.problem_name);
    ffe_sim_problem_setup(NULL, NULL);
    cow_finalize();
    return 0;
  }



  /*
   * Scan command line arguments
   * ===================================================================
   */
  for (int n=2; n<argc; ++n) {

    if (!strncmp(argv[n], "tmax=", 5)) {
      sscanf(argv[n], "tmax=%lf", &sim.time_final);
    }
    else if (!strncmp(argv[n], "cpi=", 4)) {
      sscanf(argv[n], "cpi=%lf", &sim.time_between_checkpoints);
    }
    else if (!strncmp(argv[n], "wdf=", 4)) {
      sscanf(argv[n], "wdf=%1024s", sim.write_derived_fields);
    }
    else if (!strncmp(argv[n], "N=", 2)) {
      int num = sscanf(argv[n], "N=%d,%d,%d", &sim.Ni, &sim.Nj, &sim.Nk);
      if (num != 3) {
	printf("[ffe] error: badly formatted option '%s'\n", argv[n]);
	invalid_cfg += 1;
      }
    }
    else if (!strncmp(argv[n], "outdir=", 7)) {
      sscanf(argv[n], "outdir=%1024s", sim.output_directory);
    }
    else if (!strncmp(argv[n], "cfl=", 4)) {
      sscanf(argv[n], "cfl=%lf", &sim.cfl_parameter);
    }
    else if (!strncmp(argv[n], "eps=", 4)) {
      sscanf(argv[n], "eps=%lf", &sim.eps_parameter);
    }
    else if (!strncmp(argv[n], "ko=", 3)) {
      sscanf(argv[n], "ko=%c", &sim.kreiss_oliger_mode);
    }
    else if (!strncmp(argv[n], "pfeiffer=", 9)) {
      sscanf(argv[n], "pfeiffer=%c", &sim.pfeiffer_terms);
    }
    else if (!strncmp(argv[n], "k2=", 3)) {
      sscanf(argv[n], "k2=%d", &sim.alpha_squared);
    }
    else if (!strncmp(argv[n], "helicity=", 9)) {
      sscanf(argv[n], "helicity=%lf", &sim.fractional_helicity);
    }
    else if (!strncmp(argv[n], "abc=", 2)) {
      int num = sscanf(argv[n], "abc=%lf,%lf,%lf",
		       &sim.abc_coefficients[0],
		       &sim.abc_coefficients[1],
		       &sim.abc_coefficients[2]);
      if (num != 3) {
	printf("[ffe] error: badly formatted option '%s'\n", argv[n]);
	invalid_cfg += 1;
      }
    }
    else if (!strncmp(argv[n], "tdamp=", 6)) {
      sscanf(argv[n], "tdamp=%lf", &sim.damping_timescale);
    }
    else if (!strncmp(argv[n], "post=", 5)) {
      int num = sscanf(argv[n], "post=%d,%d,%d,%d",
		       &sim.measure_cadence,
		       &sim.analyze_cadence,
		       &sim.num_pspec_bins,
		       &sim.max_pspec_bin);
      if (num != 4) {
	printf("[ffe] error: badly formatted option '%s'\n", argv[n]);
	invalid_cfg += 1;
      }
    }
    else if (!strncmp(argv[n], "io=", 3)) {
      int num = sscanf(argv[n], "io=%d,%d,%d,%d",
		       &sim.io_use_collective,
		       &sim.io_use_chunked,
		       &sim.io_align_threshold,
		       &sim.io_disk_block_size);
      if (num != 4) {
	printf("[ffe] error: badly formatted option '%s'\n", argv[n]);
	invalid_cfg += 1;
      }
    }
    else {
      printf("[ffe] error: unrecognized option '%s'\n", argv[n]);
      invalid_cfg += 1;
    }
  }





  ffe_sim_init(&sim);

  if (restarted_run) {

    cow_dfield_read(sim.electric[0], argv[1]);
    cow_dfield_read(sim.magnetic[0], argv[1]);
    cow_dfield_read(sim.psifield[0], argv[1]);

  }

  else {

    sim.status.time_last_checkpoint = -sim.time_between_checkpoints;
    sim.status.checkpoint_number = -1;
    ffe_sim_initial_data(&sim);

  }

  printf("\n-----------------------------------------\n");
  printf("resolution ................. %d %d %d\n", sim.Ni, sim.Nj, sim.Nk);
  printf("cfl_parameter .............. %12.10lf\n", sim.cfl_parameter);
  printf("eps_parameter .............. %12.10lf\n", sim.eps_parameter);
  printf("time_between_checkpoints ... %12.10lf\n", sim.time_between_checkpoints);
  printf("time_final ................. %12.10lf\n", sim.time_final);
  printf("fractional_helicity ........ %12.10lf\n", sim.fractional_helicity);
  printf("damping_timescale .......... %12.10lf\n", sim.damping_timescale);
  printf("alpha_squared .............. %d\n", sim.alpha_squared);
  printf("kriess_oliger_mode ......... %c\n", sim.kreiss_oliger_mode);
  printf("pfeiffer_terms ............. %c\n", sim.pfeiffer_terms);
  printf("output_directory ........... %s\n", sim.output_directory);
  printf("problem_name ............... %s\n", sim.problem_name);
  printf("measure_cadence ............ %d\n", sim.measure_cadence);
  printf("analyze_cadence ............ %d\n", sim.analyze_cadence);
  printf("num_pspec_bins ............. %d\n", sim.num_pspec_bins);
  printf("max_pspec_bin .............. %d\n", sim.max_pspec_bin);
  printf("io_use_collective .......... %d\n", sim.io_use_collective);
  printf("io_use_chunked ............. %d\n", sim.io_use_chunked);
  printf("io_align_threshold ......... %d\n", sim.io_align_threshold);
  printf("io_disk_block_size ......... %d\n", sim.io_disk_block_size);
  printf("-----------------------------------------\n\n");


  int local_grid_size = cow_domain_getnumlocalzonesinterior(sim.domain,
							    COW_ALL_DIMS);


  snprintf(logfile_name, 1024, "%s/ffe.dat"    , sim.output_directory);
  snprintf(anlfile_name, 1024, "%s/analysis.h5", sim.output_directory);

  if (cow_domain_getcartrank(sim.domain) == 0) {

    FILE *logf = NULL;

    mkdir(sim.output_directory, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (restarted_run) {
      logf = fopen(logfile_name, "a");
    }
    else {
      logf = fopen(logfile_name, "w");
    }

    
    if (logf == NULL) {
      printf("[ffe] error: could not open log file '%s'\n", logfile_name);
      invalid_cfg += 1;
    }
    else {
      fclose(logf);
    }
  }


  /* Propagate any errors to all procs */
  invalid_cfg += cow_domain_intsum(sim.domain, invalid_cfg);  

  if (invalid_cfg) {
    sim.time_final = 0.0;
  }


  while (sim.status.time_simulation < sim.time_final) {


    /*
     * Write a checkpoint if it's time
     * =================================================================
     */
    if (sim.status.time_simulation - sim.status.time_last_checkpoint >=
	sim.time_between_checkpoints && sim.time_between_checkpoints > 0.0) {

      sim.status.time_last_checkpoint += sim.time_between_checkpoints;
      sim.status.checkpoint_number += 1;

      ffe_sim_write_checkpoint(&sim, NULL);
    }



    /*
     * Handle post-processing and reductions
     * =================================================================
     */
    int iter = sim.status.iteration;
    if (iter % sim.measure_cadence == 0) ffe_sim_measure(&sim, &measure);
    if (iter % sim.analyze_cadence == 0) ffe_sim_analyze(&sim, anlfile_name);



    if (cow_domain_getcartrank(sim.domain) == 0) {

      FILE *logf = fopen(logfile_name, "a");

      fprintf(logf, "%+12.10e %+12.10e %+12.10e %+12.10e %+12.10e\n",
    	      sim.status.time_simulation,
    	      measure.electric_energy,
    	      measure.magnetic_energy,
    	      measure.magnetic_helicity,
    	      measure.magnetic_monopole);

      fclose(logf);
    }



    /*
     * Evolve the system
     * =================================================================
     */
    clock_t start_cycle, stop_cycle;
    start_cycle = clock();

    ffe_sim_advance(&sim);

    stop_cycle = clock();
    sim.status.kzps = 1e-3 * local_grid_size / (stop_cycle - start_cycle) *
      CLOCKS_PER_SEC;


    if (sim.status.iteration % 1 == 0) {
      printf("[ffe] n=%06d t=%6.4e %3.2f kzps\n",
	     sim.status.iteration,
	     sim.status.time_simulation,
	     sim.status.kzps);
      fflush(stdout);
    }
  }


  if (sim.time_between_checkpoints > 0.0) {
    if (invalid_cfg == 0) {
      if (sim.time_final == 0.0) {
	sim.status.checkpoint_number = 0;
	ffe_sim_write_checkpoint(&sim, NULL);
      }
      else {
	ffe_sim_write_checkpoint(&sim, "final");
      }
    }
  }


  ffe_sim_free(&sim);
  cow_finalize();
  return 0;
}




