#include <math.h>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h> /* mkdir */

#include "ffe.h"
#include "jsw_rand.h"



void nav_sim_write_checkpoint(struct nav_sim *sim, const char *base_name)
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

  cow_dfield_write(sim->velocity[0], chkpt_name);
  cow_dfield_write(sim->pressure, chkpt_name);

  if (cow_domain_getcartrank(sim->domain) == 0) {
    //read_write_status(sim, chkpt_name, 'w');
    //read_write_sim(sim, chkpt_name, 'w');
  }
}



/*
 * Evaluate initial data
 * =====================================================================
 */
void nav_sim_initial_data(struct nav_sim *sim)
{
  jsw_rand_t R;
  jsw_seed(&R, cow_domain_getcartrank(sim->domain));


  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->velocity[0], 0);
  int sj = cow_dfield_getstride(sim->velocity[0], 1);
  int sk = cow_dfield_getstride(sim->velocity[0], 2);
  double *u = cow_dfield_getdatabuffer(sim->velocity[0]);

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);

    double x[4] = {0,
    		   cow_domain_positionatindex(sim->domain, 0, i),
    		   cow_domain_positionatindex(sim->domain, 1, j),
    		   cow_domain_positionatindex(sim->domain, 2, k)};

    double w;
    sim->initial_data(sim, x, &u[m], &w);

    
    /* We add a white-noise perturbation to the velocity if pert < 0.0 */
    if (sim->perturbation < 0.0) {
      double du1 = jsw_random_double(&R, -1.0, 1.0);
      double du2 = jsw_random_double(&R, -1.0, 1.0);
      double du3 = jsw_random_double(&R, -1.0, 1.0);
      u[m+1] += du1 * sim->perturbation;
      u[m+2] += du2 * sim->perturbation;
      u[m+3] += du3 * sim->perturbation;
    }
    else {
      //u[m+1] += sim->perturbation * sin(2 * M_PI * x[2]);
      //u[m+2] += sim->perturbation * cos(2 * M_PI * x[1]);
      //u[m+3] += sim->perturbation * 0;
    }
  }

  cow_dfield_syncguard(sim->velocity[0]);
}



void nav_sim_init(struct nav_sim *sim)
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

  sim->inertial = cow_dfield_new();
  cow_dfield_setname(sim->inertial, "inertial");
  cow_dfield_setdomain(sim->inertial, sim->domain);
  cow_dfield_addmember(sim->inertial, "z1");
  cow_dfield_addmember(sim->inertial, "z2");
  cow_dfield_addmember(sim->inertial, "z3");
  cow_dfield_commit(sim->inertial);
  
  sim->divinert = cow_dfield_new();
  cow_dfield_setname(sim->divinert, "divinert");
  cow_dfield_setdomain(sim->divinert, sim->domain);
  cow_dfield_addmember(sim->divinert, "divz");
  cow_dfield_commit(sim->divinert);
  
  sim->pressure = cow_dfield_new();
  cow_dfield_setname(sim->pressure, "pressure");
  cow_dfield_setdomain(sim->pressure, sim->domain);
  cow_dfield_addmember(sim->pressure, "p");
  cow_dfield_commit(sim->pressure);
 

  for (int n=0; n<6; ++n) {
    sim->velocity[n] = cow_dfield_new();

    cow_dfield_setname(sim->velocity[n], "velocity");
    cow_dfield_setdomain(sim->velocity[n], sim->domain);

    cow_dfield_addmember(sim->velocity[n], "u1");
    cow_dfield_addmember(sim->velocity[n], "u2");
    cow_dfield_addmember(sim->velocity[n], "u3");

    cow_dfield_commit(sim->velocity[n]);
  }

}



void nav_sim_free(struct nav_sim *sim)
{
  for (int n=0; n<6; ++n) {
    cow_dfield_del(sim->velocity[n]);
  }
  cow_dfield_del(sim->inertial);
  cow_dfield_del(sim->pressure);
  cow_dfield_del(sim->divinert);
  cow_domain_del(sim->domain);
}



int nav_sim_problem_setup(struct nav_sim *sim, const char *problem_name)
{
  if (problem_name == NULL) {
    printf("1. wave\n");
    printf("2. abc\n");
    return 0;
  }
  else if (!strcmp(problem_name, "wave")) {
    sim->initial_data = nav_initial_data_wave;
    return 0;
  }
  else if (!strcmp(problem_name, "abc")) {
    sim->initial_data = nav_initial_data_abc;
    return 0;
  }
  else {
    return 1;
  }
}



/*
 * Apply Kreiss-Oliger operator to subtract high frequencies
 * =====================================================================
 */
void nav_sim_kreiss_oliger(struct nav_sim *sim)
{
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->velocity[0], 0);
  int sj = cow_dfield_getstride(sim->velocity[0], 1);
  int sk = cow_dfield_getstride(sim->velocity[0], 2);
  double *u = cow_dfield_getdatabuffer(sim->velocity[0]);
  double *du = cow_dfield_getdatabuffer(sim->velocity[1]);

  double KO_const = 0.0;

  switch (FFE_DISSIPATION_ORDER) {
  case 2: KO_const = -1./4 ; break;
  case 4: KO_const = -1./16; break;
  case 6: KO_const = -1./64; break;
  }

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);

    for (int d=1; d<=3; ++d) {
      du[m+d] = KOV(u,d);
    }
  }

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);

    double eps = sim->eps_parameter;

    for (int d=1; d<=3; ++d) {
      u[m+d] -= eps * KO_const * du[m+d];
    }
  }

  cow_dfield_syncguard(sim->velocity[0]);
}



/*
 * Advance the simulation by one Runge-Kutta substep
 * =====================================================================
 */
void nav_sim_advance_rk(struct nav_sim *sim, int RKstep)
{
  double RKparam_array[5] = {0.0, 0.5, 0.5, 1.0};
  double RKparam = RKparam_array[RKstep];

  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->velocity[0], 0);
  int sj = cow_dfield_getstride(sim->velocity[0], 1);
  int sk = cow_dfield_getstride(sim->velocity[0], 2);
  /* int ti = cow_dfield_getstride(sim->pressure, 0); */
  /* int tj = cow_dfield_getstride(sim->pressure, 1); */
  /* int tk = cow_dfield_getstride(sim->pressure, 2); */

  double *u0 = cow_dfield_getdatabuffer(sim->velocity[0]);
  double *u  = cow_dfield_getdatabuffer(sim->velocity[1]);
  /* double *p  = cow_dfield_getdatabuffer(sim->pressure); */
  /* double *z  = cow_dfield_getdatabuffer(sim->inertial); */
  /* double *g  = cow_dfield_getdatabuffer(sim->divinert); */
  
  double *dtu = cow_dfield_getdatabuffer(sim->velocity[RKstep+2]);

  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];
  double dt = sim->status.time_step;

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
    u[m+1] = u0[m+1] + dt * RKparam * dtu[m+1];
    u[m+2] = u0[m+2] + dt * RKparam * dtu[m+2];
    u[m+3] = u0[m+3] + dt * RKparam * dtu[m+3];
  }
  cow_dfield_syncguard(sim->velocity[1]);
  cow_fft_helmholtzdecomp(sim->velocity[1], COW_PROJECT_OUT_DIV);
  
  
  /* ===========================================================================
   * Fill in the n-th (n=0,1,2,3) time-derivative register, reading from the
   * [1] field register
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {
    int m = INDV(i,j,k);

    double Du1[4] = { 0, D1(u,1), D2(u,1), D3(u,1) };
    double Du2[4] = { 0, D1(u,2), D2(u,2), D3(u,2) };
    double Du3[4] = { 0, D1(u,3), D2(u,3), D3(u,3) };

    dtu[m+1] = -(u[m+1]*Du1[1] + u[m+2]*Du1[2] + u[m+3]*Du1[3]);
    dtu[m+2] = -(u[m+1]*Du2[1] + u[m+2]*Du2[2] + u[m+3]*Du2[3]);
    dtu[m+3] = -(u[m+1]*Du3[1] + u[m+2]*Du3[2] + u[m+3]*Du3[3]);
  }
  //cow_dfield_syncguard(sim->velocity[RKstep+2]);
  //cow_fft_helmholtzdecomp(sim->velocity[RKstep+2], COW_PROJECT_OUT_DIV);
}




/*
 * Average Runge-Kutta substeps to complete a full time step
 * =====================================================================
 */
void nav_sim_average_rk(struct nav_sim *sim)
{
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);

  int si = cow_dfield_getstride(sim->velocity[0], 0);
  int sj = cow_dfield_getstride(sim->velocity[0], 1);
  int sk = cow_dfield_getstride(sim->velocity[0], 2);

  double *u = cow_dfield_getdatabuffer(sim->velocity[0]);

  double *dtu0 = cow_dfield_getdatabuffer(sim->velocity[2]);
  double *dtu1 = cow_dfield_getdatabuffer(sim->velocity[3]);
  double *dtu2 = cow_dfield_getdatabuffer(sim->velocity[4]);
  double *dtu3 = cow_dfield_getdatabuffer(sim->velocity[5]);

  double dt = sim->status.time_step;


  /* ===========================================================================
   * Average the RK substeps, write result into register 0
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);

    u[m+1] += dt/6 * (dtu0[m+1] + 2*dtu1[m+1] + 2*dtu2[m+1] + dtu3[m+1]);
    u[m+2] += dt/6 * (dtu0[m+2] + 2*dtu1[m+2] + 2*dtu2[m+2] + dtu3[m+2]);
    u[m+3] += dt/6 * (dtu0[m+3] + 2*dtu1[m+3] + 2*dtu2[m+3] + dtu3[m+3]);
  }


  cow_dfield_syncguard(sim->velocity[0]);
}



/*
 * Advance the simulation by one full iteration
 * =====================================================================
 */
void nav_sim_advance(struct nav_sim *sim)
{
  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];
  double dt_light = MIN3(dx, dy, dz);

  /* !!! will need to revisit this to put in sound speed !!! */
  sim->status.time_step = sim->cfl_parameter * dt_light;

  nav_sim_advance_rk(sim, 0);
  nav_sim_advance_rk(sim, 1);
  nav_sim_advance_rk(sim, 2);
  nav_sim_advance_rk(sim, 3);
  nav_sim_average_rk(sim);
  nav_sim_kreiss_oliger(sim);
  
  sim->status.iteration += 1;
  sim->status.time_simulation += sim->status.time_step;
}


  
int main(int argc, char **argv)
{ 
  cow_init(0, NULL, 0);
  int norun_main = 0;
  int restarted_run = 0;
  
  struct nav_sim sim;
  //struct ffe_status status;
  memset(&sim, 0, sizeof(struct nav_sim));
  strcpy(sim.output_directory, ".");
  strcpy(sim.problem_name, "");


  
  sim.Ni = 128;
  sim.Nj = 128;
  sim.Nk = 1;
  sim.time_final = 1.0;
  sim.time_between_checkpoints = 1.0;
  sim.cfl_parameter = 0.10;
  sim.eps_parameter = 0.50; /* [0-1] */
  sim.perturbation = 0.0;
  sim.alpha_squared = 1.0;
  sim.abc_coefficients[0] = 1.0;
  sim.abc_coefficients[1] = 1.0;
  sim.abc_coefficients[2] = 0.0;



  /*
   * Print a help message
   * ===================================================================
   */
  printf("\nEuler flow solver\n");
  printf("Jonathan Zrake, Stanford University (2016)\n");

  if (argc == 1) {
    printf("usage: nav <problem-name> [tmax=1.0] [N=16,16,16]\n");
    printf("problems are:\n");
    nav_sim_problem_setup(NULL, NULL);
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
   
    //norun_main += read_write_sim(&sim, argv[1], 'r');
    //norun_main += read_write_status(&sim, argv[1], 'r');

    if (norun_main == 0) {
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
    else if (!strncmp(argv[n], "outdir=", 7)) {
      sscanf(argv[n], "outdir=%1024s", sim.output_directory);
    }
    else if (!strncmp(argv[n], "N=", 2)) {
      int num = sscanf(argv[n], "N=%d,%d,%d", &sim.Ni, &sim.Nj, &sim.Nk);
      if (num != 3) {
	printf("[ffe] error: badly formatted option '%s'\n", argv[n]);
	norun_main += 1;
      }
    }
    else if (!strncmp(argv[n], "cfl=", 4)) {
      sscanf(argv[n], "cfl=%lf", &sim.cfl_parameter);
    }
    else if (!strncmp(argv[n], "eps=", 4)) {
      sscanf(argv[n], "eps=%lf", &sim.eps_parameter);
    }
    else if (!strncmp(argv[n], "k2=", 3)) {
      sscanf(argv[n], "k2=%d", &sim.alpha_squared);
    }
    else if (!strncmp(argv[n], "pert=", 5)) {
      sscanf(argv[n], "pert=%lf", &sim.perturbation);
    }
    else {
      printf("[nav] error: unrecognized option '%s'\n", argv[n]);
      norun_main += 1;
    }

  }




  printf("\n-----------------------------------------\n");
  printf("resolution ................. %d %d %d\n", sim.Ni, sim.Nj, sim.Nk);
  printf("cfl_parameter .............. %12.10lf\n", sim.cfl_parameter);
  printf("eps_parameter .............. %12.10lf\n", sim.eps_parameter);
  printf("time_between_checkpoints ... %12.10lf\n", sim.time_between_checkpoints);
  printf("time_final ................. %12.10lf\n", sim.time_final);
  printf("perturbation ............... %12.10lf\n", sim.perturbation);
  printf("alpha_squared .............. %d\n", sim.alpha_squared);
  printf("-----------------------------------------\n\n");



  if (norun_main) {
    cow_finalize();
    return 0;
  }
  else if (nav_sim_problem_setup(&sim, sim.problem_name)) {
    printf("[nav] error: unkown problem name: '%s', choose one of\n",
	   sim.problem_name);
    nav_sim_problem_setup(NULL, NULL);
    cow_finalize();
    return 0;
  }


  nav_sim_init(&sim);

  /* cow_domain_setcollective(sim.domain, sim.io_use_collective); */
  /* cow_domain_setchunk(sim.domain, sim.io_use_chunked); */
  /* cow_domain_setalign(sim.domain, */
  /* 		      sim.io_align_threshold * 1024, */
  /* 		      sim.io_disk_block_size * 1024); */


  
  if (restarted_run) {
    cow_dfield_read(sim.velocity[0], argv[1]);
  }
  else {
    sim.status.time_last_checkpoint = -sim.time_between_checkpoints;
    sim.status.checkpoint_number = -1;
    nav_sim_initial_data(&sim);
  }


  int local_grid_size = cow_domain_getnumlocalzonesinterior(sim.domain,
							    COW_ALL_DIMS);


  if (cow_domain_getcartrank(sim.domain) == 0) {

    /*
     * Set up the problem directory and output log
     */


    mkdir(sim.output_directory, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
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

      nav_sim_write_checkpoint(&sim, NULL);
    }


    /*
     * Evolve the system
     * =================================================================
     */    
    void *start_cycle = cow_start_clock();

    nav_sim_advance(&sim);

    double seconds = cow_stop_clock(start_cycle);

    sim.status.kzps = 1e-3 * local_grid_size / seconds;

    if (sim.status.iteration % 1 == 0) {
      printf("[nav] n=%06d t=%6.4e dt=%6.4e %3.2f kzps\n",
	     sim.status.iteration,
	     sim.status.time_simulation,
	     sim.status.time_step,
	     sim.status.kzps);
      fflush(stdout);
    }
  }


  

  if (sim.time_between_checkpoints > 0.0) {
    if (norun_main == 0) {
      if (sim.time_final == 0.0) {
	sim.status.checkpoint_number = 0;
	nav_sim_write_checkpoint(&sim, NULL);
      }
      else {
	nav_sim_write_checkpoint(&sim, "final");
      }
    }
  }

  
  nav_sim_free(&sim);
  cow_finalize();
  return 0;
}
