#ifndef FFE_HEADER
#define FFE_HEADER

#include <stdio.h> /* FILE */
#include "cow/cow.h"



struct ffe_sim;
struct ffe_status;
struct ffe_measure;



typedef void (*InitialDataFunction)(struct ffe_sim *sim, double x[4], double E[4], double B[4]);

enum FfeSimParameter {
  FFE_OHMS_LAW_VACUUM,
  FFE_OHMS_LAW_FORCE_FREE
} ;


struct ffe_measure
{
  double time_simulation;
  double electric_energy;
  double magnetic_energy;
  double magnetic_helicity;
  double magnetic_monopole; /* |div(B)| */
} ;


struct ffe_status
{
  int iteration;
  int checkpoint_number;
  double time_simulation;
  double time_step;
  double time_last_checkpoint;
  double kzps;
} ;


struct ffe_particle
{
  double e;
  double m;
  double x[4];
  double u[4];
  double E[4];
  double B[4];
} ;


struct ffe_sim
{
  /* used by physics algorithm */
  cow_domain *domain;
  cow_dfield *electric[6]; /* 0,1: E, 4-6: dtE */
  cow_dfield *magnetic[6];
  cow_dfield *psifield[6]; /* 0,1: P, 4-6: dtP (Dedner psi field) */
  struct ffe_status status;
  struct ffe_particle *particles;
  cow_domain *particles_domain;
  cow_dfield *particles_dfield;


  /* set by problem type */
  double grid_spacing[4];
  enum FfeSimParameter flag_ohms_law;
  InitialDataFunction initial_data;


  /* set by problem type and/or user */
  int Ni, Nj, Nk;
  double cfl_parameter; /* Courant number [0.0 - 0.25]*/
  double eps_parameter; /* Kreiss-Oliger parameter [0 - 1] */
  double time_between_checkpoints;
  char write_derived_fields[1024];
  double time_final;
  double fractional_helicity;
  double abc_coefficients[3];
  double damping_timescale;
  double perturbation;
  int alpha_squared; /* wave-number (squared) of initial configuration */
  char kreiss_oliger_mode; /* 'c': coarse, 'f': fine 'n': none */
  char pfeiffer_terms; /* 't': true 'f': false */
  char output_directory[1024];
  char problem_name[1024];


  /* analysis options (e.g. post=1,100,4096,8192) */
  int measure_cadence;
  int analyze_cadence;
  int num_pspec_bins;
  int max_pspec_bin;


  /* IO options (e.g. io=1,1,1,1) */
  int io_use_collective;
  int io_use_chunked;
  int io_align_threshold; /* KB */
  int io_disk_block_size; /* KB */


  /* tracer particles */
  int num_particles;


  /* parallel options */
  int omp_num_threads;
} ;



void ffe_sim_init(struct ffe_sim *sim);
void ffe_sim_free(struct ffe_sim *sim);
void ffe_sim_analyze(struct ffe_sim *sim, struct ffe_measure *meas, char *filename);
void ffe_sim_measure(struct ffe_sim *sim, struct ffe_measure *meas);
void ffe_sim_average_rk(struct ffe_sim *sim);
void ffe_sim_advance_rk(struct ffe_sim *sim, int RKstep);
void ffe_sim_advance(struct ffe_sim *sim);
void ffe_sim_kreiss_oliger(struct ffe_sim *sim);
void ffe_sim_initial_data(struct ffe_sim *sim);
void ffe_sim_write_checkpoint(struct ffe_sim *sim, const char *base_name);
int  ffe_sim_problem_setup(struct ffe_sim *sim, const char *problem_name);

void ffe_par_move(struct ffe_sim *sim);
void ffe_par_sample(struct ffe_sim *sim);

int ffe_measure_fscanf(struct ffe_measure *meas, FILE *F);
int ffe_measure_fprintf(struct ffe_measure *meas, FILE *F);


/*
 * Initial data library
 * =====================================================================
 */
void initial_data_emwave    (struct ffe_sim *sim, double x[4], double E[4], double B[4]);
void initial_data_alfvenwave(struct ffe_sim *sim, double x[4], double E[4], double B[4]);
void initial_data_abc       (struct ffe_sim *sim, double x[4], double E[4], double B[4]);
void initial_data_beltrami  (struct ffe_sim *sim, double x[4], double E[4], double B[4]);
void initial_data_clayer    (struct ffe_sim *sim, double x[4], double E[4], double B[4]);


/*
 * Serialization library
 * =====================================================================
 */
int read_write_status(struct ffe_sim *sim, const char *chkpt_name, char mode);
int read_write_sim(struct ffe_sim *sim, const char *chkpt_name, char mode);


#define FFE_NG 3 /* number of guard zones */
#define FFE_DIFFERENCE_ORDER 4
#define FFE_DISSIPATION_ORDER 6




/*
 * Macro for a three-dimensional loop over all interior cells
 * =====================================================================
 */
#define FOR_ALL_INTERIOR(N1, N2, N3)				\
  _Pragma("omp parallel for collapse(3)")			\
  for (int i=N1==1?0:FFE_NG; i<N1+(N1==1?0:FFE_NG); ++i)	\
    for (int j=N2==1?0:FFE_NG; j<N2+(N2==1?0:FFE_NG); ++j)	\
      for (int k=N3==1?0:FFE_NG; k<N3+(N3==1?0:FFE_NG); ++k)	\


#define FOR_ALL_INTERIOR_NO_THREAD(N1, N2, N3)			\
  for (int i=N1==1?0:FFE_NG; i<N1+(N1==1?0:FFE_NG); ++i)	\
    for (int j=N2==1?0:FFE_NG; j<N2+(N2==1?0:FFE_NG); ++j)	\
      for (int k=N3==1?0:FFE_NG; k<N3+(N3==1?0:FFE_NG); ++k)	\



/*
 * Macro to calculate linear index of (i,j,k,m) ... m goes from 1, not 0
 * =====================================================================
 */

#define INDV(i,j,k) ((i)*si + (j)*sj + (k)*sk - 1)
#define INDS(i,j,k) ((i)*ti + (j)*tj + (k)*tk) /* for scalar field */



/*
 * Macros to calculate finite differences
 * =====================================================================
 */
#define DIFF1C2(F,s) ((-1*(F)[-1*s] +		\
		       +1*(F)[+1*s]) / 2.0)

#define DIFF1C4(F,s) ((+1*(F)[-2*s] +		\
		       -8*(F)[-1*s] +		\
		       +8*(F)[+1*s] +		\
		       -1*(F)[+2*s]) / 12.0)

#define DIFF2C2(F,s) ((+1*(F)[-1*s] +		\
		       -2*(F)[+0*s] +		\
		       +1*(F)[+1*s]) / 1.0)

#define DIFF2C4(F,s) ((-1 *(F)[-2*s] +		\
		       +16*(F)[-1*s] +		\
		       -30*(F)[+0*s] +		\
		       +16*(F)[+1*s] +		\
		       -1 *(F)[+2*s]) / 12.0)

#define DIFF4C2(F,s) ((+1*(F)[-2*s] +		\
		       -4*(F)[-1*s] +		\
		       +6*(F)[+0*s] +		\
		       -4*(F)[+1*s] +		\
		       +1*(F)[+2*s]) / 1.0)

#define DIFF4C4(F,s) ((-1 *(F)[-3*s] +		\
		       +12*(F)[-2*s] +		\
		       -39*(F)[-1*s] +		\
		       +56*(F)[+0*s] +		\
		       -39*(F)[+1*s] +		\
		       +12*(F)[+2*s] +		\
		       -1 *(F)[+3*s]) / 6.0)

#define DIFF6C2(F,s) (( 1 *(F)[-3*s] +		\
		       -6 *(F)[-2*s] +		\
		       +15*(F)[-1*s] +		\
		       -20*(F)[+0*s] +		\
		       +15*(F)[+1*s] +		\
		       -6 *(F)[+2*s] +		\
		        1 *(F)[+3*s]) / 1.0)


#define CROSS(E,B) {0.0,				\
      (E)[2]*(B)[3]-(E)[3]*(B)[2],			\
      (E)[3]*(B)[1]-(E)[1]*(B)[3],			\
      (E)[1]*(B)[2]-(E)[2]*(B)[1]}			\

#define DOT(E,B) ((E)[1]*(B)[1] + (E)[2]*(B)[2] + (E)[3]*(B)[3])

#define MAX3(a,b,c) (a>=b && b>=c ? a : (b >= c ? b : c))
#define MIN3(a,b,c) (a<=b && b<=c ? a : (b <= c ? b : c))


/* https://einsteintoolkit.org/documentation/ThornDoc/CactusNumerical/Dissipation */
#if (FFE_DISSIPATION_ORDER == 2)
#define KOV(F,c) ((Ni==1 ? 0.0 : DIFF2C2(F+m+c,si)) +	\
		  (Nj==1 ? 0.0 : DIFF2C2(F+m+c,sj)) +	\
		  (Nk==1 ? 0.0 : DIFF2C2(F+m+c,sk)))
#define KOS(F  ) ((Ni==1 ? 0.0 : DIFF2C2(F+n  ,ti)) +	\
		  (Nj==1 ? 0.0 : DIFF2C2(F+n  ,tj)) +	\
		  (Nk==1 ? 0.0 : DIFF2C2(F+n  ,tk)))
#elif (FFE_DISSIPATION_ORDER == 4)
#define KOV(F,c) ((Ni==1 ? 0.0 : DIFF4C2(F+m+c,si)) +	\
		  (Nj==1 ? 0.0 : DIFF4C2(F+m+c,sj)) +	\
		  (Nk==1 ? 0.0 : DIFF4C2(F+m+c,sk)))
#define KOS(F  ) ((Ni==1 ? 0.0 : DIFF4C2(F+n  ,ti)) +	\
		  (Nj==1 ? 0.0 : DIFF4C2(F+n  ,tj)) +	\
		  (Nk==1 ? 0.0 : DIFF4C2(F+n  ,tk)))
#elif (FFE_DISSIPATION_ORDER == 6)
#define KOV(F,c) ((Ni==1 ? 0.0 : DIFF6C2(F+m+c,si)) +	\
		  (Nj==1 ? 0.0 : DIFF6C2(F+m+c,sj)) +	\
		  (Nk==1 ? 0.0 : DIFF6C2(F+m+c,sk)))
#define KOS(F  ) ((Ni==1 ? 0.0 : DIFF6C2(F+n  ,ti)) +	\
		  (Nj==1 ? 0.0 : DIFF6C2(F+n  ,tj)) +	\
		  (Nk==1 ? 0.0 : DIFF6C2(F+n  ,tk)))
#endif


#if (FFE_DIFFERENCE_ORDER == 2)

#define D1(F,c)  (Ni==1 ? 0.0 : DIFF1C2(F+m+c,si)/dx)
#define D2(F,c)  (Nj==1 ? 0.0 : DIFF1C2(F+m+c,sj)/dy)
#define D3(F,c)  (Nk==1 ? 0.0 : DIFF1C2(F+m+c,sk)/dz)

#define S1(F  )  (Ni==1 ? 0.0 : DIFF1C2(F+n+0,ti)/dx)
#define S2(F  )  (Nj==1 ? 0.0 : DIFF1C2(F+n+0,tj)/dy)
#define S3(F  )  (Nk==1 ? 0.0 : DIFF1C2(F+n+0,tk)/dz)

#elif (FFE_DIFFERENCE_ORDER == 4)

#define D1(F,c)  (Ni==1 ? 0.0 : DIFF1C4(F+m+c,si)/dx)
#define D2(F,c)  (Nj==1 ? 0.0 : DIFF1C4(F+m+c,sj)/dy)
#define D3(F,c)  (Nk==1 ? 0.0 : DIFF1C4(F+m+c,sk)/dz)

#define S1(F  )  (Ni==1 ? 0.0 : DIFF1C4(F+n+0,ti)/dx)
#define S2(F  )  (Nj==1 ? 0.0 : DIFF1C4(F+n+0,tj)/dy)
#define S3(F  )  (Nk==1 ? 0.0 : DIFF1C4(F+n+0,tk)/dz)

#endif



#endif /* FFE_HEADER */
