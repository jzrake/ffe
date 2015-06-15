#include <math.h>
#include "ffe.h"



#define UNIT_INTERVAL(a) fmod(fmod(a, 1) + 1, 1)
static int move_particle(struct ffe_particle *p, double dt);




void ffe_par_move(struct ffe_sim *sim)
{
  int np = sim->num_particles / cow_domain_getcartsize(sim->domain);
  double dt = sim->status.time_step;
  for (int n=0; n<np; ++n) {
    struct ffe_particle *p = &sim->particles[n];
    move_particle(p, dt);
    /* printf("x=[%6.4e %6.4e %6.4e] E=[%6.4e %6.4e %6.4e]\n", */
    /* 	   p->x[1], p->x[2], p->x[3], */
    /* 	   p->E[1], p->E[2], p->E[3]); */
  }
}



void ffe_par_sample(struct ffe_sim *sim)
{
  int np = sim->num_particles / cow_domain_getcartsize(sim->domain);

  double *x = (double *) malloc(3 * np * sizeof(double));
  double *E = NULL;
  double *B = NULL;

  for (int n=0; n<np; ++n) {
    struct ffe_particle *p = &sim->particles[n];
    for (int d=1; d<=3; ++d) {
      x[n * 3 + d - 1] = UNIT_INTERVAL(p->x[d]);
    }
  }

  cow_dfield_setsamplecoords (sim->electric[0], x, np, 3);
  cow_dfield_setsamplecoords (sim->magnetic[0], x, np, 3);
  cow_dfield_setsamplemode   (sim->electric[0], COW_SAMPLE_NEAREST);
  cow_dfield_setsamplemode   (sim->magnetic[0], COW_SAMPLE_NEAREST);
  cow_dfield_sampleexecute   (sim->electric[0]);
  cow_dfield_sampleexecute   (sim->magnetic[0]);
  cow_dfield_getsampleresult (sim->electric[0], &E, NULL, NULL);
  cow_dfield_getsampleresult (sim->magnetic[0], &B, NULL, NULL);

  for (int n=0; n<np; ++n) {
    struct ffe_particle *p = &sim->particles[n];

    for (int d=1; d<=3; ++d) {
      p->E[d] = E[n * 3 + d - 1];
      p->B[d] = B[n * 3 + d - 1];
    }
  }

  free(x);
}













int move_particle(struct ffe_particle *p, double dt)
/* --------------------------------------------------------
 * Relativistic particle mover, attributed to Boris (1970).
 *
 * The algorithm is described in Birdsall & Langdon (1991)
 * Plasma Physics via Computer Simulation, Section 15-4.
 * It is based on a rotation of the (spatial) 4-velocity
 * components around an axis parallel to the magnetic field
 * through an angle
 *
 * theta = -2 arctan( (e B dt) / (2 gamma m c) )
 *
 * The method uses the time-centered 4-velocity, such that
 * p->E and p->B are known on integral time steps t^{n},
 * while p->u is computed for time steps t^{n+1/2}. Note
 * that the positions are known at integral time steps:
 *
 * x^{n+1} = x^{n} + v^{n+1/2} * dt
 *
 * Author: Jonathan Zrake, zrake@nyu.edu
 *
 * --------------------------------------------------------
 */
{
  double t[4], s[4];
  double u_nph[4], u_nmh[4];
  double u_plus[4], u_minus[4], u_prime[4];
  double gamma_n, gamma_nph;
  double h = (p->e/p->m) * dt;

  for (int d=1; d<=3; ++d)
    u_nmh[d] = p->u[d];

  for (int d=1; d<=3; ++d)                        /* Step 1 */
    u_minus[d] = u_nmh[d] + p->E[d] * 0.5*h;

  gamma_n = sqrt( 1.0 + DOT(u_minus, u_minus) );  /* Step 2 */

  for (int d=1; d<=3; ++d)                        /* Step 3 */
    t[d] = 0.5*h * p->B[d] / gamma_n;

  for (int d=1; d<=3; ++d)                        /* Step 4 */
    s[d] = 2.0*t[d] / ( 1.0 + DOT(t, t) );

  double u_minus_cross_t[4] = CROSS(u_minus, t);

  for (int d=1; d<=3; ++d)                        /* Step 5 */
    u_prime[d] = u_minus[d] + u_minus_cross_t[d];

  double u_prime_cross_s[4] = CROSS(u_prime, s);

  for (int d=1; d<=3; ++d)                        /* Step 6 */
    u_plus[d] = u_minus[d] + u_prime_cross_s[d];

  for (int d=1; d<=3; ++d)                        /* Step 7 */
    u_nph[d] = u_plus[d] + p->E[d] * 0.5*h;

  gamma_nph = sqrt( 1.0 + DOT(u_nph, u_nph) );

  p->u[0] = gamma_nph;
  for (int d=1; d<=3; ++d)                        /* Update the velocity */
    p->u[d] = u_nph[d];

  p->x[0] += dt;

  for (int d=1; d<=3; ++d)                        /* Update the position */
    p->x[d] += dt * u_nph[d] / gamma_nph;

  return 0;
}

