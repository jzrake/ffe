#ifndef FFE_NLE_HEADER
#define FFE_NLE_HEADER


/*
 * Data structure that describes the 2D non-linear equilibrium state, where the
 * magnetic potenial function is given by
 *
 *                    Az = (cos(x) - sin(y))^(2*n + 1)
 *
 * where the order n is an integer greater than or equal to 0, for which it
 * reduces to the linear ABC field.
 */

struct ffe_nle
{
  double *Atable;
  double *Btable;
  double *Utable;
  int order; /* integer n greater than 0 (0=linear) */
  int num_bins; /* size of the A and B tables */
  int array_size; /* resolution to use for numerical integrals */
} ;


void ffe_nle_null(struct ffe_nle *nle);
void ffe_nle_init(struct ffe_nle *nle, int order, int num_bins, int array_size);
void ffe_nle_free(struct ffe_nle *nle);
void ffe_nle_sample(struct ffe_nle *nle, double x, double y, double B[4]);
void ffe_nle_write_table(struct ffe_nle *nle, const char *fname);


#endif /* FFE_NLE_HEADER */
