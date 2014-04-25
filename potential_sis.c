/**
 * @file potential_sis.c
 * @brief Function to calculate the potential for a softened 
 * isothermal sphere
 * @author Tomas E. Tecce
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"

/**
 * @brief Calculate the potential for a softened isothermal sphere at a
 * given point in space.
 * @param x X component of position vector
 * @param y Y component of position vector
 * @param z Z component of position vector
 * @param vh2 Square of circular velocity of the potential
 * @param rh Characteristic radius of the potential
 * @param eps Gravitational softening
 */
double potential_sis(double x, double y, double z, double vh2, double rh, 
		     double eps)
{
  double r, r2;

  if (eps <= 0) {
    fprintf(stderr, "Error (potential_sis): invalid value %g for "
	    "softening\n", eps); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  r2 = x*x + y*y + z*z;
  r  = sqrt(r2);

  return vh2*(log(r/rh) + 0.5*log(1 + r2/(eps*eps)) - log(r/eps) + 
	      eps*atan(r/eps)/r);
}
