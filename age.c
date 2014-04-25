/**
 * @file age.c
 * @brief Functions to calculate time to present
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


/**
 * @brief Calculate time to present as a function of redshift.
 * @param z Redshift to calculate age to present
 * @return Time to present in code units
 */
double time_to_present(double z)
{
  double time;
  double integrand(double a);
  
  time = 1.0/Hubble*dqromb(integrand, 1/(z+1), 1);
  
  return time;
}


/**
 * @brief Integrand in function time_to_present.
 * @param a Expansion factor
 */
double integrand(double a)
{
  return 1.0/sqrt(Omega/a + (1-Omega-OmegaLambda) + OmegaLambda*a*a);
}
