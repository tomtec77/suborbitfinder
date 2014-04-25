/**
 * @file tidalradius.c
 * @brief Function to calculate the tidal radius of a subhalo
 * @author Tomas E. Tecce
 */
#include "allvars.h"
#include "readparameterfile.h"
#include "proto.h"

/**
 * @brief Function to calculate the tidal radius of the host subhalo of a
 * galaxy.
 * @param p Index of galaxy
 * @param central Index of the corresponding central galaxy
 * @return Tidal radius in main code units. If p is a central, returns the
 * virial radius
 *
 * Assumes a SIS profile for the mass density of the central halo. This
 * version includes different approaches to calculate the tidal radius:
 * either the formula from Tormen, Diaferio & Syer (1998) or the one from
 * Zentner & Bullock (2003). Both are essentially based on the King (1962)
 * approach, but the TDS98 formula neglects a term with the orbital 
 * velocity of the galaxy.
 */
float get_tidal_radius(int p, int central)
{
  int i;
  float r200, m200, v200;
  float msub, x2, v2;
  float pos[3], vel[3];


  if (p == central)
    return GalaxyA[p].Rvir;

  /* Use the coordinates of the satellite relative to the central */
  for (i=0; i<3; i++) {
    pos[i] = GalaxyA[p].Posrel[i];
    vel[i] = GalaxyA[p].Vrel[i];
  }
  
  r200 = GalaxyA[central].Rvir;
  m200 = GalaxyA[central].Mvir;

  /* Note that this is actually is the virial velocity squared */
  v200 = G*GalaxyA[central].Mvir/GalaxyA[central].Rvir;
 
  /* For a type 0 or 1 galaxy, this should be equal to Mvir */
  msub = GalaxyA[p].Mdm;
  if (msub != GalaxyA[p].Mvir) {
    fprintf(stderr, "Error (get_tidal_radius): Mdm=%g != Mvir=%g "
	    "for galaxy %d - Exit\n", msub, GalaxyA[p].Mvir, p);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  msub /= m200;

  /* Square of current radius of galaxy p within its host halo, in units 
     of the virial radius of the host */    
  x2 = (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2])/(r200*r200);

  /* Square of current velocity of galaxy p within the host halo, in 
     units of the virial velocity of the host */
  v2 = (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2])/v200;

  if (RtidalSelect == RTIDAL_TDS98)
    return r200*pow(msub*x2, 1.0/3.0);
  if (RtidalSelect == RTIDAL_ZB03)
    return r200*pow(msub*x2/(1+v2), 1.0/3.0);
}
