/**
 * @file galaxydistance.c
 * @brief Determine distance of galaxy to centre of host subhalo and FOF
 * group
 * @author Tomas E. Tecce
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"

/**
 * @brief Calculate the distance from the galaxy to the centre of its
 * host subhalo and to the centre of the corresponding main FOF halo.
 * @param p Index of selected galaxy
 * @param rsub Pointer to store the subhalocentric radius
 * @param rfof Pointer to store the FOF halocentric radius
 */
void galaxydistance(int p, float *rsub, float *rfof)
{
  int central;
  float fac;
  
  fac = 1.0/(1+ZZ[Snapshot]);

  /* If galaxy is relocated, Posrel is the distance to the FOF central.
     Otherwise it's a distance to the host subhalo centre */
  *rsub = sqrt(GalaxyA[p].Posrel[0]*GalaxyA[p].Posrel[0] +
	       GalaxyA[p].Posrel[1]*GalaxyA[p].Posrel[1] +
	       GalaxyA[p].Posrel[2]*GalaxyA[p].Posrel[2]);

  if (GalaxyA[p].Relocated == 1) {
    central = FirstGalInSubGroup_A[GalaxyA[p].ParentSubGroup];
    *rfof = *rsub;
    *rsub = sqrt(pow(GalaxyA[p].Pos[0]-GalaxyA[central].Pos[0],2) +
		 pow(GalaxyA[p].Pos[1]-GalaxyA[central].Pos[1],2) +
		 pow(GalaxyA[p].Pos[2]-GalaxyA[central].Pos[2],2))*fac;
  }
  else {
    central = FirstGalInFOFGroup_A[GalaxyA[p].ParentGroup];
    *rfof = sqrt(pow(GalaxyA[p].Pos[0]-GalaxyA[central].Pos[0],2) +
		 pow(GalaxyA[p].Pos[1]-GalaxyA[central].Pos[1],2) +
		 pow(GalaxyA[p].Pos[2]-GalaxyA[central].Pos[2],2))*fac;
  }
    
  return;
}
