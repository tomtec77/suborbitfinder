/**
 * @file init.c
 * @brief Initialisation functions
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


/**
 * @brief Determine output times either from an output list file or from
 * the specified time between snapshots.
 */
void init(void)
{
  FILE *fd;
  char buf[111];
  int snapshot;
  double a;
  double aexp[MaxSnapshot];

  if (OutputListOn == 1) { /* Reading outputs time from external file */
    sprintf(buf, "%s%s", path3, NameOutputsSelection);

    printf("Reading `%s'\n",buf);
    if (!(fd=fopen(buf,"r"))) {
      fprintf(stderr, 
	      "Error (init): cannot open output list file '%s': %s (%u)\n", 
	      buf, strerror(errno), errno); fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
    for (snapshot=0; snapshot<=MaxSnapshot; snapshot++) {
      fscanf(fd," %lf ",&aexp[snapshot]);
      /*printf("snapshot:%d, aexp: %lf \n",snapshot,aexp[snapshot]);*/
      
      ZZ[snapshot]  = 1.0/aexp[snapshot] - 1;
      Age[snapshot] = time_to_present(ZZ[snapshot]);
      /*
	printf("ZZ[%d] = %f \n",snapshot,ZZ[snapshot]); 
	printf("ZZ[%d] = %f - Age[%d] = %f \n", 
	snapshot,ZZ[snapshot],snapshot,Age[snapshot]); 
      */
    }
  }
  else {
    for (snapshot=0, a=TimeOfFirstSnapshot; snapshot<=MaxSnapshot; 
	 snapshot++, a*=TimeBetSnapshot) {
      ZZ[snapshot]  = 1.0/a - 1;
      Age[snapshot] = time_to_present(ZZ[snapshot]);
      /*   printf("ZZ[%d] = %f \n",snapshot,ZZ[snapshot]); */
      /*   printf("ZZ[%d] = %f - Age[%d] = %f \n",
	   snapshot,ZZ[snapshot],snapshot,Age[snapshot]); */
    }
  }

  Count_orphan_type1 = 0;
  
  return;
}


/**
 * @brief Set some code units.
 */
void set_units(void)
{
  
  /*
   * Default units:
   *   UnitLength_in_cm         = 3.085678e21;
   *   UnitMass_in_g            = 1.989e43; 
   *   UnitVelocity_in_cm_per_s = 1.0e5;
   */
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitTime_in_Megayears = UnitTime_in_s/SEC_PER_MEGAYEAR;
  
  printf("Unit_time_in_Megayears:%e\n", UnitTime_in_Megayears);
  
  /* Gravitational constant in code units */
  G = GRAVITY/pow(UnitLength_in_cm,3)*UnitMass_in_g*pow(UnitTime_in_s,2);
    
  /* Speed of light in code units */
  LightSpeed = C/UnitVelocity_in_cm_per_s;
  
  UnitDensity_in_cgs = UnitMass_in_g/pow(UnitLength_in_cm,3);
  UnitPressure_in_cgs = UnitMass_in_g/UnitLength_in_cm/pow(UnitTime_in_s,2);
  UnitCoolingRate_in_cgs = UnitPressure_in_cgs/UnitTime_in_s;
  
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm,2) / 
    pow(UnitTime_in_s,2);
  
  UnitMass_in_Msun = UnitMass_in_g/SOLAR_MASS;
  UnitLength_in_kpc = 1000.0*UnitLength_in_cm/CM_PER_MPC;
  UnitVelocity_in_km_per_s = UnitVelocity_in_cm_per_s/1.0e5;
  
  /* Convert some physical input parameters to internal units */
  Hubble = HUBBLE * UnitTime_in_s;
  
  /* Compute a few quantitites */
  RhoCrit = 3*Hubble*Hubble/(8*M_PI*G); /* Critical density at z=0 */
  
  return;
}
