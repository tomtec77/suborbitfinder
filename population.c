/**
 * @file population.c
 * @brief Functions to set up the initial galaxy population, and to move
 * the population from current to previous.
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


/**
 * @brief Move the galaxy population from A (current) to B (previous).
 */
void move_popA_to_popB(void)
{
  int i;
  
  
  for (i=1; i<=TotNumGalA; i++) {
    GalaxyB[i] = GalaxyA[i];
  }

  for (i=1; i<=Ngroups_A; i++) {
    NumGalInFOFGroup_B[i] = NumGalInFOFGroup_A[i];
    FirstGalInFOFGroup_B[i] = FirstGalInFOFGroup_A[i];
  }
  
  for (i=1; i<=Nsubgroups_A; i++) {
    NumGalInSubGroup_B[i] = NumGalInSubGroup_A[i];
    FirstGalInSubGroup_B[i] = FirstGalInSubGroup_A[i];
  }

  TotNumGalB = TotNumGalA;
  TotNumType0B = TotNumType0A;
  TotNumType1B = TotNumType1A;
  TotNumType2B = TotNumType2A;
  TotNumType3B = TotNumType3A;
  TotNumType4B = TotNumType4A;
  Ngroups_B = Ngroups_A;
  Nsubgroups_B = Nsubgroups_A;
  
  return;
}


/**
 * @brief Initialise the galaxy population at the first snapshot.
 */
void initialize_galaxy_population(void)
{
  char catalogue_fname[200], properties_fname[200];
  char substruc_fname[200],  volatile_fname[200], subids_fname[200];
  int  i, j, ngal, grB, check;
  int  fofgroup, magbin;
  int iele, ns;
  int central;
  long long bytes, tot_bytes=0;
  float fac;
  
  
  printf("Initializing subhaloes...\n"); fflush(stdout);

  fac = 1.0 / (1.0 + ZZ[0]);

  sprintf(catalogue_fname, "%s%s_%03d", path, NameCatalogue, 0);
  read_groups_B(catalogue_fname);
  
  sprintf(catalogue_fname,"%s%s_%03d", path, NameSublist, 0);
  sprintf(substruc_fname, "%s%s_%03d", path, NameSubstructures, 0);
  sprintf(subids_fname,   "%s%s_%03d", path, NameSubids, 0);
  sprintf(volatile_fname, "%s%s_%03d", path, NameVolatile, 0);
  read_subgroups_B(catalogue_fname, substruc_fname, volatile_fname, 
		   subids_fname);
  
  sprintf(properties_fname, "%s%s_%03d", path, NameSubproperties, 0);
  read_subgroup_properties_B(properties_fname, ZZ[0], subids_fname);
  
  sprintf(properties_fname, "%s%s_%03d", path, NameProperties, 0);
  read_group_properties_B(properties_fname, ZZ[0]);
   
  /* Allocate memory for galaxies */
  if (!(GalaxyB=malloc(bytes=sizeof(struct GALAXY)*(MAXGALAXIES+1)))) {
    fprintf(stderr,
	    "Error (initialize_galaxy_population): memory alloc failure. "
	    "(1)\n");
    fprintf(stderr,
	    "bytes=%d%09d sizeof(struct GALAXY)=%d MAXGALAXIES=%d %d \n",
	    (int) (bytes / 1000000000),
	    (int) (bytes % 1000000000),
	    (int)sizeof(struct GALAXY), MAXGALAXIES, (int)sizeof(size_t));
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  tot_bytes += bytes;
  printf("bytes=%d%09d sizeof(struct GALAXY)=%d MAXGALAXIES=%d %d \n",
         (int) (bytes / 1000000000),
         (int) (bytes % 1000000000),
         (int)sizeof(struct GALAXY), MAXGALAXIES, (int)sizeof(size_t));
  fflush(stdout);
  
  if (!(NumGalInSubGroup_B = malloc(bytes=sizeof(int)*(MAXGALAXIES+1)))) {
    fprintf(stderr,
	    "Error (initialize_galaxy_population): memory alloc failure."
	    "bytes=%d (2)\n", (int)bytes); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  tot_bytes += bytes;
  
  if (!(FirstGalInSubGroup_B = malloc(bytes=sizeof(int)*(MAXGALAXIES+1)))) {
    fprintf(stderr,
	    "Error (initialize_galaxy_population): memory alloc failure."
	    " bytes=%d  (3)\n", (int)bytes); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  tot_bytes += bytes;

  if (!(NumGalInFOFGroup_B = malloc(bytes=sizeof(int)*(MAXGROUPS+1)))) {
    fprintf(stderr,
	    "Error (initialize_galaxy_population): memory alloc failure. "
	    "bytes=%d (4)\n", (int)bytes); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  tot_bytes += bytes;
  
  if (!(FirstGalInFOFGroup_B = malloc(bytes=sizeof(int)*(MAXGROUPS+1)))) {
    fprintf(stderr,
	    "Error (initialize_galaxy_population): memory alloc failure. "
	    "bytes=%d  (5)\n", (int)bytes); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  tot_bytes += bytes;
  
  for (i=1; i<=Ngroups_B; i++)
    NumGalInFOFGroup_B[i] = FirstGalInFOFGroup_B[i] = 0;
  
  for (i=1; i<=Nsubgroups_B; i++)
    NumGalInSubGroup_B[i] = FirstGalInSubGroup_B[i] = 0;
  
  /* Libero groups por si entra al switch del calculo de Lambda */
  free_groups_B();  

  /* Initialise galaxy counters */
  TotNumType0B = TotNumType1B = TotNumType2B = 0;
  TotNumType3B = TotNumType4B = 0;

  /* Start loop over subhaloes (from 1 to Nsubgroups_B; this last value is
     read in function read_subgroups_B) */
  for (grB=1, ngal=0, TotNumGalB=0; grB<=Nsubgroups_B; grB++) {
    
    if (Volatile_B[grB] == 0) {
      TotNumGalB++;
      ngal++;
      
      GalaxyB[ngal].ParentSubGroup = grB;
      GalaxyB[ngal].ParentGroup    = SubParent_B[grB];
      
      FirstGalInSubGroup_B[grB] = ngal; /* Starts at 1 */
      NumGalInSubGroup_B[grB]   = 1;
            
      fofgroup = SubParent_B[grB];
      
      NumGalInFOFGroup_B[fofgroup]++;
      if (FirstGalInFOFGroup_B[fofgroup] == 0)
	FirstGalInFOFGroup_B[fofgroup] = ngal;
      
      /* Initialize galaxy history. For galaxy type, start by filling all
	 snapshots with -1, value updated later */
      for (j=0; j<OUTPUTS; j++) {
	GalaxyB[ngal].TypeT[j] = -1;
	GalaxyB[ngal].MvirT[j] = 0.0;
	GalaxyB[ngal].RvirT[j] = 0.0;
	GalaxyB[ngal].IndexT[j] = 0;
      }

      /* Identify if galaxy is type 0 or 1 */
      if (FirstGalInFOFGroup_B[fofgroup] == ngal) {
	GalaxyB[ngal].Type     = 0;
	GalaxyB[ngal].TypeT[0] = 0;
      }
      else {
	GalaxyB[ngal].Type     = 1;
	GalaxyB[ngal].TypeT[0] = 1;
      }
      
      GalaxyB[ngal].Id = ngal;

      GalaxyB[ngal].Progenitor = 0; /* New galaxy, no progenitor */
      GalaxyB[ngal].Descendant = 0; /* Descendant information is updated
				       at the end of the snapshot, when the
				       mergers have been determined */
      
      GalaxyB[ngal].PaIndex = SubIdMostBound_B[grB];
      GalaxyB[ngal].Len     = SubLen_B[grB];
      
      GalaxyB[ngal].Vc   = SubGroupVc_B[grB];
      GalaxyB[ngal].Vvir = SubGroupVvir_B[grB];
      GalaxyB[ngal].Mvir = SubGroupMvir_B[grB];
      GalaxyB[ngal].Rvir = SubGroupRvir_B[grB];
      
      GalaxyB[ngal].MvirT[0] = GalaxyB[ngal].Mvir;
      GalaxyB[ngal].RvirT[0] = GalaxyB[ngal].Rvir;

      /* Galaxy positions and velocities, in comoving units as in the
	 simulation snapshots. */
      for (j=0; j<3; j++) {
	GalaxyB[ngal].Pos[j] = SubGroupPos_B[grB][j];
	GalaxyB[ngal].Vel[j] = SubGroupVel_B[grB][j];
      }

      /* Initialise relative positions, velocities, orbital angular
	 momentum. Note that here we only have galaxies of types 0 
	 or 1 */
      if (GalaxyB[ngal].Type == 1) {
	central = FirstGalInFOFGroup_B[GalaxyB[ngal].ParentGroup];
	if (central <= 0) {
	  fprintf(stderr, 
		  "Error (initalize_galaxy_population): cannot find a "
		  "central galaxy for galaxy %d - Exit\n", ngal);
	  fflush(stderr);
	  exit(EXIT_FAILURE);
	}
	
	/* Relative quantities are in physical units */
	for (j=0; j<3; j++) {
	  GalaxyB[ngal].Posrel[j] = (GalaxyB[ngal].Pos[j] - 
				     GalaxyB[central].Pos[j])*fac;
	  GalaxyB[ngal].Vrel[j] = (GalaxyB[ngal].Vel[j] - 
				   GalaxyB[central].Vel[j])*sqrt(fac);
	}
	GalaxyB[ngal].Jorb[0] = 
	  (GalaxyB[ngal].Posrel[1]*GalaxyB[ngal].Vrel[2] - 
	   GalaxyB[ngal].Posrel[2]*GalaxyB[ngal].Vrel[1])*pow(fac,1.5);
	GalaxyB[ngal].Jorb[1] = 
	  (GalaxyB[ngal].Posrel[2]*GalaxyB[ngal].Vrel[0] - 
	   GalaxyB[ngal].Posrel[0]*GalaxyB[ngal].Vrel[2])*pow(fac,1.5);
	GalaxyB[ngal].Jorb[2] = 
	  (GalaxyB[ngal].Posrel[0]*GalaxyB[ngal].Vrel[1] - 
	   GalaxyB[ngal].Posrel[1]*GalaxyB[ngal].Vrel[0])*pow(fac,1.5);

	/* Save the initial value of the orbital angular momentum
	   (in physical units, squared) */
	GalaxyB[ngal].J2init = 
	  GalaxyB[ngal].Jorb[0]*GalaxyB[ngal].Jorb[0] +
	  GalaxyB[ngal].Jorb[1]*GalaxyB[ngal].Jorb[1] +
	  GalaxyB[ngal].Jorb[2]*GalaxyB[ngal].Jorb[2];
      
      } /* Close type 1 galaxy */
      else {
	/* For a central, set all relative quantities to 0 */
	for (j=0; j<3; j++) {
	  GalaxyB[ngal].Posrel[j] = 0.0;
	  GalaxyB[ngal].Vrel[j]   = 0.0;
	  GalaxyB[ngal].Jorb[j]   = 0.0;
	}
	GalaxyB[ngal].J2init = 0.0;
      }
      
      /* Initialise other quantities, this is the same for both type 0
	 and 1 galaxies */
      GalaxyB[ngal].MergerType = 0;
      GalaxyB[ngal].Relocated  = 0;

      GalaxyB[ngal].MergTime       = 0.0;
      GalaxyB[ngal].MergTimeBK     = 0.0;
      GalaxyB[ngal].MergTimeJ      = 0.0;
      GalaxyB[ngal].MergTimeActual = 0.0;
      GalaxyB[ngal].MergTimeSat    = 0.0;
      GalaxyB[ngal].Coulomb        = 0.0;
      GalaxyB[ngal].StrippedMass   = 0.0;
      GalaxyB[ngal].EnergyLoss     = 0.0;

      /* Use the virial radius as initial value for the tidal radius and 
	 the DM bounding radius */
      GalaxyB[ngal].Rtidal = GalaxyB[ngal].Rvir;
      GalaxyB[ngal].Rdm    = GalaxyB[ngal].Rvir;

      GalaxyB[ngal].Mdm    = GalaxyB[ngal].Mvir;

      /* For type 1 galaxies, update tidal radius and bounding radius for 
	 DM (not the merging time, as this is the initial situation) */
      if (GalaxyB[ngal].Type == 1) {
	GalaxyB[ngal].Rtidal = get_tidal_radius(ngal,central);
	if (GalaxyB[ngal].Rtidal < GalaxyB[ngal].Rvir)
	  GalaxyB[ngal].Rdm = GalaxyB[ngal].Rtidal;
      }

      for (j=0; j<STEPS; j++) {
	if (GalaxyB[ngal].Type == 0) {
	  GalaxyB[ngal].OrbitX[j]    = 0.0;
	  GalaxyB[ngal].OrbitY[j]    = 0.0;
	  GalaxyB[ngal].OrbitZ[j]    = 0.0;
	  GalaxyB[ngal].OrbitVx[j]   = 0.0;
	  GalaxyB[ngal].OrbitVy[j]   = 0.0;
	  GalaxyB[ngal].OrbitVz[j]   = 0.0;
	  
	  GalaxyB[ngal].Ekin = 0.0;
	  GalaxyB[ngal].Epot = 0.0;
	}
	else {
	  GalaxyB[ngal].OrbitX[j]    = GalaxyB[ngal].Posrel[0];
	  GalaxyB[ngal].OrbitY[j]    = GalaxyB[ngal].Posrel[1];
	  GalaxyB[ngal].OrbitZ[j]    = GalaxyB[ngal].Posrel[2];
	  GalaxyB[ngal].OrbitVx[j]   = GalaxyB[ngal].Vrel[0];
	  GalaxyB[ngal].OrbitVy[j]   = GalaxyB[ngal].Vrel[1];
	  GalaxyB[ngal].OrbitVz[j]   = GalaxyB[ngal].Vrel[2];

	  GalaxyB[ngal].Ekin = 0.5*
	    (GalaxyB[ngal].Vrel[0]*GalaxyB[ngal].Vrel[0] + 
	     GalaxyB[ngal].Vrel[1]*GalaxyB[ngal].Vrel[1] + 
	     GalaxyB[ngal].Vrel[2]*GalaxyB[ngal].Vrel[2]);
	  GalaxyB[ngal].Epot = potential_sis(GalaxyB[ngal].Posrel[0],
					     GalaxyB[ngal].Posrel[1],
					     GalaxyB[ngal].Posrel[2],
					     GalaxyB[ngal].Vvir*
					     GalaxyB[ngal].Vvir,
					     GalaxyB[ngal].Rvir,
					     EpsGrav*GalaxyB[ngal].Rvir);
	    
	}
	GalaxyB[ngal].OrbitMass[j] = GalaxyB[ngal].Mvir;
	GalaxyB[ngal].OrbitType[j] = GalaxyB[ngal].Type;
      }
	
    } /* Close if (Volatile_B == 0) */
  } /* Close for (grB=1, ngal=0, TotNumGalB=0... */
  
  GalaxynumberB = TotNumGalB; /* GMT */
  SubGalMerging = ivector(1,GalaxynumberB);
  for (i=1;i<=GalaxynumberB;i++) SubGalMerging[i] = 0;
  
  free_subproperties_B();
  free_subgroups_B();
  free_properties_B();
  
  printf("Allocated %g MByte for subhaloes.\n", 
	 ((double)tot_bytes)/(1024.0*1024));
  printf("Initialization done ...\n");
  printf("Have got %d subhaloes.\n", TotNumGalB); fflush(stdout);
    
  for (grB=1, check=0; grB<=Ngroups_B; grB++) {
    check += NumGalInFOFGroup_B[grB];
  }
  if (check != TotNumGalB) {
    fprintf(stderr,
	    "Error (initialize_galaxy_population): check=%d "
	    "TotNumGalB=%d\n", check, TotNumGalB ); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  /* Count types of galaxies */
  for (i=1; i<=TotNumGalB; i++) {
    if (GalaxyB[i].Type == 0) TotNumType0B++;
    if (GalaxyB[i].Type == 1) TotNumType1B++;
    if (GalaxyB[i].Type == 2) TotNumType2B++;
    if (GalaxyB[i].Type == 3) TotNumType3B++;
    if (GalaxyB[i].Type == 4) TotNumType4B++;
  }
  
  return;
}
