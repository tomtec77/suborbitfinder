/**
 * @file suborbitfinder.c
 * @brief Integrate orbits for DM subhaloes below the parent simulation's
 * resolution limit
 * @author Tomas E. Tecce
 */
#include "allvars.h"
#include "readparameterfile.h"
#include "proto.h"

int main(int argc, char **argv)
{
  FILE *fd;
  clock_t c0, c1;
  time_t t0, t1;
  struct tm *loctime;
  char hostname[256];
  char catalogue_fname[FILENAME_MAX];
  char history_fname[FILENAME_MAX];
  char input_fname0[FILENAME_MAX];
  char orderedDM_fname[FILENAME_MAX];
  char properties_fname[FILENAME_MAX];
  char substruc_fname[FILENAME_MAX];
  char subids_fname[FILENAME_MAX];
  char volatile_fname[FILENAME_MAX];
  char orbits_fname[FILENAME_MAX];
  int i, j, p, central, fof, indcen, paindex;
  int *groupCopyArray;


  t0 = time(NULL);
  c0 = clock();

  printf("#######\n"
	 "# suborbitfinder\n"
	 "#######\n\n");
  gethostname(hostname,256);
  loctime = localtime(&t0);
  printf("Running on %s - Run started %s\n", hostname, asctime(loctime));
  fflush(stdout);

  strcpy(PreTag, "");

  if (argc != 2 && argc != 3) {
    printf("  usage: suborbitfinder parameterfile [itf]\n");
    printf("  parameterfile    see readparameterfile.c\n");
    printf("  [itf]            optional string 'itf'\n\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  if (argc == 3)
    strcpy(PreTag, "itf_");
  else
    PreTag[0] = 0;

  readparameterfile(*(argv+1));

  if ((MaxSnapshot+1) > OUTPUTS) {
    fprintf(stderr, "The code is not compiled for this OUTPUTS number.\n"
	    "OUTPUTS: %d, MaxSnapshot+1: %d\n"
	    "Enter the correct value in allvars.h and recompile, or check\n"
	    "the parameters file and run again. Stop\n", 
	    OUTPUTS, MaxSnapshot+1); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  /* 
   * Read header from snapshot 0 to get number of particles 
   * and the particle mass 
   */
  sprintf(input_fname0, "%s%s%s%03d", path1, NameSnapshot, filebasename, 0);
  
#ifdef DOLAGSIMULATION
  read_header_mass(input_fname0, Files);
#else
  read_header(input_fname0, Files);
#endif
  
  if (RamPressureOn != 0)  GPMass = PartMassGas;
  
  PartMass = PartMassDM + PartMassGas;
  /* Assume that each DM particle also includes the baryon component
     to take into account the total mass.
     The baryon component is 0 if simulation is DM-only */
  
  printf("\n*******\n\n");
  printf("PartMassDM                            =%g\n", PartMassDM);
  printf("PartMassGas                           =%g\n", PartMassGas);
  printf("PartMass                              =%g\n", PartMass);
  printf("NumPart                               =%d\n", NumPart);
  fflush(stdout);
  /* PartMass takes into account the mass of gas and dark matter 
     that each particle contains */

  set_units();
  init();
  initialize_galaxy_population();

  for (Snapshot = 1; Snapshot <= MaxSnapshot; Snapshot++) {
  //for (Snapshot = 1; Snapshot <= 36; Snapshot++) {

    /* TARGET: for box #1, the first galaxies appear in snapshot 20
       at z = 10.35. First type 2 galaxy appears in snapshot 22, at
       z = 9.68 */
    
    printf("\nSnapshot %d, z = %g\n********\n", Snapshot, ZZ[Snapshot]);
    fflush(stdout);

    /*
     * Read postprocessing files
     */
    sprintf(catalogue_fname, "%s%s_%s%03d", path, NameSublist, PreTag,
	    Snapshot);
    read_groups_A(catalogue_fname);

    sprintf(substruc_fname, "%s%s_%s%03d", path, NameSubstructures, PreTag, 
	    Snapshot);
    sprintf(subids_fname, "%s%s_%s%03d", path, NameSubids, PreTag, 
	    Snapshot);
    if (Snapshot < MaxSnapshot)
      sprintf(volatile_fname, "%s%s_%s%03d", path, NameVolatile, PreTag, 
	      Snapshot);
    else
      *volatile_fname = 0;
    read_subgroups_A(catalogue_fname, substruc_fname, volatile_fname, 
		     subids_fname);
    
    sprintf(catalogue_fname, "%s%s_%s%03d", path, NameSublist, PreTag, 
	    Snapshot-1);
    read_groups_B(catalogue_fname);

    sprintf(substruc_fname, "%s%s_%s%03d", path, NameSubstructures, PreTag, 
	    Snapshot-1);
    sprintf(subids_fname, "%s%s_%s%03d", path, NameSubids, PreTag, 
	    Snapshot-1);
    sprintf(volatile_fname, "%s%s_%s%03d", path, NameVolatile, PreTag, 
	    Snapshot-1);
    read_subgroups_B(catalogue_fname, substruc_fname, volatile_fname, 
		     subids_fname);
    
    sprintf(properties_fname, "%s%s_%s%03d", path, NameSubproperties, 
	    PreTag, Snapshot);
    sprintf(subids_fname, "%s%s_%s%03d", path, NameSubids, PreTag, 
	    Snapshot);
    read_subgroup_properties_A(properties_fname, ZZ[Snapshot], 
			       subids_fname);
    
    sprintf(properties_fname, "%s%s_%s%03d", path, NameSubproperties, 
	    PreTag, Snapshot-1);
    sprintf(subids_fname, "%s%s_%s%03d", path, NameSubids, PreTag, 
	    Snapshot-1);
    read_subgroup_properties_B(properties_fname, ZZ[Snapshot-1], 
			       subids_fname);
    
    sprintf(history_fname, "%s%s_%s%03d", path, NameSubhistory, PreTag, 
	    Snapshot);
    read_subhistory(history_fname);   /* Read merger tree */
    
    sprintf(properties_fname, "%s%s_%03d", path, NameProperties, Snapshot);
    read_group_properties_A(properties_fname, ZZ[Snapshot]);
    
    sprintf(properties_fname, "%s%s_%03d", path, NameProperties, Snapshot-1);
    read_group_properties_B(properties_fname, ZZ[Snapshot-1]);
    
    sprintf(history_fname, "%s%s_%03d", path, NameHistory, Snapshot);
    read_history(history_fname);
    
    Zprev  = ZZ[Snapshot-1];         /* Start of iz'th redshift interval */
    Zcurr  = ZZ[Snapshot];             /* End of iz'th redshift interval */
    DeltaT = time_to_present(Zprev) - time_to_present(Zcurr);
    printf("DeltaT: %g Myr\n\n", DeltaT*UnitTime_in_Megayears);
    fflush(stdout);
    
    /*
     * Read DM particles data
     */
    /* TOMAS 2012-07-20: 
       This may not be necessary anymore */
    /*allocate_memory();

    strcpy(OrdTag, "DMord");
    sprintf(orderedDM_fname, "%s%s%s%s_%03d", path1, NameDMordered,
	    filebasename, OrdTag, Snapshot);
    if (!(fd=fopen(orderedDM_fname, "r"))) {
      printf("Error (main): cannot open file '%s' for input of ordered "
	     "DM particles\n", orderedDM_fname); fflush(stdout);
      exit(EXIT_FAILURE);
    }
    printf("Reading from file with %d ordered DM particles: %s\n", 
	   NumPart, orderedDM_fname); fflush(stdout);
    if (DMPosAvailableOn == 1) {
      for (p=1; p<=NumPart; p++)
	fread(&P[p].Pos[0], sizeof(float), 3, fd);
      for (p=1; p<=NumPart; p++)
	fread(&P[p].Vel[0], sizeof(float), 3, fd);
    }
    else {
      printf("Error (main): DM positions not available in file, according"
	     " to parameters file - Exit\n"); fflush(stdout);
      exit(EXIT_FAILURE);
      }*/
    
    generate_new_population();

    /* Update histories */
    for (p=1; p<=TotNumGalA; p++) {
      GalaxyA[p].TypeT[Snapshot] = GalaxyA[p].Type;
      
      GalaxyA[p].MvirT[Snapshot] = GalaxyA[p].Mvir;
      GalaxyA[p].RvirT[Snapshot] = GalaxyA[p].Rvir;

      GalaxyA[p].J2initT[Snapshot] = GalaxyA[p].J2init;

      GalaxyA[p].IndexT[Snapshot] = p;

      if (GalaxyA[p].J2init > 0) {
	GalaxyA[p].JorbT[Snapshot] = 
	  sqrt((GalaxyA[p].Jorb[0]*GalaxyA[p].Jorb[0] +
		GalaxyA[p].Jorb[1]*GalaxyA[p].Jorb[1] +
		GalaxyA[p].Jorb[2]*GalaxyA[p].Jorb[2]) / 
	       GalaxyA[p].J2init);
	if (GalaxyA[p].JorbT[Snapshot] > 1.01) {
	  fprintf(stderr, "Error: J/J_0 = %g > 1 for galaxy %d in "
		  "snapshot %d\n", GalaxyA[p].JorbT[Snapshot], p, Snapshot);
	  fprintf(stderr, "Type = %d J2init = %g\n", GalaxyA[p].Type,
		  GalaxyA[p].J2init);
	  fflush(stderr);
	  //exit(EXIT_FAILURE);
	}
      }
      else
	GalaxyA[p].JorbT[Snapshot] = 0.0;

      GalaxyA[p].MdmT[Snapshot] = GalaxyA[p].Mdm;
      GalaxyA[p].RdmT[Snapshot] = GalaxyA[p].Rdm;

      galaxydistance(p, &GalaxyA[p].SubRadiusT[Snapshot],
		     &GalaxyA[p].FOFRadiusT[Snapshot]);

      /*GalaxyA[p].RelocatedT[Snapshot] = GalaxyA[p].Relocated;*/
      if (SatelliteRelocationDisabledOn == 1 && GalaxyA[p].Relocated == 1) {
	central = FirstGalInSubGroup_A[GalaxyA[p].ParentSubGroup];
	fprintf(stderr, "Error (main): relocation flag set for galaxy %d, "
		"type %d (central %d, type %d) - Exit\n", 
		p, GalaxyA[p].Type, central, GalaxyA[central].Type);
	fflush(stderr);
	exit(EXIT_FAILURE);
      }
    } /* Close for (p=1; p<=TotNumGalA... */

    /* Save galaxies data to file */
    sprintf(orbits_fname, "%s%s_%s%03d.hdf5", path, NameSuborbits, PreTag, 
	    Snapshot);
    dumpHDF5_popA(orbits_fname);
    
    /* Now erase merged galaxies and update pointer information */
    groupCopyArray = ivector(1,Ngroups_A+1);
    for (i=1; i<=Ngroups_A; i++) groupCopyArray[i] = 0;

    for (p=1; p<=TotNumGalA; p++) {
      if (GalaxyA[p].Type > 2) {
	groupCopyArray[GalaxyA[p].ParentGroup+1] += 1;
	NumGalInFOFGroup_A[GalaxyA[p].ParentGroup]--;
      }
    }
    
    for (i=1,j=0; i<=Ngroups_A; i++) {
      j = j + groupCopyArray[i];
      if (NumGalInFOFGroup_A[i]) FirstGalInFOFGroup_A[i] -= j;
    }
    free_ivector(groupCopyArray,1,Ngroups_A+1);

    groupCopyArray = ivector(1,Nsubgroups_A+1); /* Update subgroup pointer
						   information */
    for (i=1; i<=Nsubgroups_A; i++) groupCopyArray[i] = 0;
    
    for (p=1; p<=TotNumGalA; p++) {
      if (GalaxyA[p].Type > 2) {
	groupCopyArray[GalaxyA[p].ParentSubGroup+1] += 1;
	NumGalInSubGroup_A[GalaxyA[p].ParentSubGroup]--;
      }
    }

    for (i=1,j=0; i<=Nsubgroups_A; i++) {
      j = j + groupCopyArray[i];
      if (NumGalInSubGroup_A[i]) FirstGalInSubGroup_A[i] -= j;
    }
    free_ivector(groupCopyArray,1,Nsubgroups_A+1);

    /* Copy the non-merged / non-disrupted galaxies */
    for (p=1, j=0; p<=TotNumGalA; p++) {
      if (GalaxyA[p].Type <= 2) {
	j = j+1;
	GalaxyA[j] = GalaxyA[p];
      }
    }
    TotNumGalA = j;

    /* Clear memory */
    free_subhistory();
    free_subproperties_B(); 
    free_subproperties_A(); 
    free_subgroups_B();
    free_subgroups_A();
    free_history();
    
    /*free_memory();
    Id++;
    free(Id);
    printf("Space for particle ID freed.\n"); fflush(stdout);*/

    move_popA_to_popB();

  } /* Close for (Snapshot = 1... */

  printf("Total orphan Type 1 galaxies found: %d\n", Count_orphan_type1);
  
  t1 = time(NULL);
  c1 = clock();
  printf("\nElapsed wall clock time: %ld\n", (long)(t1-t0));
  printf("Elapsed CPU time:        %f\n", (float)(c1-c0)/CLOCKS_PER_SEC);
  fflush(stdout);

  return EXIT_SUCCESS;
}


/**
 * @brief Generate the galaxy population for the current snapshot.
 *
 * When adding new data to struct GALAXY, remember that it also has to be
 * included in function initialize_galaxy_population() in file population.c.
 */
void generate_new_population(void)
{
  int grA, grB;
  int central;
  int mother_in_B, satellite_in_B;
  int i, j, g, n, numgal, vol;
  int progcount, galprogcount, newgalaxies, oldgalaxies;
  int type0to1, type1to0, type0to2, type1to2, type2to1, satflag;
  int nmergers, ndisrupt, norb;
  int status;
  float rsat;
  float feps;
  float hubble_z, fac, fac_0;
  float timetype;
  float jorb0;
  double eps, epsgrv, tstep;


  printf("Now generating new subhalo population...\n"); fflush(stdout);

  fac = 1.0/(1.0+ZZ[Snapshot]);      /* Current */
  fac_0 = 1.0/(1.0+ZZ[Snapshot-1]);  /* At previous step */

  if ((Ngroups_A + Ngroups_B) > MAXGROUPS) {
    fprintf(stderr,
	    "Already Ngroups_A + Ngroups_B = %d exceeds MAXGROUPS = %d\n",
	    Ngroups_A+Ngroups_B, MAXGROUPS);
    fprintf(stderr, "Ngroups_A =%d, Ngroups_B = %d\n", Ngroups_A,Ngroups_B);
    fprintf(stderr, "Nsubgroups_A = %d, Nsubgroups_B = %d\n",
	   Nsubgroups_A, Nsubgroups_B); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  if ((Nsubgroups_A + Nsubgroups_B) > MAXGALAXIES) {
    fprintf(stderr, 
	    "Already Nsubgroups_A + Nsubgroups_B = %d exceeds MAXGALAXIES = "
	    "%d\n", Nsubgroups_A+Nsubgroups_B, MAXGALAXIES);
    fprintf(stderr,"Ngroups_A = %d, Ngroups_B = %d\n", Ngroups_A,Ngroups_B);
    fprintf(stderr, "Nsubgroups_A = %d, Nsubgroups_B = %d\n",
	    Nsubgroups_A, Nsubgroups_B); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  printf("TotNumGalB = %d\n", TotNumGalB);
  printf("At time step B: type 0 = %d, type 1 = %d, type 2 = %d\n", 
	 TotNumType0B, TotNumType1B, TotNumType2B);
  fflush(stdout);
  
  NumGalInSubGroup_A   = NumGalInSubGroup_B + Nsubgroups_B;
  FirstGalInSubGroup_A = FirstGalInSubGroup_B + Nsubgroups_B;
  NumGalInFOFGroup_A   = NumGalInFOFGroup_B + Ngroups_B;
  FirstGalInFOFGroup_A = FirstGalInFOFGroup_B + Ngroups_B;

  GalaxyA = GalaxyB + TotNumGalB;
  
  for (i=1; i<=Ngroups_A; i++)
    NumGalInFOFGroup_A[i] = FirstGalInFOFGroup_A[i] = 0;
  
  for (i=1; i<=Nsubgroups_A; i++)
    NumGalInSubGroup_A[i] = FirstGalInSubGroup_A[i]= 0;
  
  /* Now generate the new population */
  newgalaxies  = 0;
  oldgalaxies  = 0;
  progcount    = 0;
  galprogcount = 0;
  type0to1     = 0;
  type1to0     = 0;
  type0to2     = 0;
  type1to2     = 0;
  type2to1     = 0;
  
  printf("Ngroups_A    = %d\n", Ngroups_A);
  printf("Nsubgroups_A = %d\n", Nsubgroups_A);
  
  TotNumType0A = TotNumType1A = TotNumType2A = 0;
  TotNumType3A = TotNumType4A = 0;
  for (grA=1, numgal=0, TotNumGalA=0, vol=0; grA<=Nsubgroups_A; grA++) {
    NumGalInSubGroup_A[grA] = 0;
    
    if (Volatile_A[grA] > 0) {
      vol++;
      if (CountProgenitor[grA] > 0) {
	progcount++;
	for (n=0, grB=FirstProgenitor[grA]; n<CountProgenitor[grA]; 
	     n++, grB=NextProgenitor[grB]) {
	  for (i=0; i<NumGalInSubGroup_B[grB]; i++) {
	    galprogcount++;
	  }
	}                 
      }
    } /* Close if (Volatile_A > 0) */
    
    if (Volatile_A[grA] == 0) {
      if (CountProgenitor[grA] == 0) {      /* A new galaxy has appeared */
	
	/*
	 * Note that a new galaxy must reside in a DM subhalo, so
	 * here we are dealing either with type 0 or type 1 galaxies
	 */
	
	newgalaxies++;
	numgal++;
	TotNumGalA++;
	
	NumGalInSubGroup_A[grA]   = 1;
	FirstGalInSubGroup_A[grA] = numgal;
	
	NumGalInFOFGroup_A[SubParent_A[grA]]++;
	
	if (FirstGalInFOFGroup_A[SubParent_A[grA]] == 0)
	  FirstGalInFOFGroup_A[SubParent_A[grA]] = numgal;
	
	if ((TotNumGalA+TotNumGalB) > MAXGALAXIES) {
	  fprintf(stderr, 
		  "TotNumGalA + TotNumGalB = %d exceeds MAXGALAXIES = %d\n",
		  TotNumGalA+TotNumGalB, MAXGALAXIES); fflush(stderr);
	  exit(EXIT_FAILURE);
	}
	
	/*  Here: initialize new galaxy */
	if ((FirstEntry_A[SubParent_A[grA]]+1) == grA) {
	  if (FirstGalInFOFGroup_A[SubParent_A[grA]] != numgal) {
	    fprintf(stderr, "Error (generate_new_population): grA = %d "
		    "SubParent_A[grA] = %d\n"
		    "FirstGalInFOFGroup_A[SubParent_A[grA]] = %d\n"
		    "This should match numgal = %d - Exit\n",
		    grA, SubParent_A[grA], 
		    FirstGalInFOFGroup_A[SubParent_A[grA]], numgal);
	    fflush(stderr);
	    exit(EXIT_FAILURE);
	  }
	  
	  GalaxyA[numgal].Type = 0;
	  GalaxyA[numgal].TypeT[Snapshot] = 0;
	}
	else {
	  GalaxyA[numgal].Type = 1;
	  GalaxyA[numgal].TypeT[Snapshot] = 1;
	}

	/* For new galaxies, fill the type history with -1 for 
	   snapshots previous to the formation time */
	for (j=0; j<Snapshot; j++) {
	  GalaxyA[numgal].TypeT[j]   = -1;
	  GalaxyA[numgal].J2initT[j] = 0.0;
	  GalaxyA[numgal].MvirT[j]   = 0.0;
	  GalaxyA[numgal].RvirT[j]   = 0.0;
	  GalaxyA[numgal].IndexT[j]  = 0;
	}

	/* 
	 * NOTE: in subgalaxyfinder this Id needs to be zero here, 
	 * since that is later used when constructing the galaxy history
	 * tree
	 */
	GalaxyA[numgal].Id = 0;
	
	GalaxyA[numgal].Progenitor = 0; /* New galaxy, no progenitor */
	GalaxyA[numgal].Descendant = 0;
	
	/* Set the parent halo and subhalo at the present time */
	GalaxyA[numgal].ParentSubGroup = grA;
	GalaxyA[numgal].ParentGroup    = SubParent_A[grA];
	GalaxyA[numgal].PaIndex        = SubIdMostBound_A[grA];
	GalaxyA[numgal].Len            = SubLen_A[grA];
	
	/* Set subhalo virial mass, radius and velocity, and the
	   subhalo position and velocity */
	GalaxyA[numgal].Mvir = SubGroupMvir_A[grA];
	GalaxyA[numgal].Rvir = SubGroupRvir_A[grA];
	GalaxyA[numgal].Vvir = SubGroupVvir_A[grA];
	GalaxyA[numgal].Vc   = SubGroupVc_A[grA];
	
	GalaxyA[numgal].MvirT[Snapshot] = GalaxyA[numgal].Mvir;
	GalaxyA[numgal].RvirT[Snapshot] = GalaxyA[numgal].Rvir;
	
	for (j = 0; j < 3; j++) {
	  /* Comoving coordinates */
	  GalaxyA[numgal].Pos[j] = SubGroupPos_A[grA][j];
	  GalaxyA[numgal].Vel[j] = SubGroupVel_A[grA][j];
	}
	
	/* 
	 * WARNING:
	 * It seems that in very few cases, a new group appears for 
	 * which FirstEntry_A[SubParent_A[grA]]+1 does not match grA.
	 * Then this groups ends up with a single galaxy, which is 
	 * assigned type = 1. This is apparently a bug in the 
	 * postprocessing process. As a patch, if such a galaxy is 
	 * found, it is converted here to type = 0 
	 */
	if (GalaxyA[numgal].Type == 1) {
	  
	  /* Note that here we have to use FirstGalInFOFGroup, because we
	     are dealing with type 1 galaxies, which can only be satellites
	     of type 0s. In that case, the central galaxy is always the
	     central galaxy of a FOF group */
	  central = FirstGalInFOFGroup_A[GalaxyA[numgal].ParentGroup];
	  if (central <= 0) {
	    fprintf(stderr, "Error (generate_new_population): cannot find a"
		    " central galaxy for galaxy %d - Exit\n", numgal);
	    fflush(stderr);
	    exit(EXIT_FAILURE);
	  }
	  if (GalaxyA[central].Type != 0) {
	    printf("WARNING (generate_new_population): central galaxy %d "
		   "found for galaxy %d (type %d) is of type %d != 0\n",
		   central, numgal, GalaxyA[numgal].Type,
		   GalaxyA[central].Type); 
	    printf("ParentGroup    = %d\n"
		   "ParentSubGroup = %d\n"
		   "FirstGalInFOFGroup_A[%d] = %d\n"
		   "FirstGalInSubGroup_A[%d] = %d\n"
		   "NumGalInFOFGroup_A[%d]   = %d\n", 
		   GalaxyA[numgal].ParentGroup, 
		   GalaxyA[numgal].ParentSubGroup,
		   GalaxyA[numgal].ParentGroup, 
		   FirstGalInFOFGroup_A[GalaxyA[numgal].ParentGroup],
		   GalaxyA[numgal].ParentSubGroup,
		   FirstGalInSubGroup_A[GalaxyA[numgal].ParentSubGroup],
		   GalaxyA[numgal].ParentGroup, 
		   NumGalInFOFGroup_A[GalaxyA[numgal].ParentGroup]);
	    Count_orphan_type1++;
	    printf("-------> Setting galaxy %d to type 0 "
		   "(found %d such galaxies so far)\n\n", 
		   central, Count_orphan_type1); 
	    fflush(stdout);
	    GalaxyA[central].Type = 0;
	    GalaxyA[central].TypeT[Snapshot] = 0;
	    /*exit(EXIT_FAILURE);*/
	  }
	  
	  /* For satellite galaxies, determine their positions and 
	     velocities in the reference frame of the central, and their
	     orbital angular momentum */
	  central = FirstGalInFOFGroup_A[GalaxyA[numgal].ParentGroup];
	  
	  /* Relative quantities: store them in physical units. The
	     conversion factor takes into account that position needs to
	     be multiplied by aexp and velocity by sqrt(aexp) */
	  for (j=0; j<3; j++) {
	    GalaxyA[numgal].Posrel[j] = (GalaxyA[numgal].Pos[j] - 
					 GalaxyA[central].Pos[j])*fac;
	    GalaxyA[numgal].Vrel[j] = (GalaxyA[numgal].Vel[j] - 
				       GalaxyA[central].Vel[j])*sqrt(fac);
	  }

	  /* NOTE: you've already converted Posrel and Vrel to physical -
	     no need to do it again for orbital J! */
	  GalaxyA[numgal].Jorb[0] = 
	    GalaxyA[numgal].Posrel[1]*GalaxyA[numgal].Vrel[2] - 
	    GalaxyA[numgal].Posrel[2]*GalaxyA[numgal].Vrel[1];
	  GalaxyA[numgal].Jorb[1] = 
	    GalaxyA[numgal].Posrel[2]*GalaxyA[numgal].Vrel[0] - 
	    GalaxyA[numgal].Posrel[0]*GalaxyA[numgal].Vrel[2];
	  GalaxyA[numgal].Jorb[2] = 
	    GalaxyA[numgal].Posrel[0]*GalaxyA[numgal].Vrel[1] - 
	    GalaxyA[numgal].Posrel[1]*GalaxyA[numgal].Vrel[0];
	  
	  /* Save the initial value of the orbital angular momentum
	     (in physical units, squared) */
	  GalaxyA[numgal].J2init = 
	    GalaxyA[numgal].Jorb[0]*GalaxyA[numgal].Jorb[0] +
	    GalaxyA[numgal].Jorb[1]*GalaxyA[numgal].Jorb[1] +
	    GalaxyA[numgal].Jorb[2]*GalaxyA[numgal].Jorb[2];
	} /* Close if galaxy is type 1 */
	else {
	  /* For a central, set all relative quantities to 0 */
	  for (j=0; j<3; j++) {
	    GalaxyA[numgal].Posrel[j] = 0.0;
	    GalaxyA[numgal].Vrel[j] = 0.0;
	    GalaxyA[numgal].Jorb[j] = 0.0;
	  }
	  GalaxyA[numgal].J2init  = 0.0;
	}
	
	/* Initialise other quantities */
	GalaxyA[numgal].MergerType     = 0;
	GalaxyA[numgal].Relocated      = 0;
	
	GalaxyA[numgal].MergTime       = 0.0;
	GalaxyA[numgal].MergTimeBK     = 0.0;
	GalaxyA[numgal].MergTimeJ      = 0.0;
	GalaxyA[numgal].MergTimeActual = 0.0;
	GalaxyA[numgal].MergTimeSat    = 0.0;
	GalaxyA[numgal].Coulomb        = 0.0;
	GalaxyA[numgal].StrippedMass   = 0.0;
	GalaxyA[numgal].EnergyLoss     = 0.0;
	
	GalaxyA[numgal].Rtidal = GalaxyA[numgal].Rvir;
	GalaxyA[numgal].Rdm    = GalaxyA[numgal].Rvir;
	
	GalaxyA[numgal].Mdm    = SubGroupMvir_A[grA];
	
	/* For type 1 galaxies, add the time spent as type 1 to the
	   actual merging time, and also use the tidal radius as the
	   bounding radius for DM (if it is smaller than Rvir) */
	if (GalaxyA[numgal].Type == 1) {
	  GalaxyA[numgal].MergTimeActual += DeltaT;
	  GalaxyA[numgal].MergTimeSat    += DeltaT;
	  GalaxyA[numgal].Rtidal = get_tidal_radius(numgal,central);
	  if (GalaxyA[numgal].Rtidal < GalaxyA[numgal].Rvir)
	    GalaxyA[numgal].Rdm = GalaxyA[numgal].Rtidal;
	}
	
	for (j=0; j<STEPS; j++) {
	  if (GalaxyA[numgal].Type == 0) {
	    GalaxyA[numgal].OrbitX[j]  = 0.0;
	    GalaxyA[numgal].OrbitY[j]  = 0.0;
	    GalaxyA[numgal].OrbitZ[j]  = 0.0;
	    GalaxyA[numgal].OrbitVx[j] = 0.0;
	    GalaxyA[numgal].OrbitVy[j] = 0.0;
	    GalaxyA[numgal].OrbitVz[j] = 0.0;

	    GalaxyA[numgal].Ekin = 0.0;
	    GalaxyA[numgal].Epot = 0.0;
	  }
	  else {
	    /* NOTE: Orbit is stored in physical units */
	    /* For a type 1, the position read from postprocessing files
	       is copied to every step */
	    GalaxyA[numgal].OrbitX[j]  = GalaxyA[numgal].Posrel[0];
	    GalaxyA[numgal].OrbitY[j]  = GalaxyA[numgal].Posrel[1];
	    GalaxyA[numgal].OrbitZ[j]  = GalaxyA[numgal].Posrel[2];
	    GalaxyA[numgal].OrbitVx[j] = GalaxyA[numgal].Vrel[0];
	    GalaxyA[numgal].OrbitVy[j] = GalaxyA[numgal].Vrel[1];
	    GalaxyA[numgal].OrbitVz[j] = GalaxyA[numgal].Vrel[2];
	  }
	  GalaxyA[numgal].OrbitMass[j] = GalaxyA[numgal].Mvir;
	  GalaxyA[numgal].OrbitType[j] = GalaxyA[numgal].Type;
	} /* Close STEPS */
	
	GalaxyA[numgal].Ekin = 0.5*
	  (GalaxyA[numgal].Vrel[0]*GalaxyA[numgal].Vrel[0] + 
	   GalaxyA[numgal].Vrel[1]*GalaxyA[numgal].Vrel[1] + 
	   GalaxyA[numgal].Vrel[2]*GalaxyA[numgal].Vrel[2]);
	GalaxyA[numgal].Epot = potential_sis(GalaxyA[numgal].Posrel[0],
					     GalaxyA[numgal].Posrel[1],
					     GalaxyA[numgal].Posrel[2],
					     GalaxyA[numgal].Vvir*
					     GalaxyA[numgal].Vvir,
					     GalaxyA[numgal].Rvir,
					     EpsGrav*GalaxyA[numgal].Rvir);
      } /* Close if (CountProgenitor[grA]==0) */
      /* Finish processing new galaxies for this snapshot */
      
      /* List of galaxies is such that for each subgroup, one
	 first has the central galaxy, then all satellites */
      
      if (CountProgenitor[grA] > 0) {  /* Subgroup has progenitors */
	
	for (n=0, grB=FirstProgenitor[grA]; n<CountProgenitor[grA]; 
	     n++, grB=NextProgenitor[grB]) {
	  for (i=0; i<NumGalInSubGroup_B[grB]; i++) {
	    g = FirstGalInSubGroup_B[grB] + i; 
	    oldgalaxies++;
	    numgal++;
	    TotNumGalA++;
	    
	    if ((TotNumGalA+TotNumGalB)>MAXGALAXIES) {
	      fprintf(stderr, "Error (generate_new_population): "
		      "TotNumGalA + TotNumGalB = %d  exceeds MAXGALAXIES ="
		      " %d\n", TotNumGalA+TotNumGalB, MAXGALAXIES);
	      fflush(stderr);
	      exit(EXIT_FAILURE);
	    }
	    
	    /* Match galaxies between the previous and current step */
	    GalaxyA[numgal] = GalaxyB[g];
	    
	    NumGalInSubGroup_A[grA]++;
	    
	    /* n = 0 means we are looking at the first progenitor of the
	       current subhalo, and i = 0 that we are looking at the
	       first galaxy of the first subhalo */
	    if (n==0 && i==0)
	      FirstGalInSubGroup_A[grA]= numgal;
	    
	    GalaxyA[numgal].ParentSubGroup = grA;
	    GalaxyA[numgal].ParentGroup    = SubParent_A[grA];
	    
	    NumGalInFOFGroup_A[SubParent_A[grA]]++;
	    if (FirstGalInFOFGroup_A[SubParent_A[grA]] == 0)
	      FirstGalInFOFGroup_A[SubParent_A[grA]] = numgal;
	    
	    if (FirstGalInFOFGroup_A[SubParent_A[grA]] == numgal) {
	      /* Here we have the central galaxy of the FOF halo */
	      
	      if (GalaxyA[numgal].Type == 2) {
		fprintf(stderr,
			"Error (generate_new_population): a central "
			"galaxy cannot be of Type = 2!\n"
			"  Galaxy: %d\n"
			"  ParentGroup: %d\n"
			"  ParentSubGroup: %d\n"
			"  FirstGalInFOFGroup_A[%d]: %d\n"
			"  FirstGalInSubGroup_A[%d]: %d\n"
			"  NumGalInFOFGroup[%d]: %d\n",
			numgal, GalaxyA[numgal].ParentGroup, 
			GalaxyA[numgal].ParentSubGroup, 
			GalaxyA[numgal].ParentGroup, 
			FirstGalInFOFGroup_A[GalaxyA[numgal].ParentGroup], 
			GalaxyA[numgal].ParentSubGroup, 
			FirstGalInSubGroup_A[GalaxyA[numgal].ParentSubGroup],
			GalaxyA[numgal].ParentGroup, 
			NumGalInFOFGroup_A[GalaxyA[numgal].ParentGroup]);;
		fflush(stderr);
		exit(EXIT_FAILURE);
	      }
	      
	      /* At this point, we may have a galaxy that was type 1 in 
		 the past, now reconverting to type 0 (ejected from its
		 previous host group?) */
	      if (GalaxyA[numgal].Type != 2) {
		if (GalaxyA[numgal].TypeT[Snapshot-1] == 1) type1to0++;
		GalaxyA[numgal].Type = 0;
		GalaxyA[numgal].TypeT[Snapshot] = 0;
	      }
	    } /* Close if (FirstGalInFOFGroup_A = numgal) */
	    else {
	      /* If this is not now the first galaxy in FOF group, and has
		 a subhalo, then it must be of type 1 now (could have been 
		 type 0 at B) */
	      if (GalaxyA[numgal].Type != 2) {
		if (GalaxyA[numgal].TypeT[Snapshot-1] == 0) type0to1++; 

		/* TOMAS 2013-05-08: apparently a type 1 galaxy could have
		   been of type 2 in the past! */
		if (GalaxyA[numgal].TypeT[Snapshot-1] == 2) { 
		  type2to1++;
#ifdef DEBUG 
		  printf("WARNING: galaxy %d converted from type 2 back to "
		       "type 1\n", numgal); fflush(stdout);
#endif		
		}
		GalaxyA[numgal].Type = 1;
		GalaxyA[numgal].TypeT[Snapshot] = 1;
	      }
	    } 
	    
	    /* Update new central galaxy of subhalo (type 0 or 1) */
	    if (n==0 && i==0) {
	      GalaxyA[numgal].PaIndex = SubIdMostBound_A[grA];
	      GalaxyA[numgal].Len     = SubLen_A[grA];
	      
	      GalaxyA[numgal].Mvir = SubGroupMvir_A[grA];
	      GalaxyA[numgal].Rvir = SubGroupRvir_A[grA];
	      GalaxyA[numgal].Vvir = SubGroupVvir_A[grA];
	      GalaxyA[numgal].Vc   = SubGroupVc_A[grA];
	      
	      GalaxyA[numgal].MvirT[Snapshot] = GalaxyA[numgal].Mvir;
	      GalaxyA[numgal].RvirT[Snapshot] = GalaxyA[numgal].Rvir;
	      
	      for (j=0; j<3; j++) {
		/* Comoving units */
		GalaxyA[numgal].Pos[j] = SubGroupPos_A[grA][j];
		GalaxyA[numgal].Vel[j] = SubGroupVel_A[grA][j];
	      }
	      
	      /* Relative position and velocity, 
		 and orbital angular momentum */
	      if (GalaxyA[numgal].Type == 1) {
		/* Type 1: use the central galaxy of FOF group */
		central = FirstGalInFOFGroup_A[GalaxyA[numgal].ParentGroup];
		if (central <= 0) {
		  fprintf(stderr,
			  "Error (generate_new_population): cannot find a "
			  "central galaxy for galaxy %d - Exit\n", numgal);
		  fflush(stderr);
		  exit(EXIT_FAILURE);
		}
		
		for (j=0; j<3; j++) {
		  /* In physical units. For type 1s, we read the position
		     and velocity from the postprocessing files */
		  GalaxyA[numgal].Posrel[j] = (GalaxyA[numgal].Pos[j] - 
					       GalaxyA[central].Pos[j])*fac;
		  GalaxyA[numgal].Vrel[j] = (GalaxyA[numgal].Vel[j] - 
					     GalaxyA[central].Vel[j])*
		    sqrt(fac);
		}
		GalaxyA[numgal].Jorb[0] = 
		  GalaxyA[numgal].Posrel[1]*GalaxyA[numgal].Vrel[2] - 
		  GalaxyA[numgal].Posrel[2]*GalaxyA[numgal].Vrel[1];
		GalaxyA[numgal].Jorb[1] = 
		  GalaxyA[numgal].Posrel[2]*GalaxyA[numgal].Vrel[0] - 
		  GalaxyA[numgal].Posrel[0]*GalaxyA[numgal].Vrel[2];
		GalaxyA[numgal].Jorb[2] = 
		  GalaxyA[numgal].Posrel[0]*GalaxyA[numgal].Vrel[1] - 
		  GalaxyA[numgal].Posrel[1]*GalaxyA[numgal].Vrel[0];
		
		/* Save the initial value of the orbital angular momentum
		   (in physical units, squared) */
		jorb0 = GalaxyA[numgal].Jorb[0]*GalaxyA[numgal].Jorb[0] +
		  GalaxyA[numgal].Jorb[1]*GalaxyA[numgal].Jorb[1] +
		  GalaxyA[numgal].Jorb[2]*GalaxyA[numgal].Jorb[2];
		
		if ((GalaxyA[numgal].TypeT[Snapshot-1] == 0) ||
		    (GalaxyA[numgal].Type == 1 && 
		     jorb0 > GalaxyA[numgal].J2init)) {
		  GalaxyA[numgal].J2init = jorb0;
		}
		/* TOMAS 2013-05-13:
		   For some type 1 galaxies, the angular momentum appears
		   to increase from one snapshot to the next. I don't know
		   if this is because I'm reading it wrong or if this is
		   indeed the case. If we keep the lower value as reference
		   to compute the merger, then J/J0 > 1 for most of the
		   orbit integration and the merger will never happen. To
		   avoid this, I'll keep as J0 the highest value found 
		   while the galaxy is type 1 */
		
	      } /* Close if galaxy is type 1 */
	      else {
		for (j=0; j<3; j++) {
		  GalaxyA[numgal].Posrel[j] = 0.0;
		  GalaxyA[numgal].Vrel[j] = 0.0;
		  GalaxyA[numgal].Jorb[j] = 0.0;
		}
		GalaxyA[numgal].J2init = 0.0;
	      }
	      
	      GalaxyA[numgal].Mdm  = SubGroupMvir_A[grA];
	      
	      GalaxyA[numgal].Rdm  = GalaxyA[numgal].Rvir;
	      
	      /* For type 1 galaxies, update merging time, tidal radius 
		 and bounding radius for DM */
	      if (GalaxyA[numgal].Type == 1) {
		GalaxyA[numgal].MergTimeActual += DeltaT;
		GalaxyA[numgal].MergTimeSat    += DeltaT;
		GalaxyA[numgal].Rtidal = get_tidal_radius(numgal,central);
		if (GalaxyA[numgal].Rtidal < GalaxyA[numgal].Rvir)
		  GalaxyA[numgal].Rdm = GalaxyA[numgal].Rtidal;
	      }
	      
	      for (j=0; j<STEPS; j++) {
		if (GalaxyA[numgal].Type == 0) {
		  GalaxyA[numgal].OrbitX[j]  = 0.0;
		  GalaxyA[numgal].OrbitY[j]  = 0.0;
		  GalaxyA[numgal].OrbitZ[j]  = 0.0;
		  GalaxyA[numgal].OrbitVx[j] = 0.0;
		  GalaxyA[numgal].OrbitVy[j] = 0.0;
		  GalaxyA[numgal].OrbitVz[j] = 0.0;

		  GalaxyA[numgal].Ekin = 0.0;
		  GalaxyA[numgal].Epot = 0.0;
		}
		else {
		  /* NOTE: Orbit stored in physical units */
		  GalaxyA[numgal].OrbitX[j]  = GalaxyA[numgal].Posrel[0];
		  GalaxyA[numgal].OrbitY[j]  = GalaxyA[numgal].Posrel[1];
		  GalaxyA[numgal].OrbitZ[j]  = GalaxyA[numgal].Posrel[2];
		  GalaxyA[numgal].OrbitVx[j] = GalaxyA[numgal].Vrel[0];
		  GalaxyA[numgal].OrbitVy[j] = GalaxyA[numgal].Vrel[1];
		  GalaxyA[numgal].OrbitVz[j] = GalaxyA[numgal].Vrel[2];
		}
		GalaxyA[numgal].OrbitMass[j] = GalaxyA[numgal].Mvir;
		GalaxyA[numgal].OrbitType[j] = GalaxyA[numgal].Type;
	      }

	      GalaxyA[numgal].Ekin = 0.5*
		(GalaxyA[numgal].Vrel[0]*GalaxyA[numgal].Vrel[0] + 
		 GalaxyA[numgal].Vrel[1]*GalaxyA[numgal].Vrel[1] + 
		 GalaxyA[numgal].Vrel[2]*GalaxyA[numgal].Vrel[2]);
	      GalaxyA[numgal].Epot = 
		potential_sis(GalaxyA[numgal].Posrel[0],
			      GalaxyA[numgal].Posrel[1],
			      GalaxyA[numgal].Posrel[2],
			      GalaxyA[numgal].Vvir*
			      GalaxyA[numgal].Vvir,
			      GalaxyA[numgal].Rvir,
			      EpsGrav*GalaxyA[numgal].Rvir);
	    } /* Close if (n==0 and i==0) */
	    
	    /* n > 0 (not the first progenitor of the current subhalo)
	       but i = 0 so we are looking at its first galaxy. This is
	       then an accreted galaxy that was a central galaxy in the 
	       past */
	    if (n>0 && i==0) {
	      
	      if (GalaxyA[numgal].Type == 2) {
		fprintf(stderr,
			"Error (generate_new_population): cannot have "
			"accreted a central Type 2 galaxy (%d) - Exit\n",
			numgal); fflush(stderr);
		exit(EXIT_FAILURE);
	      }
	      
	      /* Accreted central galaxy of a subhalo: the current scheme
		 does not consider subhaloes of subhaloes, so this
		 accreted galaxy has to become a type 2 satellite 
		 (remember we are going over subhaloes here) */
	      if (GalaxyA[numgal].Type == 0) type0to2++;
	      if (GalaxyA[numgal].Type == 1) type1to2++;
	      GalaxyA[numgal].Type = 2;
	      GalaxyA[numgal].TypeT[Snapshot] = 2;
	      
	      /* 
	       * Dynamical friction:
	       * Determine the Coulomb logarithm at infall. Note that
	       * FirstProgenitor points to a subhalo in B 
	       */
	      mother_in_B    = FirstProgenitor[grA];
	      satellite_in_B = grB;
	      GalaxyA[numgal].Coulomb = 
		log(1+SubLen_B[mother_in_B]/
		    ((double)SubLen_B[satellite_in_B]));
	      
	      /* eps is the eccentricity of the satellite's orbit, and 
		 feps is a function that describes the dependence of the
		 orbital decay with eps. */
	      eps  = 0.5;
	      feps = 0.5;
	      
	      /* Determine Hubble(z) */
	      hubble_z = H0*sqrt((1.0+ZZ[Snapshot])*(1.0+ZZ[Snapshot]) *
				 (Omega*(1.0+ZZ[Snapshot]) + 
				  (1-Omega-OmegaLambda)) + OmegaLambda);
	      
	      /* If the galaxy was type 0 in the previous step, then
		 the relative position and velocity will not be available.
		 We calculate them here using the positions and velocities
		 of the progenitor subhaloes (convert to physical) */
	      if (GalaxyA[numgal].TypeT[Snapshot-1] == 0) {
		for (j=0; j<3; j++) {
		  GalaxyA[numgal].Posrel[j] = 
		    (SubGroupPos_B[satellite_in_B][j] - 
		     SubGroupPos_B[mother_in_B][j])*fac_0;
		  GalaxyA[numgal].Vrel[j] = 
		    (SubGroupVel_B[satellite_in_B][j] - 
		     SubGroupVel_B[mother_in_B][j])*sqrt(fac_0);
		}
		GalaxyA[numgal].Jorb[0] = 
		  GalaxyA[numgal].Posrel[1]*GalaxyA[numgal].Vrel[2] - 
		  GalaxyA[numgal].Posrel[2]*GalaxyA[numgal].Vrel[1];
		GalaxyA[numgal].Jorb[1] = 
		  GalaxyA[numgal].Posrel[2]*GalaxyA[numgal].Vrel[0] - 
		  GalaxyA[numgal].Posrel[0]*GalaxyA[numgal].Vrel[2];
		GalaxyA[numgal].Jorb[2] = 
		  GalaxyA[numgal].Posrel[0]*GalaxyA[numgal].Vrel[1] - 
		  GalaxyA[numgal].Posrel[1]*GalaxyA[numgal].Vrel[0];
		
		/* Save the initial value of the orbital angular momentum
		   (in physical units, squared) */
		GalaxyA[numgal].J2init = 
		  GalaxyA[numgal].Jorb[0]*GalaxyA[numgal].Jorb[0] +
		  GalaxyA[numgal].Jorb[1]*GalaxyA[numgal].Jorb[1] +
		  GalaxyA[numgal].Jorb[2]*GalaxyA[numgal].Jorb[2];
		
		/* Update orbit outputs, to use them as initial 
		   conditions */
		for (j=0; j<STEPS; j++) {
		  GalaxyA[numgal].OrbitX[j]  = GalaxyA[numgal].Posrel[0];
		  GalaxyA[numgal].OrbitY[j]  = GalaxyA[numgal].Posrel[1];
		  GalaxyA[numgal].OrbitZ[j]  = GalaxyA[numgal].Posrel[2];
		  GalaxyA[numgal].OrbitVx[j] = GalaxyA[numgal].Vrel[0];
		  GalaxyA[numgal].OrbitVy[j] = GalaxyA[numgal].Vrel[1];
		  GalaxyA[numgal].OrbitVz[j] = GalaxyA[numgal].Vrel[2];
		  
		  GalaxyA[numgal].OrbitMass[j] = GalaxyA[numgal].Mvir;
		  GalaxyA[numgal].OrbitType[j] = GalaxyA[numgal].Type;
		}

		GalaxyA[numgal].Ekin = 0.5*
		  (GalaxyA[numgal].Vrel[0]*GalaxyA[numgal].Vrel[0] + 
		   GalaxyA[numgal].Vrel[1]*GalaxyA[numgal].Vrel[1] + 
		   GalaxyA[numgal].Vrel[2]*GalaxyA[numgal].Vrel[2]);
		GalaxyA[numgal].Epot = 
		  potential_sis(GalaxyA[numgal].Posrel[0],
				GalaxyA[numgal].Posrel[1],
				GalaxyA[numgal].Posrel[2],
				GalaxyA[numgal].Vvir*
				GalaxyA[numgal].Vvir,
				GalaxyA[numgal].Rvir,
				EpsGrav*GalaxyA[numgal].Rvir);
	      }
	      /*printf("numgal, Mvir: %d %g\n", numgal, GalaxyA[numgal].Mvir);*/
	      
	      /* Jiang et al. (2008) dynamical friction timescale,
		 Equation 5 */
	      GalaxyA[numgal].MergTimeJ = 0.5/0.43*(0.94*pow(eps,0.60)+0.6)/GalaxyA[numgal].Coulomb*(SubGroupMvir_B[mother_in_B]/SubGroupMvir_B[satellite_in_B])*SubGroupRvir_B[mother_in_B]/SubGroupVvir_B[mother_in_B];
	      
	      /* Boylan-Kolchin et al. (2008) dynamical friction 
		 timescale */
	      hubble_z *= Hubble_h;
	      GalaxyA[numgal].MergTimeBK = 0.1/hubble_z*0.216*pow(SubGroupMvir_B[mother_in_B]/SubGroupMvir_B[satellite_in_B],1.3)/GalaxyA[numgal].Coulomb*exp(1.9*sqrt(1-pow(eps,2)));
	      hubble_z /= Hubble_h;
	      
	      /* Chandrasekhar dynamical friction timescale */
	      GalaxyA[numgal].MergTime = 0.5/(G*0.43)*feps/GalaxyA[numgal].Coulomb*SubGroupVvir_B[mother_in_B]*pow(SubGroupRvir_B[mother_in_B],2)/SubGroupMvir_B[satellite_in_B];
	      
	    } /* Close n>0 and i=0 */ 
	    
	    /*
	     * Satellite galaxies of accreted smaller progenitor subhalos;
	     * they are already defined as type 2 galaxies, but the
	     * MergerTime must be reset and evaluated with respect to
	     * the central galaxy of the larger halo they fell in 
	     * (first progenitor of grA)
	     */  
	    if (n>0 && i>0) {
	      if (GalaxyA[numgal].Type != 2) {
		fprintf(stderr,
			"Error (generate_new_population): satellite "
			"galaxies of accreted progenitor subhaloes must "
			"be type 2\n"); fflush(stderr);
		exit(EXIT_FAILURE);
	      }
	      
	      /* Accreted satellite galaxy of a subhalo: reset MergTime
		 with respect to central galaxy of the larger subgroup */
	      mother_in_B    = FirstProgenitor[grA];
	      satellite_in_B = grB;
	      
	      /* In our subhalo scheme, these satellites are now
		 satellites of the central galaxy of the subhalo.
		 Therefore, we need to reset their relative positions,
		 velocities, and the merger time. In MergTimeSat,
		 however, we keep tracking the total time spent as a
		 satellite galaxy, regardless of central */
	      
	      /* TOMAS and IGNACIO: NEW 2012-08-29
		 Type 2 satellites of type 1 galaxies can be stripped
		 away and assigned to the corresponding FOF central
		 (Relocated = 1). In such a case, the position and 
		 velocity already are relative to the central so all one
		 has to do here is reset the Relocated flag */

	      /* TOMAS 2013-05-17:
		 Galaxy relocation does not work with the current subhalo
		 scheme. The relocation flag is disabled (it will always be
		 zero) but we do not delete the code, in case we find a 
		 satisfactory solution in the future. */
	      
	      if (SatelliteRelocationDisabledOn == 0) {
		if (GalaxyA[numgal].Relocated == 1)
		  GalaxyA[numgal].Relocated = 0;
		else {
		  for (j=0; j<3; j++) {
		    GalaxyA[numgal].Posrel[j] += 
		      (SubGroupPos_B[satellite_in_B][j] - 
		       SubGroupPos_B[mother_in_B][j])*fac_0;
		    GalaxyA[numgal].Vrel[j] += 
		      (SubGroupVel_B[satellite_in_B][j] - 
		       SubGroupVel_B[mother_in_B][j])*sqrt(fac_0);
		  }
		  GalaxyA[numgal].Jorb[0] = 
		    GalaxyA[numgal].Posrel[1]*GalaxyA[numgal].Vrel[2] - 
		    GalaxyA[numgal].Posrel[2]*GalaxyA[numgal].Vrel[1];
		  GalaxyA[numgal].Jorb[1] = 
		    GalaxyA[numgal].Posrel[2]*GalaxyA[numgal].Vrel[0] - 
		    GalaxyA[numgal].Posrel[0]*GalaxyA[numgal].Vrel[2];
		  GalaxyA[numgal].Jorb[2] = 
		    GalaxyA[numgal].Posrel[0]*GalaxyA[numgal].Vrel[1] - 
		    GalaxyA[numgal].Posrel[1]*GalaxyA[numgal].Vrel[0];
		  
		  /* Update the initial J */
		  GalaxyA[numgal].J2init = 
		    GalaxyA[numgal].Jorb[0]*GalaxyA[numgal].Jorb[0] +
		    GalaxyA[numgal].Jorb[1]*GalaxyA[numgal].Jorb[1] +
		    GalaxyA[numgal].Jorb[2]*GalaxyA[numgal].Jorb[2];
		  
		  /* Update orbit outputs, to use them as initial 
		     conditions */
		  for (j=0; j<STEPS; j++) {
		    GalaxyA[numgal].OrbitX[j]  = GalaxyA[numgal].Posrel[0];
		    GalaxyA[numgal].OrbitY[j]  = GalaxyA[numgal].Posrel[1];
		    GalaxyA[numgal].OrbitZ[j]  = GalaxyA[numgal].Posrel[2];
		    GalaxyA[numgal].OrbitVx[j] = GalaxyA[numgal].Vrel[0];
		    GalaxyA[numgal].OrbitVy[j] = GalaxyA[numgal].Vrel[1];
		    GalaxyA[numgal].OrbitVz[j] = GalaxyA[numgal].Vrel[2];
		    
		    GalaxyA[numgal].OrbitType[j] = GalaxyA[numgal].Type;
		    /* The galaxy already was a type 2, so it already 
		       experienced tidal stripping. Thus use its last 
		       recorded mass */
		    GalaxyA[numgal].OrbitMass[j] = 
		      GalaxyA[numgal].OrbitMass[STEPS-1];
		  }

		  GalaxyA[numgal].Ekin = 0.5*
		    (GalaxyA[numgal].Vrel[0]*GalaxyA[numgal].Vrel[0] + 
		     GalaxyA[numgal].Vrel[1]*GalaxyA[numgal].Vrel[1] + 
		     GalaxyA[numgal].Vrel[2]*GalaxyA[numgal].Vrel[2]);
		  GalaxyA[numgal].Epot = 
		    potential_sis(GalaxyA[numgal].Posrel[0],
				  GalaxyA[numgal].Posrel[1],
				  GalaxyA[numgal].Posrel[2],
				  GalaxyA[numgal].Vvir*
				  GalaxyA[numgal].Vvir,
				  GalaxyA[numgal].Rvir,
				  EpsGrav*GalaxyA[numgal].Rvir);
		  
		  /* Reset the merging time */
		  GalaxyA[numgal].MergTimeActual = 0.0;
		  
		  /* Also recalculate Coulomb logarithm */
		  GalaxyA[numgal].Coulomb = 
		    log(1.0 + SubLen_B[mother_in_B]/
			((double)SubLen_B[satellite_in_B]));
		  
		  eps  = 0.5;
		  feps = 0.5;
		  
		  /* Determine Hubble(z) */
		  hubble_z = H0*sqrt((1.0+ZZ[Snapshot])*(1.0+ZZ[Snapshot]) *
				     (Omega*(1.0+ZZ[Snapshot]) + 
				      (1-Omega-OmegaLambda)) + OmegaLambda);
		  
		  /* Jiang et al. (2009) DF timescale */
		  GalaxyA[numgal].MergTimeJ = 0.5/0.43*(0.94*pow(eps,0.60)+0.6)/GalaxyA[numgal].Coulomb*(SubGroupMvir_B[mother_in_B]/SubGroupMvir_B[satellite_in_B])*SubGroupRvir_B[mother_in_B]/ SubGroupVvir_B[mother_in_B];
		  
		  /* Boylan-Kolchin et al. (2008) DF timescale */
		  hubble_z *= Hubble_h;
		  GalaxyA[numgal].MergTimeBK = 0.1/hubble_z*0.216*pow(SubGroupMvir_B[mother_in_B]/SubGroupMvir_B[satellite_in_B],1.3)/GalaxyA[numgal].Coulomb*exp(1.9*sqrt(1-pow(eps,2)));
		  hubble_z /= Hubble_h;
		  
		  /* Chandrasekhar DF timescale */
		  GalaxyA[numgal].MergTime = 0.5/(G*0.43)*feps/GalaxyA[numgal].Coulomb*SubGroupVvir_B[mother_in_B]*pow(SubGroupRvir_B[mother_in_B],2)/SubGroupMvir_B[satellite_in_B];
		} /* Close if Relocated = 0 */
	      }   /* Close if relocation is enabled */

	    } /* Close n>0 and i>0 */
	  } /* Close for (i=0; i<NumGalInSubGroup_B[grB]... */
	}
      } /* Close if (CountProgenitor[grA] > 0) */
    } /* Close if (Volatile_A[grA] == 0) */
  } /* Close for (grA=1, numgal=0... */  
  
  nmergers = ndisrupt = 0;
  for (i=1; i<=TotNumGalA; i++) {
    
    GalaxyA[i].IndexT[Snapshot] = i;

    /* Process all satellite galaxies */
    if (GalaxyA[i].Type > 0) { 
#ifdef DEBUG    
      printf("======================================================\n");
      printf("Galaxy ID %d, type %d (previous %d) - index %d "
	     "(previous %d)\n", GalaxyA[i].Id, GalaxyA[i].Type, 
	     GalaxyA[i].TypeT[Snapshot-1], i, 
	     GalaxyA[i].IndexT[Snapshot-1]);
      if (GalaxyA[i].Type == 1)  testprint_type1(i,fac);
#endif
      
      if (GalaxyA[i].Type == 2)
	central = FirstGalInSubGroup_A[GalaxyA[i].ParentSubGroup];
      else
	central = FirstGalInFOFGroup_A[GalaxyA[i].ParentGroup];
      
      /* Find when it became a satellite */
      satflag = 0;
      for (j=0; j<OUTPUTS; j++) {
	if (satflag == 0) {
	  if (GalaxyA[i].TypeT[j] == 1) {
	    satflag = 1;
#ifdef DEBUG
	    printf("    Type 1 in snapshot %d\n", j);
#endif
	  }
	  if (GalaxyA[i].TypeT[j] == 2) {
	    satflag = 2;
#ifdef DEBUG
	    printf("    Type 2 in snapshot %d\n", j);
	    fflush(stdout);
#endif
	  }
	}
      }
      
      /* 
       * TOMAS AND IGNACIO 2012-08-29: NEW
       * Now search for all satellites of type 1 galaxies which are 
       * beyond the type 1's tidal radius, and reassign those galaxies
       * to the FOF central.
       * WARNING: one cannot directly change the values of ParentSubGroup
       * and NumGalInSubGroup. Galaxies are numbered in sequence, and
       * relabeling them without reordering may wreak havoc with other
       * routines. Best to just label them for separate treatment: if
       * Relocated = 1, make it orbit around the FOF central
       */

      /* TOMAS 2013-07-17
	 Satellite relocation does not work with the current subhalo
	 scheme. Disabled until we can find a better solution for stripping
	 of satellite galaxies. */
      if (SatelliteRelocationDisabledOn == 0) {
	if (GalaxyA[central].Type == 1) {
	  rsat = GalaxyA[i].Posrel[0]*GalaxyA[i].Posrel[0] + 
	    GalaxyA[i].Posrel[1]*GalaxyA[i].Posrel[1] + 
	    GalaxyA[i].Posrel[2]*GalaxyA[i].Posrel[2];
	
	  if (GalaxyA[i].Relocated == 0 &&
	      rsat > GalaxyA[central].Rvir*GalaxyA[central].Rvir) {
#ifdef DEBUG
	    printf("  ** Galaxy %d (type %d), central %d (type %d)\n"
		   "     Distance to central: %g\n"
		   "     Rvir of central:     %g\n"
		   "     Galaxy relocated to FOF central %d\n",
		   i, GalaxyA[i].Type, central, GalaxyA[central].Type,
		   sqrt(rsat), GalaxyA[central].Rvir, 
		   FirstGalInFOFGroup_A[GalaxyA[central].ParentGroup]);
	    fflush(stdout);
#endif
	    GalaxyA[i].Relocated = 1;
	    GalaxyA[i].Posrel[0] += GalaxyA[central].Posrel[0];
	    GalaxyA[i].Posrel[1] += GalaxyA[central].Posrel[1];
	    GalaxyA[i].Posrel[2] += GalaxyA[central].Posrel[2];
	    
	    GalaxyA[i].Vrel[0] += GalaxyA[central].Vrel[0];
	    GalaxyA[i].Vrel[1] += GalaxyA[central].Vrel[1];
	    GalaxyA[i].Vrel[2] += GalaxyA[central].Vrel[2];
	    
	    /* Reset merging time, but do not recalculate the timescales 
	       (to compare with the original prescription, which does not
	       consider stripping of satellites) */
	    GalaxyA[numgal].MergTimeActual = 0.0;
	    
	    /* Recalculate the orbital angular momentum */
#ifdef DEBUG
	    printf("J = %g, J_0 = %g (previous)\n", 
		   sqrt(GalaxyA[i].Jorb[0]*GalaxyA[i].Jorb[0] +
			GalaxyA[i].Jorb[1]*GalaxyA[i].Jorb[1] +
			GalaxyA[i].Jorb[2]*GalaxyA[i].Jorb[2]),
		   sqrt(GalaxyA[i].J2init));
#endif
	    GalaxyA[i].Jorb[0] = 
	      GalaxyA[i].Posrel[1]*GalaxyA[i].Vrel[2] - 
	      GalaxyA[i].Posrel[2]*GalaxyA[i].Vrel[1];
	    GalaxyA[i].Jorb[1] = 
	      GalaxyA[i].Posrel[2]*GalaxyA[i].Vrel[0] - 
	    GalaxyA[i].Posrel[0]*GalaxyA[i].Vrel[2];
	    GalaxyA[i].Jorb[2] = 
	      GalaxyA[i].Posrel[0]*GalaxyA[i].Vrel[1] - 
	      GalaxyA[i].Posrel[1]*GalaxyA[i].Vrel[0];
	    
	    /* Update the initial J */
	    GalaxyA[i].J2init = 
	      GalaxyA[i].Jorb[0]*GalaxyA[i].Jorb[0] +
	      GalaxyA[i].Jorb[1]*GalaxyA[i].Jorb[1] +
	      GalaxyA[i].Jorb[2]*GalaxyA[i].Jorb[2];
#ifdef DEBUG
	    printf("J = %g, J_0 = %g (current)\n", 
		   sqrt(GalaxyA[i].Jorb[0]*GalaxyA[i].Jorb[0] +
			GalaxyA[i].Jorb[1]*GalaxyA[i].Jorb[1] +
			GalaxyA[i].Jorb[2]*GalaxyA[i].Jorb[2]),
		   sqrt(GalaxyA[i].J2init));
#endif
	  } /* Close if Relocated = 0 and galaxy beyond tidal radius */
	}   /* Close if galaxy is type 1 */
      }     /* Close if relocation enabled */
      
      /* Integrate the orbit of type 2 satellites. If a merger happens 
	 before the time between snapshots has elapsed, flag the galaxy as 
	 merged. The orbit integration function updates the subhalo mass 
	 and the position and velocity of the galaxy, but Mvir and Rvir 
	 are kept as the last recorded values for type 1 */
      if (GalaxyA[i].Type == 2) {
#ifdef DEBUG
	testprint_0(i, central, fac);
#endif
	norb = 0;
        epsgrv = (double)EpsGrav;
        tstep  = (double)INTEGRATIONDT;
	do {      
	  status = integrate_type2_orbit(i, central, (double)DeltaT, 
                                         tstep, epsgrv);
	  if (status == -2) {
	    tstep *= 10.0; /* Note that tstep is a number of subdivisions
			      of the time interval */
	  }
	  if (status == -1)
	    epsgrv /= 10.0;
	  norb++;
	  if (norb > NORBTIMES) {
	    fprintf(stderr, "Error (generate_new_population): too many "
		    "calls to orbit integration function for galaxy %d\n"
                    "(Perhaps you have set the value of EpsGrav too "
                    " high...)\nExiting\n", i); fflush(stderr);
	    exit(EXIT_FAILURE);
	  }
        } while (status < 0);

#ifdef DEBUG	
	testprint_1(i);
#endif
	      
	/* Update the position of the galaxy in the simulation reference
	   frame (comoving coordinates) */
	for (j=0; j<3; j++) {
	  GalaxyA[i].Pos[j] = GalaxyA[central].Pos[j] + 
	    GalaxyA[i].Posrel[j]/fac;
	  GalaxyA[i].Vel[j] = GalaxyA[central].Vel[j] + 
	    GalaxyA[i].Vrel[j]/sqrt(fac);
	}

	GalaxyA[i].TypeT[Snapshot] = GalaxyA[i].Type;
      }

      if (EmergencyStop > 0) {
	printf("--- Program stopped at snapshot=%d, DeltaT=%g Myr\n", 
	       Snapshot, DeltaT*UnitTime_in_s/SEC_PER_MEGAYEAR); 
	fflush(stdout);
	exit(EXIT_FAILURE);
      }

#ifdef DEBUG
      /* Merged galaxies */
      if (GalaxyA[i].Type == 3) {
	testprint_0(i, central, fac);

	printf("    Galaxy merged with its central in this snapshot\n");
	if (GalaxyA[i].MergerType == 1)
	  printf("    (merger by angular momentum loss)\n");
	if (GalaxyA[i].MergerType == 2)
	  printf("    (merger by proximity)\n");
	printf("    Mass lost to TS: %g h^-1 M_Sun (%g per cent)\n", 
	       GalaxyA[i].StrippedMass*UnitMass_in_Msun, 
	       GalaxyA[i].StrippedMass*100.0/GalaxyA[i].Mvir);
	printf("    Estimated merging time: %g Myr (C), %g Myr (BK), "
	       "%g Myr (J)\n",
	       GalaxyA[i].MergTime*UnitTime_in_s/SEC_PER_MEGAYEAR/Hubble_h,
	       GalaxyA[i].MergTimeBK*UnitTime_in_s/SEC_PER_MEGAYEAR/
	       Hubble_h,
	       GalaxyA[i].MergTimeJ*UnitTime_in_s/SEC_PER_MEGAYEAR/
	       Hubble_h);
	printf("    Time elapsed in merger clock: %g Myr\n\n", 
	       GalaxyA[i].MergTimeActual*UnitTime_in_s/SEC_PER_MEGAYEAR/
	       Hubble_h);
      }
      fflush(stdout);
#endif
    } /* Close if galaxy type > 0 */
    
    /* Count types of galaxies */
    if (GalaxyA[i].Type == 0) TotNumType0A++;
    if (GalaxyA[i].Type == 1) TotNumType1A++;
    if (GalaxyA[i].Type == 2) TotNumType2A++;
    if (GalaxyA[i].Type == 3) TotNumType3A++;
    if (GalaxyA[i].Type == 4) TotNumType4A++;
  } /* Finish processing galaxies */

  /*printf("======================================================\n");*/
    
  printf("\nFOUND: TotNumGalA = %d, Nsubgroups_A = %d, Vol = %d\n\n",
	 TotNumGalA, Nsubgroups_A, vol); 
  printf("TotNumGalB = %d, oldgalaxies = %d, newgalaxies = %d, "
	 "progcount = %d\n", TotNumGalB, oldgalaxies, newgalaxies, progcount); 
  printf("At time step A: type 0 = %d, type 1 = %d, type 2 = %d\n", 
	 TotNumType0A, TotNumType1A, TotNumType2A);
  printf("Mergers in this snapshot: %d\n", TotNumType3A);
  printf("Disrupted galaxies in this snapshot: %d\n", TotNumType4A);
  printf("======================================================\n");
  fflush(stdout);
  
  if (TotNumGalA != TotNumType0A+TotNumType1A+TotNumType2A+TotNumType3A+
      TotNumType4A) {
    fprintf(stderr,
	    "Error (generate_new_population): numbers of galaxies don't "
	    "match - Exit\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  return;
}


/**
 * @brief Print extra information to the standard output.
 * @param numgal Index of selected galaxy
 * @param central Index of the corresponding central galaxy
 * @param fac Factor @f$1 / (1+z)@f$ to convert from comoving to physical
 * units
 */
void testprint_0(int numgal, int central, double fac)
{
  int fofcentral;

  printf("************** Galaxy %d, type %d (previous %d) *************\n",
	 numgal, GalaxyA[numgal].Type, GalaxyA[numgal].TypeT[Snapshot-1]);
  printf("Central galaxy %d (prev %d), type %d (previously %d)\n", 
	 central, GalaxyA[central].IndexT[Snapshot-1], 
	 GalaxyA[central].Type, GalaxyA[central].TypeT[Snapshot-1]);
  printf("    Mass      %g\n", GalaxyA[central].Mvir);
  printf("    Rvir      %g\n", GalaxyA[central].Rvir);
  printf("    Vvir      %g\n", GalaxyA[central].Vvir);
  printf("    Position and velocity of central (current, "
	 "physical coordinates):\n"
	 "        X     %g\n"
	 "        Y     %g\n"
	 "        Z     %g\n"
	 "        VX    %g\n"
	 "        VY    %g\n"
	 "        VZ    %g\n",
	 GalaxyA[central].Pos[0]*fac, 
	 GalaxyA[central].Pos[1]*fac,
	 GalaxyA[central].Pos[2]*fac, 
	 GalaxyA[central].Vel[0]*sqrt(fac), 
	 GalaxyA[central].Vel[1]*sqrt(fac),
	 GalaxyA[central].Vel[2]*sqrt(fac)); 
  
  if (GalaxyA[numgal].Relocated > 0)
    printf("     Galaxy relocated to the FOF central\n");

  if (GalaxyA[central].Type == 1) {
    fofcentral = FirstGalInFOFGroup_A[GalaxyA[central].ParentGroup];
    printf("FOF central galaxy %d (type %d)\n", 
	   fofcentral, GalaxyA[fofcentral].Type);
    printf("    Mass      %g\n", GalaxyA[fofcentral].Mvir);
    printf("    Rvir      %g\n", GalaxyA[fofcentral].Rvir);
    printf("    Vvir      %g\n", GalaxyA[fofcentral].Vvir);
    printf("    Position and velocity of FOF central (current, "
	   "physical coordinates):\n"
	   "        X     %g\n"
	   "        Y     %g\n"
	   "        Z     %g\n"
	   "        VX    %g\n"
	   "        VY    %g\n"
	   "        VZ    %g\n",
	   GalaxyA[fofcentral].Pos[0]*fac, 
	   GalaxyA[fofcentral].Pos[1]*fac,
	   GalaxyA[fofcentral].Pos[2]*fac, 
	   GalaxyA[fofcentral].Vel[0]*sqrt(fac), 
	   GalaxyA[fofcentral].Vel[1]*sqrt(fac),
	   GalaxyA[fofcentral].Vel[2]*sqrt(fac)); 
  }

  printf("    Mvir      %g\n", GalaxyA[numgal].Mvir);
  printf("    Rvir      %g\n", GalaxyA[numgal].Rvir);	     
  printf("    Mdm       %g\n", 
	 GalaxyA[numgal].OrbitMass[STEPS-1]);
  printf("    Rdm       %g\n", GalaxyA[numgal].Rdm);
  printf("    Relative position and velocity of galaxy "
	 "(initial, physical coordinates):\n"
	 "        x     %g\n"
	 "        y     %g\n"
	 "        z     %g\n"
	 "        vx    %g\n"
	 "        vy    %g\n"
	 "        vz    %g\n",
	 GalaxyA[numgal].Posrel[0], GalaxyA[numgal].Posrel[1], 
	 GalaxyA[numgal].Posrel[2], GalaxyA[numgal].Vrel[0], 
	 GalaxyA[numgal].Vrel[1], GalaxyA[numgal].Vrel[2]);
  printf("    Distance to centre: %g\n", 
	 sqrt(GalaxyA[numgal].Posrel[0]*GalaxyA[numgal].Posrel[0] + 
	      GalaxyA[numgal].Posrel[1]*GalaxyA[numgal].Posrel[1] + 
	      GalaxyA[numgal].Posrel[2]*GalaxyA[numgal].Posrel[2]));
  printf("    Last pericentre:    %g\n", GalaxyA[numgal].Rpericentre);
  printf("    Velocity:           %g\n", 
	 sqrt(GalaxyA[numgal].Vrel[0]*GalaxyA[numgal].Vrel[0] + 
	      GalaxyA[numgal].Vrel[1]*GalaxyA[numgal].Vrel[1] + 
	      GalaxyA[numgal].Vrel[2]*GalaxyA[numgal].Vrel[2]));
  printf("    Kinetic energy:     %g\n", GalaxyA[numgal].Ekin);
  printf("    Potential energy:   %g\n", GalaxyA[numgal].Epot);
  printf("    Total energy:       %g\n", GalaxyA[numgal].Ekin + 
	 GalaxyA[numgal].Epot);
  fflush(stdout);

  return;
}


/**
 * @brief Print additional information to the standard output.
 * @param numgal Index of selected galaxy
 */
void testprint_1(int numgal)
{
  printf("    Relative position and velocity of galaxy "
	 "(after orbit integration, physical coordinates):\n"
	 "        x     %g\n"
	 "        y     %g\n"
	 "        z     %g\n"
	 "        vx    %g\n"
	 "        vy    %g\n"
	 "        vz    %g\n",
	 GalaxyA[numgal].Posrel[0], GalaxyA[numgal].Posrel[1], 
	 GalaxyA[numgal].Posrel[2], GalaxyA[numgal].Vrel[0], 
	 GalaxyA[numgal].Vrel[1], GalaxyA[numgal].Vrel[2]);
  printf("    Mdm       %g\n", 
	 GalaxyA[numgal].OrbitMass[STEPS-1]);
  printf("    Rdm       %g\n", GalaxyA[numgal].Rdm);
  printf("    Distance to centre: %g\n", 
	 sqrt(GalaxyA[numgal].Posrel[0]*GalaxyA[numgal].Posrel[0] + 
	      GalaxyA[numgal].Posrel[1]*GalaxyA[numgal].Posrel[1] + 
	      GalaxyA[numgal].Posrel[2]*GalaxyA[numgal].Posrel[2]));
  printf("    Velocity:           %g\n", 
	 sqrt(GalaxyA[numgal].Vrel[0]*GalaxyA[numgal].Vrel[0] + 
	      GalaxyA[numgal].Vrel[1]*GalaxyA[numgal].Vrel[1] + 
	      GalaxyA[numgal].Vrel[2]*GalaxyA[numgal].Vrel[2]));
  printf("    Mass lost to TS:    %g h^-1 M_Sun\n", 
	 GalaxyA[numgal].StrippedMass*UnitMass_in_Msun);
  printf("    Estimated merging time: %g Myr (C), %g Myr (BK), "
	 "%g Myr (J)\n",
	 GalaxyA[numgal].MergTime*UnitTime_in_s/SEC_PER_MEGAYEAR/Hubble_h,
	 GalaxyA[numgal].MergTimeBK*UnitTime_in_s/SEC_PER_MEGAYEAR/Hubble_h,
	 GalaxyA[numgal].MergTimeJ*UnitTime_in_s/SEC_PER_MEGAYEAR/Hubble_h);
  printf("    Time elapsed in merger clock: %g Myr\n", 
	 GalaxyA[numgal].MergTimeActual*UnitTime_in_s/SEC_PER_MEGAYEAR/
	 Hubble_h);
  printf("    DeltaT: %g Myr\n\n", DeltaT*UnitTime_in_s/SEC_PER_MEGAYEAR/
	 Hubble_h);
  fflush(stdout);
  
  return;
}


/**
 * @brief Print some more information for type 1 galaxies.
 * @param p Index of selected galaxy
 * @param fac Factor @f$1 / (1+z)@f$ to convert from comoving to physical
 * units
 */
void testprint_type1(int p, double fac)
{
  int central;

  central = FirstGalInFOFGroup_A[GalaxyA[p].ParentGroup];

  printf("************** Galaxy %d, type %d (previous %d) *************\n",
	 p, GalaxyA[p].Type, GalaxyA[p].TypeT[Snapshot-1]);
  printf("Central galaxy %d (prev %d), type %d (previously %d)\n", 
	 central, GalaxyA[central].IndexT[Snapshot-1], 
	 GalaxyA[central].Type, GalaxyA[central].TypeT[Snapshot-1]);
  printf("    Mass      %g\n", GalaxyA[central].Mvir);
  printf("    Rvir      %g\n", GalaxyA[central].Rvir);
  printf("    Vvir      %g\n", GalaxyA[central].Vvir);
  printf("    Position and velocity of central (current, "
	 "physical coordinates):\n"
	 "        X     %g\n"
	 "        Y     %g\n"
	 "        Z     %g\n"
	 "        VX    %g\n"
	 "        VY    %g\n"
	 "        VZ    %g\n",
	 GalaxyA[central].Pos[0]*fac, 
	 GalaxyA[central].Pos[1]*fac,
	 GalaxyA[central].Pos[2]*fac, 
	 GalaxyA[central].Vel[0]*sqrt(fac), 
	 GalaxyA[central].Vel[1]*sqrt(fac),
	 GalaxyA[central].Vel[2]*sqrt(fac)); 
  
  printf("    Mass      %g\n", GalaxyA[p].Mvir);
  printf("    Rvir      %g\n", GalaxyA[p].Rvir);	     
  printf("    Mdm       %g\n", GalaxyA[p].OrbitMass[STEPS-1]);
  printf("    Rdm       %g\n", GalaxyA[p].Rdm);
  printf("    Relative position and velocity of galaxy "
	 "(physical coordinates):\n"
	 "        x     %g\n"
	 "        y     %g\n"
	 "        z     %g\n"
	 "        vx    %g\n"
	 "        vy    %g\n"
	 "        vz    %g\n",
	 GalaxyA[p].Posrel[0], GalaxyA[p].Posrel[1], 
	 GalaxyA[p].Posrel[2], GalaxyA[p].Vrel[0], 
	 GalaxyA[p].Vrel[1], GalaxyA[p].Vrel[2]);
  printf("    Distance to centre: %g\n", 
	 sqrt(GalaxyA[p].Posrel[0]*GalaxyA[p].Posrel[0] + 
	      GalaxyA[p].Posrel[1]*GalaxyA[p].Posrel[1] + 
	      GalaxyA[p].Posrel[2]*GalaxyA[p].Posrel[2]));
  printf("    Velocity: %g\n", 
	 sqrt(GalaxyA[p].Vrel[0]*GalaxyA[p].Vrel[0] + 
	      GalaxyA[p].Vrel[1]*GalaxyA[p].Vrel[1] + 
	      GalaxyA[p].Vrel[2]*GalaxyA[p].Vrel[2]));
  printf("    Angular momentum (initial) = %g\n", sqrt(GalaxyA[p].J2init));
  printf("    Angular momentum (current) = %g\n", 
	 sqrt(GalaxyA[p].Jorb[0]*GalaxyA[p].Jorb[0] + 
	      GalaxyA[p].Jorb[1]*GalaxyA[p].Jorb[1] + 
	      GalaxyA[p].Jorb[2]*GalaxyA[p].Jorb[2]));
  fflush(stdout);

  return;
}
