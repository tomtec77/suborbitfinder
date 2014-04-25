/**
 * @file io.c
 * @brief Functions to read halo/subhalo data from postprocessing files.
 */
#include <sys/types.h>
#include <sys/stat.h>

#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


/**
 * @brief Read subhaloes at the previous simulation snapshot.
 * @param catalogue_fname Name of postprocessing sublist file
 * @param substruc_fname Name of postprocessing substructures file
 * @param volatile_fname Name of postprocessing volatile file
 * @param subids_fname Name of postprocessing subids file 
 */
void read_subgroups_B(char *catalogue_fname, char *substruc_fname, 
		      char *volatile_fname, char *subids_fname)
{
  FILE *fd;
  int i;
  int nsubids_B, *idlist;
  int nsubgroubs, check;

  if (!(fd = fopen(catalogue_fname, "r"))) {
    fprintf(stderr, "Error (read_subgroups_B): can't open file '%s'\n", 
	    catalogue_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fread(&Nsubgroups_B, sizeof(int), 1, fd);

  SubLen_B = ivector(1, Nsubgroups_B);
  SubFirst_B = ivector(1, Nsubgroups_B);
  SubParent_B = ivector(1, Nsubgroups_B);
  SubIdMostBound_B = ivector(1, Nsubgroups_B);

  for(i = 1; i <= Nsubgroups_B; i++) {
    check = fread(&SubLen_B[i], sizeof(int), 1, fd);
    if (check != 1) { 
      fprintf(stderr, 
	      "Error (read_subgroups_B): Couldn't read subgroups. Stop\n"); 
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
    check = fread(&SubFirst_B[i], sizeof(int), 1, fd);
    if (check != 1) {
      fprintf(stderr, 
	      "Error (read_subgroups_B): Couldn't read subgroups. Stop\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }

    check = fread(&SubParent_B[i], sizeof(int), 1, fd);
    if (check != 1) {
      fprintf(stderr, 
	      "Error (read_subgroups_B): Couldn't read subgroups. Stop\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
  }
  fclose(fd);

  printf("Subgroup catalogue `%s' read.\n", catalogue_fname);
  printf("%d subgroups.\n\n", Nsubgroups_B);
  
  if(!(fd = fopen(substruc_fname, "r"))) {
    fprintf(stderr, 
	    "Error (read_subgroups_B): can't open file `%s`\n", 
	    substruc_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fread(&Ngroups_B, sizeof(int), 1, fd);

  GroupsSubCount_B = ivector(1, Ngroups_B);
  FirstEntry_B = ivector(1, Ngroups_B);
  
  check = fread(&GroupsSubCount_B[1], sizeof(int), Ngroups_B, fd);
  if (check != Ngroups_B) {
    fprintf(stderr, 
	    "Error (read_subgroups_B): Couldn't read substructure. Stop\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  check = fread(&FirstEntry_B[1], sizeof(int), Ngroups_B, fd);
  if (check != Ngroups_B) {
    fprintf(stderr, 
	    "Error (read_subgroups_B): Couldn't read substructure. Stop\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  fclose(fd);

  printf("Group catalogue `%s' read.\n", substruc_fname);
  printf("%d groups.\n\n", Ngroups_B);

  Volatile_B = ivector(1, Nsubgroups_B);
  if(*volatile_fname) {
    if(!(fd = fopen(volatile_fname, "r"))) {
      fprintf(stderr, "Error (read_subgroups_B): can't open file `%s`\n", 
	      volatile_fname); fflush(stderr);
      exit(EXIT_FAILURE);
    }
    fread(&nsubgroubs, sizeof(int), 1, fd);
    if(nsubgroubs != Nsubgroups_B) {
      fprintf(stderr, "Error (read_subgroups_B): nsubgroubs != "
	      "Nsubgroups_B. Stop\n"); fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
    check = fread(&Volatile_B[1], sizeof(int), Nsubgroups_B, fd);
    if(check != Nsubgroups_B)	{
      fprintf(stderr, 
	      "Error (read_subgroups_B): Couldn't read volatile. Stop\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
    fclose(fd);
  }
  else {
    for (i = 1; i <= Nsubgroups_B; i++) {
      Volatile_B[i] = 0;
    }
  }
  
  if (!(fd = fopen(subids_fname, "r"))) {
    fprintf(stderr, "Error (read_subgroups_B): can't open file `%s`\n", 
	    subids_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fread(&nsubids_B, sizeof(int), 1, fd);
  
  idlist = ivector(0, nsubids_B - 1);
  
  check = fread(&idlist[0], sizeof(int), nsubids_B, fd);
  if (check != nsubids_B) {
    fprintf(stderr, "Error (read_subgroups_B): Couldn't read subids_B. "
	    "Stop. %d %d\n", check, nsubids_B); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  fclose(fd);

  for (i = 1; i <= Nsubgroups_B; i++) {
    SubIdMostBound_B[i] = idlist[SubFirst_B[i]];
  }
  free_ivector(idlist, 0, nsubids_B - 1);

  fflush(stdout);

  return;
}


/**
 * @brief Free memory allocated to subhaloes.
 */
void free_subgroups_B(void)
{
  free_ivector(Volatile_B, 1, Nsubgroups_B);
  free_ivector(FirstEntry_B, 1, Ngroups_B);
  free_ivector(GroupsSubCount_B, 1, Ngroups_B);
  free_ivector(SubIdMostBound_B, 1, Nsubgroups_B);
  free_ivector(SubParent_B, 1, Nsubgroups_B);
  free_ivector(SubFirst_B, 1, Nsubgroups_B);
  free_ivector(SubLen_B, 1, Nsubgroups_B);

  return;
}


/**
 * @brief Read subhaloes at the current simulation snapshot.
 */
void read_subgroups_A(char *catalogue_fname, char *substruc_fname, char *volatile_fname, char *subids_fname)
{
  FILE *fd;
  int i, check;
  int nsubids_A, *idlist;
  int nsubgroubs;
  int sqa;
  float dv2,dv[3],dx[3];
  float H_of_a;
  int j, p, next, first; 
  float GroupEkin;

  if (!(fd = fopen(catalogue_fname, "r"))) {
    fprintf(stderr, "Error (read_subgroups_A): can't open file `%s`\n", 
	    catalogue_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fread(&Nsubgroups_A, sizeof(int), 1, fd);

  SubLen_A = ivector(1, Nsubgroups_A);
  SubFirst_A = ivector(1, Nsubgroups_A);
  SubParent_A = ivector(1, Nsubgroups_A);
  SubIdMostBound_A = ivector(1, Nsubgroups_A);

  for (i = 1; i <= Nsubgroups_A; i++) {
    check = fread(&SubLen_A[i], sizeof(int), 1, fd);
    if (check != 1) {
      fprintf(stderr, 
	      "Error (read_subgroups_A): Couldn't read subgroups. Stop\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }

    check = fread(&SubFirst_A[i], sizeof(int), 1, fd);
    if (check != 1) {
      fprintf(stderr, 
	      "Error (read_subgroups_A): Couldn't read subgroups. Stop\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
    check = fread(&SubParent_A[i], sizeof(int), 1, fd);
    if (check != 1) {
      fprintf(stderr, 
	      "Error (read_subgroups_A): Couldn't read subgroups. Stop\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
  }
  fclose(fd);
  
  printf("Subgroup catalogue `%s' read.\n", catalogue_fname);
  printf("%d subgroups.\n\n", Nsubgroups_A);
  
  if(!(fd = fopen(substruc_fname, "r"))) {
    fprintf(stderr, 
	    "Error (read_subgroups_A): can't open file `%s`\n", 
	    substruc_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fread(&Ngroups_A, sizeof(int), 1, fd);

  GroupsSubCount_A = ivector(1, Ngroups_A);
  FirstEntry_A = ivector(1, Ngroups_A);

  check = fread(&GroupsSubCount_A[1], sizeof(int), Ngroups_A, fd);
  if (check != Ngroups_A) {
    fprintf(stderr, 
	    "Error (read_subgroups_A): Couldn't read subcount. Stop\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  check = fread(&FirstEntry_A[1], sizeof(int), Ngroups_A, fd);
  if (check != Ngroups_A) {
    fprintf(stderr, 
	    "Error (read_subgroups_A): Couldn't read firstentry. Stop\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
  fclose(fd);
  
  printf("Group catalogue `%s' read.\n", substruc_fname);
  printf("%d groups.\n\n", Ngroups_A);
  
  Volatile_A = ivector(1, Nsubgroups_A);
  if (*volatile_fname) {
    if (!(fd = fopen(volatile_fname, "r"))) {
      fprintf(stderr, "Error (read_subgroups_A): can't open file `%s`\n", 
	      volatile_fname); fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
    fread(&nsubgroubs, sizeof(int), 1, fd);
    if (nsubgroubs != Nsubgroups_A) {
      fprintf(stderr, "Fatal error! (read_subgroups_A) nsubgroubs = %d, "
	     "Nsubgroups_A = %d\n", nsubgroubs, Nsubgroups_A);
      fflush(stderr);
      exit(EXIT_FAILURE);
    }

    check = fread(&Volatile_A[1], sizeof(int), Nsubgroups_A, fd);
    if (check != Nsubgroups_A) {
      fprintf(stderr, "Error (read_subgroups_A): Couldn't read volatile. "
	      "Stop\n"); fflush(stderr);
      exit(EXIT_FAILURE);
    }
    fclose(fd);
  }
  else {
    for(i = 1; i <= Nsubgroups_A; i++) 
      Volatile_A[i] = 0;
  }
  
  if (!(fd = fopen(subids_fname, "r"))) {
    fprintf(stderr, "Error (read_subgroups_A): can't open file `%s`\n", 
	    subids_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fread(&nsubids_A, sizeof(int), 1, fd);
  
  idlist = ivector(0, nsubids_A - 1);
  
  check = fread(&idlist[0], sizeof(int), nsubids_A, fd);
  if (check != nsubids_A) {
    fprintf(stderr, 
	    "Error (read_subgroups_A): Couldn't read subids_A. Stop\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fclose(fd);

  for (i = 1; i <= Nsubgroups_A; i++)
    SubIdMostBound_A[i] = idlist[SubFirst_A[i]];
  
  free_ivector(idlist, 0, nsubids_A - 1);
  
  return;
}


/**
 * @brief Free memory allocated to subhaloes at the current simulation 
 * snapshot.
 */
void free_subgroups_A(void)
{
  free_ivector(Volatile_A, 1, Nsubgroups_A);
  free_ivector(FirstEntry_A, 1, Ngroups_A);
  free_ivector(GroupsSubCount_A, 1, Ngroups_A);
  free_ivector(SubIdMostBound_A, 1, Nsubgroups_A);
  free_ivector(SubParent_A, 1, Nsubgroups_A);
  free_ivector(SubFirst_A, 1, Nsubgroups_A);
  free_ivector(SubLen_A, 1, Nsubgroups_A);

  return;
}


/**
 * @brief Read haloes at the previous simulation snapshot.
 * @param catalogue_fname Name of postprocessing catalogue file
 */
void read_groups_B(char *catalogue_fname)
{
  FILE *fd;

  if (!(fd = fopen(catalogue_fname, "r"))) {
    fprintf(stderr, "Error (read_groups_B): can't open file `%s`\n", 
	    catalogue_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fread(&Ngroups_B, sizeof(int), 1, fd);

  GroupLen_B = ivector(1, Ngroups_B);
  GroupTag_B = ivector(1, Ngroups_B);
  fread(&GroupLen_B[1], sizeof(int), Ngroups_B, fd);
  fread(&GroupTag_B[1], sizeof(int), Ngroups_B, fd);

  printf("Group catalogue `%s' read.\n", catalogue_fname);
  printf("%d groups.\n\n", Ngroups_B);
  fflush(stdout);

  fclose(fd);

  return;
}


/**
 * @brief Free memory allocated to haloes at the previous simulation 
 * snapshot.
 */
void free_groups_B(void)
{
  free_ivector(GroupTag_B, 1, Ngroups_B);
  free_ivector(GroupLen_B, 1, Ngroups_B);

  return;
}


/**
 * @brief Read haloes at the current simulation snapshot.
 */
void read_groups_A(char *catalogue_fname)
{
  FILE *fd;

  if (!(fd = fopen(catalogue_fname, "r"))) {
    fprintf(stderr, "Error (read_groups_A): can't open file `%s`\n", 
	    catalogue_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fread(&Ngroups_A, sizeof(int), 1, fd);

  GroupLen_A = ivector(1, Ngroups_A);
  GroupTag_A = ivector(1, Ngroups_A);
  fread(&GroupLen_A[1], sizeof(int), Ngroups_A, fd);
  fread(&GroupTag_A[1], sizeof(int), Ngroups_A, fd);

  fclose(fd);

  printf("Group catalogue `%s' read.\n", catalogue_fname);
  printf("%d groups.\n\n", Ngroups_A);
  fflush(stdout);
  
  return;
}


/**
 * @brief Free memory allocated to haloes at the current simulation 
 * snapshot.
 */
void free_groups_A(void)
{
  free_ivector(GroupTag_A, 1, Ngroups_A);
  free_ivector(GroupLen_A, 1, Ngroups_A);
  
  return;
}


/**
 * @brief Read linking list at the previous simulation snapshot.
 */
void read_linklist_B(char *linklist_fname)
{
  FILE *fd;
  struct stat buffer;
  int i, p, grB;

  if (!NumPart) {
    stat(linklist_fname, &buffer);
    NumPart = buffer.st_size / sizeof(int);
    printf("NumPart= %d\n", NumPart); fflush(stdout);
  }

  Tag_B = ivector(1, NumPart);	/* will store the group number of each particle */
  Next_B = ivector(1, NumPart);

  if (!(fd = fopen(linklist_fname, "r"))) {
    fprintf(stderr, "Error (read_linklist_B): can't open file `%s`\n", 
	    linklist_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fread(&Next_B[1], sizeof(int), NumPart, fd);
  fclose(fd);

  for (i = 1; i <= NumPart; i++)
    Tag_B[i] = 0;
  
  for (grB = 1; grB <= Ngroups_B; grB++) {
    for (i = 1, p = GroupTag_B[grB]; i <= GroupLen_B[grB]; i++) {
      Tag_B[p] = grB;
      p = Next_B[p];
    }
  }
 
  return;
}


/**
 * @brief Free memory allocated to linking list at the previous 
 * simulation snapshot.
 */
void free_linklist_B(void)
{
  free_ivector(Tag_B, 1, NumPart);
  free_ivector(Next_B, 1, NumPart);

  return;
}


/**
 * @brief Read linking list at the current simulation snapshot.
 */
void read_linklist_A(char *linklist_fname)
{
  FILE *fd;
  struct stat buffer;
  int i, p, grA;

  if (!NumPart) {
    stat(linklist_fname, &buffer);
    NumPart = buffer.st_size / sizeof(int);
    printf("NumPart= %d\n", NumPart); fflush(stdout);
  }

  Tag_A = ivector(1, NumPart);	/* will store the group number of each particle */
  Next_A = ivector(1, NumPart);

  if (!(fd = fopen(linklist_fname, "r"))) {
    fprintf(stderr, "Error (read_linklist_A): can't open file `%s`\n", 
	    linklist_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  fread(&Next_A[1], sizeof(int), NumPart, fd);
  fclose(fd);
    
  for (i = 1; i <= NumPart; i++)
    Tag_A[i] = 0;
  
  for (grA = 1; grA <= Ngroups_A; grA++) {
    for (i = 1, p = GroupTag_A[grA]; i <= GroupLen_A[grA]; i++) {
      Tag_A[p] = grA;
      p = Next_A[p];
    }
  }
 
  return;
}


/**
 * @brief Free memory allocated to linking list at the current
 * simulation snapshot.
 */
void free_linklist_A(void)
{
  free_ivector(Tag_A, 1, NumPart);
  free_ivector(Next_A, 1, NumPart);

  return;
}


void read_history(char *history_fname)
{
  FILE *fd;
  
  if (!(fd = fopen(history_fname, "r"))) {
    fprintf(stderr, "Error (read_history): can't open file `%s`\n", 
	    history_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  fread(&Ngroups_A, sizeof(int), 1, fd);
  fread(&Ngroups_B, sizeof(int), 1, fd);

  CountFOFProgenitor = ivector(1, Ngroups_A);
  FirstFOFProgenitor = ivector(1, Ngroups_A);
  NextFOFProgenitor = ivector(1, Ngroups_B);

  fread(&CountFOFProgenitor[1], sizeof(int), Ngroups_A, fd);
  fread(&FirstFOFProgenitor[1], sizeof(int), Ngroups_A, fd);
  fread(&NextFOFProgenitor[1], sizeof(int), Ngroups_B, fd);

  fclose(fd);

  printf("History-Link `%s' read.\n\n", history_fname); fflush(stdout);

  return;
}


/**
 * @brief Free memory allocated to halo history.
 */
void free_history(void)
{
  free_ivector(CountFOFProgenitor, 1, Ngroups_A);
  free_ivector(FirstFOFProgenitor, 1, Ngroups_A);
  free_ivector(NextFOFProgenitor, 1, Ngroups_B);

  return;
}


/**
 * @brief Read subhalo history from postprocessing files.
 */
void read_subhistory(char *history_fname)
{
  FILE *fd;
  int check;

  if (!(fd = fopen(history_fname, "r"))) {
    fprintf(stderr, "Error (read_subhistory): can't open file `%s`\n", 
	    history_fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  fread(&Nsubgroups_A, sizeof(int), 1, fd);
  fread(&Nsubgroups_B, sizeof(int), 1, fd);

  CountProgenitor = ivector(1, Nsubgroups_A);
  FirstProgenitor = ivector(1, Nsubgroups_A);
  NextProgenitor = ivector(1, Nsubgroups_B);

  check = fread(&CountProgenitor[1], sizeof(int), Nsubgroups_A, fd);
  checkerror(check, Nsubgroups_A);
  check = fread(&FirstProgenitor[1], sizeof(int), Nsubgroups_A, fd);
  checkerror(check, Nsubgroups_A);
  check = fread(&NextProgenitor[1], sizeof(int), Nsubgroups_B, fd);
  checkerror(check, Nsubgroups_B);

  fclose(fd);

  printf("History-Link `%s' read.\n\n", history_fname); fflush(stdout);

  return;
}


/**
 * @brief Free memory allocated to subhalo history.
 */
void free_subhistory(void)
{
  free_ivector(NextProgenitor, 1, Nsubgroups_B);
  free_ivector(FirstProgenitor, 1, Nsubgroups_A);
  free_ivector(CountProgenitor, 1, Nsubgroups_A);

  return;
}


/**
 * @brief Read subhalo properties from postprocessing files.
 * @param fname Name of postprocessing subproperties file
 * @param zcurr Current redshift
 * @param subids_fname Name of postprocessing subids file
 */
void read_subgroup_properties_B(char *fname, double zcurr, 
				char *subids_fname)
{
  FILE *fd;
  int  gr, nsubgroups, check;
  double hubble_of_z, rhocrit, fac;
  int nsubids_B, *idlist;
  int sqa;
  float dv2,dv[3],dx[3];
  float H_of_a;
  int i, j, p, next, first;
  float SubsEkin,s;
  char input_fname[256];
  int k, type;


  if (!(fd = fopen(fname, "r"))) {
    fprintf(stderr, 
	    "Error (read_subgroup_properties_B): can't open file '%s'\n", 
	    fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  fread(&nsubgroups, sizeof(int), 1, fd);

  if (Nsubgroups_B != nsubgroups) {
    fprintf(stderr, 
	    "Fatal error! (read_subgroup_properties_B): %s\n", fname);
    fprintf(stderr, 
	    "Nsubgroups_B = %d, nsubgroups = %d\n", Nsubgroups_B, 
	    nsubgroups); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  if (Nsubgroups_B) {
    SubGroupVc_B = vector(1, Nsubgroups_B);
    SubGroupVvir_B = vector(1, Nsubgroups_B);
    SubGroupRvir_B = vector(1, Nsubgroups_B);
    SubGroupMvir_B = vector(1, Nsubgroups_B);
    SubGroupPos_B = matrix(1, Nsubgroups_B, 0, 2);
    SubGroupVel_B = matrix(1, Nsubgroups_B, 0, 2);
    
    //  printf("Nsubgroup de B: %d\n",Nsubgroups_B);
    
    check = fread(&SubGroupVc_B[1], sizeof(float), Nsubgroups_B, fd);
    checkerror(check, Nsubgroups_B);
    check = fread(&SubGroupPos_B[1][0], 3 * sizeof(float), Nsubgroups_B, fd);
    checkerror(check, Nsubgroups_B);
    check = fread(&SubGroupVel_B[1][0], 3 * sizeof(float), Nsubgroups_B, fd);
    checkerror(check, Nsubgroups_B);
  }
  else {
    SubGroupVc_B = SubGroupVvir_B = SubGroupRvir_B = SubGroupMvir_B = 0;
    SubGroupPos_B = SubGroupVel_B = 0;
  }
  fclose(fd);
  
  hubble_of_z = Hubble * sqrt(Omega * pow(1+zcurr, 3) + 
			      (1-Omega-OmegaLambda)*pow(1+zcurr, 2) + 
			      OmegaLambda);
  rhocrit = 3*hubble_of_z*hubble_of_z / (8*M_PI*G);
  fac = 1.0 / (200 * 4 * M_PI / 3.0 * rhocrit);
  
  for (gr = 1; gr <= Nsubgroups_B; gr++) {
    SubGroupMvir_B[gr] = SubLen_B[gr] * PartMass;
    SubGroupRvir_B[gr] = pow(SubGroupMvir_B[gr] * fac, 1.0 / 3);
    SubGroupVvir_B[gr] = sqrt(G * SubGroupMvir_B[gr] / SubGroupRvir_B[gr]);
    SubGroupVc_B[gr] = SubGroupVvir_B[gr] * VcFactor;
    
    /*
      printf("gr=%d, SubGroupMvir_B[gr]: %g, SubGroupRvir_B[gr]:%g \n\n",
      gr,SubGroupMvir_B[gr],SubGroupRvir_B[gr]);
    */
  }
  
  return;
}


/**
 * @brief Free memory allocated to subhalo properties at the previous 
 * simulation snapshot.
 */
void free_subproperties_B(void)
{
  if (Nsubgroups_B) {
    free_matrix(SubGroupVel_B, 1, Nsubgroups_B, 0, 2);
    free_matrix(SubGroupPos_B, 1, Nsubgroups_B, 0, 2);
    free_vector(SubGroupMvir_B, 1, Nsubgroups_B);
    free_vector(SubGroupRvir_B, 1, Nsubgroups_B);
    free_vector(SubGroupVvir_B, 1, Nsubgroups_B);
    free_vector(SubGroupVc_B, 1, Nsubgroups_B);
  }
  
  return;
}


void read_subgroup_properties_A(char *fname, double zcurr, char *subids_fname)
{
  FILE *fd;
  int gr, nsubgroups, check;
  double hubble_of_z, rhocrit, fac;
  int nsubids_A, *idlist;
  int sqa;
  float dv2,dv[3],dx[3];
  float H_of_a;
  int i, j, p, next, first;
  float SubsEkin,s;
  char input_fname[256];
  int k, type;
    
  
  if (!(fd = fopen(fname, "r"))) {
    fprintf(stderr, 
	    "Error (read_subgroup_properties_A): can't open file '%s'\n", 
	    fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  fread(&nsubgroups, sizeof(int), 1, fd);

  if (Nsubgroups_A != nsubgroups) {
    fprintf(stderr,
	    "Fatal error! (read_subgroup_properties_A): %s\n", fname);
    fprintf(stderr, "Nsubgroups_A = %d, nsubgroups = %d\n", Nsubgroups_A, 
	    nsubgroups); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  SubGroupVc_A = vector(1, Nsubgroups_A);
  SubGroupVvir_A = vector(1, Nsubgroups_A);
  SubGroupRvir_A = vector(1, Nsubgroups_A);
  SubGroupMvir_A = vector(1, Nsubgroups_A);
  SubGroupPos_A = matrix(1, Nsubgroups_A, 0, 2);
  SubGroupVel_A = matrix(1, Nsubgroups_A, 0, 2);
  
  check = fread(&SubGroupVc_A[1], sizeof(float), Nsubgroups_A, fd);
  checkerror(check, Nsubgroups_A);
  check = fread(&SubGroupPos_A[1][0], 3 * sizeof(float), Nsubgroups_A, fd);
  checkerror(check, Nsubgroups_A);
  check = fread(&SubGroupVel_A[1][0], 3 * sizeof(float), Nsubgroups_A, fd);
  checkerror(check, Nsubgroups_A);
  
  fclose(fd);
  
  hubble_of_z =
    Hubble * sqrt(Omega * pow(1 + zcurr, 3) + (1 - Omega - OmegaLambda) * pow(1 + zcurr, 2) + OmegaLambda);
  rhocrit = 3 * hubble_of_z * hubble_of_z / (8 * M_PI * G);
  fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit);
  
  for (gr = 1; gr <= Nsubgroups_A; gr++) {
    SubGroupMvir_A[gr] = SubLen_A[gr] * PartMass;
    SubGroupRvir_A[gr] = pow(SubGroupMvir_A[gr] * fac, 1.0 / 3);
    SubGroupVvir_A[gr] = sqrt(G * SubGroupMvir_A[gr] / SubGroupRvir_A[gr]);
    SubGroupVc_A[gr] = SubGroupVvir_A[gr] * VcFactor;
    /*
      printf("gr=%d, PartMass: %g, SubLen_A[gr]: %d, SubGroupMvir_A[gr]: %g, SubGroupRvir_A[gr]:%g \n\n",
      gr,PartMass,SubLen_A[gr],SubGroupMvir_A[gr],SubGroupRvir_A[gr]);
    */
  }
  
  return;
}


/**
 * @brief Free memory allocated to subhalo properties at the current
 * simulation snapshot.
 */
void free_subproperties_A(void)
{
  free_matrix(SubGroupVel_A, 1, Nsubgroups_A, 0, 2);
  free_matrix(SubGroupPos_A, 1, Nsubgroups_A, 0, 2);
  free_vector(SubGroupMvir_A, 1, Nsubgroups_A);
  free_vector(SubGroupRvir_A, 1, Nsubgroups_A);
  free_vector(SubGroupVvir_A, 1, Nsubgroups_A);
  free_vector(SubGroupVc_A, 1, Nsubgroups_A);

  return;
}


/**
 * @brief Read FOF group properties from postprocessing files.
 * @param fname Name of postprocessing group properties file
 * @param zcurr Current redshift
 */
void read_group_properties_B(char *fname, double zcurr)
{
  FILE *fd;
  int gr;
  size_t check;
  double hubble_of_z, rhocrit, fac;


  if (!(fd = fopen(fname, "r"))) {
    fprintf(stderr,
	    "Error (read_group_properties_B): can't open file '%s'\n", 
	    fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  fread(&Ngroups_B, sizeof(int), 1, fd);
  GroupM200_B = vector(1, Ngroups_B);
  GroupR200_B = vector(1, Ngroups_B);
  GroupV200_B = vector(1, Ngroups_B);
  GroupMvir_B = vector(1, Ngroups_B);
  GroupRvir_B = vector(1, Ngroups_B);
  GroupVvir_B = vector(1, Ngroups_B);
  GroupVc_B = vector(1, Ngroups_B);
  GroupVc200_B = vector(1, Ngroups_B);
  GroupMostBound_B = ivector(1, Ngroups_B);
  
  GroupCenter_B = matrix(1, Ngroups_B, 0, 2);
  GroupPos_B = matrix(1, Ngroups_B, 0, 2);
  GroupVel_B = matrix(1, Ngroups_B, 0, 2);
  
  GroupL_B = matrix(1, Ngroups_B, 0, 2);
  GroupLambda_B = vector(1, Ngroups_B);
  
  GroupEkin_B = vector(1, Ngroups_B);
  GroupEpot_B = vector(1, Ngroups_B);
  
  check = fread(&GroupM200_B[1], sizeof(float), Ngroups_B, fd);
  checkerror(check, Ngroups_B);
  check = fread(&GroupR200_B[1], sizeof(float), Ngroups_B, fd);
  checkerror(check, Ngroups_B);
  
  check = fread(&GroupMvir_B[1], sizeof(float), Ngroups_B, fd);
  checkerror(check, Ngroups_B);
  check = fread(&GroupRvir_B[1], sizeof(float), Ngroups_B, fd);
  checkerror(check, Ngroups_B);
  
  fread(&GroupCenter_B[1][0], 3*sizeof(float), Ngroups_B, fd);
  fread(&GroupMostBound_B[1], sizeof(int), Ngroups_B, fd);
  
  check = fread(&GroupPos_B[1][0], 3 * sizeof(float), Ngroups_B, fd);
  checkerror(check, Ngroups_B);
  check = fread(&GroupVel_B[1][0], 3 * sizeof(float), Ngroups_B, fd);
  checkerror(check, Ngroups_B);
  
  check = fread(&GroupEkin_B[1], sizeof(float), Ngroups_B, fd); 
  checkerror(check, Ngroups_B);
  check = fread(&GroupEpot_B[1], sizeof(float), Ngroups_B, fd);
  checkerror(check, Ngroups_B);
  
  check = fread(&GroupLambda_B[1], sizeof(float), Ngroups_B, fd);
  checkerror(check, Ngroups_B);
  check = fread(&GroupL_B[1][0], 3 * sizeof(float), Ngroups_B, fd);
  checkerror(check, Ngroups_B);
  fclose(fd);
  
  hubble_of_z = Hubble * sqrt(Omega * pow(1 + zcurr, 3) + 
			      (1 - Omega - OmegaLambda) * pow(1 + zcurr, 2) 
			      + OmegaLambda);
  rhocrit = 3 * hubble_of_z * hubble_of_z / (8 * M_PI * G);
  fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit);
  
  for (gr = 1; gr <= Ngroups_B; gr++) {
    /*GroupMvir_B[gr] = GroupM200_B[gr];*/ /* here in 10^10 Msol/h */
    
    /* Recompute virial radius from group mass */
    GroupRvir_B[gr] = pow(GroupMvir_B[gr] * fac, 1.0 / 3);
    GroupR200_B[gr] = pow(GroupM200_B[gr] * fac, 1.0 / 3);
    
    GroupVvir_B[gr] = sqrt(G * GroupMvir_B[gr] / GroupRvir_B[gr]);
    GroupVc_B[gr] = VcFactor * GroupVvir_B[gr];

    GroupV200_B[gr] = sqrt(G * GroupM200_B[gr] / GroupR200_B[gr]);
    GroupVc200_B[gr] = VcFactor * GroupV200_B[gr];
  }
  
  printf("Read properties.\n"); fflush(stdout);
  
  return;
}


/**
 * @brief Free memory allocated to halo properties at the previous 
 * simulation snapshot.
 */
void free_properties_B(void)
{
  free_vector(GroupVc_B, 1, Ngroups_B);
  free_vector(GroupVvir_B, 1, Ngroups_B);
  free_vector(GroupVc200_B, 1, Ngroups_B);
  free_vector(GroupV200_B, 1, Ngroups_B);

  free_matrix(GroupL_B, 1, Ngroups_B, 0, 2);
  free_vector(GroupLambda_B,  1, Ngroups_B);

  free_vector(GroupEpot_B, 1, Ngroups_B);
  free_vector(GroupEkin_B, 1, Ngroups_B);

  free_matrix(GroupVel_B, 1, Ngroups_B, 0, 2);
  free_matrix(GroupPos_B, 1, Ngroups_B, 0, 2);
  free_matrix(GroupCenter_B, 1, Ngroups_B, 0, 2);

  free_ivector(GroupMostBound_B, 1, Ngroups_B);
  free_vector(GroupRvir_B, 1, Ngroups_B);
  free_vector(GroupMvir_B, 1, Ngroups_B);
  free_vector(GroupR200_B, 1, Ngroups_B);
  free_vector(GroupM200_B, 1, Ngroups_B);

  return;
}


void read_group_properties_A(char *fname, double zcurr)
{
  FILE *fd;
  int gr;
  size_t check;
  double hubble_of_z, rhocrit, fac;


  if (!(fd = fopen(fname, "r"))) {
    fprintf(stderr, 
	    "Error (read_group_properties_A): can't open file `%s`\n", 
	    fname); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  fread(&Ngroups_A, sizeof(int), 1, fd);
  
  GroupM200_A = vector(1, Ngroups_A);
  GroupR200_A = vector(1, Ngroups_A);
  GroupV200_A = vector(1, Ngroups_A);
  
  GroupMvir_A = vector(1, Ngroups_A);
  GroupRvir_A = vector(1, Ngroups_A);
  GroupVvir_A = vector(1, Ngroups_A);
  GroupVc_A = vector(1, Ngroups_A);
  GroupVc200_A = vector(1, Ngroups_A);
  GroupMostBound_A = ivector(1, Ngroups_A);
  
  GroupCenter_A = matrix(1, Ngroups_A, 0, 2);
  GroupPos_A = matrix(1, Ngroups_A, 0, 2);
  GroupVel_A = matrix(1, Ngroups_A, 0, 2);
  
  GroupL_A = matrix(1, Ngroups_A, 0, 2);
  GroupLambda_A = vector(1, Ngroups_A);
  
  GroupEkin_A = vector(1, Ngroups_A);
  GroupEpot_A = vector(1, Ngroups_A);
  
  check = fread(&GroupM200_A[1], sizeof(float), Ngroups_A, fd);
  checkerror(check, Ngroups_A);
  
  fread(&GroupR200_A[1], sizeof(float), Ngroups_A, fd);
  
  fread(&GroupMvir_A[1], sizeof(float), Ngroups_A, fd);
  fread(&GroupRvir_A[1], sizeof(float), Ngroups_A, fd);
    
  fread(&GroupCenter_A[1][0], 3*sizeof(float), Ngroups_A, fd);
  
  fread(&GroupMostBound_A[1], sizeof(int), Ngroups_A, fd);
  
  fread(&GroupPos_A[1][0], 3 * sizeof(float), Ngroups_A, fd);
  fread(&GroupVel_A[1][0], 3 * sizeof(float), Ngroups_A, fd);
  
  fread(&GroupEkin_A[1], sizeof(float), Ngroups_A, fd);
  fread(&GroupEpot_A[1], sizeof(float), Ngroups_A, fd);
  
  fread(&GroupLambda_A[1], sizeof(float), Ngroups_A, fd);
  fread(&GroupL_A[1][0], 3 * sizeof(float), Ngroups_A, fd);
  
  fclose(fd);
  
  hubble_of_z = Hubble * sqrt(Omega * pow(1 + zcurr, 3) + 
			      (1 - Omega - OmegaLambda) * pow(1 + zcurr, 2)
			      + OmegaLambda);
  rhocrit = 3 * hubble_of_z * hubble_of_z / (8 * M_PI * G);
  fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit);
    
  for (gr = 1; gr <= Ngroups_A; gr++) {
    /*GroupMvir_A[gr] = GroupM200_A[gr]; */ /* here in 10^10 Msol/h */

    GroupRvir_A[gr] = pow(GroupMvir_A[gr] * fac, 1.0 / 3);
    GroupR200_A[gr] = pow(GroupM200_A[gr] * fac, 1.0 / 3);

    GroupVvir_A[gr] = sqrt(G * GroupMvir_A[gr] / GroupRvir_A[gr]);
    GroupVc_A[gr] = VcFactor * GroupVvir_A[gr];

    GroupV200_A[gr] = sqrt(G * GroupM200_A[gr] / GroupR200_A[gr]);
    GroupVc200_A[gr] = VcFactor * GroupV200_A[gr];
  }

  printf("Read properties.\n"); fflush(stdout);

  return;
}


/**
 * @brief Free memory allocated to subhalo properties at the current
 * simulation snapshot.
 */
void free_properties_A(void)
{
  free_vector(GroupVc_A, 1, Ngroups_A);
  free_vector(GroupVvir_A, 1, Ngroups_A);
  free_vector(GroupVc200_A, 1, Ngroups_A);
  free_vector(GroupV200_A, 1, Ngroups_A);

  free_matrix(GroupL_A, 1, Ngroups_A, 0, 2);
  free_vector(GroupLambda_A,  1, Ngroups_A);

  free_vector(GroupEpot_A, 1, Ngroups_A);
  free_vector(GroupEkin_A, 1, Ngroups_A);

  free_matrix(GroupVel_A, 1, Ngroups_A, 0, 2);
  free_matrix(GroupPos_A, 1, Ngroups_A, 0, 2);
  free_matrix(GroupCenter_A, 1, Ngroups_A, 0, 2);

  free_ivector(GroupMostBound_A, 1, Ngroups_A);
  free_vector(GroupRvir_A, 1, Ngroups_A);
  free_vector(GroupMvir_A, 1, Ngroups_A);
  free_vector(GroupR200_A, 1, Ngroups_A);
  free_vector(GroupM200_A, 1, Ngroups_A);

  return;
}


/**
 * @brief Checks for error during reading of postprocessing data.
 */
void checkerror(int check, int length)
{
  if(check != length) {
    fprintf(stderr,
	    "ERROR (checkerror) during reading. check=%d, length=%d. Stop\n", 
	    check, length); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  return;
}
