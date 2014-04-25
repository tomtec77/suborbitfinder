/**
 * @file allvars.c
 * @brief File with global variables to link with main code.
 */
#include "allvars.h"

char OrdTag[100];
char PreTag[100];

int GalaxynumberB;
int Ngroups_A, Nsubgroups_A;
int Ngroups_B, Nsubgroups_B;
int NumPart;
int Snapshot;
int TotNumGalA;
int TotNumGalB;
int TotNumType0A, TotNumType1A, TotNumType2A, TotNumType3A;
int TotNumType0B, TotNumType1B, TotNumType2B, TotNumType3B;
int TotNumType4A;
int TotNumType4B;

int *CountFOFProgenitor, *FirstFOFProgenitor, *NextFOFProgenitor;
int *CountProgenitor, *FirstProgenitor, *NextProgenitor;
int *GroupMostBound_A, *GroupLen_A, *GroupTag_A;
int *GroupMostBound_B, *GroupLen_B, *GroupTag_B;
int *GroupsSubCount_A, *FirstEntry_A;
int *GroupsSubCount_B, *FirstEntry_B;
int *Id;
int *Next_A, *Tag_A;
int *Next_B, *Tag_B;
int *NumGalInFOFGroup_A, *FirstGalInFOFGroup_A;
int *NumGalInFOFGroup_B, *FirstGalInFOFGroup_B;
int *NumGalInSubGroup_A, *FirstGalInSubGroup_A;
int *NumGalInSubGroup_B, *FirstGalInSubGroup_B;
int *SubGalDescendant, *SubGalFirstProg, *SubGalCountProg;
int *SubGalNextProg, *SubGalMerging;
int *SubLen_A, *SubFirst_A, *SubParent_A, *SubIdMostBound_A;
int *SubLen_B, *SubFirst_B, *SubParent_B, *SubIdMostBound_B;
int *Volatile_A;
int *Volatile_B;

unsigned int Seed;
unsigned int Count_orphan_type1;

float DeltaT;
float Zprev, Zcurr;

float ZZ[OUTPUTS], Age[OUTPUTS];

float *GroupCvir_A, *GroupCvirMedian_A;
float *GroupEpot_A, *GroupEkin_A;
float *GroupEpot_B, *GroupEkin_B;
float *GroupLambda_A;
float *GroupLambda_B;
float *GroupR200_A, *GroupM200_A, *GroupV200_A, *GroupVc200_A;
float *GroupR200_B, *GroupM200_B, *GroupV200_B, *GroupVc200_B;
float *GroupRvir_A, *GroupMvir_A, *GroupVvir_A, *GroupVc_A;
float *GroupRvir_B, *GroupMvir_B, *GroupVvir_B, *GroupVc_B;
float *SubGroupVc_A, *SubGroupVvir_A, *SubGroupRvir_A, *SubGroupMvir_A;
float *SubGroupVc_B, *SubGroupVvir_B, *SubGroupRvir_B, *SubGroupMvir_B;

float **GroupL_A;
float **GroupL_B;
float **GroupPos_A, **GroupVel_A, **GroupCenter_A;
float **GroupPos_B, **GroupVel_B, **GroupCenter_B;
float **SubGroupPos_A, **SubGroupVel_A;
float **SubGroupPos_B, **SubGroupVel_B;

double GPMass;
double Hubble;
double LightSpeed;
double PartMass, PartMassGas, PartMassDM;
double RhoCrit;
double UnitTime_in_s, UnitTime_in_Megayears, UnitMass_in_Msun;
double UnitLength_in_kpc, UnitVelocity_in_km_per_s;
double UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
double UnitCoolingRate_in_cgs;

unsigned int TSTimescaleSelect;
unsigned int RtidalSelect;

struct particle_data *P;
struct GALAXY *GalaxyA, *GalaxyB;

unsigned int EmergencyStop;

#ifdef DEBUG
unsigned int Nprints;
#endif

