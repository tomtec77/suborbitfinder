/**
 * @file allvars.h
 * @brief Defines and global variable declarations.
 * @author Tomas E. Tecce
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "lib/nrsrc/nrutil.h"
#include "lib/nrsrc/nrsag.h"

/* Set the total number of snapshots */
#ifdef DOLAGSIMULATION
#define OUTPUTS 93
#else
#define OUTPUTS 100
#endif

/** Output steps to store orbit information. */
#define STEPS 50    

/** Number of intervals to subdivide the orbital timescale, to use for 
    orbits integration */
#define INTEGRATIONDT 5000.0 

/** Maximum number of allowed consecutive calls to the orbit integration
    function. */
#define NORBTIMES 10

/** Maximum number of groups allowed. Used for memory allocation. */
#define MAXGROUPS   120000

/** Maximum number of galaxies allowed. Used for memory allocation. */
#define MAXGALAXIES 120000

/** Speed of light in cm/s. */
#define C                2.9979e10

/** 1 Mpc in cm. */
#define CM_PER_MPC       3.085678e24

/** Gravitational constant in SI units. */
#define GRAVITY          6.672e-8

/** Hubble constant in h/sec. */
#define HUBBLE           3.2407789e-18

/** 1 Myr in seconds. */ 
#define SEC_PER_MEGAYEAR 3.155e13

/** Mass of the Sun in g. */
#define SOLAR_MASS       1.989e33

/** Gravitational constant in kpc M_Sun^-1 (km/s)^2 */
#define GRAVITY_KPCMSUN  4.3e-6

/** Satellite disruption threshold (fraction of initial mass). */
#define SAT_DISRUPT      0.01 

/** Minimum ratio of current to initial angular momentum, used to determine
    when a merger occurs. */
#define SAT_MINJORB      0.01

/** Numerical ID for the approximation for the tidal radius from
    Tormen, Diaferio & Syer (1998). */
#define RTIDAL_TDS98 1

/** Numerical ID for the approximation for the tidal radius from 
    Zentner & Bullock (2003). */
#define RTIDAL_ZB03 2

/** Numerical ID for the choice of subhalo dynamical time at infall as tidal 
    stripping timescale. */
#define TSTIME_SUBINFALL 1

/** Numerical ID for the choice of current subhalo dynamical time as tidal
    stripping timescale. */
#define TSTIME_SUBHALO 2


/* Unit length in kpc for orbit calculation */
/* TOMAS 2013-04-24: from now on, use the virial radius of the host halo
   as the unit length. This has the advantage that it sets both G and the
   virial velocity to unity */
/*#define ULENGTH_ORBIT    1000.0 */

extern char OrdTag[100];
extern char PreTag[100];

extern int GalaxynumberB;
extern int Ngroups_A, Nsubgroups_A;
extern int Ngroups_B, Nsubgroups_B;
extern int NumPart;
extern int Snapshot;

/** Counter for the total number of galaxies at the current snapshot. */
extern int TotNumGalA;

/** Counter for the total number of galaxies at the previous snapshot. */
extern int TotNumGalB;

/** Counter for the total number of type 0 galaxies at the current 
    snapshot. */
extern int TotNumType0A;

/** Counter for the total number of type 1 galaxies at the current 
    snapshot. */
extern int TotNumType1A;

/** Counter for the total number of type 2 galaxies at the current 
    snapshot. */
extern int TotNumType2A;

/** Counter for the total number of merged (type 3) galaxies at the current 
    snapshot. */
extern int TotNumType3A;

/** Counter for the total number of disrupted (type 4) galaxies at the 
    current snapshot. */
extern int TotNumType4A;

/** Counter for the total number of type 0 galaxies at the previous 
    snapshot. */
extern int TotNumType0B;

/** Counter for the total number of type 1 galaxies at the previous 
    snapshot. */
extern int TotNumType1B;

/** Counter for the total number of type 2 galaxies at the previous 
    snapshot. */
extern int TotNumType2B;

/** Counter for the total number of merged (type 3) galaxies at the 
    previous snapshot. */
extern int TotNumType3B;

/** Counter for the total number of disrupted (type 4) galaxies at the 
    previous snapshot. */
extern int TotNumType4B;

extern int *CountFOFProgenitor, *FirstFOFProgenitor, *NextFOFProgenitor;
extern int *CountProgenitor, *FirstProgenitor, *NextProgenitor;
extern int *GroupMostBound_A, *GroupLen_A, *GroupTag_A;
extern int *GroupMostBound_B, *GroupLen_B, *GroupTag_B;
extern int *GroupsSubCount_A, *FirstEntry_A;
extern int *GroupsSubCount_B, *FirstEntry_B;
extern int *Id;
extern int *Next_A, *Tag_A;
extern int *Next_B, *Tag_B;

/** Total number of galaxies in FOF group at current snapshot. */
extern int *NumGalInFOFGroup_A;

/** Index of first galaxy in FOF group at current snapshot. */
extern int *FirstGalInFOFGroup_A;

/** Total number of galaxies in FOF group at previous snapshot. */
extern int *NumGalInFOFGroup_B;

/** Index of first galaxy in FOF group at previous snapshot. */
extern int *FirstGalInFOFGroup_B;

/** Total number of galaxies in subhalo at current snapshot. */
extern int *NumGalInSubGroup_A;

/** Index of first galaxy in subhalo at current snapshot. */
extern int *FirstGalInSubGroup_A;

/** Total number of galaxies in subhalo at previous snapshot. */
extern int *NumGalInSubGroup_B;

/** Index of first galaxy in FOF group at previous snapshot. */
extern int *FirstGalInSubGroup_B;

extern int *SubGalDescendant, *SubGalFirstProg, *SubGalCountProg;
extern int *SubGalNextProg, *SubGalMerging;

/** Number of particles bound to the subhalo at the current snapshot. */
extern int *SubLen_A;

extern int *SubFirst_A, *SubParent_A, *SubIdMostBound_A;

/** Number of particles bound to the subhalo at the previous snapshot. */
extern int *SubLen_B;

extern int *SubFirst_B, *SubIdMostBound_B;

/** Index of parent FOF halo of a subhalo at previous snapshot. */
extern int *SubParent_B;

extern int *Volatile_A;
extern int *Volatile_B;

extern unsigned int Seed;
extern unsigned int Count_orphan_type1;

extern float DeltaT;
extern float Zprev, Zcurr;

extern float ZZ[OUTPUTS], Age[OUTPUTS];

extern float *GroupCvir_A, *GroupCvirMedian_A;
extern float *GroupEpot_A, *GroupEkin_A;
extern float *GroupEpot_B, *GroupEkin_B;
extern float *GroupLambda_A;
extern float *GroupLambda_B;
extern float *GroupR200_A, *GroupM200_A, *GroupV200_A, *GroupVc200_A;
extern float *GroupR200_B, *GroupM200_B, *GroupV200_B, *GroupVc200_B;
extern float *GroupRvir_A, *GroupMvir_A, *GroupVvir_A, *GroupVc_A;
extern float *GroupRvir_B, *GroupMvir_B, *GroupVvir_B, *GroupVc_B;

/** Virial velocity times VcFactor, at the current snapshot. */
extern float *SubGroupVc_A;

/** Circular velocity at the virial radius in the current snapshot, 
    calculated from the virial radius and mass read from the postprocessing 
    files: 
    @f$V_\mathrm{vir} = \sqrt{G M_\mathrm{vir} / r_\mathrm{vir}}@f$ */
extern float *SubGroupVvir_A;

/** Subhalo virial radius at the current snapshot, calculated from the
    value of the virial mass. */
extern float *SubGroupRvir_A;

/** Subhalo virial mass at the current snapshot, determined from the 
    postprocessing files as the number of bound particles times the
    particle mass. */
extern float  *SubGroupMvir_A;

/** Virial velocity times VcFactor, at the previous snapshot. */
extern float *SubGroupVc_B;

/** Circular velocity at the virial radius in the previous snapshot, 
    calculated from the virial radius and mass read from the postprocessing 
    files (see SubGroupVvir_A). */
extern float *SubGroupVvir_B;

/** Subhalo virial radius at the previous snapshot, calculated from the
    value of the virial mass. */
extern float  *SubGroupRvir_B;

/** Subhalo virial mass at the previous snapshot, determined from the 
    postprocessing files as the number of bound particles times the
    particle mass. */
extern float *SubGroupMvir_B;

extern float **GroupL_A;
extern float **GroupL_B;
extern float **GroupPos_A, **GroupVel_A, **GroupCenter_A;
extern float **GroupPos_B, **GroupVel_B, **GroupCenter_B;

/** Position of the subhalo at the current snapshot. */
extern float **SubGroupPos_A;

/** Velocity of the subhalo at the current snapshot. */
extern float  **SubGroupVel_A;

/** Position of the subhalo at the previous snapshot. */
extern float **SubGroupPos_B;

/** Velocity of the subhalo at the previous snapshot. */
extern float **SubGroupVel_B;

extern double Dt, Dtout, Tmax, Tdyn;
extern double E_in;
extern double GPMass;
extern double Hubble;
extern double LightSpeed;
extern double PartMass, PartMassGas, PartMassDM;
extern double RhoCrit;
extern double Rperi, Rt_peri;
extern double UnitTime_in_s, UnitTime_in_Megayears, UnitMass_in_Msun;
extern double UnitLength_in_kpc, UnitVelocity_in_km_per_s;
extern double UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
extern double UnitCoolingRate_in_cgs;
extern double UnitLength_orbit, UnitVelocity_orbit, UnitTime_orbit;
extern double UnitMass_orbit;

/** Numerical ID for ease of identification of the chosen tidal stripping
    timescale. */
extern unsigned int TSTimescaleSelect;

/** Numerical ID for ease of identification of the chosen approximation
    for the tidal radius. */
extern unsigned int RtidalSelect;

struct GALAXY
{
  /** Index of parent subhalo. */
  int ParentSubGroup;

  /** Index of parent FOF halo. */
  int ParentGroup;

  /** Index of most bound DM particle in the base simulation. */
  int PaIndex;

  /** Number of simulation particles bound to the subhalo. */
  int Len;

  /** Galaxy type: 0 FOF central, 1 central galaxy of a subhalo which is 
      not the main structure in the FOF, 2 satellite that had a subhalo in
      the past, but now is lost */
  int Type;

  int Id;

  int MergerType;

  /** Flag to identify satellites of type 1 galaxies which have been 
      tidally stripped and now orbit the FOF central */
  int Relocated;

  /** Index of progenitor galaxy in the previous step. A value of 0 
      indicates no progenitor, i.e. galaxy appears in the current 
      snapshot. */
  int Progenitor;

  /** Index of descendant galaxy. A value of 0 indicates no descendant,
      either because the galaxy has merged or because current snapshot
      is the final one. */
  int Descendant;

  float Mvir;

  float Rvir;

  float Vvir;

  float Vc;

  /** Current DM bound mass. For a type 0 or 1 galaxy, this should be
      equal to Mvir. For type 2 galaxies, this is equal to OrbitMass at
      the final STEP. */
  float Mdm;

  /** Current tidal radius. This is the value determined from the chosen
      approximation to the tidal radius, and if the TS timescale is longer
      than the time elapsed between snapshots this may be smaller than 
      Rdm. */
  float Rtidal;
 
  /** Current outer radius for the DM. This should be initially equal to 
      the virial radius at infall. */
  float Rdm; 

  /** Coulomb logarithm at infall. */
  float Coulomb;
       
  /** Merger time from the default analytical estimation. */
  float MergTime;
  
  /** Merger time from the analytical estimation by Boylan-Kolchin et al.
      (2008). */
  float MergTimeBK;

  /** Merger time from the analytical estimation by Jiang et al. (2008). */
  float MergTimeJ;

  /** Actual time elapsed from infall until the merger, or until the present
     if no merger happened. */
  float MergTimeActual;

  /** Total time spent as satellite. */
  float MergTimeSat;

  /** Kinetic energy of the galaxy. */
  float Ekin;

  /** Potential energy of the galaxy. */
  float Epot;

  /** Position in simulation coordinates (comoving). */
  float Pos[3];          

  /** Velocity in simulation coordinates (comoving). */
  float Vel[3];

  /** Position relative to central galaxy (physical coordinates). */
  float Posrel[3];  

  /** Velocity relative to central galaxy (physical coordinates). */
  float Vrel[3];

  /** Orbital angular momentum relative to central galaxy (physical 
      coordinates). */
  float Jorb[3];

  /** Initial value of the orbital angular momentum (squared) - used to 
      compare with the current @f$J^2@f$ to determine when a merger 
      happens. */
  float J2init;

  /** Energy lost by the galaxy via dynamical friction. */
  float EnergyLoss;

  /** Last known pericentric radius of the satellite. */
  float Rpericentre;

  /** Last known apocentric radius of the satellite. */
  float Rapocentre;

  /** Stores the X coordinate of the satellite's orbit. */
  float OrbitX[STEPS];

  /** Stores the Y coordinate of the satellite's orbit. */
  float OrbitY[STEPS];

  /** Stores the Z coordinate of the satellite's orbit. */
  float OrbitZ[STEPS];

  /** Stores the X coordinate of the satellite's orbital velocity. */
  float OrbitVx[STEPS];

  /** Stores the Y coordinate of the satellite's orbital velocity. */
  float OrbitVy[STEPS];

  /** Stores the Z coordinate of the satellite's orbital velocity. */
  float OrbitVz[STEPS];

  /** Stores the satellite mass along the orbit. */
  float OrbitMass[STEPS];

  /** Stores the satellite's type along the orbit. */
  int OrbitType[STEPS];

  /** Total mass lost to tidal stripping. */
  float StrippedMass;

  /** Stores the type history of the galaxy. */
  int TypeT[OUTPUTS];
  
  /** Stores the virial mass history of the galaxy. */
  float MvirT[OUTPUTS];
  
  /** Stores the virial radius history of the galaxy. */
  float RvirT[OUTPUTS];

  /** Store the initial angular momentum, for control */
  float J2initT[OUTPUTS];

  /** Stores the evolution of orbital angular momentum. */
  float JorbT[OUTPUTS];

  /** Stores the evolution of DM bound mass. */
  float MdmT[OUTPUTS];

  /** Stores the evolution of DM bounding radius. */
  float RdmT[OUTPUTS];

  /** Stores the halocentric radius. */
  float FOFRadiusT[OUTPUTS];

  /** Stores the subhalocentric radius. This is always relative to the 
      original central galaxy, even for relocated galaxies. */
  float SubRadiusT[OUTPUTS];

  /** Stores the history of galaxy relocation. */
  int RelocatedT[OUTPUTS];

  /** Stores the galaxy index. */
  int IndexT[OUTPUTS];
};

struct particle_data
{
  float Pos[3];
  float Vel[3];
  int Nextnode;
};

/** Stores the particle data from the simulation. */
extern struct particle_data *P;

/** Stores the galaxy population at the current snapshot. */
extern struct GALAXY *GalaxyA;

/** Stores the galaxy population at the previous snapshot. */
extern struct GALAXY *GalaxyB;

/** Flag to signal the code to stop because an error happened. */
extern unsigned int EmergencyStop;

#ifdef DEBUG
/** Used to determine whether it is the first time that orbit information
    is printed to the standard output. If it is, print a table header. */
extern unsigned int Nprints;
#endif
