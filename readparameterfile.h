/**
 * @file readparameterfile.h
 * @brief This is (almost) the same file from SAG, to use here the same 
 * format of parameters file. Currently the parameter EpsGrav is not 
 * present in SAG.
 */
extern char     path1[256];
extern char     path[256];
extern char     path2[256];
extern char     path3[256];
extern char     filebasename[256];
extern int      Files;
extern int 	MaxSnapshot;
extern double 	Omega; 		
extern double   OmegaLambda;
extern double 	Hubble_h;
extern double	H0;
extern double 	G;  		
extern int    	GroupMinLen;
extern double   BaryonFrac;
extern double	EnergySN;
extern double 	EtaSN;
extern double	EnergySNgamma;
extern double   BinFrac;
extern double	VcFactor;
extern double	Alpha;
extern double 	Epsilon;
extern double   ThreshMerger;
extern double   ThreshGasBurst;
extern double   LowerMerger;
extern double   ThreshDisk;
extern double   PertDist;
extern double   FracEj;
extern float    FracH;
extern float    FracHe;
extern float    FracBH;
extern float    EtaBH;
extern float    K_AGN;
extern float    SNIa_TimeToExplode; 
extern float    SNIa_RateCte; 
extern float    SpinMeanLambda;
extern float    SpinSigmaLambda;
extern float    GasRscale;
extern float    AlphaRP;
extern float    MaxCvirRedshift;
extern float    MajorWetParameter;

/** Gravitational softening to use in the host halo potential for orbit
    calculations. */
extern float    EpsGrav;

/** IMF to use in SAG. Options: Salpeter, Kroupa, Chabrier, WKTopheavy */
extern char     IMFModel[16];

/** Approximation of the tidal radius to use. Options are: ZB03 (use the 
    approximation of King 1962 as used by Zentner & Bullock 2003), TDS98
    (use the variant in Tormen, Diaferio & Syer 1998). */
extern char RtidalModel[8];

/** Select which timescale to use for tidal stripping. Options are:
    subhalo (use the current subhalo dynamical time), subinfall (use the
    subhalo dynamical time at infall). */
extern char TSTimescaleModel[16];

extern char	identifier[256];
extern double	TimeBetSnapshot;
extern double	TimeOfFirstSnapshot;
extern double	UnitLength_in_cm;
extern double	UnitMass_in_g;
extern double	UnitVelocity_in_cm_per_s;
extern char     NameOutputsSelection[256];
extern char     NameCatalogue[256];
extern char     NameContamination[256];
extern char     NameHistory[256];
extern char     NameProperties[256];
extern char     NameSnapshot[256];
extern char     NameSubgalaxies[256];
extern char     NameSubgalfuture[256];
extern char     NameSubgalhistory[256];
extern char     NameSubhistory[256];
extern char     NameSubids[256];
extern char     NameSublist[256];
extern char     NameSubproperties[256];
extern char     NameSubstructures[256];
extern char     NameVolatile[256];
extern char     NameDMordered[256];
extern char     NameConcentration[256];
extern char     NameRscaleTables[256];
extern char     NameVcTables[256];
extern char     NameICMProperties[256];
extern char     NameSuborbits[256];

/*
 * Code options
 */
extern unsigned int EjectionOn;
extern unsigned int ReIncorporateOn;
extern unsigned int SubGalaxyTreeOn;
extern unsigned int DustOn;
extern unsigned int IRAOn;
extern unsigned int AlphaConstantOn;
extern unsigned int FormatGadgetHDM5On;
extern unsigned int DMPartOn;
extern unsigned int AGNOn;
extern unsigned int AGNFeedbackInQSOModeOn;
extern unsigned int DiskInstabilitiesOn;
extern unsigned int GradualInstabilitiesOn;
extern unsigned int EddingtonOn;
extern unsigned int GrasilOn;
extern unsigned int DumpNoPositions;
extern unsigned int DumpNoVelocities;
extern unsigned int OutputListOn;
extern unsigned int RamPressureOn;
extern unsigned int RamPressureHaloOn;
extern unsigned int IncreaseYieldsOn;
extern unsigned int JiangTimeFrictionOn;
extern unsigned int BoylanTimeFrictionOn;
extern unsigned int ComputeBHSpinOn;
extern unsigned int RandomBHSpinOn;
extern unsigned int FollowingJGal;
extern unsigned int SpecialDumpOn;
extern unsigned int UseInclinationRPOn;
extern unsigned int EllipticalViaMinorMergerOn;
extern unsigned int RPFitOn;
extern unsigned int CrotonSFOn;
extern unsigned int CoolingFeHOn;
extern unsigned int DMPosAvailableOn;
extern unsigned int StripColdGasToStarsOn;
extern unsigned int BC2003OldTablesOn;
extern unsigned int QeffectiveOn;
extern unsigned int UseRscaleTablesOn;

/* New implementation of IMF selection 1/12/11 
   Modified the selection method 2012-01-02 
   IMFSelect = 0 for Salpeter, 1 for Kroupa, 2 for Chabrier, 3 for
   WKTopheavy */
#define IMFSALPETER 0
#define IMFKROUPA 1
#define IMFCHABRIER 2
#define IMFTOPHEAVY 3
extern unsigned int IMFSelect;

/** Enable dynamical friction on the orbiting satellite. */
extern unsigned int DynamicalFrictionOn;

/** Enable tidal stripping of the satellite's mass. */
extern unsigned int TidalStrippingOn;

/** Enable mergers by proximity to the central galaxy. */
extern unsigned int ProximityMergerOn;

/** Disable relocation of type 2 satellites of type 1. */
extern unsigned int SatelliteRelocationDisabledOn;
