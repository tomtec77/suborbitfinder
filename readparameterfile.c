/**
 * @file readparameterfile.c
 * @brief This is (almost) the same file from SAG, to use the same format of 
 * parameters file, plus modifications necessary for @a suborbitfinder
 * @author Tomas E. Tecce
 */
#include "allvars.h"
#include "proto.h"
#include "readparameterfile.h"

char    path1[256];
char    path[256];
char    path2[256];
char    path3[256];     /* Path of files with ordered gas particles */
char	filebasename[256];
int	Files;
int 	MaxSnapshot;
double 	Omega; 		/* Used together with particle mass (from header 
			   of snapshot file) */
double  OmegaLambda;
double 	Hubble_h;	/* Hubble constant in units of 100 km/s/Mpc */
double	H0;		/* Hubble constant in code units (see gadget output)*/
double 	G;  		/* G to compute interparticle separation */
int    	GroupMinLen; 	/* Store only groups in the catalogue 
                           with at least this number of particles */
double 	BaryonFrac;
double	EnergySN;
double 	EtaSN;
double	EnergySNgamma;
double  BinFrac;
double	VcFactor;
double	Alpha;
double 	Epsilon;
double	ThreshMerger;
double  ThreshGasBurst;
double  LowerMerger;
double  ThreshDisk;
double  PertDist;
double  FracEj;
float   FracH;
float   FracHe;
float   FracBH;
float   EtaBH;
float   K_AGN;
float   SNIa_TimeToExplode;
float   SNIa_RateCte;
float   SpinMeanLambda;      /* For distribution of halo spin parameter */
float   SpinSigmaLambda; 
float   GasRscale;
float   AlphaRP;
float   MaxCvirRedshift;
float   MajorWetParameter;   /* To decide if a major merger is wet or dry */
float   EpsGrav;
char    IMFModel[16];        
char    RtidalModel[8];
char    TSTimescaleModel[16];
char 	identifier[256];
double	TimeBetSnapshot;
double	TimeOfFirstSnapshot;
double	UnitLength_in_cm;
double	UnitMass_in_g;
double	UnitVelocity_in_cm_per_s;

char    NameOutputsSelection[256];
char    NameCatalogue[256];
char    NameContamination[256];
char    NameHistory[256];
char    NameProperties[256];
char    NameSnapshot[256];
char    NameSubgalaxies[256];
char    NameSubgalfuture[256];
char    NameSubgalhistory[256];
char    NameSubhistory[256];
char    NameSubids[256];
char    NameSublist[256];
char    NameSubproperties[256];
char    NameSubstructures[256];
char    NameVolatile[256];
char    NameDMordered[256];
char    NameConcentration[256];
char    NameRscaleTables[256];
char    NameVcTables[256];
char    NameICMProperties[256];
char    NameSuborbits[256];

/*
 * Code options
 */
unsigned int EjectionOn;
unsigned int ReIncorporateOn;
unsigned int SubGalaxyTreeOn;
unsigned int DustOn;
unsigned int IRAOn;
unsigned int AlphaConstantOn;
unsigned int FormatGadgetHDM5On;
unsigned int DMPartOn;
unsigned int AGNOn;
unsigned int AGNFeedbackInQSOModeOn;
unsigned int DiskInstabilitiesOn;
unsigned int GradualInstabilitiesOn;
unsigned int EddingtonOn;
unsigned int GrasilOn;
unsigned int DumpNoPositions;
unsigned int DumpNoVelocities;
unsigned int OutputListOn;
unsigned int RamPressureOn;
unsigned int RamPressureHaloOn;
unsigned int IncreaseYieldsOn;
unsigned int JiangTimeFrictionOn;
unsigned int BoylanTimeFrictionOn;
unsigned int ComputeBHSpinOn;
unsigned int RandomBHSpinOn;
unsigned int FollowingJGal;
unsigned int SpecialDumpOn;
unsigned int UseInclinationRPOn;
unsigned int EllipticalViaMinorMergerOn;
unsigned int RPFitOn;
unsigned int CrotonSFOn;
unsigned int CoolingFeHOn;
unsigned int DMPosAvailableOn;
unsigned int StripColdGasToStarsOn;
unsigned int IGIMFTopheavyOn;
unsigned int BC2003OldTablesOn;
unsigned int QeffectiveOn;
unsigned int UseRscaleTablesOn;

unsigned int IMFSelect;

unsigned int DynamicalFrictionOn;
unsigned int TidalStrippingOn;
unsigned int ProximityMergerOn;
unsigned int SatelliteRelocationDisabledOn;

/*=======================================================================*/
void readparameterfile(char *filename)
/*=======================================================================
  This routine reads the parameterfile that is valid for ALL postprocessing 
  programs. The parameters are declared in readparameterfile.h

  The codes (especially 'galaxies') work only for the Choice of standard 
  units.
   
  Seed                          Seed for the random number generator.
                                Use any POSITIVE integer between 1 and 100.
  Path 				Path of the directory containing a Snapshot 
                                subdirectory with the snapshot files of one 
                                simulation.
  FileBasename			Snapshotnames without the 3-digit number.
  MaxSnapshot			Snapshot number for z=0.
  Omega				Omega Matter
  OmegaLambda			Omega Lambda
  Hubble_h			Hubble constant in units of 100 km/s/Mpc
  H0				Hubble constant in code units of Gadget, as 
                                written out during the simulation. 0.1 for 
                                standard choice of units.
  G				Gravitational constant in Gadget code units. 
                                As above.
  GroupMinLen   		Minimum particle number for the FOF 
                                identification.
  BaryonFrac			Original baryon fraction.
  EnergySN			Energy released by one supernova explosion. 
                                (In erg)
  EtaSN				Relative number of SN per amount of stars.
  Yield				Mass fraction of metals produced during 
                                star formation.
  VcFactor			Assumption: Vmax = Vcircular * VcFactor
  Alpha				Star formation efficiency.
  Epsilon			Feedback efficiency.	
  ThreshMerger			Threshold for the mass ratio when a merging 
                                event will be treated as "major merger".
  ThreshGasBurst		Threshold for the Cold gas-disc ratio which 
                                will defined if would be star burst in a minor 
                                merger
  LowerMerger			Limit of mass ratio used to determined if 
                                would be a star burst in a minor merger
  Identifier			Text string written into the final galaxy 
                                output filename.
  TimeBetSnapshot		As entered in Gadget.
  TimeOfFirstSnapshot 		As entered in Gadget.
  UnitLength_in_cm 		As entered in Gadget. (3.085678e21 for 1.0 
                                kpc)
  UnitMass_in_g       		As entered in Gadget. (1.989e43 for 1.0e10 
                                solar masses)
  UnitVelocity_in_cm_per_s 	As entered in Gadget. (1e5 for 1 km/sec)
  =======================================================================*/
{
  FILE *fp;
  int elements, check;
  char *error;
  char buffer[256], buffer2[256];
#ifndef HDF5OUTPUT
  FILE *fw;
  char saveparam_fname[FILENAME_MAX];
#endif
  void checkforerror(char *error, int elements, char *buf);

  
  fp = fopen(filename, "r");
  if (fp==NULL) {
    fprintf(stdout, 
	    "Error (readparameterfile): file '%s' not found - Stop.\n",
	    filename); fflush(stdout);
    exit(EXIT_FAILURE);     
  }
  fprintf(stdout,"Reading from parameters file %s...\n", filename);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Seed%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  Seed = atoi(buffer2);
  printf("  Seed            %d\n", Seed);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Path1%s",path1);
  checkforerror(error,elements,&buffer[0]);
  printf("  Path1           %s\n", path1);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Path%s",path);
  checkforerror(error,elements,&buffer[0]);
  printf("  Path            %s\n", path);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Path2%s",path2);
  checkforerror(error,elements,&buffer[0]);
  printf("  Path2           %s\n", path2);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Path3%s",path3);
  checkforerror(error,elements,&buffer[0]);
  printf("  Path3           %s\n", path3);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"FileBasename%s",filebasename);
  checkforerror(error,elements,&buffer[0]);
  printf("  FileBasename    %s\n", filebasename);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Files%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  Files = atoi(buffer2);
  printf("  Files           %d\n", Files);
  
  printf("\n");

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"MaxSnapshot%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  MaxSnapshot = atoi(buffer2);
  printf("  MaxSnapshot                  %d\n", MaxSnapshot);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Omega%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  Omega = atof(buffer2);
  printf("  Omega                        %g\n", Omega);
  
  error = fgets(buffer,sizeof(buffer),fp);  
  elements = sscanf(buffer,"OmegaLambda%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  OmegaLambda = atof(buffer2);
  printf("  OmegaLambda                  %g\n", OmegaLambda);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Hubble_h%s",buffer2); 
  checkforerror(error,elements,&buffer[0]);    
  Hubble_h = atof(buffer2);
  printf("  Hubble_h                     %g\n", Hubble_h);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"H0%s",buffer2); 
  checkforerror(error,elements,&buffer[0]);    
  H0 = atof(buffer2);
  printf("  H0                           %g\n", H0);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"G%s",buffer2); 
  checkforerror(error,elements,&buffer[0]);    
  G = atof(buffer2);
  printf("  G                            %g\n", G);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"GroupMinLen%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  GroupMinLen = atoi(buffer2);
  printf("  GroupMinLen                  %d\n", GroupMinLen);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"BaryonFrac%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  BaryonFrac = atof(buffer2);
  printf("  BaryonFrac                   %g\n", BaryonFrac);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"EnergySN%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  EnergySN = atof(buffer2);
  printf("  EnergySN                     %g\n", EnergySN);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"EtaSN%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  EtaSN = atof(buffer2);
  printf("  EtaSN                        %g\n", EtaSN);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"EnergySNgamma%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  EnergySNgamma = atof(buffer2);
  printf("  EnergySNgamma                %g\n", EnergySNgamma);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"BinFrac%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  BinFrac = atof(buffer2);
  printf("  BinFrac                      %g\n", BinFrac);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"VcFactor%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  VcFactor = atof(buffer2);
  printf("  VcFactor                     %g\n", VcFactor);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Alpha%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  Alpha = atof(buffer2);
  printf("  Alpha                        %g\n", Alpha);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Epsilon%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  Epsilon = atof(buffer2);
  printf("  Epsilon                      %g\n", Epsilon);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"ThreshMerger%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  ThreshMerger = atof(buffer2);
  printf("  ThreshMerger                 %g\n", ThreshMerger);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"ThreshGasBurst%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  ThreshGasBurst = atof(buffer2);
  printf("  ThreshGasBurst               %g\n", ThreshGasBurst);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"LowerMerger%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  LowerMerger = atof(buffer2);
  printf("  LowerMerger                  %g\n", LowerMerger);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"ThreshDisk%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  ThreshDisk = atof(buffer2);
  printf("  ThreshDisk                   %g\n", ThreshDisk);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"PertDist%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  PertDist = atof(buffer2);
  printf("  PertDist                     %g\n", PertDist);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"FracEj%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  FracEj = atof(buffer2);
  printf("  FracEj                       %g\n", FracEj);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"FracH%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  FracH = atof(buffer2);
  printf("  FracH                        %g\n", FracH);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"FracHe%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  FracHe = atof(buffer2);
  printf("  FracHe                       %g\n", FracHe);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"FracBH%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  FracBH = atof(buffer2); 
  printf("  FracBH                       %g\n", FracBH);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"EtaBH%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  EtaBH = atof(buffer2);
  printf("  EtaBH                        %g\n", EtaBH);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"K_AGN%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  K_AGN = atof(buffer2);
  printf("  K_AGN                        %g\n", K_AGN);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"SNIa_TimeToExplode%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  SNIa_TimeToExplode = atof(buffer2);
  printf("  SNIa_TimeToExplode           %g\n", SNIa_TimeToExplode);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"SNIa_RateCte%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  SNIa_RateCte = atof(buffer2);
  printf("  SNIa_RateCte                 %g\n", SNIa_RateCte);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"SpinMeanLambda%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  SpinMeanLambda = atof(buffer2);
  printf("  SpinMeanLambda               %g\n", SpinMeanLambda);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"SpinSigmaLambda%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  SpinSigmaLambda = fabsf(log(atof(buffer2)));
  printf("  SpinSigmaLambda              %g\n", atof(buffer2));
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"GasRscale%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  GasRscale = atof(buffer2);
  printf("  GasRscale                    %g\n", GasRscale);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"AlphaRP%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  AlphaRP = atof(buffer2);
  printf("  AlphaRP                      %g\n", AlphaRP);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"MaxCvirRedshift%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  MaxCvirRedshift = atof(buffer2);
  printf("  MaxCvirRedshift              %g\n", MaxCvirRedshift);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"MajorWetParameter%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  MajorWetParameter = atof(buffer2);
  printf("  MajorWetParameter            %g\n", MajorWetParameter);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"EpsGrav%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  EpsGrav = atof(buffer2);
  printf("  EpsGrav                      %g\n", EpsGrav);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"IMFModel%s",IMFModel);
  checkforerror(error,elements,&buffer[0]);
  printf("  IMFModel                     %s\n", IMFModel);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"RtidalModel%s",RtidalModel);
  checkforerror(error,elements,&buffer[0]);
  printf("  RtidalModel                  %s\n", RtidalModel);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"TSTimescaleModel%s",TSTimescaleModel);
  checkforerror(error,elements,&buffer[0]);
  printf("  TSTimescaleModel             %s\n", TSTimescaleModel);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Identifier%s",identifier);
  checkforerror(error,elements,&buffer[0]);
  printf("  Identifier                   %s\n", identifier);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"TimeBetSnapshot%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  TimeBetSnapshot = atof(buffer2);
  printf("  TimeBetSnapshot              %g\n", TimeBetSnapshot);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"TimeOfFirstSnapshot%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  TimeOfFirstSnapshot = atof(buffer2);
  printf("  TimeOfFirstSnapshot          %g\n", TimeOfFirstSnapshot);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UnitLength_in_cm%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  UnitLength_in_cm = atof(buffer2);
  printf("  UnitLength_in_cm             %g\n", UnitLength_in_cm);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UnitMass_in_g%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  UnitMass_in_g = atof(buffer2);
  printf("  UnitMass_in_g                %g\n", UnitMass_in_g);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UnitVelocity_in_cm_per_s%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  UnitVelocity_in_cm_per_s = atof(buffer2);
  printf("  UnitVelocity_in_cm_per_s     %g\n", UnitVelocity_in_cm_per_s);
  
  printf("\n");
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameOutputsSelection%s",NameOutputsSelection);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameOutputsSelection    %s\n", NameOutputsSelection);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameCatalogue%s",NameCatalogue);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameCatalogue           %s\n", NameCatalogue);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameContamination%s",NameContamination);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameContamination       %s\n", NameContamination);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameHistory%s",NameHistory);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameHistory             %s\n", NameHistory);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameProperties%s",NameProperties);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameProperties          %s\n", NameProperties);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSnapshot%s",NameSnapshot);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameSnapshot            %s\n", NameSnapshot);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSubgalaxies%s",NameSubgalaxies);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameSubgalaxies         %s\n", NameSubgalaxies);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSubgalfuture%s",NameSubgalfuture);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameSubgalfuture        %s\n", NameSubgalfuture);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSubgalhistory%s",NameSubgalhistory);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameSubgalhistory       %s\n", NameSubgalhistory);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSubhistory%s",NameSubhistory);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameSubhistory          %s\n", NameSubhistory);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSubids%s",NameSubids);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameSubids              %s\n", NameSubids);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSublist%s",NameSublist);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameSublist             %s\n", NameSublist);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSuborbits%s",NameSuborbits);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameSuborbits           %s\n", NameSuborbits);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSubproperties%s",NameSubproperties);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameSubproperties       %s\n", NameSubproperties);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameSubstructures%s",NameSubstructures);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameSubstructures       %s\n", NameSubstructures);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameVolatile%s",NameVolatile);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameVolatile            %s\n", NameVolatile);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameDMordered%s",NameDMordered);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameDMordered           %s\n", NameDMordered);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameConcentration%s",NameConcentration);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameConcentration       %s\n", NameConcentration);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameRscaleTables%s",NameRscaleTables);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameRscaleTables        %s\n", NameRscaleTables);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameVcTables%s",NameVcTables);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameVcTables            %s\n", NameVcTables);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameICMProperties%s",NameICMProperties);
  checkforerror(error,elements,&buffer[0]);
  printf("  NameICMProperties       %s\n", NameICMProperties);
  
  /*
   * Code options
   */
  printf("\n");
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"EjectionOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  EjectionOn = atoi(buffer2);
  printf("  EjectionOn                    %d\n", EjectionOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"ReIncorporateOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  ReIncorporateOn = atoi(buffer2);
  printf("  ReIncorporateOn               %d\n", ReIncorporateOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"SubGalaxyTreeOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  SubGalaxyTreeOn = atoi(buffer2);
  printf("  SubGalaxyTreeOn               %d\n", SubGalaxyTreeOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"DustOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  DustOn = atoi(buffer2);
  printf("  DustOn                        %d\n", DustOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"IRAOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  IRAOn = atoi(buffer2);
  printf("  IRAOn                         %d\n", IRAOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"AlphaConstantOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  AlphaConstantOn = atoi(buffer2);
  printf("  AlphaConstantOn               %d\n", AlphaConstantOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"FormatGadgetHDM5On%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  FormatGadgetHDM5On = atoi(buffer2);
  printf("  FormatGadgetHDM5On            %d\n", FormatGadgetHDM5On);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"DMPartOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  DMPartOn = atoi(buffer2);
  printf("  DMPartOn                      %d\n", DMPartOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"AGNOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  AGNOn = atoi(buffer2);
  printf("  AGNOn                         %d\n", AGNOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"AGNFeedbackInQSOModeOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  AGNFeedbackInQSOModeOn = atoi(buffer2);
  printf("  AGNFeedbackInQSOModeOn        %d\n", AGNFeedbackInQSOModeOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"DiskInstabilitiesOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  DiskInstabilitiesOn = atoi(buffer2);
  printf("  DiskInstabilitiesOn           %d\n", DiskInstabilitiesOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"GradualInstabilitiesOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  GradualInstabilitiesOn = atoi(buffer2);
  printf("  GradualInstabilitiesOn        %d\n", GradualInstabilitiesOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"EddingtonOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  EddingtonOn = atoi(buffer2);
  printf("  EddingtonOn                   %d\n", EddingtonOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"GrasilOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  GrasilOn = atoi(buffer2);
  printf("  GrasilOn                      %d\n", GrasilOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"DumpNoPositions%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  DumpNoPositions = atoi(buffer2);
  printf("  DumpNoPositions               %d\n", DumpNoPositions);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"DumpNoVelocities%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  DumpNoVelocities = atoi(buffer2);
  printf("  DumpNoVelocities              %d\n", DumpNoVelocities);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"OutputListOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  OutputListOn = atoi(buffer2);
  printf("  OutputListOn                  %d\n", OutputListOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"RamPressureOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  RamPressureOn = atoi(buffer2);
  printf("  RamPressureOn                 %d\n", RamPressureOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"RamPressureHaloOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  RamPressureHaloOn = atoi(buffer2);
  printf("  RamPressureHaloOn             %d\n", RamPressureHaloOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"IncreaseYieldsOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  IncreaseYieldsOn = atoi(buffer2);
  printf("  IncreaseYieldsOn              %d\n", IncreaseYieldsOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"JiangTimeFrictionOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  JiangTimeFrictionOn = atoi(buffer2);
  printf("  JiangTimeFrictionOn           %d\n", JiangTimeFrictionOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"BoylanTimeFrictionOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  BoylanTimeFrictionOn = atoi(buffer2);
  printf("  BoylanTimeFrictionOn          %d\n", BoylanTimeFrictionOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"ComputeBHSpinOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  ComputeBHSpinOn = atoi(buffer2);
  printf("  ComputeBHSpinOn               %d\n", ComputeBHSpinOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"RandomBHSpinOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  RandomBHSpinOn = atof(buffer2);
  printf("  RandomBHSpinOn                %d\n", RandomBHSpinOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"FollowingJGal%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  FollowingJGal = atof(buffer2);
  printf("  FollowingJGal                 %d\n", FollowingJGal);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"SpecialDumpOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  SpecialDumpOn = atoi(buffer2);
  printf("  SpecialDumpOn                 %d\n", SpecialDumpOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UseInclinationRPOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  UseInclinationRPOn = atoi(buffer2);
  printf("  UseInclinationRPOn            %d\n", UseInclinationRPOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"EllipticalViaMinorMergerOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  EllipticalViaMinorMergerOn = atoi(buffer2);
  printf("  EllipticalViaMinorMergerOn    %d\n", EllipticalViaMinorMergerOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"RPFitOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  RPFitOn = atoi(buffer2);
  printf("  RPFitOn                       %d\n", RPFitOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"CrotonSFOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  CrotonSFOn = atoi(buffer2);
  printf("  CrotonSFOn                    %d\n", CrotonSFOn);
   
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"CoolingFeHOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  CoolingFeHOn = atoi(buffer2);
  printf("  CoolingFeHOn                  %d\n", CoolingFeHOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"DMPosAvailableOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  DMPosAvailableOn = atof(buffer2);
  printf("  DMPosAvailableOn              %d\n", DMPosAvailableOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"StripColdGasToStarsOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  StripColdGasToStarsOn = atoi(buffer2);
  printf("  StripColdGasToStarsOn         %d\n", StripColdGasToStarsOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"BC2003OldTablesOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  BC2003OldTablesOn = atoi(buffer2);
  printf("  BC2003OldTablesOn             %d\n", BC2003OldTablesOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"QeffectiveOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  QeffectiveOn = atoi(buffer2);
  printf("  QeffectiveOn                  %d\n", QeffectiveOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UseRscaleTablesOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  UseRscaleTablesOn = atoi(buffer2);
  printf("  UseRscaleTablesOn             %d\n", UseRscaleTablesOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"DynamicalFrictionOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  DynamicalFrictionOn = atoi(buffer2);
  printf("  DynamicalFrictionOn           %d\n", DynamicalFrictionOn);
  
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"TidalStrippingOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  TidalStrippingOn = atoi(buffer2);
  printf("  TidalStrippingOn              %d\n", TidalStrippingOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"ProximityMergerOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  ProximityMergerOn = atoi(buffer2);
  printf("  ProximityMergerOn             %d\n", ProximityMergerOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"SatelliteRelocationDisabledOn%s",buffer2);
  checkforerror(error,elements,&buffer[0]);
  SatelliteRelocationDisabledOn = atoi(buffer2);
  printf("  SatelliteRelocationDisabledOn %d\n", 
	 SatelliteRelocationDisabledOn);

  printf("\n");
  fclose(fp);

  /*
   * Safety checks
   */
  if (SpecialDumpOn == 1 && RamPressureOn == 1) {
    printf("WARNING: RamPressureOn = 1 but special dump "
	   "has been selected\n");
    printf("Check values in parameters file and try again - Exit\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  IMFSelect = 10;
  if (strcmp(IMFModel, "Salpeter") == 0) IMFSelect = IMFSALPETER;
  if (strcmp(IMFModel, "Kroupa") == 0) IMFSelect = IMFKROUPA;
  if (strcmp(IMFModel, "Chabrier") == 0) IMFSelect = IMFCHABRIER;
  if (strcmp(IMFModel, "WKTopheavy") == 0) IMFSelect = IMFTOPHEAVY;
  if (IMFSelect >= 10) {
    fprintf(stderr, 
	    "Error (readparameterfile): invalid IMF model '%s' selected\n",
	    IMFModel);
    fprintf(stderr, "Check parameters file and try again - Exit\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  RtidalSelect = 10;
  if (strcmp(RtidalModel, "TDS98") == 0) RtidalSelect = RTIDAL_TDS98;
  if (strcmp(RtidalModel, "ZB03") == 0) RtidalSelect = RTIDAL_ZB03; 
  if (RtidalSelect >= 10) {
    fprintf(stderr, 
	    "Error (readparameterfile): invalid tidal radius approximation"
	    " '%s' selected\n", RtidalModel);
    fprintf(stderr, "Check parameters file and try again - Exit\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  TSTimescaleSelect = 10;
  if (strcmp(TSTimescaleModel, "subhalo") == 0) 
    TSTimescaleSelect = TSTIME_SUBHALO;
  if (strcmp(TSTimescaleModel, "subinfall") == 0) 
    TSTimescaleSelect = TSTIME_SUBINFALL; 
  if (TSTimescaleSelect >= 10) {
    fprintf(stderr, 
	    "Error (readparameterfile): invalid TS timescale option "
	    "'%s' selected\n", TSTimescaleModel);
    fprintf(stderr, "Check parameters file and try again - Exit\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
 
  /* 
   * TOMAS 2013-05-17:
   * Satellite relocation, with the current subhalo scheme, does not work.
   * Some type 1 galaxies are later ejected from their host halo, so if any
   * of its satellites are relocated they will no longer see the correct 
   * FOF central. In other cases, the subhalo goes from type 1 to 0 because
   * the original main subhalo disappears in the current snapshot. There
   * are also cases where the type 0 and the type 1 trade places. The 
   * simple scheme cannot handle these cases, especially in the case where
   * the original main subhalo disappears. Therefore we suspend the 
   * relocation of galaxies until we can find a better treatment for the
   * stripping of satellites.
   */
  if (SatelliteRelocationDisabledOn != 1) {
    fprintf(stderr, "Error (readparameterfile): in this version of the "
	    "code, satellite relocation must be disabled.\n"
	    "If you are sure you need to use that, edit file "
	    "readparameterfile.c to remove this control.\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  return;
}

 
void checkforerror(char *error, int elements, char *buf)
{
  char buffer2[256]; 

  
  if (error==NULL) {
    fprintf(stderr, "Error: couldn't read from parameter file. Stop.\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  if (elements==0) {
    strncpy(buffer2,buf,strlen(buf));
    /* strncpy does not add a final \0 terminator. It has to be added by hand*/
    buffer2[strlen(buf)-1]='\0';
    fprintf(stderr, "Error: couldn't convert '%s'. Stop.\n", buffer2);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
}


