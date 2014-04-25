/**
 * @file dumpHDF5_param.c
 * @brief Output code parameters in HDF5 format. 
 * @author Tomas E. Tecce
 * @warning Requires HDF5 version 1.8 or greater. 
 *
 * This file contains a function which adds, as attributes to the root 
 * group in the HDF5 output file, the parameters used in the run.
 */
#include <hdf5.h>
#include <hdf5_hl.h>
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


/**
 * @brief Add the parameters used as attributes to the root group in the
 * output file.
 * @param file_id HDF5 file ID
 */
herr_t dumpHDF5_file_attr(hid_t file_id)
{
  herr_t status;
  char hostname[256];
  time_t currtime;
  struct tm *loctime;
  float buf;


  status = H5LTset_attribute_string(file_id, "/", "Format", "SAG-HDF5");
  checkwrite_HDF5(status);

  gethostname(hostname, 256);
  status = H5LTset_attribute_string(file_id, "/", "System", hostname);
  checkwrite_HDF5(status);
  
  currtime = time(NULL);
  loctime = localtime(&currtime);
  status = H5LTset_attribute_string(file_id, "/", "Created",
				    asctime(loctime));

  status = H5LTset_attribute_int(file_id, "/", "Snapshot", &Snapshot, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "Redshift", 
				   &ZZ[Snapshot], 1);
  checkwrite_HDF5(status);

  /*
   * Write the parameters used in the run
   */
  status = H5LTset_attribute_uint(file_id, "/", "Seed", &Seed, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "Path1", path1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "Path", path);
  checkwrite_HDF5(status);
  
  status = H5LTset_attribute_string(file_id, "/", "Path2", path2);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "Path3", path3);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "FileBasename", 
				    filebasename);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_int(file_id, "/", "Files", &Files, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_int(file_id, "/", "MaxSnapshot", 
				 &MaxSnapshot, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "Omega", &Omega, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "OmegaLambda", 
				    &OmegaLambda, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "Hubble_h", 
				    &Hubble_h, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "H0", &H0, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "G", &G, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_int(file_id, "/", "GroupMinLen", 
				 &GroupMinLen, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "BaryonFrac", 
				    &BaryonFrac, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "EnergySN", 
				    &EnergySN, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "EtaSN", &EtaSN, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "EnergySNgamma", 
				    &EnergySNgamma, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "VcFactor", 
				    &VcFactor, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "Alpha", &Alpha, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "Epsilon", &Epsilon, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "ThreshMerger", 
				    &ThreshMerger, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "ThreshGasBurst", 
				    &ThreshGasBurst, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "LowerMerger", 
				    &LowerMerger, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "EtaSN", &EtaSN, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "ThreshDisk", 
				    &ThreshDisk, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "PertDist", &PertDist, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "FracEj", &FracEj, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "FracH", &FracH, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "FracHe", &FracHe, 1);
  checkwrite_HDF5(status);
  
  status = H5LTset_attribute_float(file_id, "/", "FracBH", &FracBH, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "EtaBH", &EtaBH, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "K_AGN", &K_AGN, 1);
  checkwrite_HDF5(status);

  /* Convert back from code units */
  buf = SNIa_TimeToExplode * UnitTime_in_Megayears / Hubble_h / 1.0e3;
  status = H5LTset_attribute_float(file_id, "/", "SNIa_TimeToExplode", 
				   &buf, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "SNIa_RateCte", 
				    &SNIa_RateCte, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "SpinMeanLambda", 
				   &SpinMeanLambda, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "SpinSigmaLambda", 
				   &SpinSigmaLambda, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "GasRscale", 
				   &GasRscale, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "AlphaRP", &AlphaRP, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "MaxCvirRedshift", 
				   &MaxCvirRedshift, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_float(file_id, "/", "MajorWetParameter", 
				   &MajorWetParameter, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "IMFModel", IMFModel);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_string(file_id, "/", "Identifier", 
				   identifier);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "TimeBetSnapshot", 
				    &TimeBetSnapshot, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "TimeOfFirstSnapshot", 
				    &TimeOfFirstSnapshot, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "UnitLength_in_cm", 
				    &UnitLength_in_cm, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", "UnitMass_in_g", 
				    &UnitMass_in_g, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_double(file_id, "/", 
				    "UnitVelocity_in_cm_per_s", 
				    &UnitVelocity_in_cm_per_s, 1);
  checkwrite_HDF5(status);


  /*
   * Code options
   */
  status = H5LTset_attribute_uint(file_id, "/", "EjectionOn", 
				 &EjectionOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "ReIncorporateOn", 
				 &ReIncorporateOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "SubGalaxyTreeOn", 
				 &SubGalaxyTreeOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "DustOn", &DustOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "IRAOn", &IRAOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "AlphaConstantOn", 
				 &AlphaConstantOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "FormatGadgetHDM5On", 
				 &FormatGadgetHDM5On, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "DMPartOn", &DMPartOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "AGNOn", &AGNOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "AGNFeedbackInQSOModeOn", 
				 &AGNFeedbackInQSOModeOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "DiskInstabilitiesOn", 
				 &DiskInstabilitiesOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "GradualInstabilitiesOn", 
				 &GradualInstabilitiesOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "EddingtonOn", 
				  &EddingtonOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "GrasilOn", &GrasilOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "DumpNoPositions", 
				 &DumpNoPositions, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "DumpNoVelocities", 
				 &DumpNoVelocities, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "OutputListOn", 
				 &OutputListOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "RamPressureOn", 
				  &RamPressureOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "RamPressureHaloOn", 
				  &RamPressureHaloOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "IncreaseYieldsOn", 
				 &IncreaseYieldsOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "JiangTimeFrictionOn", 
				 &JiangTimeFrictionOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "BoylanTimeFrictionOn", 
				 &BoylanTimeFrictionOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "ComputeBHSpinOn", 
				 &ComputeBHSpinOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "RandomBHSpinOn", 
				 &RandomBHSpinOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "FollowingJGal", 
				 &FollowingJGal, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "SpecialDumpOn", 
				 &SpecialDumpOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "UseInclinationRPOn", 
				  &UseInclinationRPOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", 
				  "EllipticalViaMinorMergerOn", 
				  &EllipticalViaMinorMergerOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "RPFitOn", &RPFitOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "CrotonSFOn", 
				 &CrotonSFOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "CoolingFeHOn", 
				 &CoolingFeHOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "DMPosAvailableOn", 
				  &DMPosAvailableOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "StripColdGasToStarsOn", 
				 &StripColdGasToStarsOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "BC2003OldTablesOn", 
				 &BC2003OldTablesOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "QeffectiveOn", 
				 &QeffectiveOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "UseRscaleTablesOn", 
				  &UseRscaleTablesOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "DynamicalFrictionOn", 
				  &DynamicalFrictionOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "TidalStrippingOn", 
				 &TidalStrippingOn, 1);
  checkwrite_HDF5(status);

  status = H5LTset_attribute_uint(file_id, "/", "ProximityMergerOn", 
				 &ProximityMergerOn, 1);
  checkwrite_HDF5(status);

#ifdef DOLAGSIMULATION
  status = H5LTset_attribute_string(file_id, "/", "DOLAGSIMULATION", 
				    "YES");
  checkwrite_HDF5(status);
#else
  status = H5LTset_attribute_string(file_id, "/", "DOLAGSIMULATION", "NO");
  checkwrite_HDF5(status);
#endif
  
  return status;
}
