/**
 * @file dumpHDF5_orbits.c
 * @brief Output orbits data in HDF5 format. 
 * @warning Requires HDF5 version 1.8 or greater. 
 */
#include <hdf5.h>
#include <hdf5_hl.h>
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


void dumpHDF5_popA(char *fname)
{
  int i, index, p;
  int *ibuffer, *ibufferT, *ibufferN;
  float buf;
  float *buffer, *bufferT, *buffer3d, *bufferN;
  double dbuf;
  hid_t file_id, subgroup_id;
  herr_t status;
  hsize_t dims[2]={TotNumGalA,1};
  hsize_t dim3d[2]={TotNumGalA,3};
  hsize_t dimT[2]={TotNumGalA,OUTPUTS};
  hsize_t dimN[2]={TotNumGalA,STEPS};
  char label[256], desc[256], path[256];
  

  /* Create a new HDF5 output file using default properties */
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /*
   * Write file attributes:
   * Store the parameters used and the code and compilation options
   */
  status = dumpHDF5_file_attr(file_id);

  ibuffer  = malloc(TotNumGalA*sizeof(int));
  ibufferT = malloc(TotNumGalA*OUTPUTS*sizeof(int));
  ibufferN = malloc(TotNumGalA*STEPS*sizeof(int));
  buffer   = malloc(TotNumGalA*sizeof(float));
  buffer3d = malloc(3*TotNumGalA*sizeof(float));
  bufferT  = malloc(TotNumGalA*OUTPUTS*sizeof(float));
  bufferN  = malloc(TotNumGalA*STEPS*sizeof(float));
 
  /*
   * Save data to file
   */
  sprintf(label, "/ParentGroup");
  sprintf(desc, "Index of parent FOF group");
  for (p=1; p<=TotNumGalA; p++) {
    ibuffer[p-1] = GalaxyA[p].ParentGroup;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, 
			    ibuffer);

  sprintf(label, "/ParentSubGroup");
  sprintf(desc, "Index of parent subgroup");
  for (p=1; p<=TotNumGalA; p++) {
    ibuffer[p-1] = GalaxyA[p].ParentSubGroup;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, 
			    ibuffer);

  sprintf(label, "/Type");
  sprintf(desc, "Galaxy type");
  for (p=1; p<=TotNumGalA; p++) {
    ibuffer[p-1] = GalaxyA[p].Type;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, 
			    ibuffer);

  sprintf(label, "/PaIndex");
  sprintf(desc, "Index of most bound DM particle");
  for (p=1; p<=TotNumGalA; p++) {
    ibuffer[p-1] = GalaxyA[p].PaIndex;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, 
			    ibuffer);

  /*sprintf(label, "/Len");
  sprintf(desc, "Number of DM particles in subgroup");
  for (p=1; p<=TotNumGalA; p++) {
    ibuffer[p-1] = GalaxyA[p].Len;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, 
  ibuffer);*/

  sprintf(label, "/CentralGal");
  sprintf(desc, "Index of central galaxy");
  for (p=1; p<=TotNumGalA; p++) {
    if (GalaxyA[p].Type == 0) ibuffer[p-1] = p;
    if (GalaxyA[p].Type == 1)
      ibuffer[p-1] = FirstGalInFOFGroup_A[GalaxyA[p].ParentGroup];
    if (GalaxyA[p].Type >= 2)
      ibuffer[p-1] = FirstGalInSubGroup_A[GalaxyA[p].ParentSubGroup];
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, 
			    ibuffer);

  sprintf(label, "/Relocated");
  sprintf(desc, "Flag to signal that the galaxy has been stripped from "
	  "its parent subgroup and now orbits the FOF central");
  for (p=1; p<=TotNumGalA; p++)
    ibuffer[p-1] = GalaxyA[p].Relocated;
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, 
			    ibuffer);

  sprintf(label, "/Mvir");
  sprintf(desc, "Virial mass of parent subhalo. For satellites it is "
	  "the value at infall");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].Mvir;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT,
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitMass_in_g, 
				   "Unit of mass given in h^-1 g.");

  sprintf(label, "/Rvir");
  sprintf(desc, "Virial radius of parent subhalo. For satellites it is "
	  "the value at infall");	  
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].Rvir;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units. Unit length "
				   "given in h^-1 cm.");

  sprintf(label, "/Vc");
  sprintf(desc, "Circular velocity (of halo?)");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].Vc;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitVelocity_in_cm_per_s, 
				   "Unit velocity given in cm s^-1.");

  sprintf(label, "/Rtidal");
  sprintf(desc, "Current tidal radius");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].Rtidal;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units. Unit length "
				   "given in h^-1 cm.");

  sprintf(label, "/Rdm");
  sprintf(desc, "Current DM bounding radius");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].Rdm;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units. Unit length "
				   "given in h^-1 cm.");

  sprintf(label, "/Mdm");
  sprintf(desc, "Current DM bound mass");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].Mdm;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units. Unit mass "
				   "given in h^-1 g.");

  sprintf(label, "/J2init");
  sprintf(desc, "Initial orbital angular momentum (squared)");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].J2init;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  dbuf = pow(UnitLength_in_cm*UnitVelocity_in_cm_per_s,2);
  status = set_float_variable_attr(file_id, label, desc, &dbuf,
				   "In physical units. Unit angular "
				   "momentum (squared) given in "
				   "h^-2 cm^4 s^-2.");

  sprintf(label, "/EnergyLoss");
  sprintf(desc, "Total energy lost by the subhalo via dynamical friction");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].EnergyLoss;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc, &UnitEnergy_in_cgs,
				   "In physical units. Unit energy "
				   "given in h^-1 erg. Only stores the"
				   "energy loss computed along the orbit "
				   "(i.e. the energy lost while type 1 is "
				   "not recorded).");

  sprintf(label, "/Coulomb");
  sprintf(desc, "Coulomb logarithm ln(1 + M/M_host) at infall");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].Coulomb;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  dbuf = 0.0;
  status = set_float_variable_attr(file_id, label, desc, &dbuf,
				   "Adimensional quantity.");

  sprintf(label, "/MergTime");
  sprintf(desc, "Time to merger from Chandrasekhar's formula");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].MergTime;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc, 
				   &UnitTime_in_s,
				   "Unit time given in h^-1 sec.");

  sprintf(label, "/MergTimeBK");
  sprintf(desc, "Time to merger from Boylan-Kolchin et al. (2008)");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].MergTimeBK;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc, 
				   &UnitTime_in_s,
				   "Unit time given in h^-1 sec.");

  sprintf(label, "/MergTimeJ");
  sprintf(desc, "Time to merger from Jiang et al. (2009)");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].MergTimeJ;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc, 
				   &UnitTime_in_s,
				   "Unit time given in h^-1 sec.");
  
  sprintf(label, "/MergTimeActual");
  sprintf(desc, "Actual time elapsed since infall");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].MergTimeActual;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			    buffer);
  status = set_float_variable_attr(file_id, label, desc, 
				   &UnitTime_in_s,
				   "Unit time given in sec.");

  sprintf(label, "/StrippedMass");
  sprintf(desc, "Total mass lost to tidal stripping");
  for (p=1; p<=TotNumGalA; p++) {
    buffer[p-1] = GalaxyA[p].StrippedMass;
  }
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT,
			    buffer);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitMass_in_g, 
				   "Unit of mass given in h^-1 g.");

  sprintf(label, "/Pos");
  sprintf(desc, "Position of galaxy");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<3; i++) {
      index = indexrm(p-1, i, TotNumGalA, 3) + 3; 
      buffer3d[index] = GalaxyA[p].Pos[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dim3d, H5T_NATIVE_FLOAT, 
			    buffer3d);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In comoving units. Unit length "
				   "given in h^-1 cm.");
  
  sprintf(label, "/Vel");
  sprintf(desc, "Velocity of galaxy");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<3; i++) {
      index = indexrm(p-1, i, TotNumGalA, 3) + 3;
      buffer3d[index] = GalaxyA[p].Vel[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dim3d, H5T_NATIVE_FLOAT, 
			    buffer3d);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitVelocity_in_cm_per_s,
				   "In funny units (look up "
				   "conversion). Unit velocity "
				   "given in cm s^-1.");

  sprintf(label, "/Posrel");
  sprintf(desc, "Position of galaxy relative to its central");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<3; i++) {
      index = indexrm(p-1, i, TotNumGalA, 3) + 3; 
      buffer3d[index] = GalaxyA[p].Posrel[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dim3d, H5T_NATIVE_FLOAT, 
			    buffer3d);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units. Unit length "
				   "given in h^-1 cm.");
  
  sprintf(label, "/Vrel");
  sprintf(desc, "Velocity of galaxy relative to its central");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<3; i++) {
      index = indexrm(p-1, i, TotNumGalA, 3) + 3;
      buffer3d[index] = GalaxyA[p].Vel[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dim3d, H5T_NATIVE_FLOAT, 
			    buffer3d);
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitVelocity_in_cm_per_s,
				   "In physical units. Unit velocity "
				   "given in cm s^-1.");
  
  sprintf(label, "/Jorb");
  sprintf(desc, "Orbital angular momentum of galaxy relative to its "
	  "central");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<3; i++) {
      index = indexrm(p-1, i, TotNumGalA, 3) + 3;
      buffer3d[index] = GalaxyA[p].Jorb[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dim3d, H5T_NATIVE_FLOAT, 
			    buffer3d);
  dbuf = UnitLength_in_cm*UnitVelocity_in_cm_per_s;
  status = set_float_variable_attr(file_id, label, desc, &dbuf,
				   "In physical units. Unit angular "
				   "momentum given in h^-1 cm^2 s^-1.");
  
  /* Create subgroup named "/Histories" */
  sprintf(path, "/Histories");
  subgroup_id = H5Gcreate(file_id, path, H5P_DEFAULT, H5P_DEFAULT,
			  H5P_DEFAULT);
  
  sprintf(label, "%s/TypeT", path);
  sprintf(desc, "History of galaxy type");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<OUTPUTS; i++) {
      index = indexrm(p-1, i, TotNumGalA, OUTPUTS) + OUTPUTS; 
      ibufferT[index] = GalaxyA[p].TypeT[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimT, H5T_NATIVE_INT, 
			    ibufferT);

  sprintf(label, "%s/RelocatedT", path);
  sprintf(desc, "History of galaxy relocation");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<OUTPUTS; i++) {
      index = indexrm(p-1, i, TotNumGalA, OUTPUTS) + OUTPUTS; 
      ibufferT[index] = GalaxyA[p].RelocatedT[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimT, H5T_NATIVE_INT, 
			    ibufferT);

  sprintf(label, "%s/JorbT", path);
  sprintf(desc, "Angular momentum history");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<OUTPUTS; i++) {
      index = indexrm(p-1, i, TotNumGalA, OUTPUTS) + OUTPUTS; 
      bufferT[index] = GalaxyA[p].JorbT[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimT, H5T_NATIVE_FLOAT, 
			    bufferT);
  dbuf = UnitLength_in_cm*UnitVelocity_in_cm_per_s;
  status = set_float_variable_attr(file_id, label, desc, &dbuf,
				   "In physical units. Unit angular "
				   "momentum given in h^-1 cm^2 s^-1.");

  sprintf(label, "%s/J2initT", path);
  sprintf(desc, "Initial angular momentum history");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<OUTPUTS; i++) {
      index = indexrm(p-1, i, TotNumGalA, OUTPUTS) + OUTPUTS; 
      bufferT[index] = GalaxyA[p].J2initT[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimT, H5T_NATIVE_FLOAT, 
			    bufferT);
  dbuf = UnitLength_in_cm*UnitVelocity_in_cm_per_s;
  status = set_float_variable_attr(file_id, label, desc, &dbuf,
				   "In physical units. Unit angular "
				   "momentum given in h^-1 cm^2 s^-1.");

  sprintf(label, "%s/RdmT", path);
  sprintf(desc, "DM bounding radius history");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<OUTPUTS; i++) {
      index = indexrm(p-1, i, TotNumGalA, OUTPUTS) + OUTPUTS; 
      bufferT[index] = GalaxyA[p].RdmT[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimT, H5T_NATIVE_FLOAT, 
			    bufferT);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units. Unit length "
				   "given in h^-1 cm.");

  sprintf(label, "%s/MdmT", path);
  sprintf(desc, "DM bound mass history");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<OUTPUTS; i++) {
      index = indexrm(p-1, i, TotNumGalA, OUTPUTS) + OUTPUTS; 
      bufferT[index] = GalaxyA[p].MdmT[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimT, H5T_NATIVE_FLOAT, 
			    bufferT);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitMass_in_g,
				   "In physical units. Unit mass "
				   "given in h^-1 g.");

  sprintf(label, "%s/RvirT", path);
  sprintf(desc, "Galaxy virial radius history");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<OUTPUTS; i++) {
      index = indexrm(p-1, i, TotNumGalA, OUTPUTS) + OUTPUTS; 
      bufferT[index] = GalaxyA[p].RvirT[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimT, H5T_NATIVE_FLOAT, 
			    bufferT);	
  status = set_float_variable_attr(file_id, label, desc, &UnitLength_in_cm,
				   "In physical units. Unit length "
				   "given in h^-1 cm.");

  sprintf(label, "%s/MvirT", path);
  sprintf(desc, "Galaxy virial mass history");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<OUTPUTS; i++) {
      index = indexrm(p-1, i, TotNumGalA, OUTPUTS) + OUTPUTS; 
      bufferT[index] = GalaxyA[p].MvirT[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimT, H5T_NATIVE_FLOAT, 
			    bufferT);	
  status = set_float_variable_attr(file_id, label, desc, &UnitMass_in_g,
				   "In physical units. Unit mass "
				   "given in h^-1 g.");
  
  sprintf(label, "%s/FOFRadiusT", path);
  sprintf(desc, "Halocentric distance history");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<OUTPUTS; i++) {
      index = indexrm(p-1, i, TotNumGalA, OUTPUTS) + OUTPUTS; 
      bufferT[index] = GalaxyA[p].FOFRadiusT[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimT, H5T_NATIVE_FLOAT, 
			    bufferT);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units. Unit length "
				   "given in h^-1 cm.");

  sprintf(label, "%s/SubRadiusT", path);
  sprintf(desc, "Subhalocentric distance history");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<OUTPUTS; i++) {
      index = indexrm(p-1, i, TotNumGalA, OUTPUTS) + OUTPUTS; 
      bufferT[index] = GalaxyA[p].SubRadiusT[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimT, H5T_NATIVE_FLOAT, 
			    bufferT);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units. Unit length "
				   "given in h^-1 cm.");
  
  /* Close subgroup Histories */
  status = H5Gclose(subgroup_id);

  /* Create subgroup named "/Orbits" */
  sprintf(path, "/Orbits");
  subgroup_id = H5Gcreate(file_id, path, H5P_DEFAULT, H5P_DEFAULT,
			  H5P_DEFAULT);

  sprintf(label, "%s/OrbitX", path);
  sprintf(desc, " ");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<STEPS; i++) {
      index = indexrm(p-1, i, TotNumGalA, STEPS) + STEPS; 
      bufferN[index] = GalaxyA[p].OrbitX[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimN, H5T_NATIVE_FLOAT, 
			    bufferN);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units, relative to the "
				   "central galaxy. Unit length "
				   "given in h^-1 cm.");

  sprintf(label, "%s/OrbitY", path);
  sprintf(desc, " ");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<STEPS; i++) {
      index = indexrm(p-1, i, TotNumGalA, STEPS) + STEPS; 
      bufferN[index] = GalaxyA[p].OrbitY[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimN, H5T_NATIVE_FLOAT, 
			    bufferN);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units, relative to the "
				   "central galaxy. Unit length "
				   "given in h^-1 cm.");

  sprintf(label, "%s/OrbitZ", path);
  sprintf(desc, " ");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<STEPS; i++) {
      index = indexrm(p-1, i, TotNumGalA, STEPS) + STEPS; 
      bufferN[index] = GalaxyA[p].OrbitZ[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimN, H5T_NATIVE_FLOAT, 
			    bufferN);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitLength_in_cm,
				   "In physical units, relative to the "
				   "central galaxy. Unit length "
				   "given in h^-1 cm.");

  sprintf(label, "%s/OrbitVx", path);
  sprintf(desc, " ");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<STEPS; i++) {
      index = indexrm(p-1, i, TotNumGalA, STEPS) + STEPS; 
      bufferN[index] = GalaxyA[p].OrbitVx[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimN, H5T_NATIVE_FLOAT, 
			    bufferN);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitVelocity_in_cm_per_s,
				   "In physical units, relative to the "
				   "central galaxy. Unit velocity "
				   "given in cm/s.");

  sprintf(label, "%s/OrbitVy", path);
  sprintf(desc, " ");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<STEPS; i++) {
      index = indexrm(p-1, i, TotNumGalA, STEPS) + STEPS; 
      bufferN[index] = GalaxyA[p].OrbitVy[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimN, H5T_NATIVE_FLOAT, 
			    bufferN);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitVelocity_in_cm_per_s,
				   "In physical units, relative to the "
				   "central galaxy. Unit velocity "
				   "given in cm/s.");

  sprintf(label, "%s/OrbitVz", path);
  sprintf(desc, " ");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<STEPS; i++) {
      index = indexrm(p-1, i, TotNumGalA, STEPS) + STEPS; 
      bufferN[index] = GalaxyA[p].OrbitVz[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimN, H5T_NATIVE_FLOAT, 
			    bufferN);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitVelocity_in_cm_per_s,
				   "In physical units, relative to the "
				   "central galaxy. Unit velocity "
				   "given in cm/s.");

  sprintf(label, "%s/OrbitMass", path);
  sprintf(desc, "Subhalo DM mass along the orbit");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<STEPS; i++) {
      index = indexrm(p-1, i, TotNumGalA, STEPS) + STEPS; 
      bufferN[index] = GalaxyA[p].OrbitMass[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimN, H5T_NATIVE_FLOAT, 
			    bufferN);	
  status = set_float_variable_attr(file_id, label, desc,
				   &UnitMass_in_g,
				   "Unit mass given in h^-1 g.");
  
  sprintf(label, "%s/OrbitType", path);
  sprintf(desc, "Galaxy type along the orbit: 2 orbiting satellite, "
	  "3 merged, 4 disrupted");
  for (p=1; p<=TotNumGalA; p++) {
    for (i=0; i<STEPS; i++) {
      index = indexrm(p-1, i, TotNumGalA, STEPS) + STEPS; 
      ibufferN[index] = GalaxyA[p].OrbitType[i];
    }
  }
  status = H5LTmake_dataset(file_id, label, 2, dimN, H5T_NATIVE_INT, 
			    ibufferN);	
  
  sprintf(label, "%s/Rpericentre", path);
  sprintf(desc, "Last known pericentric radius of the galaxy");
  for (p=1; p<=TotNumGalA; p++)
     buffer[p-1] = GalaxyA[p].Rpericentre;
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT,
                            buffer);
  status = set_float_variable_attr(file_id, label, desc,
                                   &UnitLength_in_cm,
                                   "Unit length given in h^-1 cm.");
  
  sprintf(label, "%s/Rapocentre", path);
  sprintf(desc, "Last known apocentric radius of the galaxy");
  for (p=1; p<=TotNumGalA; p++) 
     buffer[p-1] = GalaxyA[p].Rapocentre;
  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT,
                            buffer);
  status = set_float_variable_attr(file_id, label, desc,
                                   &UnitLength_in_cm,
                                   "Unit length given in h^-1 cm.");

  /* Close subgroup Orbits */
  status = H5Gclose(subgroup_id);

  /* Close HDF5 file */
  status = H5Fclose(file_id);

  free(buffer);
  free(bufferT);
  free(buffer3d);
  free(ibuffer);
  free(ibufferT);
  free(ibufferN);
  
  return;
}
