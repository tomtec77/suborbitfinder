/**
 * @file dumpHDF5_aux.c
 * @brief Auxiliary functions for HDF5 output.
 */
#include <hdf5.h>
#include <hdf5_hl.h>
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


/**
 * @brief Set attributes of float-type output variables. 
 * @param file_id HDF5 file ID
 * @param label Label of the output dataset to which the attribute will be
 * attached
 * @param desc Description of the attribute
 * @param unit The corresponding unit of measure for the attribute
 * @param comment Additional commentary or notes
 *
 * For the selected variable, store a description, the value of the 
 * corresponding unit in the code, and a comment. 
 */
int set_float_variable_attr(hid_t file_id, char *label, char *desc, 
			    double *unit, char *comment)
{
  herr_t status;

  status = H5LTset_attribute_string(file_id, label, "Description", desc);
  if (status < 0) {
    fprintf(stderr, "Error (HDF5): cannot create attribute 'Description'"
	    " - Exit\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  status = H5LTset_attribute_double(file_id, label, "Unit", unit, 1);
  if (status < 0) {
    fprintf(stderr, "Error (HDF5): cannot create attribute 'Unit'"
	    " - Exit\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  status = H5LTset_attribute_string(file_id, label, "Comment", comment);
  if (status < 0) {
    fprintf(stderr, "Error (HDF5): cannot create attribute 'Comment'"
	    " - Exit\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}


/**
 * @brief Check for errors in creating HDF5 attributes for data.
 * @param status HDF5 error status variable
 */
void checkwrite_HDF5(herr_t status)
{
  if (status < 0) {
    fprintf(stderr, "Error (HDF5): cannot create attribute - Exit|\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  return;
}
