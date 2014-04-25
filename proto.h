/**
 * @file proto.h
 * @brief Function prototypes for use with suborbitfinder
 */
#include <hdf5.h>
#include <hdf5_hl.h>

double accel(double, double, double, double, double, double, double);
void   allocate_memory();
void   checkerror(int, int);
void   checkwrite_HDF5(herr_t);
herr_t dumpHDF5_file_attr(hid_t);
void   dumpHDF5_popA(char *);
void   free_concentration_A();
void   free_groups_A();
void   free_groups_B();
void   free_history();
void   free_linklist_A();
void   free_linklist_B();
void   free_memory();
void   free_properties_A();
void   free_properties_B();
void   free_subgroups_A();
void   free_subgroups_B();
void   free_subhistory();
void   free_subproperties_A();
void   free_subproperties_B();
void   galaxydistance(int, float *, float *);
void   generate_new_population();
float  get_tidal_radius(int, int);
int    getNumPart(char *, int, int);
int    indexrm(int, int, int, int);
void   init();
void   initialize_galaxy_population();
int    integrate_type2_orbit(int, int, double, double, double);
void   move_popA_to_popB();
void   output(double, double, double);
double potential_sis(double, double, double, double, double, double);
void   read_concentration_A(char *);
void   read_group_properties_A(char *, double);
void   read_group_properties_B(char *, double);
void   read_groups_A(char *);
void   read_groups_B(char *);
void   read_header(char *, int);
void   read_header_mass(char *, int);
void   read_history(char *);
void   read_linklist_A(char *);
void   read_linklist_B(char *);
void   read_subgroup_properties_A(char *, double, char *);
void   read_subgroup_properties_B(char *, double, char *);
void   read_subgroups_A(char *, char *, char *, char*);
void   read_subgroups_B(char *, char *, char *, char*);
void   read_subhistory(char *);
void   readparameterfile(char *);
int    set_float_variable_attr(hid_t, char *, char *, double *, char *);
void   set_units();
void   testprint_0(int, int, double);
void   testprint_1(int);
void   testprint_type1(int, double);
double time_to_present(double);
