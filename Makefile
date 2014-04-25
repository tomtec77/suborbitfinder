####### Select target computer
#SYSTYPE = 'FCAGLP-horus'
#SYSTYPE = 'FCAGLP-seminare'
SYSTYPE = 'PUC-geryon'
#SYSTYPE = 'i686'

####### Select simulation to use
# DOLAGSIM = 'ON' corresponds in general to cluster resimulations, due
# to differences in the headers of the output files.
DOLAGSIM = 'ON'

####### Print extra information for debugging
#DEBUGMODE = 'ON'

####### Path to the HDF5 library
# For some of the systems, the location of the HDF5 library is known
# (horus, geryon). In those systems the actual path is specified below,
# when setting the optimizations.
# If you need to specify a different path, do it here and then comment
# all other instances of HDF5PATH.
#HDF5PATH = /home/ttecce/libs/hdf5-1.8.9-linux-x86_64-static/
#HDF5PATH = /usr/local/hdf5-1.8.8-linux-static/
#HDF5PATH = /home/tomas/work/hdf5-1.8.8-linux-static/
HDF5PATH = /usr/local/lib/hdf5-1.8.10-linux-static/

##########################################################################
# Code configuration complete. In routine use, it should not be necessary 
# to modify anything beyond this point.
##########################################################################

####### Code compilation options
# These are activated by the selections above.
OPT =
ifeq ($(DOLAGSIM),'ON')
  OPT   = -DDOLAGSIMULATION
endif
ifeq ($(DEBUGMODE),'ON')
  OPT  += -DDEBUG
endif

####### Files
#EXEC = suborbitfinder
#EXEC = t-subof
EXEC = g-subof

OBJS = suborbitfinder.o readparameterfile.o allvars.o io.o age.o init.o \
       io_input_G2.o population.o integrate_orbit.o tidalradius.o \
       indexrm.o dumpHDF5_orbits.o dumpHDF5_param.o dumpHDF5_aux.o \
       potential_sis.o galaxydistance.o \
       lib/nrsrc/nrutil.o lib/nrsrc/dqromb.o lib/nrsrc/dtrapzd.o \
       lib/nrsrc/dpolint.o

INCL = proto.h allvars.h readparameterfile.h \
       lib/nrsrc/nrutil.h lib/nrsrc/nrsag.h

####### Compiler settings
.KEEP_STATE:

OPTIMIZE = -O2

LIBS = -lm -L$(HDF5PATH)/lib

ifeq ($(SYSTYPE),'i686')
  CC = $(HDF5PATH)/bin/h5cc
  CFLAGS = $(OPTIMIZE) $(OPT) -I$(HDF5PATH)/include -mtune=native -Wno-unused-result
endif

ifeq ($(SYSTYPE),'PUC-geryon')
  HDF5PATH = /usr/local/hdf5-1.8.7-linux-x86_64-static/
  CC = $(HDF5PATH)/bin/h5cc
  CFLAGS = $(OPTIMIZE) $(OPT) -I$(HDF5PATH)/include -march=nocona -pipe -I/usr/local/include
endif

ifeq ($(SYSTYPE),'FCAGLP-seminare')
  CC = $(HDF5PATH)/bin/h5cc
  CFLAGS = $(OPTIMIZE) $(OPT) -I$(HDF5PATH)/include -Wall -pedantic -std=c99 -D_BSD_SOURCE -Wno-unused-variable -march=native -I/opt/apps/gsl/gsl-1.15/include
  LIBS     += -L/opt/apps/gsl/gsl-1.15/lib
endif

ifeq ($(SYSTYPE),'FCAGLP-horus')
  HDF5PATH = /usr/local/hdf5-1.8.9-linux-x86_64-static/
  CC = $(HDF5PATH)/bin/h5cc
  CFLAGS = $(OPTIMIZE) $(OPT) -I$(HDF5PATH)/include -mtune=native
endif


$(EXEC): $(OBJS) Makefile
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LIBS)

$(OBJS): $(INCL) Makefile

clean:
	rm -f $(OBJS) $(EXEC) *~ lib/nrsrc/*~

tidy:
	rm -f $(OBJS) *~ lib/nrsrc/*~
