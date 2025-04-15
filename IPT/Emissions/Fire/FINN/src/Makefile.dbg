.SUFFIXES: .o .f90

F90    = $(FC)
C_C    = $(CC)
LIBS   = -L$(NETCDF_DIR)/lib $(AR_FILES)
INCLUDE_MODULES = -I$(NETCDF_DIR)/include

ifeq ($(C_C),pgcc)
  CCFLAGS += -g -O2 -Ktrap=fp
endif

ifeq ($(C_C),gcc)
# CCFLAGS += -O -std=c99
  CCFLAGS += -O2 -std=c99 -g
endif

ifeq ($(FC),lf95)
  FFLAGS = --g
else
  FFLAGS = -g
endif

ifeq ($(FC),pgf90)
  FFLAGS += -O0 -Mnosave -Ktrap=fp
endif

ifeq ($(FC),gfortran)
  FFLAGS += -O0 -ffree-line-length-none
endif

ifeq ($(FC),lf95)
  FFLAGS += --nap --chk --pca --nsav --trace --trap --wide
endif

ifeq ($(FC),ifort)
  FFLAGS += -fpe0 -ftrapuv
endif

EXEC = fire_emis 

OBJS = attr_types.o \
       fire_types.o\
       netcdf_utils.o\
       hires_timer.o\
       module_timing.o\
       camse_utils.o\
       PolyUtils.o\
       fire_utils.o\
       wrf_utils.o\
       glb_utils.o\
       fire_srf_types.o \
       fire_file.o \
       fire_emis.o

##dependencies
#$(OBJECTS) : makefile

${EXEC} :       ${OBJS}
		${F90} -o $@ ${OBJS} ${LIBS} 

.f90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

.c.o:
		${C_C} ${CCFLAGS} -c $<

cleanup:
		rm -f ${OBJS} *.mod

clean:
		rm -f core ${EXEC} ${OBJS} *.mod
