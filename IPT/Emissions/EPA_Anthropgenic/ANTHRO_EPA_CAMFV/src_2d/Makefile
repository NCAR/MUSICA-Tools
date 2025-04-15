.SUFFIXES: .o .f90

F90    = $(FC)
C_C     = $(CC)
LIBS   = -L$(NETCDF_DIR)/lib $(AR_LIBS)
INCLUDE_MODULES = -I$(NETCDF_DIR)/include

ifeq ($(C_C),pgcc)
  CCFLAGS += -O
endif
ifeq ($(FC),pgf90)
# FFLAGS += -Mnosave -O -Mbounds
# FFLAGS += -Mnosave -O  -Ktrap=fp
  FFLAGS += -Mnosave -O0 -g
# FFLAGS += -Mnosave -O
endif

ifeq ($(FC),pgf95)
# FFLAGS += -Mnosave -O -Mbounds -Ktrap=fp
# FFLAGS += -Mnosave -O0 -Ktrap=fp
# FFLAGS += -Mnosave -O0
  FFLAGS += -Mnosave -O
endif

ifeq ($(C_C),gcc)
  CCFLAGS += -O -std=c99
endif
ifeq ($(FC),gfortran)
# FFLAGS += -ffree-line-length-none -O
  FFLAGS += -ffree-line-length-none -O -g -fbacktrace -fcheck=bounds,do,mem,pointer -ffpe-trap=invalid,zero,overflow
endif

ifeq ($(FC),lf95)
  FFLAGS += --nap --chk --pca --nsav --trace --trap --wide
endif

ifeq ($(C_C),icc)
  CCFLAGS += -O -std=99
endif
ifeq ($(FC),ifort)
  FFLAGS += -fpe0 -ftrapuv -O
endif

EXEC = anthro_emis 

OBJS = QuadXsectArea.v4.o misc_definitions_module.o constants_module.o mo_calendar.o anthro_types.o mapper_types.o netcdf_utils.o area_mapper.o anthro_utils.o epa.o stack.o data_file_utils.o anthro_emis.o

##dependencies
#$(OBJECTS) : makefile

${EXEC} :       ${OBJS}
		${F90} -o $@ ${OBJS} ${LIBS} 

.f90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

.c.o:
		${C_C} ${CCFLAGS} -c ${INCLUDE_MODULES} $<

cleanup:
		rm -f ${OBJS} *.mod

clean:
		rm -f core ${EXEC} ${OBJS} *.mod
