.SUFFIXES: .o .f90

F90    = $(FC)
C_C     = $(CC)
LIBS   = -L$(NETCDF_DIR)/lib $(AR_LIBS)
INCLUDE_MODULES = -I$(NETCDF_DIR)/include

ifeq ($(C_C),pgcc)
  CCFLAGS += -g -O2 -Ktrap=fp
# CCFLAGS += -g -O0 -Ktrap=fp
endif
ifeq ($(FC),pgf90)
# FFLAGS += -Mnosave -O -Mbounds
  FFLAGS += -g -Mnosave -O2 -Ktrap=fp
  FFLAGS1 := $(FFLAGS)
# FFLAGS1 += -Mnosave -O0 -g
# FFLAGS += -Mnosave -O0 -g
# FFLAGS += -Mnosave -O
endif

ifeq ($(FC),pgf95)
# FFLAGS += -Mnosave -O -Mbounds -Ktrap=fp
# FFLAGS += -Mnosave -O0 -Ktrap=fp
# FFLAGS += -Mnosave -O0
  FFLAGS += -Mnosave -O
endif

ifeq ($(C_C),gcc)
# CCFLAGS += -O -std=c99
  CCFLAGS += -O2 -std=c99 -g
endif
ifeq ($(FC),gfortran)
# FFLAGS += -ffree-line-length-none -O -g
  FFLAGS += -ffree-line-length-none -O2 -g -fbacktrace -ffpe-trap=invalid,zero,overflow
  FFLAGS1 += -ffree-line-length-none -O2 -g -fbacktrace -ffpe-trap=invalid,zero,overflow
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

#OBJS = QuadXsectArea.v4.o misc_definitions_module.o constants_module.o mo_calendar.o mapper_types.o anthro_types.o netcdf_utils.o area_mapper.v2.o camse_utils.v1.o anthro_utils.o epa.o stack.v2.o data_file_utils.o anthro_emis.v1.o
OBJS = PolyXsectArea.o misc_definitions_module.o constants_module.o mo_calendar.o mapper_types.o anthro_types.o netcdf_utils.o area_mapper.o camse_utils.o anthro_utils.o epa.o stack.o data_file_utils.o anthro_emis.o

##dependencies
#$(OBJECTS) : makefile

${EXEC} :       ${OBJS}
		${F90} -o $@ ${OBJS} ${LIBS} 

.f90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

.c.o:
		${C_C} ${CCFLAGS} -c $<

anthro_emis.o: anthro_emis.f90
		${F90} ${FFLAGS1} -c ${INCLUDE_MODULES} $<

camse_utils.o: camse_utils.f90
		${F90} ${FFLAGS1} -c ${INCLUDE_MODULES} $<

cleanup:
		rm -f ${OBJS} *.mod

clean:
		rm -f core ${EXEC} ${OBJS} *.mod
