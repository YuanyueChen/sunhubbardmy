# first delete all default suffixes
.SUFFIXES:
# then add your suffixes
.SUFFIXES: .f90 .F90 .o

include ./make.sys

default: all

modu1 = lattice.o dqmc_basic_data.o model_para.o 
modu2 = dqmc_config_u.o dqmc_config_plqu.o dqmc_config_h0.o dqmc_config_v.o dqmc_para.o dqmc_util.o dqmc_measure.o
modu3 = dqmc_ft_core.o  dqmc_ft_ctrl.o
modu4 = dqmc_proj_core.o dqmc_proj_ctrl.o
main = dqmc_main.o

objects1 = $(modu1) $(modu2) $(modu3) $(main)
objects2 = $(modu1) $(modu2) $(modu4) $(main)

#LIBS+= ../lib/libMM.a

all: cpmod dqmc_ft dqmc_proj
ft: cpmod dqmc_ft
proj: cpmod dqmc_proj

dqmc_ft: $(objects1)
	$(LINKER) $(objects1) -o dqmc_ft $(LFLAGS) $(LIBS)

dqmc_proj: $(objects2)
	$(LINKER) $(objects2) -o dqmc_proj $(LFLAGS) $(LIBS)

cpmod:
	cp ../lib/*.mod .

.f90.o:
	$(F90) $(FFLAGS) $*.f90

.F90.o:
	$(F90) $(FFLAGS) $*.F90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f dqmc_ft
	rm -f dqmc_proj

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
