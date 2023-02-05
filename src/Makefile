# first delete all default suffixes
.SUFFIXES:
# then add your suffixes
.SUFFIXES: .f90 .F90 .o

include ./make.sys

default: all

modu1 = lattice.o model_para.o dqmc_ft_basic_data.o
modu2 = dqmc_config_u.o dqmc_config_plqu.o dqmc_config_h0.o dqmc_measure.o
modu3 = dqmc_ft_core.o dqmc_ctrl.o
modu4 = dqmc_proj_core.o dqmc_proj_ctrl.o
modu5 = dqmc_ft_core.o dqmc_ctrl3d.o
modu6 = dqmc_proj_core.o dqmc_proj_ctrl3d.o
#dqmcf1 = sli.o salph.o npbc.o seth.o thop_mag.o \
#         sthop.o dqmc_ft_initial.o inconfc.o outconfc.o
#dqmcf2 = mmthl.o mmthlm1.o mmthr.o mmthrm1.o mmuul.o mmuulm1.o mmuur.o \
#         mmuurm1.o upgradeuf_delay.o upgradejk.o upgradejh.o preq.o \
#         prtau.o \
#         set_phase.o
main = dqmc_main.o

objects = $(modu1) $(modu2) $(modu3) $(dqmcf1) $(dqmcf2) $(main)
objects2 = $(modu1) $(modu2) $(modu4) $(dqmcf1) $(dqmcf2) $(main)
objects3 = $(modu1) $(modu2) $(modu5) $(dqmcf1) $(dqmcf2) $(main)
objects4 = $(modu1) $(modu2) $(modu6) $(dqmcf1) $(dqmcf2) $(main)

#LIBS+= ../lib/libMM.a

all: cpmod dqmc_ft dqmc_proj dqmc_ft3d dqmc_proj3d
ft: cpmod dqmc_ft
proj: cpmod dqmc_proj
ft3d: cpmod dqmc_ft3d
proj3d: cpmod dqmc_proj3d

dqmc_ft: $(objects)
	$(LINKER) $(objects) -o dqmc_ft $(LFLAGS) $(LIBS)

dqmc_proj: $(objects2)
	$(LINKER) $(objects2) -o dqmc_proj $(LFLAGS) $(LIBS)

dqmc_ft3d: $(objects3)
	$(LINKER) $(objects3) -o dqmc_ft3d $(LFLAGS) $(LIBS)

dqmc_proj3d: $(objects4)
	$(LINKER) $(objects4) -o dqmc_proj3d $(LFLAGS) $(LIBS)

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
	rm -f dqmc_ft3d
	rm -f dqmc_proj3d

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
