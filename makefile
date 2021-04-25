###################################################################################
#Define all the compilers and their options
FC= ifort
MPIFC= mpiifort
FFLAGS= -warn -g -O2 -i8 -r8 -w -no-wrap-margin -module ../build/modules
FLIBS_PAR=-L${MKLROOT}/lib/intel64 -mkl=sequential -lmkl_rt -lpthread -lm -ldl
FLIBS_SEQ=-L${MKLROOT}/lib/intel64 -mkl=sequential -lmkl_rt -lpthread -lm -ldl
INCLUDE= -I${MKLROOT}/include
BUILDDIR= ../build
OBJDIR= ../build/objects
MODDIR= ../build/modules
###################################################################################
#List the object files and their dependencies
#Define some variables with common files
COMFILES= $(OBJDIR)/blas.o $(OBJDIR)/linpack.o $(OBJDIR)/timer.o \
	$(OBJDIR)/lbfgsb.o $(OBJDIR)/nr_fft.o \
	$(OBJDIR)/instantonmod.o $(OBJDIR)/utils.o $(OBJDIR)/verletmodule.o

#1D double well potential:
1DFILES= $(OBJDIR)/mcmod_1d.o

#2D C_6 double well potential:
2DFILES= $(OBJDIR)/mcmod_2dtest.o

#SO(2) potential:
SO2FILES= $(OBJDIR)/mcmod_so2.o

#Malonaldehyde:
MALONFILES= $(OBJDIR)/pes_malonaldehyde.o $(OBJDIR)/mcmod_malon.o

#Malonaldehyde Wang:
WMALONFILES= $(OBJDIR)/mcmod_wmalon.o
pimd_wmalon_par: FFLAGS+= -O -I./mod_malon
pimd_wmalon_ser: FFLAGS+= -O -I./mod_malon
WMALONLIBS= -L. -lpes_malon

#Formic acid dimer:
FORMICFILES= $(OBJDIR)/mcmod_formic.o
pimd_formic_par: FFLAGS+= -O -I./mod_formic
pimd_formic_ser: FFLAGS+= -O -I./mod_formic
FORMICLIBS= -L. -lpes_formic

#Water dimer (MB-pol):
WATDIMFILES= $(OBJDIR)/mcmod_waterdimer.o
WATDIMLIBS= -L. -lmbpol -cxxlib

#Water trimer (MB-pol):
WATTRIMFILES= $(OBJDIR)/mcmod_watertrimer.o

#Water hexamer (MB-pol):
WATHEXFILES= $(OBJDIR)/mcmod_waterhexamer.o

#MOLPRO interface
MOLPROFILES= $(OBJDIR)/mcmod_molpro.o
rpi_molpro_par: FFLAGS+= -fpp -save-temps -DMOLPRO

#Water dimer (HBB2):
HBB2FILES= $(OBJDIR)/mcmod_waterdimer_hbb2.o
pimd_hbb2_par: FFLAGS+= -O -I./mod_water
pimd_hbb2_ser: FFLAGS+= -O -I./mod_water
rpi_hbb2_ser: FFLAGS+= -O -I./mod_water
HBB2LIBS= -L. -lpes_water -lpesd_water

#Water dimer (CCPol):
CCPOLFILES= $(OBJDIR)/H2O.pjt2.o $(OBJDIR)/main_CCpol-8sf.o $(OBJDIR)/proc_ccpol8s-dimer_xyz_ncd.o $(OBJDIR)/proc_sapt5sf_new_ncd.o $(OBJDIR)/ang_interface.o $(OBJDIR)/mcmod_waterdimer_ccpol.o


#Water-Methane Dimer
WATMETHFILES= $(OBJDIR)/watermethane.o $(OBJDIR)/mcmod_watmeth.o

#Methane Clathrate:
CLATHFILES= $(OBJDIR)/watermethane.o $(OBJDIR)/mcmod_clathrate.o

#Graphene:
GRAPHENEFILES= $(OBJDIR)/graphene.o $(OBJDIR)/mcmod_graphene.o

#Water dimer and halide:
WATHALFILES= $(OBJDIR)/mcmod_waterhalide.o
WATHALLIBS= -L. -lpot_halides -ltools_halides -cxxlib
pimd_wathal_par: FFLAGS+= -lstdc++
pimd_wathal_ser: FFLAGS+= -lstdc++ 

###################################################################################
#Compilation commands for object files
$(OBJDIR)/%.o: %.f
	$(FC) $(INCLUDE) -c $(FFLAGS) $< -o $@

$(OBJDIR)/%.o: %.f90
	$(FC)  $(INCLUDE) -c $(FFLAGS) $< -o $@

$(OBJDIR)/pimd_par.o: pimd_par.f90
	$(MPIFC) $(INCLUDE) -c $(FFLAGS) $< -o $@

$(OBJDIR)/rpi_par.o: rpi_par.f90
	$(MPIFC) $(INCLUDE) -c $(FFLAGS) $< -o $@

$(OBJDIR)/parallelmod.o: parallelmod.f90
	$(MPIFC) $(INCLUDE) -c $(FFLAGS) $< -o $@

$(OBJDIR)/mcmod_molpro.o: mcmod_molpro.f90
	$(MPIFC) $(INCLUDE) -c $(FFLAGS) $< -o $@

###################################################################################
#Rules for the final executables

################################
#1D double well potential:
pimd_1d_par: $(1DFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(1DFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR)

pimd_1d_ser:  $(1DFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(1DFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#2D C_6 potential:
pimd_2dtest_par: $(2DFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(2DFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR)

pimd_2dtest_ser: $(2DFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(2DFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpi_2dtest_ser: $(2DFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(2DFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#SO(2) potential:
pimd_so2_par:  $(SO2FILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(SO2FILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR)

pimd_so2_ser: $(SO2FILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(SO2FILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpi_so2_ser: $(SO2FILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(SO2FILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#Malonaldehyde:
pimd_malon_par: $(MALONFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(MALONFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR)

pimd_malon_ser: $(MALONFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(MALONFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpi_malon_ser: $(MALONFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(MALONFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#Wang Malonaldehyde:
pimd_wmalon_par: $(WMALONFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(WMALONFLAGS) $(WMALONFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@  $(WMALONLIBS) $(FLIBS_PAR)

pimd_wmalon_ser: $(WMALONFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(WMALONFLAGS) $(WMALONFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(WMALONLIBS) $(FLIBS_SEQ)

################################
#Formic Acid Dimer
pimd_formic_par: $(FORMICFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(FORMICFLAGS) $(FORMICFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@  $(FORMICLIBS) $(FLIBS_PAR)

pimd_formic_ser: $(FORMICFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(FORMICFLAGS) $(FORMICFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FORMICLIBS) $(FLIBS_SEQ)

################################
#Water Dimer HBB2
pimd_hbb2_par: $(HBB2FILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(HBB2FLAGS) $(HBB2FILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@  $(HBB2LIBS) $(FLIBS_PAR)

pimd_hbb2_ser: $(HBB2FILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(HBB2FLAGS) $(HBB2FILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(HBB2LIBS) $(FLIBS_SEQ)

rpi_hbb2_ser: $(HBB2FILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(HBB2FLAGS) $(HBB2FILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(HBB2LIBS) $(FLIBS_SEQ)

################################
#Water Dimer CCPol
pimd_ccpol_par: $(CCPOLFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(CCPOLFLAGS) $(CCPOLFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@  $(FLIBS_PAR)

pimd_ccpol_ser: $(CCPOLFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(CCPOLFLAGS) $(CCPOLFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpi_ccpol_ser: $(CCPOLFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(CCPOLFLAGS) $(CCPOLFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#Water Dimer MB-pol
pimd_watdim_par: $(WATDIMFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(WATDIMFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATDIMLIBS)

pimd_watdim_ser: $(WATDIMFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(WATDIMFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

rpi_watdim_ser: $(WATDIMFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(WATDIMFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

rpi_watdim_par: $(WATDIMFILES) $(COMFILES)  $(OBJDIR)/parallelmod.o $(OBJDIR)/rpi_par.o
	$(MPIFC) $(FFLAGS) $(WATDIMFILES) $(COMFILES)  $(OBJDIR)/parallelmod.o $(OBJDIR)/rpi_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATDIMLIBS)

crossover_watdim: $(WATDIMFILES) $(COMFILES)  $(OBJDIR)/crossover.o
	$(FC) $(FFLAGS) $(WATDIMFILES) $(COMFILES)  $(OBJDIR)/crossover.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

################################
#Water Trimer MB-pol
pimd_wattrim_par: $(WATTRIMFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(WATTRIMFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATDIMLIBS)

pimd_wattrim_ser: $(WATTRIMFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(WATTRIMFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

rpi_wattrim_ser: $(WATTRIMFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(WATTRIMFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

rpi_wattrim_par: $(WATTRIMFILES) $(COMFILES)  $(OBJDIR)/parallelmod.o $(OBJDIR)/rpi_par.o
	$(MPIFC) $(FFLAGS) $(WATTRIMFILES) $(COMFILES)  $(OBJDIR)/parallelmod.o $(OBJDIR)/rpi_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATDIMLIBS)

crossover_wattrim: $(WATTRIMFILES) $(COMFILES)  $(OBJDIR)/crossover.o
	$(FC) $(FFLAGS) $(WATTRIMFILES) $(COMFILES)  $(OBJDIR)/crossover.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

################################
#Water Hexamer MB-pol
pimd_wathex_par: $(WATHEXFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(WATHEXFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATDIMLIBS)

pimd_wathex_ser: $(WATHEXFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(WATHEXFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

rpi_wathex_ser: $(WATHEXFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(WATHEXFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

rpi_wathex_par: $(WATHEXFILES) $(COMFILES)  $(OBJDIR)/parallelmod.o $(OBJDIR)/rpi_par.o
	$(MPIFC) $(FFLAGS) $(WATHEXFILES) $(COMFILES)  $(OBJDIR)/parallelmod.o $(OBJDIR)/rpi_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATDIMLIBS)

################################
#Clathrate
pimd_clath_par: $(CLATHFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(CLATHFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATDIMLIBS)

pimd_clath_ser: $(CLATHFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(CLATHFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

rpi_clath_ser: $(CLATHFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(CLATHFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

################################
#Graphene
pimd_graphene_par: $(GRAPHENEFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(GRAPHENEFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR)

pimd_graphene_ser: $(GRAPHENEFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(GRAPHENEFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpi_graphene_ser: $(GRAPHENEFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(GRAPHENEFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#Water-methane (for Eszter!)

rpi_watmeth_ser: $(WATMETHFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(WATMETHFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

################################
#Water Dimer and Halide
pimd_wathal_par: $(WATHALFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(WATHALFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATHALLIBS)

pimd_wathal_ser: $(WATHALFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(WATHALFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATHALLIBS)

rpi_wathal_ser: $(WATHALFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(WATHALFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATHALLIBS)

rpi_wathal_par: $(WATHALFILES) $(COMFILES)  $(OBJDIR)/parallelmod.o $(OBJDIR)/rpi_par.o
	$(MPIFC) $(FFLAGS) $(WATHALFILES) $(COMFILES)  $(OBJDIR)/parallelmod.o $(OBJDIR)/rpi_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATHALLIBS)

crossover_wathal: $(WATHALFILES) $(COMFILES)  $(OBJDIR)/crossover.o
	$(FC) $(FFLAGS) $(WATHALFILES) $(COMFILES)  $(OBJDIR)/crossover.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATHALLIBS)

################################
#Water Dimer and Halide
pimd_molpro_par: $(MOLPROFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(MOLPROFILES) $(COMFILES)  $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR)

pimd_molpro_ser: $(MOLPROFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(MOLPROFILES) $(COMFILES)  $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpi_molpro_ser: $(MOLPROFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(MOLPROFILES) $(COMFILES)  $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpi_molpro_par: $(MOLPROFILES) $(COMFILES) $(OBJDIR)/parallelmod.o $(OBJDIR)/rpi_par.o
	$(MPIFC) $(FFLAGS) $(MOLPROFILES) $(COMFILES) $(OBJDIR)/parallelmod.o $(OBJDIR)/rpi_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR)

crossover_molpro: $(MOLPROFILES) $(COMFILES)  $(OBJDIR)/crossover.o
	$(FC) $(FFLAGS) $(MOLPROFILES) $(COMFILES)  $(OBJDIR)/crossover.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)


###################################################################################
#Rules for cleanup

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(MODDIR)/*.f90

all: pimd_1d_par pimd_1d_ser pimd_2dtest_par pimd_2dtest_ser rpi_2dtest_ser pimd_so2_par pimd_so2_ser rpi_so2_ser pimd_malon_par pimd_malon_ser rpi_malon_ser pimd_wmalon_par pimd_wmalon_ser pimd_formic_par pimd_formic_ser pimd_ccpol_par pimd_ccpol_ser rpi_ccpol_ser pimd_watdim_par pimd_watdim_ser rpi_watdim_ser rpi_watdim_par crossover_watdim pimd_wathex_par pimd_wathex_ser rpi_wathex_ser pimd_clath_par pimd_clath_ser rpi_clath_ser pimd_graphene_par pimd_graphene_ser rpi_graphene_ser rpi_watmeth_ser pimd_wathal_par pimd_wathal_ser rpi_wathal_ser rpi_wathal_par crossover_wathal
