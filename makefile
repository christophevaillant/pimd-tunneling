###################################################################################
#Define all the compilers and their options
FC= ifort
MPIFC= mpif90
FFLAGS= -g -O2 -i8 -r8 -w -no-wrap-margin -module ../build/modules
FLIBS_PAR=-L${MKLROOT}/lib/intel64 -mkl=sequential -lmkl_rt -lpthread -lm -ldl
FLIBS_SEQ=-L${MKLROOT}/lib/intel64 -mkl=sequential -lmkl_rt -lpthread -lm -ldl
INCLUDE= -I${MKLROOT}/include
BUILDDIR= ../build
OBJDIR= ../build/objects

###################################################################################
#List the object files and their dependencies
#Define some variables with common files
COMFILES= $(OBJDIR)/blas.o $(OBJDIR)/linpack.o $(OBJDIR)/timer.o \
	$(OBJDIR)/lbfgsb.o $(OBJDIR)/nr_fft.o \
	$(OBJDIR)/verletmodule.o

#1D double well potential:
1DFILES= $(OBJDIR)/mcmod_1d.o

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

#Water dimer (HBB2):
HBB2FILES= $(OBJDIR)/mcmod_waterdimer_hbb2.o
pimd_hbb2_par: FFLAGS+= -O -I./mod_water
pimd_hbb2_ser: FFLAGS+= -O -I./mod_water
rpi_hbb2_ser: FFLAGS+= -O -I./mod_water
HBB2LIBS= -L. -lpes_water -lpesd_water

#Water dimer (CCPol):
# CCPOLFILES= $(OBJDIR)/H2O.pjt2.o $(OBJDIR)/main_CCpol-8sf.o $(OBJDIR)/proc_ccpol8s-dimer_xyz_ncd.o $(OBJDIR)/proc_sapt5sf_new_ncd.o $(OBJDIR)/ang_interface.o $(OBJDIR)/mcmod_waterdimer_ccpol.o
CCPOLFILES= $(OBJDIR)/h2o_ps_pes.o $(OBJDIR)/main_CCpol-8sf.o $(OBJDIR)/proc_ccpol8s-dimer_xyz_ncd.o $(OBJDIR)/proc_sapt5sf_new_ncd.o $(OBJDIR)/ang_interface.o $(OBJDIR)/mcmod_waterdimer_ccpol.o

#Methane Clathrate:
CLATHFILES= $(OBJDIR)/watermethane.o $(OBJDIR)/mcmod_clathrate.o

###################################################################################
#Compilation commands for object files
$(OBJDIR)/%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJDIR)/%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJDIR)/pimd_par.o: pimd_par.f90
	$(MPIFC) -c $(FFLAGS) $< -o $@

$(OBJDIR)/rpi_par.o: rpi_par.f90
	$(MPIFC) -c $(FFLAGS) $< -o $@


###################################################################################
#Rules for the final executables

################################
#1D double well potential:
pimd_1d_par: $(1DFILES) $(COMFILES) $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(COMFILES) $(1DFILES) $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR)

pimd_1d_ser: $(1DFILES) $(COMFILES) $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(COMFILES) $(1DFILES) $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#SO(2) potential:
pimd_so2_par: $(SO2FILES) $(COMFILES) $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(COMFILES) $(SO2FILES) $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR)

pimd_so2_ser: $(SO2FILES) $(COMFILES) $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(COMFILES) $(SO2FILES) $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpi_so2_ser: $(SO2FILES) $(COMFILES) $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(COMFILES) $(SO2FILES) $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#Malonaldehyde:
pimd_malon_par: $(MALONFILES) $(COMFILES) $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(COMFILES) $(MALONFILES) $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR)

pimd_malon_ser: $(MALONFILES) $(COMFILES) $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(COMFILES) $(MALONFILES) $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpi_malon_ser: $(MALONFILES) $(COMFILES) $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(COMFILES) $(MALONFILES) $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#Wang Malonaldehyde:
pimd_wmalon_par: $(WMALONFILES) $(COMFILES) $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(WMALONFLAGS) $(COMFILES) $(WMALONFILES) $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@  $(WMALONLIBS) $(FLIBS_PAR)

pimd_wmalon_ser: $(WMALONFILES) $(COMFILES) $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(WMALONFLAGS) $(COMFILES) $(WMALONFILES) $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(WMALONLIBS) $(FLIBS_SEQ)

################################
#Formic Acid Dimer
pimd_formic_par: $(FORMICFILES) $(COMFILES) $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(FORMICFLAGS) $(COMFILES) $(FORMICFILES) $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@  $(FORMICLIBS) $(FLIBS_PAR)

pimd_formic_ser: $(FORMICFILES) $(COMFILES) $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(FORMICFLAGS) $(COMFILES) $(FORMICFILES) $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FORMICLIBS) $(FLIBS_SEQ)

################################
#Water Dimer HBB2
pimd_hbb2_par: $(HBB2FILES) $(COMFILES) $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(HBB2FLAGS) $(COMFILES) $(HBB2FILES) $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@  $(HBB2LIBS) $(FLIBS_PAR)

pimd_hbb2_ser: $(HBB2FILES) $(COMFILES) $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(HBB2FLAGS) $(COMFILES) $(HBB2FILES) $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(HBB2LIBS) $(FLIBS_SEQ)

rpi_hbb2_ser: $(HBB2FILES) $(COMFILES) $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(HBB2FLAGS) $(COMFILES) $(HBB2FILES) $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(HBB2LIBS) $(FLIBS_SEQ)

################################
#Water Dimer CCPol
pimd_ccpol_par: $(CCPOLFILES) $(COMFILES) $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(CCPOLFLAGS) $(COMFILES) $(CCPOLFILES) $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@  $(FLIBS_PAR)

pimd_ccpol_ser: $(CCPOLFILES) $(COMFILES) $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(CCPOLFLAGS) $(COMFILES) $(CCPOLFILES) $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpi_ccpol_ser: $(CCPOLFILES) $(COMFILES) $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(CCPOLFLAGS) $(COMFILES) $(CCPOLFILES) $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#Water Dimer MB-pol
pimd_watdim_par: $(WATDIMFILES) $(COMFILES) $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(COMFILES) $(WATDIMFILES) $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATDIMLIBS)

pimd_watdim_ser: $(WATDIMFILES) $(COMFILES) $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(COMFILES) $(WATDIMFILES) $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

rpi_watdim_ser: $(WATDIMFILES) $(COMFILES) $(OBJDIR)/rpi_ser.o
	$(FC) $(FFLAGS) $(COMFILES) $(WATDIMFILES) $(OBJDIR)/rpi_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

rpi_watdim_par: $(WATDIMFILES) $(COMFILES) $(OBJDIR)/rpi_par.o
	$(MPIFC) $(FFLAGS) $(COMFILES) $(WATDIMFILES) $(OBJDIR)/rpi_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATDIMLIBS)

################################
#Clathrate
pimd_clath_par: $(CLATHFILES) $(COMFILES) $(OBJDIR)/pimd_par.o
	$(MPIFC) $(FFLAGS) $(COMFILES) $(CLATHFILES) $(OBJDIR)/pimd_par.o -o $(BUILDDIR)/$@ $(FLIBS_PAR) $(WATDIMLIBS)

pimd_clath_ser: $(CLATHFILES) $(COMFILES) $(OBJDIR)/pimd_ser.o
	$(FC) $(FFLAGS) $(COMFILES) $(CLATHFILES) $(OBJDIR)/pimd_ser.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ) $(WATDIMLIBS)

###################################################################################
#Rules for cleanup

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/modules/*.mod
