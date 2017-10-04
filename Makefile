#_a makefile



#FC= gfortran-mp-4.5 -O2   -Warray-bounds -O0 -g -fbounds-check
FC= gfortran -std=f2008  -O2 #  -Warray-bounds -p -O0 -g -fbounds-check

INC=/Users/evenhuis/Dropbox/fortran_programs/lib/
BLAS= -L/usr/lib -lblas -llapack

LIB_DIR = /Users/evenhuis/Dropbox/fortran_programs/lib/
LIB_O   = -I$(LIB_DIR) -L$(LIB_DIR) -lmylib 

%: %.f90 util_mod.o
	$(FC) $(BLAS) -o $* $*.f90   $(LIB_O) util_mod.o


# - - - - - - - Biomass - - - -  - --  - - --  - - - -
biomass_model_mod.o : biomass_model_mod.f90  							      ext_driver_mod.o
	$(FC) $(BLAS) -c               			   biomass_model_mod.f90 ext_driver_mod.o $(LIB_O)

biomass_model.so: biomass_model_mod.o 
	touch biomass_model.so
	rm    biomass_model.so
	$(FC) $(BLAS) -c     py_biomass_mod.f90 biomass_model_mod.o ext_driver_mod.o $(LIB_O)
	f2py -c -m biomass_model py_biomass_mod.f90 biomass_model_mod.o ext_driver_mod.o $(LIB_O) $(BLAS)

biomass_model:                     biomass_model.f90     biomass_model_mod.o  ext_driver_mod.o
	$(FC) $(BLAS) -o biomass_model  biomass_model.f90     biomass_model_mod.o  ext_driver_mod.o $(LIB_O)



util_mod.o: util_mod.f90
	$(FC) $(BLAS) -c util_mod.f90 $(LIB_O)

ext_driver_mod.o : ext_driver_mod.f90
	 $(FC) $(BLAS) -c ext_driver_mod.f90 $(LIB_O)





