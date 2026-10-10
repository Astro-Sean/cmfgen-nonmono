# Include local system definitions

include ../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries for library
#
OBJ = $(LIB)(electron_density_calc_v2.o)\
      $(LIB)(eta_scat_v6.o)\
      $(LIB)(gamma_emiss_eta_v9.o)\
      $(LIB)(gamma_energy_dep_v7.o)\
      $(LIB)(gamma_nu_grid_v16.o)\
      $(LIB)(gamma_xray_opac_v4_test.o)\
      $(LIB)(gamray_sub_v3.o)\
      $(LIB)(klein_freq.o)\
      $(LIB)(read_nuc_decay_data_v3.o)\
      $(LIB)(rd_gamray_cntrl_v2.o)\
      $(LIB)(write_array_v2.o)\
      $(LIB)(write_array_v2_ISO.o)

all : $(LIB) 
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o

$(LIB) : $(OBJ) 

STEQ=$(LIB_MOD)(steq_data_mod.o)
MOD=$(LIB_MOD)(mod_cmfgen.o)
NUC=$(LIB_MOD)(nuc_iso_mod.o)
GAM_NUC=$(LIB_MOD)(mod_read_nuc_decay_data_v2.o)
GAM=$(LIB_MOD)(mod_gamma_trans_v3.o)
GAM_CNTRL=$(LIB_MOD)(mod_rd_gamray_cntrl_variables.o)
GAM_MU=$(LIB_MOD)(gam_mu_mod.o)

$(LIB)(%.o) : %.f $(STEQ) $(MOD) $(NUC)
	$(F90) -c $(FG) $<
	ar ruv $(LIB) $*.o

$(LIB)(electron_density_calc_v2.o) : electron_density_calc_v2.f90 $(MOD)
	$(F90) -c $(FFREE) -l$(MOD_DIR) electron_density_calc_v2.f90
	ar ruv $(LIB) electron_density_calc_v2.o

$(LIB)(eta_scat_v6.o) : eta_scat_v6.f90 $(GAM_MU)
	$(F90) -c $(FFREE) -l$(MOD_DIR)  eta_scat_v6.f90
	ar ruv $(LIB) eta_scat_v6.o

$(LIB)(gamma_emiss_eta_v9.o) : gamma_emiss_eta_v9.f90 $(NUC) $(GAM_NUC)
	$(F90) -c $(FFREE) -l$(MOD_DIR)  gamma_emiss_eta_v9.f90
	ar ruv $(LIB) gamma_emiss_eta_v9.o

$(LIB)(gamma_energy_dep_v7.o) : gamma_energy_dep_v7.f90 $(GAM_MU)
	$(F90) -c $(FFREE) -l$(MOD_DIR)  gamma_energy_dep_v7.f90
	ar ruv $(LIB) gamma_energy_dep_v7.o

$(LIB)(gamma_nu_grid_v16.o) : gamma_nu_grid_v16.f90 $(NUC) $(GAM_NUC) $(GAM_CNTRL)
	$(F90) -c $(FFREE) -l$(MOD_DIR)  gamma_nu_grid_v16.f90
	ar ruv $(LIB) gamma_nu_grid_v16.o

$(LIB)(gamma_xray_opac_v4_test.o) : gamma_xray_opac_v4_test.f $(MOD)
	$(F90) -c $(FG) -l$(MOD_DIR) gamma_xray_opac_v4_test.f
	ar ruv $(LIB) gamma_xray_opac_v4_test.o

$(LIB)(gamray_sub_v3.o) : gamray_sub_v3.f90 $(GAM) $(GAM_CNTRL) $(NUC)
	$(F90) -c $(FFREE) -l$(MOD_DIR)  gamray_sub_v3.f90
	ar ruv $(LIB) gamray_sub_v3.o

$(LIB)(klein_freq.o) : klein_freq.f90
	$(F90) -c $(FFREE)  klein_freq.f90
	ar ruv $(LIB) klein_freq.o

$(LIB)(read_nuc_decay_data_v3.o) : read_nuc_decay_data_v3.f90 $(GAM_NUC)
	$(F90) -c $(FFREE) -l$(MOD_DIR)  read_nuc_decay_data_v3.f90
	ar ruv $(LIB) read_nuc_decay_data_v3.o

$(LIB)(rd_gamray_cntrl_v2.o) : rd_gamray_cntrl_v2.f90 $(GAM_CNTRL)
	$(F90) -c $(FFREE) -l$(MOD_DIR)  rd_gamray_cntrl_v2.f90
	ar ruv $(LIB) rd_gamray_cntrl_v2.o

$(LIB)(write_array_v2.o) : write_array_v2.f90 $(GAM_MU)
	$(F90) -c $(FFREE) -l$(MOD_DIR)  write_array_v2.f90
	ar ruv $(LIB) write_array_v2.o

$(LIB)(write_array_v2_ISO.o) : write_array_v2_ISO.f90
	$(F90) -c $(FFREE) -l$(MOD_DIR)  write_array_v2_ISO.f90
	ar ruv $(LIB) write_array_v2_ISO.o

$(STEQ):
	(cd ../mod_subs ; make)

$(MOD):
	(cd ../mod_subs ; make)

$(NUC):
	(cd ../mod_subs ; make)

$(GAM_NUC):
	(cd ../mod_subs ; make)

$(GAM):
	(cd ../mod_subs ; make)

$(GAM_CNTRL):
	(cd ../mod_subs ; make)

$(GAM_MU):
	(cd ../mod_subs ; make)

clean :
	rm -f *.o 

clean_lib :
	rm $(LIB)
