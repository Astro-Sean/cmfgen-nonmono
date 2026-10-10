# Include local system definitions

include ../../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libsubs.a
LIB_MOD=$(LIB_DIR)libmod.a

OPTION=$(FG)

# Entries for library

MOD = $(LIB_MOD)(mod_cmfgen.o)\
      $(LIB_MOD)(ang_qw_mod.o)\
      $(LIB_MOD)(chg_exch_mod_v3.o)\
      $(LIB_MOD)(cmf_sob_mod.o)\
      $(LIB_MOD)(control_variable_mod.o)\
      $(LIB_MOD)(eddfac_rec_defs_mod.o)\
      $(LIB_MOD)(gam_mu_mod.o)\
      $(LIB_MOD)(line_mod.o)\
      $(LIB_MOD)(line_vec_mod.o)\
      $(LIB_MOD)(mod_gamma_trans_v3.o)\
      $(LIB_MOD)(mod_rd_gamray_cntrl_variables.o)\
      $(LIB_MOD)(mod_read_nuc_decay_data_v2.o)\
      $(LIB_MOD)(mod_lev_dis_blk.o)\
      $(LIB_MOD)(nuc_iso_mod.o)\
      $(LIB_MOD)(opac_mod.o)\
      $(LIB_MOD)(radiation_mod.o)\
      $(LIB_MOD)(shock_power_mod.o)\
      $(LIB_MOD)(steq_data_mod.o)\
      $(LIB_MOD)(two_phot_mod.o)\
      $(LIB_MOD)(var_rad_mod.o)

OBJ = $(LIB)(adjust_r_grid_v4.o)\
      $(LIB)(adjust_den_vecs.o)\
      $(LIB)(check_param_consistency.o)\
      $(LIB)(comp_j_blank.o)\
      $(LIB)(comp_grey_v2.o)\
      $(LIB)(comp_grey_v4.o)\
      $(LIB)(comp_opac.o)\
      $(LIB)(comp_var_opac.o)\
      $(LIB)(do_var_cont.o)\
      $(LIB)(eval_lte_v5.o)\
      $(LIB)(grey_t_iterate.o)\
      $(LIB)(init_opac_mod.o)\
      $(LIB)(init_line_opac_var_v2.o)\
      $(LIB)(rd_control_variables.o)\
      $(LIB)(rd_cont_j.o)\
      $(LIB)(set_line_opac.o)\
      $(LIB)(set_ang_qw.o)\
      $(LIB)(set_ang_qw_v2.o)\
      $(LIB)(set_frequency_grid_v2.o)\
      $(LIB)(set_new_model_estimates.o)\
      $(LIB)(set_abund_clump.o)\
      $(LIB)(solve_for_pops.o)\
      $(LIB)(jfnk_cmfgen_int_mod.o)\
      $(LIB_MOD)(jfnk_solver_v1.o)\
      $(LIB)(jfnk_cmfgen_interface_v1.o)\
      $(LIB)(jfnk_solve_corr_v1.o)\
      $(LIB)(sub_cmf_line.o)\
      $(LIB)(sub_sob_line_v3.o)\
      $(LIB)(wind_scale_pops_v1.o)\
      $(LIB)(sup_to_full_v4.o)

all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o

$(LIB) : $(MOD) $(OBJ)	

UK=$(LIB)(update_keyword.o)

# General rules to make obects in a library and object files

$(LIB_MOD)(%.o) : %.f 
	$(F90) -c $(OPTION) $< 
	ar $(AR_OPTION) $(LIB_MOD) $*.o

$(LIB)(%.o) : %.f 
	$(F90) -c $(OPTION) $< 
	ar $(AR_OPTION) $(LIB) $*.o

$(MOD_LIB)(%.o) : %.f90
	$(F90) -c $(FFREE) $<
	ar $(AR_OPTION) $(LIB_MOD) $*.o

# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)

$(LIB)(adjust_r_grid_v4.o) : adjust_r_grid_v4.f $(LIB_MOD)(mod_cmfgen.o) $(UK)
	$(F90) -c $(OPTION) adjust_r_grid_v4.f
	ar $(AR_OPTION) $(LIB) adjust_r_grid_v4.o

$(LIB)(check_param_consistency.o) : check_param_consistency.f $(LIB_MOD)(control_variable_mod.o)
	$(F90) -c $(OPTION) check_param_consistency.f
	ar $(AR_OPTION) $(LIB) check_param_consistency.o

$(LIB)(comp_grey_v2.o) : comp_grey_v2.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
                                  $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(OPTION) comp_grey_v2.f
	ar $(AR_OPTION) $(LIB) comp_grey_v2.o

$(LIB)(comp_grey_v4.o) : comp_grey_v4.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
                                  $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(OPTION) comp_grey_v4.f
	ar $(AR_OPTION) $(LIB) comp_grey_v4.o

# MOD_RAY_MOM_STORE lives in plane/; comp_j_blank.f imports HNU_AT_OB for the
# DO_NONMONO_KNOP dispatch.  Build it here so a clean build is not order-dependent.
$(LIB)(mod_ray_mom_store.o) : ../../plane/mod_ray_mom_store.f
	$(F90) -c $(OPTION) ../../plane/mod_ray_mom_store.f
	ar $(AR_OPTION) $(LIB) mod_ray_mom_store.o

$(LIB)(comp_j_blank.o) : comp_j_blank.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(radiation_mod.o) $(LIB_MOD)(mod_cmfgen.o) $(LIB_MOD)(ang_qw_mod.o) $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) comp_j_blank.f
	ar $(AR_OPTION) $(LIB) comp_j_blank.o

# FG_J_CMF_NONMONO_MOD_V1 lives in plane/; jfnk_solve_corr_v1.f imports
# FREEZE_I_NEXT (freezes the lagged Knop I_NEXT store during JFNK
# residual evals).  Built here for the .mod file, same as
# mod_ray_mom_store above; plane/Makefile archives the same member.
$(LIB)(fg_j_cmf_nonmono_v1.o) : ../../plane/fg_j_cmf_nonmono_v1.f $(LIB)(mod_ray_mom_store.o)
	$(F90) -c $(OPTION) ../../plane/fg_j_cmf_nonmono_v1.f
	ar $(AR_OPTION) $(LIB) fg_j_cmf_nonmono_v1.o

$(LIB)(comp_opac.o) : comp_opac.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(OPTION) comp_opac.f
	ar $(AR_OPTION) $(LIB) comp_opac.o

$(LIB)(comp_var_opac.o) : comp_var_opac.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(OPTION) comp_var_opac.f
	ar $(AR_OPTION) $(LIB) comp_var_opac.o

$(LIB)(do_var_cont.o) : do_var_cont.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
             $(LIB_MOD)(line_mod.o) $(LIB_MOD)(var_rad_mod.o) $(LIB_MOD)(radiation_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
             $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(OPTION) do_var_cont.f
	ar $(AR_OPTION) $(LIB) do_var_cont.o
                                                                                                                      
$(LIB)(eval_lte_v5.o) : eval_lte_v5.f $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(OPTION) eval_lte_v5.f
	ar $(AR_OPTION) $(LIB) eval_lte_v5.o

$(LIB)(grey_t_iterate.o) : grey_t_iterate.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
             $(LIB_MOD)(line_mod.o) $(LIB_MOD)(line_vec_mod.o) $(LIB_MOD)(mod_cmfgen.o) $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(OPTION) grey_t_iterate.f
	ar $(AR_OPTION) $(LIB) grey_t_iterate.o
                                                                                                                      
$(LIB)(init_line_opac_var_v2.o) : init_line_opac_var_v2.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
             $(LIB_MOD)(line_mod.o) $(LIB_MOD)(mod_cmfgen.o) $(LIB_MOD)(var_rad_mod.o)
	$(F90) -c $(OPTION) init_line_opac_var_v2.f
	ar $(AR_OPTION) $(LIB) init_line_opac_var_v2.o
                                                                                                                      
$(LIB)(rd_cont_j.o) : rd_cont_j.f $(LIB)(opac_mod.o) $(LIB_MOD)(radiation_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
             $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(OPTION) rd_cont_j.f
	ar $(AR_OPTION) $(LIB) rd_cont_j.o
                                                                                                                      
$(LIB)(rd_control_variables.o) : rd_control_variables.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o) \
              $(LIB_MOD)(jfnk_solver_v1.o)
	$(F90) -c $(OPTION) rd_control_variables.f
	ar $(AR_OPTION) $(LIB) rd_control_variables.o

$(LIB)(set_abund_clump.o) : set_abund_clump.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(OPTION) set_abund_clump.f
	ar $(AR_OPTION) $(LIB) set_abund_clump.o
                                                                                                                      
$(LIB)(set_ang_qw.o) : set_ang_qw.f $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(OPTION) set_ang_qw.f
	ar $(AR_OPTION) $(LIB) set_ang_qw.o

$(LIB)(set_frequency_grid_v2.o) : set_frequency_grid_v2.f $(LIB_MOD)(control_variable_mod.o)\
              $(LIB_MOD)(line_vec_mod.o) $(LIB_MOD)(mod_cmfgen.o) 
	$(F90) -c $(OPTION) set_frequency_grid_v2.f
	ar $(AR_OPTION) $(LIB) set_frequency_grid_v2.o
                                                                                                                      
$(LIB)(set_line_opac.o) : set_line_opac.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
              $(LIB_MOD)(line_mod.o) $(LIB_MOD)(mod_cmfgen.o) $(LIB_MOD)(line_vec_mod.o)
	$(F90) -c $(OPTION) set_line_opac.f
	ar $(AR_OPTION) $(LIB) set_line_opac.o
                                                                                                                      
$(LIB)(set_new_model_estimates.o) : set_new_model_estimates.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
              $(LIB_MOD)(line_mod.o) $(LIB_MOD)(line_vec_mod.o) $(LIB_MOD)(mod_cmfgen.o) $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(OPTION) set_new_model_estimates.f
	ar $(AR_OPTION) $(LIB) set_new_model_estimates.o
                                                                                                                      
$(LIB)(solve_for_pops.o) : solve_for_pops.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
                                  $(LIB_MOD)(ang_qw_mod.o)
	$(F90) -c $(OPTION) solve_for_pops.f
	ar $(AR_OPTION) $(LIB) solve_for_pops.o

# JFNK population solver (V7 port of V4/V6 machinery, 02-Oct-2026):
# jfnk_solver_v1.f90 is a free-format self-contained GMRES/Newton
# module; jfnk_cmfgen_* provide the CMFGEN callbacks; the coupled-
# residual radiation sweep is the internal subroutine
# JFNK_RAD_SWEEP_V1 in cmfgen_sub.f (JFNK_COUPLED);
# jfnk_solve_corr_v1 is the driver called in place of SOLVEBA_V13.
# editted by Devin (V7-JFNK, 02-Oct-2026)

$(LIB_MOD)(jfnk_solver_v1.o) : jfnk_solver_v1.f90 $(LIB)(jfnk_cmfgen_int_mod.o)
	$(F90) -c $(FFREE) jfnk_solver_v1.f90
	ar $(AR_OPTION) $(LIB_MOD) jfnk_solver_v1.o

$(LIB)(jfnk_cmfgen_int_mod.o) : jfnk_cmfgen_int_mod.f $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(OPTION) jfnk_cmfgen_int_mod.f
	ar $(AR_OPTION) $(LIB) jfnk_cmfgen_int_mod.o

$(LIB)(jfnk_cmfgen_interface_v1.o) : jfnk_cmfgen_interface_v1.f $(LIB_MOD)(mod_cmfgen.o) \
                                  $(LIB_MOD)(steq_data_mod.o) $(LIB_MOD)(control_variable_mod.o) \
                                  $(LIB_MOD)(ang_qw_mod.o) $(LIB_MOD)(jfnk_solver_v1.o) \
                                  $(LIB)(jfnk_cmfgen_int_mod.o)
	$(F90) -c $(OPTION) jfnk_cmfgen_interface_v1.f
	ar $(AR_OPTION) $(LIB) jfnk_cmfgen_interface_v1.o

$(LIB)(jfnk_solve_corr_v1.o) : jfnk_solve_corr_v1.f $(LIB_MOD)(control_variable_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
                                  $(LIB_MOD)(steq_data_mod.o) $(LIB_MOD)(jfnk_solver_v1.o) \
                                  $(LIB)(jfnk_cmfgen_int_mod.o) $(LIB)(jfnk_cmfgen_interface_v1.o) \
                                  $(LIB)(fg_j_cmf_nonmono_v1.o)
	$(F90) -c $(OPTION) jfnk_solve_corr_v1.f
	ar $(AR_OPTION) $(LIB) jfnk_solve_corr_v1.o

$(LIB)(sub_cmf_line.o) : sub_cmf_line.f LINEGEN.INC $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
              $(LIB_MOD)(line_mod.o) $(LIB_MOD)(var_rad_mod.o) $(LIB_MOD)(radiation_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
              $(LIB_MOD)(cmf_sob_mod.o) $(LIB_MOD)(steq_data_mod.o) $(LIB_MOD)(ang_qw_mod.o) $(LIB_MOD)(line_vec_mod.o)
	$(F90) -c $(OPTION) sub_cmf_line.f
	ar $(AR_OPTION) $(LIB) sub_cmf_line.o
                                                                                                                      
$(LIB)(sub_sob_line_v3.o) : sub_sob_line_v3.f $(LIB)(opac_mod.o) $(LIB_MOD)(control_variable_mod.o)\
              $(LIB_MOD)(line_mod.o) $(LIB_MOD)(var_rad_mod.o) $(LIB_MOD)(radiation_mod.o) $(LIB_MOD)(mod_cmfgen.o)\
              $(LIB_MOD)(cmf_sob_mod.o) $(LIB_MOD)(steq_data_mod.o) $(LIB_MOD)(ang_qw_mod.o) $(LIB_MOD)(line_vec_mod.o)
	$(F90) -c $(OPTION) sub_sob_line_v3.f
	ar $(AR_OPTION) $(LIB) sub_sob_line_v3.o
                                                                                                                      
$(LIB)(sup_to_full_v4.o) : sup_to_full_v4.f $(LIB_MOD)(mod_cmfgen.o)
	$(F90) -c $(OPTION) sup_to_full_v4.f
	ar $(AR_OPTION) $(LIB) sup_to_full_v4.o

$(UK):
	(cd ../../tools ; make)

clean :
	rm -f *.o 
