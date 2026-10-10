#        Include file with system dependent definitions

include ../Makefile_definitions

OPTION=$(FG)


# These names are for the fortran linker

LN1=subs
LN2=pgplt
LN3=disp
LN4=spec


# These names refer to the full file name

FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a


# We only make dispgen.exe if there has been a change to the main *.f files,
# or if the libraries have changed since the last link.

 
OBJ1  =  set_hydro_hd_loc.o update_pg_filename.o adjust_atm_r_grid.o rd_f_to_s_ids.o rd_ba_mat.o rd_f_to_s_ids_v2.o write_osc_v2.o
 
EXE   =  $(EXE_DIR)wr_f_to_s.exe\
         $(EXE_DIR)ave_dc.exe\
         $(EXE_DIR)check_ba_estimates.exe\
         $(EXE_DIR)create_batobs_ins.exe\
         $(EXE_DIR)create_batch_ins.exe\
         $(EXE_DIR)create_kappa_file.exe\
         $(EXE_DIR)create_if_file.exe\
         $(EXE_DIR)create_mod_sum.exe\
         $(EXE_DIR)create_pp_mod_sum.exe\
         $(EXE_DIR)count_phot_data.exe\
         $(EXE_DIR)do_recom_sum.exe\
         $(EXE_DIR)ed_vadat.exe\
         $(EXE_DIR)check_energy_cons.exe\
         $(EXE_DIR)chk_ba_est.exe\
         $(EXE_DIR)chk_phot.exe\
         $(EXE_DIR)pack_osc.exe\
         $(EXE_DIR)do_ng_v2.exe\
         $(EXE_DIR)extrap_scr.exe\
         $(EXE_DIR)guess_dc.exe\
         $(EXE_DIR)land_col_merge.exe\
         $(EXE_DIR)land_multi_merge.exe\
         $(EXE_DIR)mod_cool.exe\
         $(EXE_DIR)mod_prrr.exe\
         $(EXE_DIR)n_col_merge.exe\
         $(EXE_DIR)n_multi_merge.exe\
         $(EXE_DIR)rewrite_scr.exe\
         $(EXE_DIR)rewrite_dc_v2.exe\
         $(EXE_DIR)rewrite_old_model_data.exe\
         $(EXE_DIR)rewrite_jh.exe\
         $(EXE_DIR)rev_rdinr.exe\
         $(EXE_DIR)rev_rvsig.exe\
         $(EXE_DIR)recomz.exe\
         $(EXE_DIR)scale_omega_v1.exe\
         $(EXE_DIR)show_float_kind.exe\
         $(EXE_DIR)append_dc.exe\
         $(EXE_DIR)append_r_dc.exe\
         $(EXE_DIR)rev_hydro_file.exe\
         $(EXE_DIR)rev_hydro_turb.exe\
         $(EXE_DIR)mass_sc.exe\
         $(EXE_DIR)plt_alpha.exe\
         $(EXE_DIR)plt_cmf_lum.exe\
         $(EXE_DIR)plt_ddt_work.exe\
         $(EXE_DIR)plt_gf.exe\
         $(EXE_DIR)plt_hydro.exe\
         $(EXE_DIR)plt_ion_mf.exe\
         $(EXE_DIR)plt_lte_ross.exe\
         $(EXE_DIR)plt_raw_xray.exe\
         $(EXE_DIR)plt_cool.exe\
         $(EXE_DIR)plt_cool_sort.exe\
         $(EXE_DIR)plt_ln_heat.exe\
         $(EXE_DIR)plt_j_comp.exe\
         $(EXE_DIR)plt_q_col.exe\
         $(EXE_DIR)plt_net.exe\
         $(EXE_DIR)plt_non_therm.exe\
         $(EXE_DIR)plt_rates.exe\
         $(EXE_DIR)plt_steq.exe\
         $(EXE_DIR)set_new_sn_mod.exe\
         $(EXE_DIR)solve_ba_mat.exe\
         $(EXE_DIR)tst_ba_mat_sol.exe\
         $(EXE_DIR)comp_ba_mat_sol.exe\
         $(EXE_DIR)tst_char_v2.exe\
         $(EXE_DIR)tlusty_vel.exe
 
all: OBJ ALLEXE 
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)

ALLEXE : $(EXE)
OBJ    :  $(OBJ1) $(FN1) $(FN2) $(FN3)

OANDF = $(OBJ1) $(FN1) $(FN2) $(FN3)

# Give the dependencies of ?.exe and how it is created.

$(EXE_DIR)%.exe : %.o $(OANDF) 
	$(F90) -g $*.o \
              $(OBJ1) -L$(LIB_DIR) -l$(LN3) -l$(LN2) -l$(LN1) -l$(LN4) -l$(LN1) \
              -lmod \
              $(PGLIB) $(LAPACK) $(BLAS) $(LOCLIB) $(X11LIB) $(G2CLIB) \
              -g -o $(EXE_DIR)$*.exe

$(EXE_DIR)plt_rates.exe : plt_rates.o 
	$(F90) -g plt_rates.o \
              -L$(LIB_DIR) -l$(LN3) -l$(LN2) -l$(LN1) -l$(LN4) -l$(LN1) \
              $(PGLIB) $(LAPACK) $(BLAS) $(LOCLIB) $(X11LIB) $(G2CLIB) \
              -g -o $(EXE_DIR)plt_rates.exe

$(EXE_DIR)chk_phot.exe : chk_phot.o recom_chk_v2.o get_exc_freq.o $(OANDF)
	$(F90) chk_phot.o recom_chk_v2.o get_exc_freq.o \
              $(OBJ1) -L$(LIB_DIR) -l$(LN2) -l$(LN1) \
              $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)chk_phot.exe

$(EXE_DIR)guess_dc.exe : guess_dc.o $(OANDF) 
	$(F90) guess_dc.o \
              $(OBJ1) -L$(LIB_DIR) -l$(LN2) -l$(LN1) \
              $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)guess_dc.exe

$(EXE_DIR)pack_osc.exe : pack_osc.o $(OANDF)
	$(F90) pack_osc.o \
              $(OBJ1) -L$(LIB_DIR)  -l$(LN2) -l$(LN1) \
              $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)pack_osc.exe

$(EXE_DIR)rev_rdinr.exe : rev_rdinr.o 
	$(F90) rev_rdinr.o \
              $(OBJ1) -L$(LIB_DIR) -l$(LN2) -l$(LN1) \
              $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)rev_rdinr.exe

$(EXE_DIR)tlusty_vel.exe : tlusty_vel.o $(OANDF) 
	$(F90) tlusty_vel.o \
              $(OBJ1) -L$(LIB_DIR)  -l$(LN2) -l$(LN1) \
              $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)tlusty_vel.exe

$(EXE_DIR)solve_ba_mat.exe : solve_ba_mat.o
	$(F90) solve_ba_mat.o -L$(LIB_DIR) $(LAPACK) $(BLAS) -l$(LN2) -l$(LN1) \
              -o $(EXE_DIR)solve_ba_mat.exe

$(EXE_DIR)wr_f_to_s.exe : wr_f_to_s.o output_term_seq.o $(OANDF)
	$(F90) wr_f_to_s.o output_term_seq.o \
              $(OBJ1) -L$(LIB_DIR) -l$(LN2) -l$(LN1) -l$(LN2)  \
              $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -o $(EXE_DIR)wr_f_to_s.exe

# General rule to compile a fortran ".f" file.

%.o : %.f 
	$(F90)  -c $(OPTION)  $<

USR = $(LIB_DIR)libsubs.a(usr_option.o)
GEN = $(LIB_DIR)libsubs.a(gen_in.o)

wr_f_to_s.o     : wr_f_to_s.f  $(USR)

$(USR):
	(cd ../tools; make )
$(GEN):
	(cd ../tools; make )

 
$(FN1) : lib1
$(FN2) : lib2
$(FN3) : lib3
$(FN4) : lib4
$(FN5) : lib5

# We now MAKE the required libraries.

lib1:
	(cd ../subs; make )

lib2:
	(cd ../pgplt; make )

clean :
	rm -f *.o
