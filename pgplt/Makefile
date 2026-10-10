# Include local system definitions

include ../Makefile_definitions

# Library where object files will/are stored

LIB=$(LIB_DIR)libpgplt.a
LIB_MOD=$(LIB_DIR)libmod.a

#
# Entries in a library
#
OBJ = $(LIB)(mod_curve_data.o)\
      $(LIB)(line_id_mod.o)\
      $(LIB)(adjust_id_waves.o)\
      $(LIB)(amoeba.o)\
      $(LIB)(amotry.o)\
      $(LIB)(arrange_pg_curve_ids.o)\
      $(LIB)(change_color.o)\
      $(LIB)(change_pen.o)\
      $(LIB)(change_xaxis_gridding.o)\
      $(LIB)(curve.o)\
      $(LIB)(curve_lab.o)\
      $(LIB)(cut_points_from_plot.o)\
      $(LIB)(define_more_pens.o)\
      $(LIB)(do_cursor_balmer.o)\
      $(LIB)(do_gauss_smooth.o)\
      $(LIB)(do_pg_deriv.o)\
      $(LIB)(do_file_balmer_v1.o)\
      $(LIB)(do_file_ew_v1.o)\
      $(LIB)(do_fill.o)\
      $(LIB)(do_many_ew.o)\
      $(LIB)(do_mult_gf.o)\
      $(LIB)(do_pg_regrid.o)\
      $(LIB)(do_vec_op.o)\
      $(LIB)(do_cursor_ew_v3.o)\
      $(LIB)(draw_line_ids_v2.o)\
      $(LIB)(draw_ew_ids.o)\
      $(LIB)(draw_ew_lines.o)\
      $(LIB)(draw_right_y_axis.o)\
      $(LIB)(ed_gauss_fit.o)\
      $(LIB)(ew_group.o)\
      $(LIB)(find_lines.o)\
      $(LIB)(fntick.o)\
      $(LIB)(gauss_fit.o)\
      $(LIB)(gauss_fit_data.o)\
      $(LIB)(gauss_fit_func.o)\
      $(LIB)(gauss_fit_er.o)\
      $(LIB)(gauss_romb_v2.o)\
      $(LIB)(get_gramon_min_max.o)\
      $(LIB)(get_line_id_pg.o)\
      $(LIB)(get_pgi_file_info.o)\
      $(LIB)(get_titles.o)\
      $(LIB)(gramon_pgplot.o)\
      $(LIB)(hist_adj.o)\
      $(LIB)(justify_convert_v2.o)\
      $(LIB)(modify_pgi_ps.o)\
      $(LIB)(mod_ew_variables.o)\
      $(LIB)(mon_num.o)\
      $(LIB)(monbord_v4.o)\
      $(LIB)(moninside.o)\
      $(LIB)(mon_interp_sp.o)\
      $(LIB)(pg_add_noise.o)\
      $(LIB)(pg_def_continuum.o)\
      $(LIB)(pg_replace_data_int.o)\
      $(LIB)(pg_mod_cont_nodes.o)\
      $(LIB)(pg_set_curve_tit_length.o)\
      $(LIB)(put_text.o)\
      $(LIB)(rd_ew_id.o)\
      $(LIB)(rd_line_ids.o)\
      $(LIB)(modify_curve.o)\
      $(LIB)(mod_smear_pg.o)\
      $(LIB)(set_upper_axis.o)\
      $(LIB)(simp_ew.o)\
      $(LIB)(strip_slash_p.o)\
      $(LIB)(shrink_vectors.o)\
      $(LIB)(smooth_plt.o)\
      $(LIB)(spacing.o)\
      $(LIB)(step_plot_v1.o)\
      $(LIB)(vectorpen.o)\
      $(LIB)(write_line_lambdas.o)

all : $(LIB)
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
	rm -f *.o

$(LIB) : $(OBJ)	

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG)  $<
	ar $(AR_OPTION) $(LIB) $*.o

# 
# Specific dependency of source files which uses external 
# fortran moudles. (The order of compilation matters.)
#
# NB: change_folor.f, change_pen.f, gramon_pgplot.f, vectorpen.f
#     depend on ../subs/gen_in

GEN = ../lib/libsubs.a(gen_in.o)
COL=  ../lib/libsubs.a(mod_color_pen_def.o)

$(LIB)(change_color.o)    : change_color.f  $(GEN)
$(LIB)(change_pen.o)      : change_pen.f    $(GEN)

$(LIB)(draw_line_ids.o)      : draw_line_ids_v2.f $(LIB)(line_id_mod.o)  $(COL) 
$(LIB)(draw_line_ids_v2.o)   : draw_line_ids.f    $(LIB)(line_id_mod.o)  $(COL) 
$(LIB)(draw_ew_ids.o)        : draw_ew_ids.f      $(LIB)(line_id_mod.o)  $(COL) 
$(LIB)(draw_ew_lines.o)      : draw_ew_lines.f    $(LIB)(line_id_mod.o)  $(COL) 
$(LIB)(rd_ew_id.o)           : rd_ew_id.f         $(LIB)(line_id_mod.o)  $(COL) 
$(LIB)(get_line_id.o)        : get_line_id.f      $(LIB)(line_id_mod.o)  $(COL) 
$(LIB)(rd_line_ids.o)        : rd_line_ids.f      $(LIB)(line_id_mod.o)  $(COL) 

$(LIB)(gramon_pgplot.o)         : gramon_pgplot.f       $(LIB)(mod_curve_data.o) $(LIB)(line_id_mod.o) $(GEN) $(COL) 
$(LIB)(curve.o)                 : curve.f                  $(LIB)(mod_curve_data.o) $(GEN)
$(LIB)(curve_lab.o)             : curve_lab.f              $(LIB)(mod_curve_data.o) $(GEN)
$(LIB)(vectorpen.o)             : vectorpen.f              $(LIB)(mod_curve_data.o) $(GEN)
$(LIB)(arrange_pg_curve_ids.o)  : arrange_pg_curve_ids.f   $(LIB)(mod_curve_data.o) $(GEN)

$(LIB)(do_cursor_ew_v3.o) : do_cursor_ew_v3.f   $(LIB)(mod_curve_data.o) $(GEN) $(COL)  $(LIB)(mod_ew_variables.o)
$(LIB)(do_file_ew_v1.o)   : do_file_ew_v1.f     $(LIB)(mod_curve_data.o) $(GEN) $(COL)  $(LIB)(mod_ew_variables.o)

$(LIB)(gauss_fit.o)        : gauss_fit.f         $(LIB)(gauss_fit_data.o)
$(LIB)(gauss_fit_func.o)   : gauss_fit_func.f    $(LIB)(gauss_fit_data.o)
$(LIB)(gauss_fit_er.o)     : gauss_fit_er.f      $(LIB)(gauss_fit_data.o)
$(LIB)(gauss_romb_v2.o)    : gauss_romb_v2.f     $(LIB)(gauss_fit_data.o)
$(LIB)(ed_gauss_fit.o)     : ed_gaus_fit.f       $(LIB)(gauss_fit_data.o)

$(LIB)(do_mult_gf.o)          : do_mult_gf.f            $(LIB)(gaus_fit_data.o)
$(LIB)(do_gauss_smooth.o)     : do_gauss_smooth.f       $(LIB)(mod_smear_pg.o)

$(GEN):
	(cd ../tools; make )

$(COL):
	(cd ../tools; make )

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
