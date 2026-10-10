# TEST MAKEFILE

#Include file with system dependent definitions

include ../Makefile_definitions

LIB=$(LIB_DIR)libsubs.a

#
# Entries for library
#

OBJ = $(LIB)(det_r_grid_v1.o)\
      $(LIB)(det_r_grid_v2.o)\
      $(LIB)(get_lte_ross_v3.o)\
      $(LIB)(get_lte_ross_v2.o)\
      $(LIB)(wind_vel_law_v1.o)\
      $(LIB)(wind_vel_law_v2.o)\
      $(LIB)(wind_vel_law_v3.o)
#
# These names are for the fortran linker
#
LN1=subs
LN2=pgplt
LN6=mod
LN7=disp

#
# These names refer to the full file name
#

FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN6=$(LIB_DIR)lib$(LN6).a
FN7=$(LIB_DIR)lib$(LN7).a

#
# We only make dispgen.exe if there has been a change to the main *.f files,
# or if the libraries have changed since the last link.
#
  
all : $(OBJ)   $(EXE_DIR)main_lte.exe $(EXE_DIR)wind_hyd.exe
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR)
$(LIB) : $(OBJ)

OBJ_FILES : $(OBJ)

# General rules to make obects in a library and object files

$(LIB)(%.o) : %.f 
	$(F90) -c $(FG)  $<
	ar $(AR_OPTION) $(LIB) $*.o

#
# Give the dependencies of dispgen.exe and how it is created.
#

$(EXE_DIR)main_lte.exe : lte.o lte_sub.o $(OBJ_FILES) set_lte_est.o det_ed_v2.o lte_pop_sl.o \
                         ltepop_wld_v1.o genopaeta_v10.o set_line_opac.o writedc_v3.o \
                          $(FN1) $(FN2) $(FN6)
	$(F90) lte.o lte_sub.o set_lte_est.o det_ed_v2.o lte_pop_sl.o ltepop_wld_v1.o \
                     genopaeta_v10.o set_line_opac.o writedc_v3.o \
             -L$(LIB_DIR) -l$(LN1) -l$(LN2) -l$(LN6)  \
               $(LOCLIB) $(LAPACK) $(BLAS) \
               -o $(EXE_DIR)main_lte.exe

$(EXE_DIR)wind_hyd.exe : wind_hyd.o $(OBJ_FILES) \
                          $(FN1) $(FN6) $(FN7)
	$(F90)  wind_hyd.o  \
             -L$(LIB_DIR) -l$(LN7) -l$(LN2) -l$(LN1) -l$(LN6) \
               $(LOCLIB) $(LAPACK) $(BLAS) $(PGLIB) $(G2CLIB) $(X11LIB) \
               -o $(EXE_DIR)wind_hyd.exe


USR=$(LIB_DIR)libsubs.a(usr_option.o)\
    $(LIB_DIR)libsubs.a(gen_in.o)

MODS= $(LIB_DIR)libmod.a(mod_cmfgen.o)\
      $(LIB_DIR)libmod.a(cmf_sob_mod.o)\
      $(LIB_DIR)libmod.a(control_variable_mod.o)\
      $(LIB_DIR)libmod.a(ang_qw_mod.o)\
      $(LIB_DIR)libmod.a(line_mod.o)\
      $(LIB_DIR)libmod.a(line_vec_mod.o)\
      $(LIB_DIR)libmod.a(mod_lev_dis_blk.o)\
      $(LIB_DIR)libmod.a(opac_mod.o)\
      $(LIB_DIR)libmod.a(radiation_mod.o)\
      $(LIB_DIR)libmod.a(steq_data_mod.o)\
      $(LIB_DIR)libmod.a(var_rad_mod.o)

lte.o : lte.f $(MODS) $(USR)
	$(F90) -c $(FG) lte.f 


lte_sub.o : lte_sub.f $(MODS)
	$(F90) -c $(FG) lte_sub.f 

wind_hyd.o : wind_hyd.f $(USR)
	$(F90) -c $(FG) wind_hyd.f 

det_ed_v2.o : det_ed_v2.f $(MODS)
	$(F90) -c $(FG) det_ed_v2.f 

genopaeta_v10.o : genopaeta_v10.f $(MODS)
	$(F90) -c $(FG) genopaeta_v10.f 

writedc_v3.o : writedc_v3.f $(MODS)
	$(F90) -c $(FG) writedc_v3.f 

set_line_opac.o : set_line_opac.f $(MODS)
	$(F90) -c $(FG) set_line_opac.f 

set_lte_est.o : set_lte_est.f $(MODS)
	$(F90) -c $(FG) set_lte_est.f 

%.o : %.f 
	$(F90) -c $(FG) $< 

$(USR):
	(cd ../tools ; make )

$(MODS): mod

mod:
	(cd ../new_main/mod_subs ; make )

# We now MAKE the required libraries. Note that this double
# dependence is rquired to ensure that the libraries are updated.
# Otherwise Makfile simple checks date of library.

$(FN1): lib1
$(FN2): lib2
$(FN6): lib6
$(FN7): lib7

lib1:
	(cd ../tools; make )
	(cd ../subs; make )
	(cd ../newsubs; make )
	(cd ../stark; make )
	(cd ../new_main/subs; make )
lib2:
	(cd ../pgplt; make )
lib6:
	(cd ../new_main/mod_subs; make )
lib7:
	(cd ../disp; make )

clean :
	rm -f *.o

clean_lib :
	rm -f $(LIB)
