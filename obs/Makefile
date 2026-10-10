#Include file with system dependent definitions

include ../Makefile_definitions

OPTION=$(FG)

#
# These names are for the fortran linker
#

LN1=subs
LN2=disp
LN5=pgplt
LN6=mod

#
# These names refer to the full file name
#

FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN5=$(LIB_DIR)lib$(LN5).a
FN6=$(LIB_DIR)lib$(LN6).a

#
# We only make the executables if there has been a change to the main *.f files,
# or if the libraries have changed since the last link.
#
 
OBJ1  =  cmf_flux_v5.o cmf_flux_sub_v5.o mod_cmf_obs.o mod_freq_obs.o \
        cmf_flux_cntrl_var_mod.o rd_cmf_flux_controls.o check_cmf_flux_param_consis.o \
        rd_model_file.o revise_obs_p.o  ins_line_v5.o obs_frame_sub_v9.o\
        ins_line_obs_v4.o rd_cont_j_obs.o set_forbid_zero.o form_ew.o \
        comp_ew_abs_ew.o
 
OBJ2  =  obs_frame_v2.o obs_frame_sub_v9.o 

INC = ABS_EW.INC CMF_ABS_EW.INC COMP_JCONT_V4.INC EVAL_LTE_INC_V5.INC  OPACITIES_V5.INC XOPAC_V4.INC  

all: $(EXE_DIR)cmf_flux.exe $(EXE_DIR)obs_frame.exe

#
# Give the dependencies of cmf_flux.exe and how it is created.
# 

$(EXE_DIR)cmf_flux.exe :  $(OBJ1) $(FN1) $(FN2) $(FN6) $(FN7) $(INC)
	ar -xo ../lib/libsubs.a griem_stark_mod.o
	# Publish atomically: KEEP_EXE=1 leaves the previous binary live while
	# build.sh rebuilds, so a queued OBS job must not see a truncated exe.
	$(F90) $(OBJ1) griem_stark_mod.o \
              -L$(LIB_DIR) -l$(LN1) -l$(LN5) -l$(LN2) -l$(LN6) \
               $(LOCLIB) \
               -o $(EXE_DIR)cmf_flux.exe.tmp
	mv -f $(EXE_DIR)cmf_flux.exe.tmp $(EXE_DIR)cmf_flux.exe

$(EXE_DIR)obs_frame.exe :  $(OBJ2) $(FN1) $(FN2)
	$(F90) $(OBJ2) \
              -L$(LIB_DIR) -l$(LN2) -l$(LN1) \
               $(LOCLIB) \
              -o  $(EXE_DIR)obs_frame.exe.tmp
	mv -f $(EXE_DIR)obs_frame.exe.tmp $(EXE_DIR)obs_frame.exe

# General rule to compile a fortran ".f" file.

%.o : %.f
	$(F90) -c $(OPTION)  $<

#
# I use F90noamp because of an issue with the PGF95 compiler on the Mac PRO.
#

cmf_flux_v5.o : cmf_flux_v5.f  mod_cmf_obs.o  cmf_flux_cntrl_var_mod.o $(USR)
	$(F90) -c $(OPTION) cmf_flux_v5.f 

#	$(F90nomp) -c $(OPTION) cmf_flux_v5.f 



cmf_flux_sub_v5.o : cmf_flux_sub_v5.f  mod_cmf_obs.o  mod_freq_obs.o cmf_flux_cntrl_var_mod.o $(GEN) \
   COMP_JCONT_V4.INC  EVAL_LTE_INC_V5.INC\
   OPACITIES_V5.INC   XOPAC_V4.INC
	$(F90) -c $(OPTION) cmf_flux_sub_v5.f 

#	$(F90) -c $(OPTION) cmf_flux_sub_v5.f 




USR = ../lib/libtools.a(usr_option.o)
GEN = ../lib/libmod.a(mod_lev_dis_blk.o)

rd_cmf_flux_controls.o : cmf_flux_cntrl_var_mod.o
check_cmf_flux_param_consis.o : cmf_flux_cntrl_var_mod.o

$(GEN):
	(cd ../subs; make )
$(USR):
	(cd ../tools; make )

 
$(FN1) : lib1
$(FN2) : lib2
$(FN5) : lib5
$(FN6) : lib6

# We now MAKE the required libraries.

lib1:
	(cd ../tools; make)
	(cd ../subs; make )
	(cd ../unix; make )
	(cd ../stark; make )
	(cd ../newsubs; make)
	(cd ../new_main/subs; make)

lib2:
	(cd ../disp/subs; make)
lib5:
	(cd ../pgplt; make)
lib6:
	(cd ../new_main/mod_subs; make)

clean :
	rm -f *.o

mod_freq_obs.o : mod_freq_obs.f
mod_freq:
	$(F90) -c $(OPTION) mod_freq_obs.f
	cp -f $(INSTALL_DIR)perm_mod_file jnk.mod
	mv -f *.mod $(MOD_DIR) 
