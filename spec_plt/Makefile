#Include file with system dependent definitions

include ../Makefile_definitions

# These names are for the fortran linker

LN1=subs
LN2=pgplt
LN4=disp
LN5=spec

# These names refer to the full file name

FN1=$(LIB_DIR)lib$(LN1).a
FN2=$(LIB_DIR)lib$(LN2).a
FN4=$(LIB_DIR)lib$(LN4).a
FN5=$(LIB_DIR)lib$(LN5).a
     
#  $(EXE_DIR)plt_jh_shell.exe\

EXE  = $(EXE_DIR)plt_jh.exe\
       $(EXE_DIR)plt_jh_cur.exe\
       $(EXE_DIR)plt_ip.exe\
       $(EXE_DIR)plt_delta.exe\
       $(EXE_DIR)plt_dfr.exe\
       $(EXE_DIR)plt_many_mod.exe\
       $(EXE_DIR)plt_many_sn_spec.exe\
       $(EXE_DIR)plt_phot_raw.exe\
       $(EXE_DIR)plt_rj.exe\
       $(EXE_DIR)plt_scr.exe\
       $(EXE_DIR)plt_steq.exe\
       $(EXE_DIR)plt_spec.exe\
       $(EXE_DIR)new_plt_rj.exe

all : OBJ ALLEXE


GEN =  $(LIB_DIR)libsubs.a(gen_in.o)\
		$(LIB_DIR)libsubs.a(usr_option.o)\
		$(LIB_DIR)libsubs.a(usr_hidden.o)\
		$(LIB_DIR)libmod.a(mod_lev_dis_blk.o)\
		$(LIB_DIR)libsubs.a(mod_color_pen_def.o)\
		$(LIB_DIR)libsubs.a(zm_data_module.o)\
		$(LIB_DIR)libsubs.a(wr_string.o)

# Give the dependencies of dispgen.exe and how it is created.


OBJ       : $(FN1) $(FN2) $(FN4) $(FN5) 
ALLEXE    : $(EXE)


# A MACRO deinition could probably be used here, but it is difficult to get it to work on all
# UNIX systems, hence we take the easy but tediuos method.
#	$(F90) $*.o /Users/jdh/pgplot/pgplotsrc/grsy00.o  \

$(EXE_DIR)%.exe : %.o $(OBJ) $(FN1) $(FN2) $(FN4) $(FN5)
	$(F90) $*.o \
		$(OBJ1) -L$(LIB_DIR) -l$(LN5) -l$(LN4) -l$(LN2) -l$(LN1) -lmod \
		$(LOCLIB) $(PGLIB) $(LAPACK) $(BLAS) $(X11LIB) $(G2CLIB) \
		-o $(EXE_DIR)$*.exe

$(EXE_DIR)plt_phot_raw.exe : plt_phot_raw.o $(FN1) $(FN2) $(FN4) $(FN5) 
	$(F90) plt_phot_raw.o \
             -L$(LIB_DIR) -l$(LN5)  -l$(LN4) -l$(LN2) -l$(LN1) \
               $(LOCLIB) $(PGLIB) $(X11LIB) $(G2CLIB) \
              -g -o $(EXE_DIR)plt_phot_raw.exe

%.o : %.f $(GEN)
	$(F90) -c $(FG) $<

plt_phot_raw.o : plt_phot_raw.f $(GEN) $(hyd_bf)
	$(F90) -c $(FG) plt_phot_raw.f

HYD_BF = ../lib/libnewsubs.a(rd_hyd_bf_data.o)

$(HYD_BF):
	(cd ../newsubs; make)

$(SUBS):
	(cd subs; make)

$(GEN):
	(cd ../tools; make)

$(FN1) : lib1
$(FN2) : lib2
$(FN4) : lib4
$(FN5) : lib5

# We now MAKE the required libraries.

lib1:
	(cd ../tools; make )
	(cd ../subs; make ) 
lib2:
	(cd ../pgplt; make )
lib4:
	(cd ../disp/subs; make ) 
lib5:
	(cd subs; make ) 

#clean :
#	rm -f *.o

clean_lib :
	rm -f $(LIB)
