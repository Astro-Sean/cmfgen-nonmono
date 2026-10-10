# Makefile for PGPLOT stub library
#
# Builds a static library of empty PGPLOT routine stubs so that
# PGPLOT-dependent executables can link on systems without PGPLOT.
# The stubs do nothing - plotting tools will not produce output,
# but the CMFGEN solver is completely unaffected.
#
# Added 2026-09-07 by Devin for SN2024pba V4 build on Raven.

include ../Makefile_definitions

LIB = $(LIB_DIR)libpgplot_stubs.a

all : $(LIB)

$(LIB) : pgplot_stubs.o
	ar $(AR_OPTION) $(LIB) pgplot_stubs.o

pgplot_stubs.o : pgplot_stubs.f
	$(F90nomp) -c -extend-source -O0 -I$(MOD_DIR) pgplot_stubs.f

clean :
	rm -f *.o $(LIB)
