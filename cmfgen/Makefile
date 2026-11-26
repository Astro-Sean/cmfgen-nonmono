# TEST MAKEFILE


# Set parameters for arechiv .o files into libraries.
# On UBUNTU systems msy nrrf Uruv

AR_OPTION = ruv

#Include local system dependent definitions

include Makefile_definitions


# We will access the makefile in each local directory tp create the
# libraries and executables.

all : info d_tools d_blas d_dust d_lpack d_pgplot_stubs d_nm_subs d_stark d_unix d_subs d_plane d_newsubs d_pgplt \
         d_lte d_new_main d_disp d_spec_plt d_obs d_misc

# Convenience target: build only cmfgen_dev.exe and its dependencies.
# This skips all plotting tools and is much faster for development.
cmfgen : info d_tools d_blas d_lpack d_pgplot_stubs d_mod_subs d_nm_subs d_stark d_unix d_subs d_plane d_newsubs d_pgplt d_lte d_new_main

# Build PGPLOT stub library (needed for linking on systems without PGPLOT).
d_pgplot_stubs:
	(cd pgplot_stubs; make)

# We now MAKE the required libraries and executables.

info:
	@echo " "
	@echo " Host         = "$(HOST)
	@echo " Compile opt. = "$(COMPILE_OPT)
	@echo " Install dir. = "$(INSTALL_DIR)
	@echo " (F90)        = "$(F90)
	@echo " (FG)         = "$(FG)
	@echo "ar option     = "$(AR_OPTION)
	@echo " "

d_blas:
	$(MAKE) -C blas
d_lpack:
	$(MAKE) -C lpack
d_dust:
	$(MAKE) -C dust
d_tools:
	$(MAKE) -C tools
d_mod_subs:
	$(MAKE) -C new_main/mod_subs
d_nm_subs: d_mod_subs
	$(MAKE) -C new_main/subs
d_plane:
	$(MAKE) -C plane
d_stark:
	$(MAKE) -C stark
d_subs:
	$(MAKE) -C subs
d_unix:
	$(MAKE) -C unix
d_newsubs:
	$(MAKE) -C newsubs
d_pgplt:
	$(MAKE) -C pgplt

# The following will create the executables

# Only the archive objects are needed: do_cmf_hydro_v2.o in libsubs.a
# references det_r_grid_v1/v2, wind_vel_law_v3 and get_lte_ross_v2.
# The utility exes are not part of the cmfgen link; main_lte.exe cannot
# link under -O0 anyway (lte_sub.f calls WR_ABUND_INFO_V2, which has no
# implementation in this tree -- only V3 exists in subs/wr_x_info.f).
d_lte:
	$(MAKE) -C lte_hydro OBJ_FILES
d_new_main:
	$(MAKE) -C new_main
d_disp:
	$(MAKE) -C disp
d_spec_plt:
	$(MAKE) -C spec_plt
d_obs:
	$(MAKE) -C obs
d_misc:
	$(MAKE) -C misc
#
# To use the following command enter"
#         make -i clean
# The -i is necessary in case some files don't exist.

clean:
	rm -f lib/*.a
	rm -f mod/*.mod
	# build.sh sets KEEP_EXE=1 so the previous executables stay installed
	# while a production rebuild runs; each production link then publishes
	# through a temporary name and mv.
	if [ -z "$(KEEP_EXE)" ]; then rm -f exe/*.exe; fi
	rm -f */*.o
	rm -f */*/*.o

dbg_clean:
	rm -f lib_dbg/*.a
	rm -f mod_dbg/*.mod
	rm -f exe_dbg/*.exe
	rm -f */*.o
	rm -f */*/*.o
