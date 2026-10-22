# Toy NLTE convergence bake-off harness.  Links the production
# newton_hook_v1 (new_main/mod_subs/jfnk_solver_v1.f90, compiled
# unchanged) against the toy model.
#   make           build bakeoff (ifort if available, else gfortran)
#   make test      build + sanity + easy-instance smoke run
# FC may be a full path; the basename selects the flag set.
ifeq ($(origin FC),default)
ifneq ($(shell command -v ifort 2>/dev/null),)
FC := ifort
else
FC := gfortran
endif
endif
FC_BASE := $(notdir $(FC))

ifeq ($(FC_BASE),ifort)
FFLAGS = -O2 -g -traceback -fpp
LDLIBS = -qmkl -Wl,-rpath,$(MKLROOT)/lib/intel64
else ifeq ($(FC_BASE),gfortran)
FFLAGS = -O2 -g -fcheck=all -fbacktrace -cpp -ffree-line-length-none -fno-range-check
LDLIBS = -Wl,-rpath,/usr/lib64 -l:libopenblas.so.0
endif

SOLVER = ../../new_main/mod_subs/jfnk_solver_v1.f90

bakeoff: jfnk_int_stub.o toy_se_mod.o toy_ctx_mod.o jfnk_solver_v1.o bakeoff.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

jfnk_int_stub.o: jfnk_int_stub.f90 checkfc
	$(FC) $(FFLAGS) -c $<

toy_se_mod.o: toy_se_mod.f90 checkfc
	$(FC) $(FFLAGS) -c $<

toy_ctx_mod.o: toy_ctx_mod.f90 checkfc toy_se_mod.o
	$(FC) $(FFLAGS) -c $<

jfnk_solver_v1.o: $(SOLVER) checkfc jfnk_int_stub.o
	$(FC) $(FFLAGS) -c $(SOLVER)

bakeoff.o: bakeoff.f90 checkfc jfnk_int_stub.o toy_se_mod.o toy_ctx_mod.o jfnk_solver_v1.o
	$(FC) $(FFLAGS) -c $<

checkfc:
	@if [ "$(FC_BASE)" != ifort ] && [ "$(FC_BASE)" != gfortran ]; then \
	  echo "unsupported FC '$(FC)' - use ifort or gfortran"; exit 1; fi

test: bakeoff
	./bakeoff sanity params_easy.in
	./bakeoff s0 params_easy.in
	./bakeoff s2f params_easy.in
	./bakeoff s2c params_easy.in
	./bakeoff noise params_easy.in

clean:
	rm -f *.o *.mod bakeoff

.PHONY: test clean checkfc
