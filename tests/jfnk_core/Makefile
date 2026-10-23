# Standalone unit test of the Newton core (newton_hook_v1 in
# new_main/mod_subs/jfnk_solver_v1.f90, compiled unchanged).
# Run: module load intel/2024.0 mkl; make test
#      make FC=gfortran test   (no modules needed; links openblas)
# FC may also be a full path; the basename selects the flag set.
# origin check, not '?=': GNU make predefines FC=f77 as a builtin
# default; this sets ifort only when FC comes from that builtin, while
# a command-line or environment FC still wins.
ifeq ($(origin FC),default)
FC := ifort
endif
FC_BASE := $(notdir $(FC))

ifeq ($(FC_BASE),ifort)
FFLAGS = -O2 -g -traceback -fpp
LDLIBS = -qmkl -Wl,-rpath,$(MKLROOT)/lib/intel64
else ifeq ($(FC_BASE),gfortran)
FFLAGS = -O2 -g -fcheck=all -fbacktrace -cpp -ffree-line-length-none -fno-range-check -fallow-argument-mismatch
# System openblas carries the LAPACK symbols; there is no unversioned
# libopenblas.so symlink, so -l:libopenblas.so.0 is required.
LDLIBS = -Wl,-rpath,/usr/lib64 -l:libopenblas.so.0
endif

SOLVER = ../../new_main/mod_subs/jfnk_solver_v1.f90

# .mod files are compiler-specific; refuse to mix toolchains and purge
# any stale set before compiling with the other compiler.
checkfc:
	@if [ "$(FC_BASE)" != ifort ] && [ "$(FC_BASE)" != gfortran ]; then \
	  echo "unsupported FC '$(FC)' - use ifort or gfortran"; exit 1; fi

test: test_newton_core test_guards test_precond_blocks test_ptc_sign test_trace_state
	./test_newton_core
	./test_guards
	./test_precond_blocks
	./test_ptc_sign
	./test_trace_state

test_newton_core: jfnk_int_stub.o jfnk_solver_v1.o test_newton_core.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_guards: jfnk_int_stub.o jfnk_solver_v1.o test_guards.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_guards.o: test_guards.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o
	$(FC) $(FFLAGS) -c $<

test_precond_blocks: jfnk_int_stub.o jfnk_solver_v1.o test_precond_blocks.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_precond_blocks.o: test_precond_blocks.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o
	$(FC) $(FFLAGS) -c $<

test_ptc_sign: jfnk_int_stub.o jfnk_solver_v1.o test_ptc_sign.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_ptc_sign.o: test_ptc_sign.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o
	$(FC) $(FFLAGS) -c $<

# Stage-B state machine: links the production control_variable_mod
# (fixed-form, needs -extend-source) so BA_TRACE_EVAL is tested, not
# a mirror.
set_kind_module.o: ../../tools/set_kind_module.f checkfc
	$(FC) $(FFLAGS) -c $< -o set_kind_module.o

control_variable_mod.o: ../../new_main/mod_subs/control_variable_mod.f set_kind_module.o checkfc
	$(FC) $(FFLAGS) -extend-source -c $< -o control_variable_mod.o

test_trace_state: control_variable_mod.o test_trace_state.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_trace_state.o: test_trace_state.f90 checkfc control_variable_mod.o
	$(FC) $(FFLAGS) -c $<

jfnk_int_stub.o: jfnk_int_stub.f90 checkfc
	$(FC) $(FFLAGS) -c $<

jfnk_solver_v1.o: $(SOLVER) checkfc jfnk_int_stub.o
	$(FC) $(FFLAGS) -c $(SOLVER)

test_newton_core.o: test_newton_core.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod test_newton_core test_guards test_precond_blocks test_ptc_sign test_trace_state

.PHONY: test clean checkfc
