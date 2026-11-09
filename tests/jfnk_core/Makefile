# Standalone unit test of the Newton core (newton_hook_v1 in
# new_main/mod_subs/jfnk_solver_v1.f90).
# Run: module load intel/2024.0 mkl; make test
# The solver and control_variable_mod under test are the production
# objects extracted verbatim from lib/libmod.a; only the test drivers
# and the jfnk_cmfgen_int_mod stub compile here.  ifort is required -
# the production .o/.mod files are ifort artifacts gfortran cannot use.
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
$(error gfortran cannot link the production ifort objects in lib/libmod.a - build with ifort)
endif

# The solver under test is the production object, extracted from
# lib/libmod.a (audit directive: link the shipped .o, not a recompile).
# Its matching .mod is copied into pmod/ so the stub's
# jfnk_cmfgen_int_mod.mod stays the only int-module visible here.
LIBMOD = ../../lib/libmod.a
PMOD = pmod

# Member-level archive freshness: fails when any lib/*.a member is
# older than its production source (run make DEBUG=nodebug all first).
.PHONY: fresh
fresh:
	@python3 -c "import sys; sys.path.insert(0,'../../unit_tests'); \
	  from common import require_fresh_archives; \
	  require_fresh_archives()"

# .mod files are compiler-specific; refuse to mix toolchains and purge
# any stale set before compiling with the other compiler.
checkfc:
	@if [ "$(FC_BASE)" != ifort ] && [ "$(FC_BASE)" != gfortran ]; then \
	  echo "unsupported FC '$(FC)' - use ifort or gfortran"; exit 1; fi

test: test_newton_core test_guards test_precond_blocks test_ptc_sign test_trace_state test_nm_window test_pc_quality test_schur test_gmres_edge
	./test_newton_core
	./test_guards
	./test_precond_blocks
	./test_ptc_sign
	./test_trace_state
	./test_nm_window
	./test_pc_quality
	./test_schur
	./test_gmres_edge

test_newton_core: jfnk_int_stub.o jfnk_solver_v1.o test_newton_core.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_guards: jfnk_int_stub.o jfnk_solver_v1.o test_guards.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_guards.o: test_guards.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o $(PMOD)/jfnk_solver_v1_mod.mod
	$(FC) $(FFLAGS) -I$(PMOD) -c $<

test_precond_blocks: jfnk_int_stub.o jfnk_solver_v1.o test_precond_blocks.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_precond_blocks.o: test_precond_blocks.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o $(PMOD)/jfnk_solver_v1_mod.mod
	$(FC) $(FFLAGS) -I$(PMOD) -c $<

test_ptc_sign: jfnk_int_stub.o jfnk_solver_v1.o test_ptc_sign.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_ptc_sign.o: test_ptc_sign.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o $(PMOD)/jfnk_solver_v1_mod.mod
	$(FC) $(FFLAGS) -I$(PMOD) -c $<

test_nm_window: jfnk_int_stub.o jfnk_solver_v1.o test_nm_window.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_nm_window.o: test_nm_window.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o $(PMOD)/jfnk_solver_v1_mod.mod
	$(FC) $(FFLAGS) -I$(PMOD) -c $<

test_pc_quality: jfnk_int_stub.o jfnk_solver_v1.o test_pc_quality.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_pc_quality.o: test_pc_quality.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o $(PMOD)/jfnk_solver_v1_mod.mod
	$(FC) $(FFLAGS) -I$(PMOD) -c $<

test_schur: jfnk_int_stub.o jfnk_solver_v1.o test_schur.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_schur.o: test_schur.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o $(PMOD)/jfnk_solver_v1_mod.mod
	$(FC) $(FFLAGS) -I$(PMOD) -c $<

test_gmres_edge: jfnk_int_stub.o jfnk_solver_v1.o test_gmres_edge.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_gmres_edge.o: test_gmres_edge.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o $(PMOD)/jfnk_solver_v1_mod.mod
	$(FC) $(FFLAGS) -I$(PMOD) -c $<

# Stage-B state machine: links the production control_variable_mod.o
# from lib/libmod.a so BA_TRACE_EVAL is tested, not a mirror.  Its
# .mod is staged into pmod/ like the solver's.
$(PMOD)/control_variable_mod.mod: ../../mod/control_variable_mod.mod
	mkdir -p $(PMOD) && cp $< $@

$(PMOD)/set_kind_module.mod: ../../mod/set_kind_module.mod
	mkdir -p $(PMOD) && cp $< $@

control_variable_mod.o: $(LIBMOD) | fresh
	ar x $(LIBMOD) control_variable_mod.o

test_trace_state: control_variable_mod.o test_trace_state.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_trace_state.o: test_trace_state.f90 checkfc control_variable_mod.o $(PMOD)/control_variable_mod.mod $(PMOD)/set_kind_module.mod
	$(FC) $(FFLAGS) -I$(PMOD) -c $<

jfnk_int_stub.o: jfnk_int_stub.f90 checkfc
	$(FC) $(FFLAGS) -c $<

# Production object, verbatim from the archive the exe links.
jfnk_solver_v1.o: $(LIBMOD) | fresh
	ar x $(LIBMOD) jfnk_solver_v1.o

$(PMOD)/jfnk_solver_v1_mod.mod: ../../mod/jfnk_solver_v1_mod.mod
	mkdir -p $(PMOD) && cp $< $@

test_newton_core.o: test_newton_core.f90 checkfc jfnk_int_stub.o jfnk_solver_v1.o $(PMOD)/jfnk_solver_v1_mod.mod
	$(FC) $(FFLAGS) -I$(PMOD) -c $<

clean:
	rm -f *.o *.mod test_newton_core test_guards test_precond_blocks test_ptc_sign test_trace_state test_nm_window test_pc_quality test_schur test_gmres_edge
	rm -rf $(PMOD)

.PHONY: test clean checkfc
