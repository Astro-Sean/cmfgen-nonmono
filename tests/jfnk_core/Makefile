# Standalone unit test of the Newton core (newton_hook_v1 in
# new_main/mod_subs/jfnk_solver_v1.f90, compiled unchanged).
# Run: module load intel/2024.0 mkl; make test
FC = ifort
FFLAGS = -O2 -g -traceback -fpp
LDLIBS = -qmkl -Wl,-rpath,$(MKLROOT)/lib/intel64

SOLVER = ../../new_main/mod_subs/jfnk_solver_v1.f90

test: test_newton_core test_guards
	./test_newton_core
	./test_guards

test_newton_core: jfnk_int_stub.o jfnk_solver_v1.o test_newton_core.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_guards: jfnk_int_stub.o jfnk_solver_v1.o test_guards.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

test_guards.o: test_guards.f90 jfnk_int_stub.o jfnk_solver_v1.o
	$(FC) $(FFLAGS) -c $<

jfnk_int_stub.o: jfnk_int_stub.f90
	$(FC) $(FFLAGS) -c $<

jfnk_solver_v1.o: $(SOLVER) jfnk_int_stub.o
	$(FC) $(FFLAGS) -c $(SOLVER)

test_newton_core.o: test_newton_core.f90 jfnk_int_stub.o jfnk_solver_v1.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod test_newton_core

.PHONY: test clean
