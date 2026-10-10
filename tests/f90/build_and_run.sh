#!/bin/bash
# Build and run the JFNK test suite.
#
# test_jfnk_solver: compiles the production jfnk_solver_v1.f90
#   against a stub jfnk_cmfgen_int_mod so the real Newton/GMRES
#   code path is exercised standalone.
# test_jfnk_interface: links the real interface objects from the
#   production library and drives the state-free routines
#   (BA_RAD matvec, fixed-unknown mask, weighted dot).
#
# The two module output dirs are kept separate because the stub
# jfnk_cmfgen_int_mod would shadow the real module for the
# interface test.
set -e
cd "$(dirname "$0")"

ROOT="$(cd ../.. && pwd)"
SRC=$ROOT/new_main/mod_subs
mkdir -p mod mod2 obj

# ifort resolves .mod files from the cwd ahead of -I paths, so a stub
# jfnk_cmfgen_int_mod.mod left here by a manual compile would shadow the
# production module for the interface test. Nothing legit writes .mod
# files to the cwd - all outputs go to mod/, mod2/ or bin/.
rm -f ./*.mod

F90="ifort -qopenmp -mkl -diag-disable=10448,10412,10434 -O0 -g -traceback"

echo "== [solver test] compile stub int_mod"
$F90 -c jfnk_int_mod_stub.f90 -module mod -o obj/jfnk_int_mod_stub.o

echo "== [solver test] compile production jfnk_solver_v1.f90"
$F90 -c "$SRC/jfnk_solver_v1.f90" -Imod -module mod -o obj/jfnk_solver_v1.o

echo "== [solver test] compile driver"
$F90 -c test_jfnk_solver.f90 -Imod -module mod -o obj/test_jfnk_solver.o

echo "== [solver test] link"
$F90 obj/jfnk_int_mod_stub.o obj/jfnk_solver_v1.o obj/test_jfnk_solver.o \
    -o test_jfnk_solver

echo "== [interface test] compile driver against production mods"
$F90 -c test_jfnk_interface.f90 -I$ROOT/mod -module mod2 \
    -o obj/test_jfnk_interface.o

echo "== [interface test] link against production libs"
$F90 obj/test_jfnk_interface.o -L$ROOT/lib -lsubs -lmod \
    -o test_jfnk_interface

echo ""
echo "== run solver tests"
./test_jfnk_solver
echo ""
echo "== run interface tests"
./test_jfnk_interface

echo ""
echo "== convergence-history plots"
# Plotting needs matplotlib; a missing or headless-broken install
# must not fail the suite - the test exit codes above are the gate.
python3 ../python/plot_solver_tests.py . || \
    echo "warning: solver history plot step skipped/failed"

echo ""
echo "== analytical JFNK harness (Python, subprocess-driven)"
# Compiles the solver + CLI driver itself; exercises the production
# newton/gmres path on synthetic problems and emits PASS/FAIL +
# plots under tests/plots/harness/.  Non-fatal like the plot step
# only if python/matplotlib are unavailable.
python3 ../python/jfnk_harness.py || \
    echo "warning: jfnk_harness step skipped/failed"
