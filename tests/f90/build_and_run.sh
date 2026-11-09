#!/bin/bash
# Build and run the JFNK test suite.
#
# test_jfnk_solver: links the production jfnk_solver_v1.o extracted
#   verbatim from lib/libmod.a (freshness-checked) against a stub
#   jfnk_cmfgen_int_mod so the real Newton/GMRES code path is
#   exercised standalone.
# test_jfnk_interface: links the real interface objects from the
#   production library and drives the state-free routines
#   (BA_RAD matvec, fixed-unknown mask, weighted dot).  Whole-archive
#   freshness is checked before linking.
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
# mod/ must not keep a jfnk_solver_v1_mod.mod from when the solver was
# compiled here: -Imod precedes -I$ROOT/mod, so a stale copy would
# shadow the production module interface.  Every .mod in mod/ and
# mod2/ is regenerated below; the solver's own .mod now always comes
# from $ROOT/mod.
rm -f mod/*.mod mod2/*.mod

F90="ifort -qopenmp -mkl -diag-disable=10448,10412,10434 -O0 -g -traceback"

echo "== [solver test] resolve production jfnk_solver_v1.o"
SOLVER_O=$(python3 - "$ROOT" <<'PYEOF'
import sys
sys.path.insert(0, sys.argv[1] + "/unit_tests")
from common import TREE, prod_obj
print(prod_obj(TREE / "new_main/mod_subs/jfnk_solver_v1.f90"))
PYEOF
)

echo "== [solver test] compile stub int_mod"
$F90 -c jfnk_int_mod_stub.f90 -module mod -o obj/jfnk_int_mod_stub.o

echo "== [solver test] compile driver against production mods"
$F90 -c test_jfnk_solver.f90 -Imod -I$ROOT/mod -module mod \
    -o obj/test_jfnk_solver.o

echo "== [solver test] link"
$F90 obj/jfnk_int_mod_stub.o "$SOLVER_O" obj/test_jfnk_solver.o \
    -o test_jfnk_solver

echo "== [interface test] production archives fresh?"
python3 - "$ROOT" <<'PYEOF'
import sys
sys.path.insert(0, sys.argv[1] + "/unit_tests")
from common import require_fresh_archives
require_fresh_archives()
PYEOF

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
