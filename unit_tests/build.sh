#!/bin/bash
# Build the surrogate-SE CLI driver against the production solver.
set -e
TD="$(cd "$(dirname "$0")" && pwd)"
TREE="$(dirname "$TD")"
mkdir -p "$TD/bin"
if [ -z "$FC" ] && command -v ifort >/dev/null 2>&1; then
   FC=ifort
fi
if [ "$(basename "$FC" 2>/dev/null)" = "ifort" ]; then
   FLAGS="-mkl=sequential -O0 -g -traceback -diag-disable=10448,10412,10434 -module $TD/bin"
   LIBS=""
else
   FC="${FC:-gfortran}"
   FLAGS="-O0 -g -fcheck=all -fbacktrace -ffree-line-length-none -fno-range-check -J$TD/bin"
   if ldconfig -p 2>/dev/null | grep -q 'liblapack'; then
      LIBS="-llapack -lblas"
   else
      LIBS="-l:libopenblas.so.0"
   fi
fi
$FC $FLAGS \
   "$TD/jfnk_int_mod_stub.f90" \
   "$TREE/new_main/mod_subs/jfnk_solver_v1.f90" \
   "$TD/se_prob_mod.f90" \
   "$TD/se_driver.f90" \
   $LIBS \
   -o "$TD/bin/se_driver"
echo "built $TD/bin/se_driver"
