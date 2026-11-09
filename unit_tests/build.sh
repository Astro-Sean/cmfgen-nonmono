#!/bin/bash
# Build the surrogate-SE CLI driver against the exact production
# solver objects: jfnk_solver_v1.o and jfnk_op_probe_v1.o are the
# lib/libmod.a members the production executable links (resolved and
# freshness-checked by unit_tests/common.py).  Test-side sources
# (stub interface module, SE surrogate, driver) compile fresh and
# link with the production link line.  Requires ifort and a
# completed production build (make DEBUG=nodebug all / ./build.sh).
set -e
TD="$(cd "$(dirname "$0")" && pwd)"
python3 - "$TD" <<'EOF'
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from common import TREE, UT, compile_and_link_prod

compile_and_link_prod(
    "se_driver",
    [
        UT / "jfnk_int_mod_stub.f90",
        UT / "se_prob_mod.f90",
        UT / "se_driver.f90",
        TREE / "new_main" / "mod_subs" / "jfnk_solver_v1.f90",
        TREE / "new_main" / "mod_subs" / "jfnk_op_probe_v1.f90",
    ],
)
print("built unit_tests/bin/se_driver")
EOF
