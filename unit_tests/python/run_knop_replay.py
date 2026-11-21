#!/usr/bin/env python3
"""Phase 2.22 capped-ray replay runner.

Builds bin/knray_replay (drv_knop_replay.f90 + the same production
helper set as the knray driver) and runs it on a KNRDUMP deck
written by the armed production solve.  The deck is produced by a
model run with VADAT keys KN_RDUMP_LS / KN_RDUMP_LS2 (default -1
= off) -- see docs/audits/CMFGEN_KNOP_PHASE222_REFRESH_STALL_AUDIT.md.

Usage:
  python3 unit_tests/python/run_knop_replay.py <deck.bin> [mode] [value]

Modes (passed through to the driver): base | loop <n> | om <v> |
btol <x> | stall <x> | maxp <n> | seed cold|scale | verify
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from common import (  # noqa: E402
    F90DIR,
    TREE,
    WORKDIR,
    compile_and_link_prod,
    run,
)

HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "unix/gen_asci_open.f",
    "unix/verbose_info.f",
    "subs/set_case.f",
    "subs/phys_con.f",
    "subs/mon_int_funs_v2.f",
    "subs/mon_interp.f",
    "subs/mod_trap_derivatives.f",
    "subs/nordtau.f",
    "subs/qfromf.f",
    "subs/integrate.f",
    "subs/write_vec.f",
    "unix/dir_acc_pars_gen.f",
    "plane/mod_ray_mom_store.f",
    "new_main/mod_subs/knop_status_mod.f90",
    "plane/get_moms_non_rel.f",
    "new_main/mod_subs/mod_cmfgen.f",
    "new_main/mod_subs/control_variable_mod.f",
    "new_main/mod_subs/jfnk_cmfgen_int_mod.f",
    "new_main/mod_subs/jfnk_solver_v1.f90",
    "new_main/mod_subs/check_param_consistency.f",
    "newsubs/rd_store_log.f",
    "new_main/mod_subs/iter_ledger_mod.f90",
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
    "unit_tests/f90/stub_vadat.f",
]


def build():
    return compile_and_link_prod(
        "knray_replay",
        [TREE / h for h in HELPERS] + [F90DIR / "drv_knop_replay.f90"],
    )


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: run_knop_replay.py <deck.bin> [mode] [value]")
    deck = Path(sys.argv[1]).resolve()
    args = [str(deck)] + sys.argv[2:]
    bin_ = build()
    r = run(bin_, *args, cwd=WORKDIR, check=False)
    sys.stdout.write(r.stdout)
    sys.stderr.write(r.stderr)


if __name__ == "__main__":
    main()
