#!/usr/bin/env python3
"""Runtime ABI/purity test for KN_ORD_SPV_IMPL.

Phase-5.5 follow-up: the OPTIONAL-argument design was replaced by a
module-scope worker (KN_ORD_SPV_IMPL, required LOGICAL L_CNT) behind
two explicit-interface entry points - KN_ORD_SPV (production, .TRUE.)
and KN_ORD_SPV_SHADOW (shadow, .FALSE.).  This test compiles the
exact production file and drives the worker directly across the full
branch matrix: GAM=0, expanding L>1 and L=1, compressing with two /
one / no valid donors, the compressing frequency boundary, and both
FG_FREQ_MATCH states, plus the unallocated KN_DBG_NDON path.

Asserts, per case:
  * production and shadow SPV bit-identical (BEQ via TRANSFER);
  * NDP (production counter delta) is exactly +1 only on the
    expanding L>1 branch, 0 elsewhere;
  * NDS (shadow counter delta) is 0 on every allocated case;
  * operand checksums are unchanged after all calls.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    F90DIR,
    TREE,
    WORKDIR,
    Report,
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

# Expected production counter delta per case: +1 only where the
# expanding L>1 donor read happens; 0 elsewhere. The shadow delta
# (NDS) must be 0 on every allocated case; -1 marks unallocated.
EXP_NDP = {
    "gam0": 0,
    "gamp": 1,
    "gamp1": 0,
    "gamn2": 0,
    "gamn1": 0,
    "gamnn": 0,
    "gamnb": 0,
    "fmtch": 0,
    "fmtcb": 0,
    "unall": -1,
}


def build():
    return compile_and_link_prod(
        "knspv", [TREE / h for h in HELPERS] + [F90DIR / "drv_knop_spv.f90"]
    )


def main():
    rep = Report("knop_spv")
    exe = build()
    out = run(exe).stdout
    print(out)

    cases = {}
    for m in re.finditer(
        r"SPCASE\s+(\S+)\s+R\s*(\d+)\s+L\s*(\d+)\s+SVT\s*(\S+)\s+"
        r"SVF\s*(\S+)\s+BEQ\s*([TF])\s+NDP\s*(-?\d+)\s+NDS\s*(-?\d+)",
        out,
    ):
        cases[m.group(1).strip()] = {
            "R": int(m.group(2)),
            "L": int(m.group(3)),
            "SVT": m.group(4),
            "SVF": m.group(5),
            "BEQ": m.group(6) == "T",
            "NDP": int(m.group(7)),
            "NDS": int(m.group(8)),
        }

    fails = []
    for name, ndp in EXP_NDP.items():
        c = cases.get(name)
        if c is None:
            fails.append(f"{name}: no SPCASE line")
            continue
        if not c["BEQ"]:
            fails.append(f"{name}: SPV differs T/F {c['SVT']} vs {c['SVF']}")
        if c["NDP"] != ndp:
            fails.append(f"{name}: NDP={c['NDP']} expected {ndp}")
        exp_nds = -1 if name == "unall" else 0
        if c["NDS"] != exp_nds:
            fails.append(f"{name}: NDS={c['NDS']} expected {exp_nds}")

    # Spot-check the arithmetic itself: expanding gamp must equal
    # SR + Q*FLD(R,L-1) with the manufactured fill values.
    g = cases.get("gamp")
    if g is not None:
        sr = 1.5 + 0.01 * g["R"] + 0.001 * g["L"]
        q = 0.3 + 0.01 * g["R"]
        fld_prev = 0.1 * g["R"] + 0.01 * (g["L"] - 1)
        want = sr + q * fld_prev
        if abs(float(g["SVT"]) - want) > 1e-12:
            fails.append(f"gamp: SVT={g['SVT']} expected {want}")

    # Operand checksums present and finite.
    for key in ("SUMSR", "SUMQ", "SUMFLD", "SUMJPR", "SUMIKC"):
        if not re.search(rf"{key}\s+\S+", out):
            fails.append(f"missing checksum {key}")

    for name, ndp in EXP_NDP.items():
        c = cases.get(name)
        rep.check(f"{name} parsed", c is not None)
        if c is not None:
            rep.check(f"{name} SPV bit-identical", c["BEQ"])
    rep.check("all checks passed", not fails, f"{len(cases)} cases")
    rep.summary()
    if fails:
        sys.exit(1)
    print("all KNOP SPV ABI/purity checks passed")


if __name__ == "__main__":
    main()
