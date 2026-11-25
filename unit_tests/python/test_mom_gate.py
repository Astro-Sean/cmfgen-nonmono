#!/usr/bin/env python3
"""Phase-2 FIX 1 / FIX 2 units for the moment refresh.

FIX 1 (audit finding F): MOM_SAME_DIR_RHO must estimate the
contraction factor from same-direction residual pairs only. The
manufactured scalar-map sequences have known per-direction radii,
so the alternating schedule must return the composite factor
rho1*rho2, the single-direction and SKIP_D schedules the running
map's own factor, and a cross-map beat sequence must still read as
contracting where the old pass-ordered estimator read garbage.

FIX 2 (audit finding G): MOM_RF_COMMIT_GATE applies one acceptance
policy to the exact field about to be committed - finite entries,
no new non-positive J, the 50x jump guard vs the retained store,
and the KEEP_BEST lexicographic bulk bound referenced to the pre-AA
field. drv_mom_gate feeds it manufactured seed / pre-AA / candidate
triples and this battery checks each verdict and the warn text
naming the failing check.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    F90DIR,
    TREE,
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
    "new_main/mod_subs/iter_log_fmt_mod.f90",
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
    "subs/thomas.f",
    "plane/mom_j_cmf_v11.f",
]

RHO_CASES = ["rhoalt", "rhosingle", "rhodiv", "rhobeat",
             "rhoskip", "rhofew"]
GATE_CASES = ["clean", "aafix", "aaimp", "neg", "nonfin", "jump",
              "bulk", "seedgate"]


def build():
    return compile_and_link_prod(
        "mom_gate",
        [TREE / h for h in HELPERS]
        + [F90DIR / "stub_vadat.f", F90DIR / "drv_mom_gate.f"],
    )


def val(out, key):
    m = re.search(rf"^{key}=\s*(\S+)", out, re.M)
    return m.group(1) if m else None


def main():
    rep = Report("mom_gate")
    binary = build()
    rep.check("driver built", binary is not None and binary.exists())

    outs = {}
    for case in RHO_CASES + GATE_CASES:
        r = run(binary, case)
        outs[case] = r.stdout
        rep.check(f"{case}: driver ran", "STOP" not in r.stderr or
                  r.returncode == 0, detail=r.stdout[-400:])

    # FIX 1: same-direction rho estimates.
    rho = float(val(outs["rhoalt"], "RHO"))
    rep.check("rhoalt: composite rho = 0.5*0.9",
              abs(rho - 0.45) < 1e-12, detail=outs["rhoalt"])
    rep.check("rhoalt: enough ratios", int(val(outs["rhoalt"], "NRAT")) >= 3)
    rep.check("rhoalt: LOK", val(outs["rhoalt"], "LOK") == "T")

    rho = float(val(outs["rhosingle"], "RHO"))
    rep.check("rhosingle: single-map rho = 0.5",
              abs(rho - 0.5) < 1e-12, detail=outs["rhosingle"])

    rho = float(val(outs["rhodiv"], "RHO"))
    rep.check("rhodiv: divergent composite rho = 1.69 unclamped",
              abs(rho - 1.69) < 1e-12, detail=outs["rhodiv"])

    rho = float(val(outs["rhobeat"], "RHO"))
    rep.check("rhobeat: cross-map beat still reads contraction 0.5",
              abs(rho - 0.5) < 1e-12, detail=outs["rhobeat"])

    rho = float(val(outs["rhoskip"], "RHO"))
    rep.check("rhoskip: pooled median = 0.8 of surviving direction",
              abs(rho - 0.8) < 1e-12, detail=outs["rhoskip"])
    rep.check("rhoskip: 3 pooled ratios",
              int(val(outs["rhoskip"], "NRAT")) == 3)

    rep.check("rhofew: under-sampled -> LOK=F",
              val(outs["rhofew"], "LOK") == "F",
              detail=outs["rhofew"])

    # FIX 2: uniform commit gate.
    for case in ["clean", "aafix", "aaimp"]:
        rep.check(f"{case}: committed", val(outs[case], "KEEP") == "T",
                  detail=outs[case])
        rep.check(f"{case}: no gate warn",
                  "commit gate rejected" not in outs[case])
    for case, frag in [
        ("neg", "cells driven to J<=0"),
        ("nonfin", "non-finite or >1E40"),
        ("jump", "jumped more than 50x"),
        ("bulk", "bulk regression"),
        ("seedgate", "bulk regression"),
    ]:
        rep.check(f"{case}: rejected", val(outs[case], "KEEP") == "F",
                  detail=outs[case])
        rep.check(f"{case}: warn names {frag}",
                  frag in outs[case], detail=outs[case])
        rep.check(f"{case}: restore announced",
                  "pre-AA field restored" in outs[case],
                  detail=outs[case])

    rep.summary()


if __name__ == "__main__":
    main()
