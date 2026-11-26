"""WP-5 signed-Aitken store relaxation (KN_OMEGA_SGN).

drv_knop_osgn.f90 exercises the estimator at two levels:

  scalar  KN_OSGN_STEP on a one-dimensional lagged map
          x_{k+1} = x_k + omega_k * (a + lam*(x_k - a) - x_k).
          A linear map hands the kernel the fixed-point omega
          1/(1-lam) from one (dx,dr) pair.

  vector  the module arrays: FG_OSGN_ALLOC, column-buffer pokes
          replicating the sweep write-site captures,
          FG_OSGN_COL_FIN at column close, FG_OSGN_UPDATE at the
          refresh tail.

Audit: docs/audits/KNOP_REMEDIATION_WP5_2026-11-26.md.
"""

import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import SRC, F90DIR, TREE, WORKDIR, Report, compile_and_link_prod, run

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
]


def _f(x):
    if "E" not in x:
        x = x[:-4] + "E" + x[-4:]
    return float(x)


def parse_out(text):
    out = {}
    for key in (
        "OSGN_OMEGA",
        "OSGN_RESID",
        "OSGN_RC2",
        "OSGN_RP2",
        "OSGN_OMEGA_HI",
        "OSGN_OMEGA_LO",
        "OSGN_OMEGA2",
    ):
        m = re.search(rf"{key}=\s*(\S+)", text)
        out[key.lower()] = _f(m.group(1)) if m else None
    for key in ("OSGN_ICODE", "OSGN_ICODE_HI", "OSGN_ICODE_LO",
                "OSGN_ICODE2", "OSGN_NSTEP", "OSGN_NRET"):
        m = re.search(rf"{key}=\s*(\S+)", text)
        out[key.lower()] = int(m.group(1)) if m else None
    m = re.search(r"OSGN_LFALL2=\s*(\S+)", text)
    out["osgn_lfall2"] = m.group(1).strip() if m else None
    out["steps"] = [
        (int(m[0]), int(m[1]), _f(m[2]))
        for m in re.findall(
            r"OSGN_STEP\s+(\d+)\s+ICODE=\s*(-?\d+)\s+OMEGA=\s*(\S+)",
            text,
        )
    ]
    out["sweeps"] = [
        (int(m[0]), int(m[1]), _f(m[2]), m[3])
        for m in re.findall(
            r"OSGN_SWEEP\s+(\d+)\s+ICODE=\s*(-?\d+)\s+OMEGA=\s*(\S+)"
            r"\s+LFALL=\s*(\S+)",
            text,
        )
    ]
    out["ok"] = "RESULT_OK" in text
    out["text"] = text
    return out


def main():
    rep = Report("KNOP signed-Aitken store relaxation (WP-5)")
    bin_ = compile_and_link_prod(
        "knosgn",
        [SRC / h for h in HELPERS]
        + [F90DIR / "stub_vadat.f", F90DIR / "drv_knop_osgn.f90"],
    )

    def case(name):
        return parse_out(run(bin_, name).stdout)

    # Contracting map lam=0.3: the relaxed iteration converges to
    # the root and omega lands on the fixed-point value 10/7.
    d = case("sc")
    rep.check(
        "sc: contracting map reaches the root",
        d["ok"]
        and d["osgn_resid"] is not None
        and d["osgn_resid"] < 1.0e-12
        and d["osgn_icode"] == 1,
        d["text"][-300:] if not d["ok"] else "",
    )
    rep.check(
        "sc: omega lands on 1/(1-lam) = 10/7",
        d["osgn_omega"] is not None
        and abs(d["osgn_omega"] - 10.0 / 7.0) < 1.0e-8,
    )

    # Oscillating map lam=-0.8: sign-flipped residual, same
    # fixed-point omega 1/(1-lam) = 5/9.
    d = case("so")
    rep.check(
        "so: oscillating map reaches the root",
        d["ok"]
        and d["osgn_resid"] is not None
        and d["osgn_resid"] < 1.0e-12
        and d["osgn_icode"] == 1,
        d["text"][-300:] if not d["ok"] else "",
    )
    rep.check(
        "so: omega lands on 5/9",
        d["osgn_omega"] is not None
        and abs(d["osgn_omega"] - 5.0 / 9.0) < 1.0e-8,
    )

    # Expanding map lam=1.5 under omega0=0.5: residual grows, the
    # kernel halves once (ICODE=2) then spends the retry budget
    # (ICODE=3) so the caller drops to the unaccelerated path.
    d = case("sd")
    rep.check(
        "sd: growing residual halves once then falls back",
        d["ok"]
        and len(d["steps"]) == 2
        and d["steps"][0][1] == 2
        and d["steps"][1][1] == 3
        and abs(d["steps"][0][2] - 0.25) < 1.0e-12,
        d["text"][-300:] if not d["ok"] else "",
    )

    # Collapsed denominator: dr ~ 0 must hold omega, not divide.
    d = case("sc0")
    rep.check(
        "sc0: collapsed denominator holds omega",
        d["ok"]
        and d["osgn_icode"] == 0
        and abs(d["osgn_omega"] - 0.5) < 1.0e-12,
        d["text"][-300:] if not d["ok"] else "",
    )

    # Clamps: lam=0.55 wants 20/9 > KN_OSGN_MAX; lam=-1.5 wants
    # 0.4 < the 0.5 floor the driver passes.
    d = case("scl")
    rep.check(
        "scl: over-omega clamps at KN_OSGN_MAX",
        d["ok"]
        and d["osgn_icode_hi"] == 1
        and abs(d["osgn_omega_hi"] - 2.0) < 1.0e-12,
        d["text"][-300:] if not d["ok"] else "",
    )
    rep.check(
        "scl: under-omega clamps at the floor",
        d["ok"]
        and d["osgn_icode_lo"] == 1
        and abs(d["osgn_omega_lo"] - 0.5) < 1.0e-12,
        d["text"][-300:] if not d["ok"] else "",
    )

    # Vector machinery: first sweep has no history (ICODE=0
    # hold), the second prices omega to 10/7, the third confirms.
    d = case("vec")
    rep.check(
        "vec: sweep 1 holds omega (no history)",
        d["ok"]
        and len(d["sweeps"]) == 3
        and d["sweeps"][0][1] == 0,
        d["text"][-300:] if not d["ok"] else "",
    )
    rep.check(
        "vec: sweep 2 prices omega to 10/7",
        d["ok"]
        and len(d["sweeps"]) == 3
        and d["sweeps"][1][1] == 1
        and abs(d["sweeps"][1][2] - 10.0 / 7.0) < 1.0e-8,
        d["text"][-300:] if not d["ok"] else "",
    )
    rep.check(
        "vec: two-slot map converges to both roots",
        d["ok"]
        and d["osgn_resid"] is not None
        and d["osgn_resid"] < 1.0e-12,
        d["text"][-300:] if not d["ok"] else "",
    )

    # Fingerprint move: a coefficient change resets the slot
    # history, so the same-sweep update sees an empty residual
    # and holds omega instead of consuming a cross-map dx.
    d = case("vfp")
    rep.check(
        "vfp: fingerprint move holds omega on the changed sweep",
        d["ok"]
        and len(d["sweeps"]) == 3
        and d["sweeps"][2][1] == 0,
        d["text"][-300:] if not d["ok"] else "",
    )
    rep.check(
        "vfp: omega survives the reset unchanged",
        d["ok"]
        and d["osgn_omega"] is not None
        and abs(d["osgn_omega"] - 10.0 / 7.0) < 1.0e-8,
        d["text"][-300:] if not d["ok"] else "",
    )

    # Generation guard: a second update inside one sweep must
    # fall back -- the residual pair was already consumed.
    d = case("vgen")
    rep.check(
        "vgen: second update in one generation falls back",
        d["ok"]
        and d["osgn_icode2"] == -1
        and d["osgn_lfall2"] == "T",
        d["text"][-300:] if not d["ok"] else "",
    )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
