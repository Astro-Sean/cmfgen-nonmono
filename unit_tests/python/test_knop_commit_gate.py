#!/usr/bin/env python3
"""Committed-state gate for the Knop non-monotonic solver (WP-3).

The maturity gate in SOLVE_FOR_POPS used to read the solve's own
best-pass residual (KN_RF_RESID -> GS_MAXCHG) and the MIN of the
last/best-pass bad censuses. A certified commit can still leave a
clipped or otherwise off-fixed-point field in the store -- the
population gate kept reading the field the solver *proposed*, not
the field it *kept*.

WP-3 [KN_GATE_COMMIT_RES] runs one non-committing L_QP pass on the
store each commit wrote and exports its max-norm as
FG_RF_COMMIT_RES plus a per-ray committed census
(FG_NBAD_CM_RAY/FG_NTOT_CM_RAY). The gate then folds COMMIT_RES and
the committed census into GS_MAX_ALL / GS_NBAD_EFF in place of the
pre-commit metrics.

Cases (drv_knop_ordering.f90, NI=140, NCF=200, compressing ray):
  capok   flag on, frequency-major march: the certified commit
          stores the converged field -> CRES ~ tol, clean census.
  nccm    flag on, ordered march: the stalled exit is not
          certified, so FG_COMMIT_NC writes a bounded provisional
          step from the 1e-6 seed. The store parks far off the
          fixed point and the committed verify reads >>tol with
          a nonzero committed census -- the WP-3 hold signature.
  nclg    same provisional commit with the flag off: no verify
          runs, so CRES and the census slots stay -1 -- the legacy
          gate cannot see what the store actually holds.
  lsok    WP-4 [KN_STORE_LS] on the nccm stall: the search's
          baseline (1e-6 seed) reads ~1e5, alpha=1 hands the store
          the fresh field in one step, and its committed residual
          drops to the solve's own level.
  lsrj    WP-4 rejection: solve 1 stalls and the LS parks fresh1;
          solve 2 over-relaxes (KN_RAY_OMEGA=1.36) and is capped
          after one march, so its iterate exits non-certified away
          from the fixed point. No blend descends below fresh1's
          residual -> the search rejects and keeps fresh1.

Assertions:
  * nccm: CRES > KN_RAY_TOL and CM_NB > 0 -- the committed-state
    gate holds on the field the store actually kept.
  * capok: CRES < KN_RAY_TOL and CM_NB == 0 -- an honest full
    commit still certifies.
  * nclg: CRES == -1 and CM_NB == -1 -- verify disarmed under
    the flag.
  * lsok: VC=3, LSA=1.0, LSN=1, LSR=0 -- full-step accept; CRES
    lands orders below the nccm baseline.
  * lsrj: VC=2, LSA=-2, LSR=1 -- every rung rejected; CRES still
    reads fresh1's committed residual (store retained).
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import SRC, F90DIR, WORKDIR, Report, compile_and_link_prod, run

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

KN_RAY_TOL = 2.0e-3


def _f(x):
    if "E" not in x:
        x = x[:-4] + "E" + x[-4:]
    return float(x)


def parse_out(text):
    out = {"text": text}
    out["ok"] = "RESULT_OK" in text
    for key in ("RESID", "CRES", "LSA"):
        m = re.search(rf"{key}=\s*(\S+)", text)
        out[key.lower()] = _f(m.group(1)) if m else None
    for key in ("CM_NB", "CM_NT", "VC", "NXIT", "LSN", "LSR"):
        m = re.search(rf"{key}=\s*(\S+)", text)
        out[key.lower()] = int(m.group(1)) if m else None
    return out


def main():
    rep = Report("KNOP committed-state gate (WP-3) + store line search (WP-4)")
    bin_ = compile_and_link_prod(
        "kncg",
        [SRC / h for h in HELPERS]
        + [F90DIR / "stub_vadat.f", F90DIR / "drv_knop_ordering.f90"],
    )

    data = {}
    for case in ("capok", "nccm", "nclg", "lsok", "lsrj"):
        r = run(bin_, case, "", cwd=WORKDIR, check=False)
        d = parse_out(r.stdout)
        data[case] = d
        rep.check(
            f"{case}: driver completed",
            r.returncode == 0 and d["ok"],
            f"rc={r.returncode}",
        )
        sys.stdout.write(
            f"    {case:<7s} resid={d['resid']:.3e} nxit={d['nxit']}"
            f" vc={d['vc']} cres={d['cres']:.3e}"
            f" cm_nb={d['cm_nb']} cm_nt={d['cm_nt']}\n"
        )

    # -- T1: the provisional NC commit parks the store far off the
    #    fixed point. The committed-state gate reads that directly:
    #    CRES far above tolerance and a nonzero committed census.
    d = data["nccm"]
    rep.check(
        "nccm: committed residual above tol (gate holds)",
        d["cres"] is not None and d["cres"] > KN_RAY_TOL,
        f"cres={d['cres']:.3e}",
    )
    rep.check(
        "nccm: committed census nonzero (cells written bad)",
        d["cm_nb"] is not None and d["cm_nb"] > 0
        and d["cm_nt"] is not None and d["cm_nt"] > 0,
        f"cm_nb={d['cm_nb']} cm_nt={d['cm_nt']}",
    )

    # -- T2: an honest commit still certifies. With the cap above the
    #    fixed point the store keeps the converged field, so the
    #    committed verify reads ~tol and the census stays clean.
    d = data["capok"]
    rep.check(
        "capok: committed residual below tol",
        d["cres"] is not None and 0.0 <= d["cres"] < KN_RAY_TOL,
        f"cres={d['cres']:.3e}",
    )
    rep.check(
        "capok: committed census empty (cells written good)",
        d["cm_nb"] is not None and d["cm_nb"] == 0
        and d["cm_nt"] is not None and d["cm_nt"] > 0,
        f"cm_nb={d['cm_nb']} cm_nt={d['cm_nt']}",
    )

    # -- T3: the flag disarms the verify entirely. The same
    #    provisional commit leaves CRES and the census slots at the
    #    -1 sentinel, while the write itself still happened (VC=3).
    d = data["nclg"]
    rep.check(
        "nclg: no committed residual (verify skipped)",
        d["cres"] is not None and d["cres"] < 0.0,
        f"cres={d['cres']:.3e}",
    )
    rep.check(
        "nclg: census slots never filled",
        d["cm_nb"] is not None and d["cm_nb"] < 0
        and d["cm_nt"] is not None and d["cm_nt"] < 0,
        f"cm_nb={d['cm_nb']} cm_nt={d['cm_nt']}",
    )
    rep.check(
        "nclg: provisional commit still wrote (VC=3)",
        d["vc"] is not None and d["vc"] == 3,
        f"vc={d['vc']}",
    )

    # -- T4: the LS hands the store the fresh field in one step
    #    where the fixed NC bound would ratchet over O(100)
    #    refreshes. alpha=1 descends, VC=3, one accept, no rejects.
    d = data["lsok"]
    rep.check(
        "lsok: alpha=1 accepted (LSA=1, VC=3)",
        d["lsa"] is not None and d["lsa"] == 1.0
        and d["vc"] is not None and d["vc"] == 3,
        f"lsa={d['lsa']} vc={d['vc']}",
    )
    rep.check(
        "lsok: one LS accept, no rejections",
        d["lsn"] is not None and d["lsn"] == 1
        and d["lsr"] is not None and d["lsr"] == 0,
        f"lsn={d['lsn']} lsr={d['lsr']}",
    )
    rep.check(
        "lsok: committed residual descends far below nccm baseline",
        d["cres"] is not None and data["nccm"]["cres"] is not None
        and 0.0 < d["cres"] < data["nccm"]["cres"] / 100.0,
        f"cres={d['cres']:.3e} nccm={data['nccm']['cres']:.3e}",
    )

    # -- T5: an expansive step cannot ratchet. Solve 1 parks fresh1;
    #    solve 2's over-relaxed march kicks the iterate off the fixed
    #    point, every rung reads worse than the retained store, the
    #    search rejects (LSA=-2, VC=2), and the committed residual
    #    still reports fresh1.
    d = data["lsrj"]
    rep.check(
        "lsrj: all rungs rejected (LSA=-2, VC=2)",
        d["lsa"] is not None and d["lsa"] == -2.0
        and d["vc"] is not None and d["vc"] == 2,
        f"lsa={d['lsa']} vc={d['vc']}",
    )
    rep.check(
        "lsrj: one accept (solve 1) then one rejection",
        d["lsn"] is not None and d["lsn"] == 1
        and d["lsr"] is not None and d["lsr"] == 1,
        f"lsn={d['lsn']} lsr={d['lsr']}",
    )
    rep.check(
        "lsrj: retained store still reads fresh1 residual",
        d["cres"] is not None and data["lsok"]["cres"] is not None
        and abs(d["cres"] - data["lsok"]["cres"]) < 1.0e-10,
        f"cres={d['cres']:.6e} lsok={data['lsok']['cres']:.6e}",
    )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
