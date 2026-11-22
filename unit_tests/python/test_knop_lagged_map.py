#!/usr/bin/env python3
"""Lagged-field map tests for the Knop refresh (fg_j_cmf_nonmono_v1.f).

The retained I_*_DEPTH_NEXT store is iterated by two nested loops:
the in-sweep relaxed writes and the post-sweep FG_KNOP_REFRESH_V1
re-solve that commits candidate fields back into the store.  Running
the refresh N times on a fixed captured problem steps the outer
lagged map directly -- the committed-store trajectory measures the
map's contraction, and a committed FULL field is (by contract) a
fixed point of the per-ray solver regardless of the seed.

Cases (drv_knop_ray.f90):
  lagc  cold store (all zero) x 3 refreshes, one compressing band
  lagw  warm smooth seed x 3 refreshes, same band (seed A/B)
  lagr  cold store x 3 refreshes, two-band reversal profile
  lags  cold store x 3 refreshes, strong-compression band
  sdab  warm seed refresh, then a garbage seed refresh (seed B)
  sdp   poisoned store seed (NaN/Inf/negative/huge) x 1 refresh
  qres  KN_Q_RESID=T x 1 refresh -> independent fixed-point residual

The committed store must always be finite, certifications must be
honest (no FULL claim while the independent residual fails), and
fields committed as converged from different seeds must agree --
the lagged map is only a solver if its fixed point is unique.
"""

import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, TREE, WORKDIR, Report, compile_and_link_prod, run

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
        "knray", [TREE / h for h in HELPERS] + [F90DIR / "drv_knop_ray.f90"]
    )


def _f(x):
    # ES16.8 drops the E on 3-digit exponents: '1.79+308'.
    if "E" not in x:
        x = x[:-4] + "E" + x[-4:]
    return float(x)


def parse_refreshes(out):
    """Split driver output into per-refresh records.

    Each refresh prints RF_PASS i, then the refresh's own LUER
    summary, RF_* metrics, VC/CLS/SEED lines and the RSTORE store
    dump.  The KN_Q_RESID diagnostic reports through the LUER line
    'independent q residual max=' (KN_Q_LS is deallocated before the
    caller returns).  Returns one dict per refresh.
    """
    refreshes = []
    cur = {"rstore": {}, "cls": {}, "vc": {}, "qmax": None, "qn": 0}
    saw_pass = False
    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "RF_PASS":
            refreshes.append(cur)
            cur = {"rstore": {}, "cls": {}, "vc": {}, "qmax": None, "qn": 0}
            saw_pass = True
            continue
        mq = re.search(r"independent q residual max=\s*(\S+)", ln)
        if mq:
            cur["qmax"] = _f(mq.group(1))
        mq = re.search(r"rays reported=\s*(\d+)", ln)
        if mq:
            cur["qn"] = int(mq.group(1))
        m = re.match(r"^(\w+)=\s*(\S+)", ln)
        if m:
            try:
                cur[m.group(1)] = float(m.group(2))
            except ValueError:
                try:
                    cur[m.group(1)] = _f(m.group(2))
                except (ValueError, IndexError):
                    pass
            continue
        if t[0] == "RSTORE":
            cur["rstore"][(int(t[1]), int(t[2]), int(t[3]))] = (
                float(t[4]),
                float(t[5]),
            )
        elif t[0] == "CLS":
            cur["cls"][int(t[1])] = (int(t[2]), int(t[3]))
        elif t[0] == "VC":
            cur["vc"][int(t[1])] = (int(t[2]), _f(t[3]))
    refreshes.append(cur)
    if saw_pass and not refreshes[0]["rstore"] and not refreshes[0]["cls"]:
        refreshes = refreshes[1:]
    return refreshes


def _finite(v):
    return math.isfinite(v)


def store_finite(rec):
    return all(_finite(p) and _finite(m) for p, m in rec["rstore"].values())


def store_rel_diff(a, b):
    """Max |a-b|/max(|b|, tiny) over shared committed cells."""
    keys = set(a["rstore"]) & set(b["rstore"])
    worst = 0.0
    for k in keys:
        for j in (0, 1):
            av, bv = a["rstore"][k][j], b["rstore"][k][j]
            worst = max(worst, abs(av - bv) / max(abs(bv), 1e-30))
    return worst


def main():
    rep = Report("KNOP lagged-field map")
    bin_ = build()

    runs = {}
    for case in ("lagc", "lagw", "lagr", "lags", "sdab", "sdp", "qres"):
        r = run(bin_, case, cwd=WORKDIR, check=False)
        runs[case] = r
        rep.check(
            f"{case}: driver completed",
            r.returncode == 0 and "RESULT_OK" in r.stdout,
            f"rc={r.returncode}",
        )

    refr = {c: parse_refreshes(r.stdout) for c, r in runs.items()}

    # -- T1: every committed store is finite, cold or warm, on the
    # single-band profile and on the complementary reversal and
    # strong-compression profiles. Poison never passes through.
    # sdab refresh 2 is exempt: it exposes the retained-poison defect
    # documented under T4 and is checked separately there.
    for case in ("lagc", "lagw", "lagr", "lags", "sdp"):
        allfin = all(store_finite(rec) for rec in refr[case])
        rep.check(
            f"{case}: committed store finite on every refresh",
            allfin,
            f"nrf={len(refr[case])}",
        )

    # -- T2: certification classes in range. CLS describes the solve,
    # VC the commit decision -- a certified solve can still be
    # commit-vetoed (VC=1/2/3), so no cross-consistency is claimed.
    for case in ("lagc", "lagw", "lagr", "lags", "sdab", "sdp", "qres"):
        ok = True
        for rec in refr[case]:
            for ls, (cls, _rea) in rec["cls"].items():
                if cls < 0 or cls > 5:
                    ok = False
        rep.check(f"{case}: certification classes in range", ok)

    # -- T3: cold vs warm seeds converge to the same committed field
    # on the last refresh when both certify. The lagged map's fixed
    # point must not depend on the initial generation-0 state.
    last_c, last_w = refr["lagc"][-1], refr["lagw"][-1]
    full_c = all(v[0] == 0 for v in last_c["cls"].values())
    full_w = all(v[0] == 0 for v in last_w["cls"].values())
    d_cw = store_rel_diff(last_c, last_w)
    rep.check(
        "lagc-vs-lagw: final committed stores agree to 20%",
        (full_c and full_w and d_cw < 0.2) or not (full_c and full_w),
        f"full_c={full_c} full_w={full_w} d={d_cw:.3e}",
    )

    # -- T4: sdab seed A/B. Refresh 1 commits a converged field from
    # the warm seed. The store is then re-poisoned (garbage + NaN +
    # Inf) and refresh 2 converges the ray again -- but the commit
    # jump veto compares against the poisoned retained generation and
    # rejects, so the garbage stays in the store. HEAL_STORE only runs
    # on diverged rays, so this door has no recovery path: the slot is
    # stuck on every later refresh. Documented defect, kept failing:
    # a converged certified field must be able to replace a poisoned
    # retained generation.
    ra, rb = refr["sdab"][0], refr["sdab"][1]
    full_a = all(v[0] == 0 for v in ra["cls"].values())
    full_b = all(v[0] == 0 for v in rb["cls"].values())
    d_ab = store_rel_diff(ra, rb)
    rep.check(
        "sdab: converged field replaces poisoned retained store "
        "[KNOWN DEFECT: jump veto keeps garbage]",
        store_finite(rb) and full_b and d_ab < 0.05,
        f"full_a={full_a} full_b={full_b} "
        f"vc_b={rb['vc']} d={d_ab:.3e}",
    )

    # -- T5: qres independent fixed-point residual -- reported for
    # every ray, finite, non-negative, and <= KN_RAY_TOL (2e-3) when
    # the ray certified FULL.
    qrec = refr["qres"][-1]
    all_full = all(v[0] == 0 for v in qrec["cls"].values())
    q_ok = (
        qrec["qn"] >= 1
        and qrec["qmax"] is not None
        and _finite(qrec["qmax"])
        and qrec["qmax"] >= 0.0
        and (not all_full or qrec["qmax"] <= 2.0e-3)
    )
    rep.check(
        "qres: independent residual reported, finite, consistent",
        q_ok,
        f"qmax={qrec['qmax']} n={qrec['qn']}",
    )

    # -- T6: outer-map trajectory is bounded -- the committed-store
    # change between consecutive refreshes never grows without bound
    # (divergence would mean the lagged map amplifies the seed error).
    for case in ("lagc", "lagw", "lagr", "lags"):
        diffs = [
            store_rel_diff(refr[case][i - 1], refr[case][i])
            for i in range(1, len(refr[case]))
        ]
        bounded = all(math.isfinite(d) and d < 1.0e6 for d in diffs)
        rep.check(
            f"{case}: outer-map committed-store trajectory bounded",
            bounded,
            " ".join(f"{d:.2e}" for d in diffs),
        )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
