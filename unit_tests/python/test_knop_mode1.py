#!/usr/bin/env python3
"""MODE 1 diagnostic-probe isolation tests (KN_NOCP_PROBE).

Drives drv_knop_ray cases that arm ray 3 as the non-committing
diagnostic probe on a three-ray problem identical to 'rf3' and
asserts that every controller-facing export stays on the
production rays while the probe keeps its own telemetry:

  punit  unit checks on KN_PROD_RESID_V1 / KN_OMEGA_POOL_V1 EXCL /
         KN_NOCP_PROBE_RAY with synthetic arrays
  prb    probe store x2.5, converges: production exports = rf3
  prbv   probe store x1E-4, real jump veto: counters stay zero
  prbs   probe + KN_SEED_ARM=1: production exports = rf3

Baseline is 'rf3' (three identical rays, no probe). Isolation is
proven on the underlying values, not the printed summary alone.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import WORKDIR, Report, run

from test_knop_ray import build, parse


def parse_extra(out):
    """Pull the probe-only lines the generic parse() skips."""
    x = {}
    m = re.search(r"probe ray=\s*(\d+)\s+resid=\s*(\S+)\s+div=\s*(\d+)\s+nbad=\s*(\d+)", out)
    if m:
        x["PROBE_LS"] = int(m.group(1))
        x["PROBE_RESID"] = float(m.group(2))
        x["PROBE_NDIV"] = int(m.group(3))
        x["PROBE_NBAD"] = int(m.group(4))
    m = re.search(r"would-be verdict=\s*(\d+)", out)
    if m:
        x["PROBE_WBV"] = int(m.group(1))
    m = re.search(r"KN-NOP store ray=.*delta M/P=\s*(\S+)\s+(\S+)", out)
    if m:
        x["NOP_DM"] = float(m.group(1))
        x["NOP_DP"] = float(m.group(2))
    m = re.search(r"KN-SEED ray=\s*(\d+)\s+arm=\s*(\d+)\s+used=\s*(\d+)", out)
    if m:
        x["SEED_ARM"] = int(m.group(2))
        x["SEED_USED"] = int(m.group(3))
    return x


def case(bin_, name):
    r = run(bin_, name, cwd=WORKDIR, check=False)
    d, store = parse(r.stdout)
    d.update(parse_extra(r.stdout))
    return d, store


PROD_KEYS = [
    "RF_RESID", "RF_DCHG", "RF_PEND", "RF_NDIV", "RF_MAXCHG",
    "RF_STAL", "RF_PLAT", "RF_INEXT_RELAX", "RF_HEAL", "RF_JREJ",
    "RF_NCREJ", "RF_NCCAP", "RF_NBAD", "RF_CAPPED", "RF_ICAP",
    "RF_FAIL",
]


def prod_equal(d, base):
    return all(d.get(k) == base.get(k) for k in PROD_KEYS)


def main():
    rep = Report("KNOP_MODE1_PROBE_ISOLATION")
    bin_ = build()

    # -- M1: helper unit checks on synthetic arrays ------------------
    d, _ = case(bin_, "punit")
    ok = "RESULT_OK" in d
    exp = {
        "M1U_RES_PX": 0.3,       # probe excluded from max
        "M1U_RES_P0": 5.0,       # XRAY=0 reproduces old pool
        "M1U_RES_LOW": 0.3,
        "M1U_RES_ALLDV": 1.0e30,  # all-diverged fallback
        "M1U_RES_ALLDVP": 8.0,   # diverged production max, not probe
        "M1U_RES_DVIST": 1.0e30,
        # degenerate all-probe run: fixed fail-closed sentinel,
        # independent of the probe residual value entirely
        "M1U_RES_SOLO": 1.0,
        "M1U_RES_SOLO2": 1.0,
        "M1U_RES_SOLO3": 1.0,
        "M1U_IST_SOLO": 2,
        "M1U_IST_SOLO2": 2,
        "M1U_IST_ALLDV": 1,
        # empty production omega pool retains OMEGA_PREV whatever
        # the probe rho looks like
        "M1U_OM_NOPROD": 0.75,
        "M1U_OM_NOPROD2": 0.75,
        "M1U_NSAM_NOPROD": 0,
        "M1U_OM_BASE": None,
        "M1U_RHOP_DIAG": 1.5,    # diagnostic rho preserved
    }
    for k, v in exp.items():
        if v is not None:
            ok = ok and abs(d.get(k, -1e9) - v) < 1e-12 * max(1.0, abs(v))
    ok = (
        ok
        and abs(d.get("M1U_OM_EXHI", 0) - d["M1U_OM_BASE"]) < 1e-15
        and abs(d.get("M1U_OM_EXLO", 0) - d["M1U_OM_BASE"]) < 1e-15
        and d.get("M1U_NSAM_EXHI") == d.get("M1U_NSAM_BASE") == 40
        and d.get("M1U_NSAM_INCL") == 42
        and abs(d.get("M1U_OM_FEW", 0) - 0.75) < 1e-15
        and d.get("M1U_NSAM_FEW") == 8
    )
    rep.check(
        "M1 helper unit checks: residual pool, omega EXCL, fallbacks",
        ok,
        f"punit={ {k: d.get(k) for k in sorted(d) if k.startswith('M1U')} }",
    )

    # -- M2: converging probe, production exports identical to rf3 ---
    db, sb = case(bin_, "rf3")
    d, s = case(bin_, "prb")
    ok = "RESULT_OK" in d and prod_equal(d, db)
    # probe's own cells leave the denominator (90 per ray, 2 prod rays)
    ok = (
        ok
        and d.get("RF_NTOT") == db["RF_NTOT"] - 90
        and d.get("RF_NTB") == db["RF_NTB"] - 90
        and d.get("PROBE_LS") == 3
        and d.get("PROBE_RESID") is not None
        and s.get(("vc", 3), (0,))[0] == 7
        and abs(d.get("NOP_DM", 1)) == 0.0
        and abs(d.get("NOP_DP", 1)) == 0.0
    )
    rep.check(
        "M2 prb: production exports equal rf3, probe telemetry kept",
        ok,
        f"RF_RESID={d.get('RF_RESID'):.3e} vs rf3 {db.get('RF_RESID'):.3e} "
        f"NTOT={d.get('RF_NTOT')}/{db.get('RF_NTOT')} "
        f"probe resid={d.get('PROBE_RESID')}",
    )
    # committed store slices for the production rays match rf3
    dd = 0.0
    for ls in (1, 2):
        for k in range(1, 4):
            for l in range(1, 6):
                a = s.get(("r", ls, k, l))
                b = sb.get(("r", ls, k, l))
                if a and b:
                    dd = max(dd, abs(a[0] - b[0]), abs(a[1] - b[1]))
    rep.check(
        "M3 prb: production ray stores bitwise-unchanged vs rf3",
        dd == 0.0,
        f"max abs store diff {dd:.3e}",
    )
    # probe slice untouched by the refresh: it still holds the
    # perturbed SEED fill (2.5 x the generic fill), not a commit
    # of its solved field (rf3's slice was overwritten by a solve).
    def _fill(ls_p, ls_m, k, l):
        fp = 2.5 * (2.0e-4 * (1 + 0.1 * l) * (1 + 0.1 * k))
        fm = 2.5 * (1.0e-4 * (1 + 0.1 * l) * (1 + 0.1 * k))
        return abs(ls_p - fp) < 1e-8 * fp and abs(ls_m - fm) < 1e-8 * fm

    ok = all(
        _fill(*s.get(("r", 3, k, l), (0.0, 0.0)), k, l)
        for k in range(1, 4)
        for l in range(1, 6)
    )
    rep.check(
        "M4 prb: probe store slice retains perturbed seed (no commit)",
        ok,
        "ray3 RSTORE = 2.5 x seed fill",
    )

    # -- M5: probe jump veto stays off the production counters -------
    d, s = case(bin_, "prbv")
    ok = (
        "RESULT_OK" in d
        and d.get("PROBE_WBV") == 1
        and d.get("RF_JREJ") == 0
        and d.get("RF_NCREJ") == 0
        and d.get("RF_NCCAP") == 0
        and d.get("RF_HEAL") == 0
        and d.get("RF_PEND") == 0.0
        and prod_equal(d, db)
        and s.get(("vc", 3), (0,))[0] == 7
    )
    rep.check(
        "M5 prbv: real probe veto, production counters stay zero",
        ok,
        f"wbv={d.get('PROBE_WBV')} JREJ={d.get('RF_JREJ')} "
        f"PEND={d.get('RF_PEND')} vc3={s.get(('vc', 3))}",
    )

    # -- M6: KN_SEED_ARM=1 leaves production exports identical -------
    d, s = case(bin_, "prbs")
    ok = (
        "RESULT_OK" in d
        and d.get("SEED_ARM") == 1
        and d.get("SEED_USED") == 1
        and prod_equal(d, db)
        and d.get("RF_NTOT") == db["RF_NTOT"] - 90
        and s.get(("vc", 3), (0,))[0] == 7
    )
    rep.check(
        "M6 prbs: sibling-median probe seed, production unchanged",
        ok,
        f"arm={d.get('SEED_ARM')} RF_RESID={d.get('RF_RESID'):.3e} "
        f"vs rf3 {db.get('RF_RESID'):.3e}",
    )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
