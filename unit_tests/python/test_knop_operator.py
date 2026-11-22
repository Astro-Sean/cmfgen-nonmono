#!/usr/bin/env python3
"""Discrete-operator checks for the Knop formal solve.

The A/B parity tests in test_knop_ray.py prove the sweep and refresh
build the same operator. This battery checks the operator itself
against the mathematics it discretises:

  * A0..A4 march coefficients (refresh-side BUILD_FREQ dump, DBRC):
    A0 = exp(-DTAU) exactly, and A0+A1+A2 = 1 -- a field in local
    source equilibrium I = S must march to itself, so the constant-
    source weights must sum to 1 - A0. Both identities hold exactly
    in the production coefficient build (E0 = 1-EEL for DTAU>0.5,
    series evaluation below, closed form at DTAU>=40).
  * chi_hat construction: CHIH = chi_base + |GAM|/DLN at coupled
    cells of either sign and CHIH = chi_base at uncoupled cells --
    the "coupling or nothing" rule (KN_COUPLING_OK drops the
    implicit outflow together with the explicit inflow).
  * Q coefficient: Q = GAM/(DLN*CHIH) at coupled cells (both signs
    couple; only the donor field differs -- live previous-frequency
    for expansion, lagged store for compression), zero elsewhere ->
    Q*CHIH*DLNM reproduces the signed GAM.
  * SR = ETA/CHIH against the analytic ETA profile.
  * DTAU >= 0 and, on the pure-expansion profile, its depth sum
    matches the analytic integral of chi_base (R1/R)^2 dZ.
  * Solved field: finite, non-negative, and thermalises onto the
    modified source where the cells are optically thick.

Cases:
  gposo  all-expanding profile (NCASE=16, ordered sweep). No
         compressing cell exists, so chi_hat = chi_base everywhere
         and DTAU is analytically integrable.
  gcr    linear sign-reversal GAM(i) = GAMP*0.3*(i-4.5): cells 1-4
         compressing, 5-8 expanding, continuous through zero.
  gnego  all-compressing (NCASE=17, ordered) - the pure lagged-store
         regime: every donor read goes through I_*_DEPTH_NEXT.
  gzb    discrete combination (NCASE=31): expansion at cell 3, exact
         zero at cell 4, compression at cell 5.
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

# Manufactured model (must mirror drv_knop_ray.f90):
RR = (8.0e4, 6.4e4, 5.0e4, 3.9e4, 3.0e4, 2.2e4, 1.5e4, 1.0e4)
ND = 8
NCFX = 6
GAMP = 4.0e-6


def chi_base(i, l):
    return 1.0e-3 * (1.0 + 0.05 * l) * (RR[0] / RR[i - 1]) ** 2


def eta_base(i, l):
    return 2.0e-6 * (1.0 + 0.02 * l) * (RR[0] / RR[i - 1])


def gam_of(case, i):
    if case == "gposo":
        return GAMP * (0.5 + 0.1 * i)
    if case == "gcr":
        return GAMP * 0.3 * (i - 4.5)
    if case == "gnego":
        return -GAMP * (0.5 + 0.1 * i)
    if case == "gzb":
        if i == 3:
            return GAMP * (1.0 + 0.3 * i)
        if i == 4:
            return 0.0
        if i == 5:
            return -GAMP * (1.0 + 0.3 * i)
        return GAMP * (0.5 + 0.1 * i)
    raise ValueError(case)


# All-expanding (monotonic outflow), all-compressing (monotonic
# inflow - every donor read hits the lagged store), continuous sign
# reversal, and a discrete expansion/zero/compression combination.
CASES = ("gposo", "gcr", "gnego", "gzb")


def _f(x):
    # ES16.8 drops the E on 3-digit exponents: '1.79+308'.
    if "E" not in x:
        x = x[:-4] + "E" + x[-4:]
    return float(x)


_FLOAT = re.compile(r"[+-]?\d+\.\d+(?:E[+-]\d+|[+-]\d{3,})")


def parse_dbr(text):
    """{(i,l): dict} for DBR/DBRC rows; sentinel slots kept as -1."""
    dbr, dbrc = {}, {}
    for ln in text.splitlines():
        t = ln.split()
        if len(t) < 3 or t[0] not in ("DBR", "DBRC"):
            continue
        key = (int(t[1]), int(t[2]))
        val = [_f(x) for x in _FLOAT.findall(ln)]
        if t[0] == "DBR" and len(val) == 11:
            dbr[key] = dict(
                zip(
                    ("chih", "q", "dchr", "dtau", "spm", "spp",
                     "wx", "dlnm", "fldm", "fldp", "sr"),
                    val,
                )
            )
        elif t[0] == "DBRC" and len(val) == 6:
            dbrc[key] = dict(
                zip(("a0", "a1", "a2", "a3", "a4", "dtau"), val)
            )
    return dbr, dbrc


def main():
    rep = Report("KNOP discrete operator")
    bin_ = compile_and_link_prod(
        "knray", [TREE / h for h in HELPERS] + [F90DIR / "drv_knop_ray.f90"]
    )

    data = {}
    for case in CASES:
        r = run(bin_, case, cwd=WORKDIR, check=False)
        rep.check(
            f"{case}: driver completed and converged",
            r.returncode == 0 and "RESULT_OK" in r.stdout,
            f"rc={r.returncode}",
        )
        data[case] = parse_dbr(r.stdout)

    # -- T1: dump coverage and march-coefficient structure. A slots
    #    are filled on intervals 1..NI-1; i=NI keeps the -1 sentinel.
    for case in CASES:
        dbr, dbrc = data[case]
        ok = len(dbr) == ND * NCFX and len(dbrc) == ND * NCFX
        interior = True
        for (i, l), a in dbrc.items():
            if i < ND:
                interior &= all(
                    math.isfinite(a[k]) for k in ("a0", "a1", "a2", "a3", "a4")
                ) and a["dtau"] >= 0.0
            else:
                interior &= a["a0"] == -1.0
        rep.check(
            f"{case}: DBR/DBRC coverage, finite coefficients, DTAU>=0",
            ok and interior,
            f"rows={len(dbr)}/{len(dbrc)}",
        )

    # -- T2: A0 = exp(-DTAU) at every interval (attenuation of the
    #    upwind field). Exact in the production build.
    for case in CASES:
        worst = 0.0
        for (i, l), a in data[case][1].items():
            if i >= ND:
                continue
            worst = max(worst, abs(a["a0"] - math.exp(-a["dtau"])))
        rep.check(
            f"{case}: A0 = exp(-DTAU) at every interval",
            worst < 1.0e-10,
            f"max diff {worst:.2e}",
        )

    # -- T3: A0+A1+A2 = 1 -- a constant-source equilibrium field
    #    must march to itself. The E-series build preserves this
    #    exactly at all DTAU; the DBRC print is ES16.8, so the
    #    residual is bounded by print precision (~5e-9), not fp.
    for case in CASES:
        worst = 0.0
        for (i, l), a in data[case][1].items():
            if i >= ND:
                continue
            worst = max(worst, abs(a["a0"] + a["a1"] + a["a2"] - 1.0))
        rep.check(
            f"{case}: A0+A1+A2 = 1 (constant-source equilibrium)",
            worst < 5.0e-9,
            f"max diff {worst:.2e} (ES16.8 print floor ~5e-9)",
        )

    # -- T4: generic coupling identities, both profiles and both GAM
    #    signs. Where Q != 0 the cell is coupled: Q*CHIH*DLNM must
    #    reproduce the signed GAM and CHIH must carry |GAM|/DLNM over
    #    chi_base. Where Q = 0 the whole coupling is dropped and
    #    chi_hat reduces to chi_base. On this uniform log grid the
    #    expanding span equals DLNM, so one identity covers both.
    for case in CASES:
        dbr, _ = data[case]
        q_w, c_w, ncoupled = 0.0, 0.0, 0
        uncoupled_neg = []
        for (i, l), v in dbr.items():
            g = gam_of(case, i)
            if abs(v["q"]) > 0.0 and v["dlnm"] > 0.0:
                ncoupled += 1
                q_w = max(
                    q_w,
                    abs(v["q"] * v["chih"] * v["dlnm"] - g)
                    / max(abs(g), 1e-30),
                )
                c_w = max(
                    c_w,
                    abs(v["chih"] - (chi_base(i, l) + abs(g) / v["dlnm"]))
                    / chi_base(i, l),
                )
            else:
                if g < 0.0:
                    uncoupled_neg.append((i, l))
                c_w = max(
                    c_w,
                    abs(v["chih"] - chi_base(i, l)) / chi_base(i, l),
                )
        rep.check(
            f"{case}: coupled cells satisfy Q*CHIH*DLNM = GAM and "
            f"CHIH = chi + |GAM|/DLN",
            ncoupled > 0 and q_w < 0.05 and c_w < 0.05,
            f"ncoupled={ncoupled} q rel dev {q_w:.2e} "
            f"chih rel dev {c_w:.2e} uncoupled_neg={uncoupled_neg[:4]}",
        )

    # -- T5: DTAU interval construction. The ZERO-method build is
    #    the plain trapezoid 0.5*(CHIH_i+CHIH_{i+1})*dZ, so the
    #    dumped DTAU must reproduce it exactly from the dumped CHIH
    #    (checked cell by cell -- no discretization ambiguity, and
    #    this folds the chi_hat coupling into the opacity correctly).
    #    As a sanity cross-check the trapezoid sum must sit near the
    #    analytic int chi dz (few % at this grid resolution).
    for case in CASES:
        dbr, _ = data[case]
        d_w, s_w = 0.0, 0.0
        for (i, l), v in dbr.items():
            if i >= ND:
                continue
            expect = (
                0.5
                * (v["chih"] + dbr[(i + 1, l)]["chih"])
                * (RR[i - 1] - RR[i])
            )
            d_w = max(d_w, abs(v["dtau"] - expect) / expect)
        for l in range(1, NCFX + 1):
            numer = sum(
                v["dtau"]
                for (i, ll), v in dbr.items()
                if ll == l and i < ND
            )
            analy = (
                1.0e-3
                * (1.0 + 0.05 * l)
                * RR[0] ** 2
                * (1.0 / RR[-1] - 1.0 / RR[0])
            )
            s_w = max(s_w, abs(numer - analy) / analy)
        rep.check(
            f"{case}: DTAU = trapezoid of dumped CHIH (exact) and "
            f"consistent with analytic int chi dz",
            d_w < 1.0e-6 and s_w < 0.10,
            f"trapezoid rel dev {d_w:.2e}  analytic-sum {s_w:.2e}",
        )

    # -- T5b: SR*CHIH = ETA against the analytic profile (both cases).
    for case in CASES:
        worst = 0.0
        for (i, l), v in data[case][0].items():
            e = eta_base(i, l)
            worst = max(worst, abs(v["sr"] * v["chih"] - e) / e)
        rep.check(
            f"{case}: SR*CHIH = analytic ETA",
            worst < 0.05,
            f"max rel dev {worst:.3e}",
        )

    # -- T6: solved field is finite, non-negative, vacuum outer BC,
    #    and thermalises onto the modified source at optically thick
    #    interior cells (DTAU per cell ~ 20-600 on this model, so the
    #    field must sit close to the local SP).
    for case in CASES:
        dbr, _ = data[case]
        finite = all(
            math.isfinite(v["fldm"]) and math.isfinite(v["fldp"])
            for v in dbr.values()
        )
        nonneg = all(
            v["fldm"] >= 0.0 and v["fldp"] >= 0.0 for v in dbr.values()
        )
        obc = all(v["fldm"] == 0.0 for (i, l), v in dbr.items() if i == 1)
        th_w = 0.0
        nthick = 0
        for (i, l), v in dbr.items():
            if 2 <= i <= ND - 1 and v["spm"] > 0.0 and v["dtau"] > 3.0:
                nthick += 1
                th_w = max(
                    th_w, abs(v["fldm"] - v["spm"]) / v["spm"]
                )
        rep.check(
            f"{case}: field finite, >=0, vacuum BC, interior ~ SP",
            finite and nonneg and obc and nthick > 0 and th_w < 0.5,
            f"nthick={nthick} max |FLDM-SPM|/SPM {th_w:.2e}",
        )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
