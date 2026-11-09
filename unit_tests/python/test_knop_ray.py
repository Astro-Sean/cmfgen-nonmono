#!/usr/bin/env python3
"""Runtime test for KNOP_RAY_SOLVE_V1 (plane/fg_j_cmf_nonmono_v1.f).

Drives the production per-ray refresh solver directly through the
FG_J_CMF_NONMONO_MOD_V1 module state on a manufactured core-ray
problem (NI=8, NCF=6, one compressing band at depths 3-5).  This is
the first runtime coverage of the refresh fixed-point solver - the
machinery previously had static checks only.

Compiles the exact production file (HELPERS recipe identical to
test_mom_sweep_table.py) plus f90/drv_knop_ray.f90 and asserts on
the printed census fields: exit class (NXIT), residual, store
writes, matching census (NMATCH/NMAX), and store finiteness.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, TREE, WORKDIR, Report, compile_objs, link, run

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
    objs = compile_objs(
        "knray", [TREE / h for h in HELPERS] + [F90DIR / "drv_knop_ray.f90"]
    )
    return link("knray", objs)


def parse(out):
    d = {}
    store = {}
    import re

    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "STORE":
            store[(int(t[1]), int(t[2]))] = (float(t[3]), float(t[4]))
            continue
        m = re.match(r"^(\w+)=\s*(\S+)", ln)
        if m:
            k, v = m.group(1), m.group(2)
            try:
                d[k] = float(v)
            except ValueError:
                if v in ("T", "F"):
                    d[k] = v == "T"
                else:
                    d[k] = v
        elif t[0] == "RESULT_OK":
            d["RESULT_OK"] = True
    return d, store


def main():
    rep = Report("KNOP_RAY_SOLVE_V1")
    bin_ = build()

    # -- T1: alternating sweep on the manufactured compressing band.
    r = run(bin_, "conv", cwd=WORKDIR, check=False)
    d, store = parse(r.stdout)
    ok = (
        "RESULT_OK" in d
        and d.get("NXIT", -1) in (0, 1, 2)
        and d.get("RESID", 1e9) < 1.0e-1
        and d.get("NCOMP", 0) > 0
        and all(abs(v) < 1e30 for v in store.values() for v in v)
    )
    rep.check(
        "T1 alternating sweep converges on manufactured band",
        ok,
        f"NXIT={d.get('NXIT')} RESID={d.get('RESID'):.3e} "
        f"NPASS={d.get('NPASS')} NCOMP={d.get('NCOMP')} rc={r.returncode}",
    )

    # -- T2: ordered sweep (KN_RAY_ORDERED) solves the same field.
    r = run(bin_, "ord", cwd=WORKDIR, check=False)
    d2, store2 = parse(r.stdout)
    ok = "RESULT_OK" in d2 and d2.get("RESID", 1e9) < 1.0e-1
    rep.check(
        "T2 ordered sweep converges",
        ok,
        f"NXIT={d2.get('NXIT')} RESID={d2.get('RESID'):.3e} "
        f"NPASS={d2.get('NPASS')}",
    )

    # -- T3: ordered vs alternating committed stores agree to the
    #    solve tolerance (both march the same discrete operator).
    common = set(store) & set(store2)
    if common:
        dd = max(
            abs(store[k][j] - store2[k][j]) / max(store2[k][j], 1e-30)
            for k in common
            for j in (0, 1)
        )
    else:
        dd = 9.9e9
    rep.check(
        "T3 ord/alt committed stores agree",
        "RESULT_OK" in d2 and dd < 5.0e-2,
        f"max rel store diff {dd:.2e} over {len(common)} slots",
    )

    # -- T4: DO_NONMONO_FREQ_MATCH with drift > 1 bin must activate
    #    multi-bin matched slots (NMATCH>0 and NMAX>1).
    r = run(bin_, "match", cwd=WORKDIR, check=False)
    d3, _ = parse(r.stdout)
    ok = "RESULT_OK" in d3 and d3.get("NMATCH", 0) > 0 and d3.get(
        "NMAX", 1
    ) > 1
    rep.check(
        "T4 FREQ_MATCH activates multi-bin slots under strong drift",
        ok,
        f"NMATCH={d3.get('NMATCH')} NMAX={d3.get('NMAX')} "
        f"RESID={d3.get('RESID'):.3e}",
    )

    # -- T5: a non-converged solve must exit visibly, never as
    #    'converged': forced stall/cap via KN_RAY_TOL=1e-14, MAXPASS=4
    #    must produce NXIT=4 (cap) - or another failure class - and
    #    a residual that stays above the production tolerance.
    r = run(bin_, "cap", cwd=WORKDIR, check=False)
    d4, _ = parse(r.stdout)
    xit = d4.get("NXIT", -1)
    ok = "RESULT_OK" in d4 and xit in (3, 4, 5, 6) and d4.get(
        "RESID", 0.0
    ) > 1.0e-14
    rep.check(
        "T5 forced non-convergence exits visibly (not as converged)",
        ok,
        f"NXIT={xit} RESID={d4.get('RESID'):.3e} NPASS={d4.get('NPASS')}",
    )

    # -- T6: warm start near the fixed point must show a small first-
    #    pass residual and converge in very few passes.
    r = run(bin_, "fixed", cwd=WORKDIR, check=False)
    d5, _ = parse(r.stdout)
    ok = "RESULT_OK" in d5 and d5.get("NXIT", -1) == 1 and d5.get(
        "NPASS", 99
    ) <= 4
    rep.check(
        "T6 near-fixed-point seed converges immediately",
        ok,
        f"NXIT={d5.get('NXIT')} RESID={d5.get('RESID'):.3e} "
        f"NPASS={d5.get('NPASS')}",
    )

    # -- T7: checkpoint round-trip: dump/restore identical, and the
    #    failure classes are all visible and distinguishable.
    r = run(bin_, "ckpt", cwd=WORKDIR, check=False)
    d6, _ = parse(r.stdout)
    ok = (
        "RESULT_OK" in d6
        and d6.get("CKPT_OK") is not False
        and d6.get("CKPT_MAXDIFF", 1.0) == 0.0
        and d6.get("CKPT_FP_OK") is False
        and d6.get("CKPT_FP_STORE", 1.0) == 0.0
        and d6.get("CKPT_TAG_OK") is False
        and d6.get("CKPT_FOR_OK") is False
        and d6.get("CKPT_FOR_STORE", 0.0) > 0.0
    )
    rep.check(
        "T7 FGINEXT_GS round-trip + reject classes (fp, tag, foreign)",
        ok,
        f"diff={d6.get('CKPT_MAXDIFF')} fp_ok={d6.get('CKPT_FP_OK')} "
        f"tag_ok={d6.get('CKPT_TAG_OK')} "
        f"for_ok={d6.get('CKPT_FOR_OK')} "
        f"for_store={d6.get('CKPT_FOR_STORE')}",
    )

    # -- T8: end-to-end FG_KNOP_REFRESH_V1 over the OMP ray loop:
    #    residual lands in KN_RF_RESID, the exit census is populated,
    #    and multi-thread output is identical to one thread.
    d7, s1 = parse(run(bin_, "rfsh", cwd=WORKDIR, check=False).stdout)
    import os
    import subprocess

    env = dict(os.environ, OMP_NUM_THREADS="4")
    r8 = subprocess.run(
        [str(bin_), "rfsh"], cwd=WORKDIR, env=env,
        capture_output=True, text=True, timeout=300,
    )
    d8, s2 = parse(r8.stdout)
    omp_diff = max(
        abs(s1[k][j] - s2[k][j]) / max(s1[k][j], 1e-30)
        for k in set(s1) & set(s2)
        for j in (0, 1)
    ) if (s1 and s2) else 9.9e9
    ok = (
        "RESULT_OK" in d7
        and "RESULT_OK" in d8
        and d7.get("RF_RESID", 1.0) < 1.0e-2
        and d7.get("RF_NDIV", 9) == 0
        and omp_diff < 1.0e-12
    )
    rep.check(
        "T8 FG_KNOP_REFRESH_V1 end-to-end + 1 vs 4 threads identical",
        ok,
        f"RF_RESID={d7.get('RF_RESID'):.3e} NDIV={d7.get('RF_NDIV')} "
        f"omp_maxdiff={omp_diff:.2e}",
    )

    # -- T9: two compressing bands separated by an expanding cell -
    #    the gradient sign reverses twice along the ray, so the
    #    sign-dependent upwind pick runs in both directions on one
    #    profile. Solver must still converge with finite stores.
    r = run(bin_, "rev", cwd=WORKDIR, check=False)
    d9, s9 = parse(r.stdout)
    ok = (
        "RESULT_OK" in d9
        and d9.get("NXIT", -1) in (0, 1, 2)
        and d9.get("RESID", 1e9) < 1.0e-1
        and all(abs(v) < 1e30 for v in s9.values() for v in v)
    )
    rep.check(
        "T9 multi-reversal profile converges (two GAM bands)",
        ok,
        f"NXIT={d9.get('NXIT')} RESID={d9.get('RESID'):.3e} "
        f"NPASS={d9.get('NPASS')} NCOMP={d9.get('NCOMP')}",
    )

    # -- T10: strong compression - same band geometry, GAM ~150x
    #    deeper negative (sigma << -1 regime). Converged exit,
    #    finite stores, no divergence census.
    r = run(bin_, "strg", cwd=WORKDIR, check=False)
    d10, s10 = parse(r.stdout)
    ok = (
        "RESULT_OK" in d10
        and d10.get("NXIT", -1) in (0, 1, 2)
        and d10.get("RESID", 1e9) < 1.0e-1
        and d10.get("NDIV", 9) == 0
        and all(abs(v) < 1e30 for v in s10.values() for v in v)
    )
    rep.check(
        "T10 strong-compression band converges (GAM ~ 150x)",
        ok,
        f"NXIT={d10.get('NXIT')} RESID={d10.get('RESID'):.3e} "
        f"NPASS={d10.get('NPASS')} NDIV={d10.get('NDIV')}",
    )

    # -- T11: donor-mask drop equivalence (audit P2). 'gap' puts a
    #    compressing GAM cell at ray point 5 with no sparse donor slot
    #    in its bracket, so the whole frequency coupling is dropped
    #    there. 'gapr' gives the same cell GAM=0, which produces the
    #    identical discrete operator: plain chi, zero Q, no upwind
    #    term. Both must converge and commit the same store.
    r = run(bin_, "gap", cwd=WORKDIR, check=False)
    dg, sg = parse(r.stdout)
    r = run(bin_, "gapr", cwd=WORKDIR, check=False)
    dr, sr = parse(r.stdout)
    common = set(sg) & set(sr)
    if common:
        dd = max(
            abs(sg[k][j] - sr[k][j]) / max(sr[k][j], 1e-30)
            for k in common
            for j in (0, 1)
        )
    else:
        dd = 9.9e9
    ok = (
        "RESULT_OK" in dg
        and "RESULT_OK" in dr
        and dg.get("NXIT", -1) in (0, 1, 2)
        and dr.get("NXIT", -1) in (0, 1, 2)
        and dd < 5.0e-2
    )
    rep.check(
        "T11 dropped coupling == never-coupled operator (gap vs gapr)",
        ok,
        f"NXIT={dg.get('NXIT')}/{dr.get('NXIT')} "
        f"RESID={dg.get('RESID'):.3e}/{dr.get('RESID'):.3e} "
        f"max rel store diff {dd:.2e}",
    )

    # -- T12: NaN opacity in the captured store must fail the ray
    #    (diverged exit, residual pinned at 1), never commit a
    #    repaired field.
    r = run(bin_, "nanchi", cwd=WORKDIR, check=False)
    dn, sn = parse(r.stdout)
    ok = (
        "RESULT_OK" in dn
        and dn.get("NXIT", -1) == 5
        and dn.get("NDIV", 0) == 1
        and dn.get("RESID", 0.0) == 1.0
    )
    rep.check(
        "T12 NaN CHI fails the ray closed (no commit)",
        ok,
        f"NXIT={dn.get('NXIT')} RESID={dn.get('RESID'):.3e} "
        f"NDIV={dn.get('NDIV')}",
    )

    # -- T13: infinite GAM likewise: Inf slipped past the old
    #    NaN-only guards into DTAU and the march coefficients.
    r = run(bin_, "infgam", cwd=WORKDIR, check=False)
    di, si = parse(r.stdout)
    ok = (
        "RESULT_OK" in di
        and di.get("NXIT", -1) == 5
        and di.get("NDIV", 0) == 1
        and di.get("RESID", 0.0) == 1.0
    )
    rep.check(
        "T13 Inf GAM fails the ray closed (no commit)",
        ok,
        f"NXIT={di.get('NXIT')} RESID={di.get('RESID'):.3e} "
        f"NDIV={di.get('NDIV')}",
    )

    # -- T14: formal-solver/refresh operator A/B (audit phase 2, task
    #    5). 'ab' drives FG_J_CMF_NONMONO_V1 itself for six outer
    #    iterations on the manufactured problem, dumps the discrete
    #    operator per (cell, frequency), then runs the refresh over
    #    the captured store. 'abx' is the same comparison with the
    #    ordered refresh path disabled (alternating SWEEP_FLD).
    #    Field-independent coefficients must agree to the span
    #    evaluation (~1e-6); index maps and interpolation weights
    #    must be identical; the red-edge cell (compressing, last
    #    frequency) must show the same pure-outflow closure.
    import re as _re

    _F = _re.compile(r"[+-]?\d+\.\d+(?:E[+-]\d+|[+-]\d{3,})")

    def _f(x):
        if "E" not in x:
            x = x[:-4] + "E" + x[-4:]
        return float(x)

    def _db(text):
        d = {}
        for ln in text.splitlines():
            t = ln.split()
            if len(t) < 4 or t[0] not in ("DBS", "DBR", "DBSA", "DBRA",
                                        "DBSI", "DBRI"):
                continue
            if t[0] in ("DBSI", "DBRI"):
                vals = [float(x) for x in t[3:]]
            else:
                vals = [_f(x) for x in _F.findall(ln)]
            d.setdefault(t[0], {})[(int(t[1]), int(t[2]))] = vals
        return d

    def _abdiffs(text):
        d = _db(text)
        coef = [0.0] * 4
        nmap = 0
        for i in range(1, 8):
            for l in range(1, 7):
                s = d["DBS"].get((i, l))
                r = d["DBR"].get((i, l))
                if not s or not r:
                    continue
                for j in range(4):
                    coef[j] = max(
                        coef[j],
                        abs(s[j] - r[j]) / max(abs(r[j]), 1e-30),
                    )
                a = d["DBSI"].get((i, l))
                b = d["DBRI"].get((i, l))
                if a and b:
                    if int(a[0]) != int(b[0]):
                        nmap += 1
                    # sweep KC are sparse-store slots, refresh IKC the
                    # availability flags for the same bracket ends --
                    # the raw indices differ by convention (JPR is the
                    # ray-row of the bracketing depth via J_PNT, which
                    # is not identity on the production grid), so the
                    # meaningful equality is donor availability.
                    if (a[1] >= 1) != (int(b[1]) == 1):
                        nmap += 1
                    if (a[2] >= 1) != (int(b[2]) == 1):
                        nmap += 1
        amax = 0.0
        for i in range(1, 8):
            for l in range(1, 7):
                s = d["DBSA"].get((i, l))
                r = d["DBRA"].get((i, l))
                if not s or not r:
                    continue
                for j in range(5):
                    amax = max(
                        amax,
                        abs(s[j] - r[j]) / max(abs(r[j]), 1e-30),
                    )
        red = None
        s = d["DBS"].get((5, 6))
        r = d["DBR"].get((5, 6))
        a = d["DBSI"].get((5, 6))
        b = d["DBRI"].get((5, 6))
        if s and r and a and b:
            red = (
                max(
                    abs(s[j] - r[j]) / max(abs(r[j]), 1e-30)
                    for j in range(4)
                ),
                int(a[0]),
                int(b[0]),
            )
        return coef, amax, nmap, red

    for tag in ("ab", "abx"):
        rr = run(bin_, tag, cwd=WORKDIR, check=False)
        coef, amax, nmap, red = _abdiffs(rr.stdout)
        ok = (
            "RESULT_OK" in rr.stdout
            and max(coef) < 1.0e-6
            and amax < 1.0e-5
            and nmap == 0
            and red is not None
            and red[0] < 1.0e-6
            and red[1] == 1
            and red[2] == 1
        )
        rep.check(
            f"T14 {tag}: refresh == in-sweep operator "
            f"({'ordered' if tag == 'ab' else 'alternating'})",
            ok,
            f"chih/q/dchr/dtau max={max(coef):.2e} A0-A4 max={amax:.2e} "
            f"index mismatches={nmap} red-edge={red}",
        )

    # -- T15: GAM zero-crossing continuity (audit phase 2, task 2).
    #    gpp/gpn0 put GAM(4) at +/-1e-9*GAMP: the discrete operator is
    #    continuous in |GAM| with the sign flipping only the upwind
    #    pick, so the committed stores must be bit-identical. gzz is
    #    the exact-zero neighbour; its field must sit within the
    #    fixed-point tolerance of gpp (no jump). All six profiles
    #    must converge with no divergence, floor, or cap census.
    store_map = {}
    dgam = {}
    for tag in ("gpos", "gneg", "gpp", "gpn0", "gzz", "gcr"):
        rr = run(bin_, tag, cwd=WORKDIR, check=False)
        dgam[tag], store_map[tag] = parse(rr.stdout)

    def _sdiff(a, b):
        common = set(a) & set(b)
        if not common:
            return 9.9e9
        return max(
            abs(a[k][j] - b[k][j]) / max(b[k][j], 1e-30)
            for k in common
            for j in (0, 1)
        )

    clean = all(
        "RESULT_OK" in dgam[t]
        and dgam[t].get("NDIV", 9) == 0
        and dgam[t].get("NNEG", 9) == 0
        and dgam[t].get("NCAP", 9) == 0
        and dgam[t].get("NXIT", -1) in (0, 1, 2)
        for t in dgam
    )
    d_micro = _sdiff(store_map["gpp"], store_map["gpn0"])
    d_zero = _sdiff(store_map["gpp"], store_map["gzz"])
    ok = clean and d_micro == 0.0 and d_zero < 2.0e-2
    rep.check(
        "T15 GAM sign flip is continuous (microstep identical, "
        "exact-zero within tol)",
        ok,
        f"clean={clean} microstep={d_micro:.2e} exact0={d_zero:.2e}",
    )

    rep.summary()


if __name__ == "__main__":
    main()
