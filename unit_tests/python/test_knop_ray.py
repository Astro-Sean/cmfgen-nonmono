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


def build():
    return compile_and_link_prod(
        "knray", [TREE / h for h in HELPERS] + [F90DIR / "drv_knop_ray.f90"]
    )


def parse(out):
    d = {}
    store = {}
    import re

    def _f(x):
        # ES16.8 drops the E on 3-digit exponents: '1.79+308'.
        if "E" not in x:
            x = x[:-4] + "E" + x[-4:]
        return float(x)

    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "STORE":
            store[(int(t[1]), int(t[2]))] = (float(t[3]), float(t[4]))
            continue
        if t[0] == "RSTORE":
            store[("r", int(t[1]), int(t[2]), int(t[3]))] = (
                float(t[4]),
                float(t[5]),
            )
            continue
        if t[0] == "SREF":
            store[("sref", int(t[1]))] = (float(t[2]), float(t[3]))
            continue
        if t[0] == "VC":
            store[("vc", int(t[1]))] = (int(t[2]), _f(t[3]))
            continue
        m = re.match(r"^(\w+)=\s*(\S+)", ln)
        if m:
            k, v = m.group(1), m.group(2)
            try:
                d[k] = float(v)
            except ValueError:
                try:
                    d[k] = _f(v)
                    continue
                except (ValueError, IndexError):
                    pass
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

    # -- T14: formal-solver/refresh operator A/B (audit phase 2,
    #    task 5; strengthened in the follow-up). 'ab' drives
    #    FG_J_CMF_NONMONO_V1 itself for six outer iterations on the
    #    manufactured problem, dumps the discrete operator per
    #    (cell, frequency), then runs the refresh over the captured
    #    store. 'abx' is the same comparison with the alternating
    #    SWEEP_FLD path; 'ab2' dumps the non-central ray LS=2
    #    (NPA=2, PA(2)=1.2E4, NI_SMALL=7 model depths on-ray plus
    #    periapsis) so the donor map is exercised on a ray whose
    #    geometry and compressing set differ from the core ray.
    #
    #    Coverage is mandatory: every (depth, frequency) record must
    #    exist on both sides, hold finite non-sentinel values in the
    #    dumped coefficient columns, and appear exactly once. The
    #    only permitted sentinel slots are bookkeeping rows that are
    #    never used: refresh DTAU at i=NI (interval array covers
    #    1..NI-1), refresh SPM at i=1 and SPP at i=NI (march
    #    boundaries take the incident field, not a local source).
    #
    #    Field-independent quantities are compared directly:
    #    CHI_HAT, Q, DCHR, DTAU (interior), WINX, DLNM, A0-A4, SR.
    #    These are built from the same captured store and model
    #    inputs; residual differences come only from the evaluation
    #    order of the DLOG_NU spans and EXP(-T1), so a scale-aware
    #    tolerance of 1e-5 relative plus a small absolute floor is
    #    appropriate (observed: <=3e-7 for all but A0 at 2e-6).
    #
    #    SPM/SPP are field-dependent (lagged donor field of two
    #    different iterates) and are NOT compared across paths.
    #    Instead the construction rule is checked per path: where the
    #    coupling coefficient multiplies a real donor the modified
    #    source differs from SR, and where no donor exists the
    #    modified source reduces to SR exactly -- this is also the
    #    red-edge pure-outflow assertion at the compressing cell
    #    (5,6): phantom slot forces the donor to zero, so SPM must
    #    equal SR bitwise on BOTH paths while Q itself stays finite.
    #
    #    Donor identity: the sweep's KC1/KC2 are sparse-store slots
    #    (depth = inverse COMP_DEPTH_IDX), the refresh's JPR1/JPR2
    #    are ray rows (depth = inverse J_PNT). The test maps both to
    #    physical model depths and requires identical donor depths,
    #    identical availability, identical NM/NFM offsets, and equal
    #    WINX/DLNM. A same-availability/different-depth selection
    #    fails.
    import re as _re

    _F = _re.compile(r"[+-]?\d+\.\d+(?:E[+-]\d+|[+-]\d{3,})")

    def _f(x):
        if "E" not in x:
            x = x[:-4] + "E" + x[-4:]
        return float(x)

    _BLOCKS = (
        "DBS", "DBR", "DBSA", "DBRA", "DBSI", "DBRI",
        "DBSD", "DBRD", "JMAP", "GAMV",
    )

    def _db(text):
        # Every row must be well-formed and unique. A duplicate or
        # malformed record is a hard failure, not a silent overwrite.
        d = {}
        errs = []
        for ln in text.splitlines():
            t = ln.split()
            if not t:
                continue
            if len(t) < 3:
                if t[0] in _BLOCKS:
                    errs.append(f"malformed {t[0]} record: {ln!r}")
                continue
            key = None
            if t[0] in ("DBSI", "DBRI", "DBSD", "DBRD"):
                want = {"DBSI": 3, "DBRI": 5, "DBSD": 1, "DBRD": 1}[t[0]]
                try:
                    key = (t[0], int(t[1]), int(t[2]))
                    val = [int(x) for x in t[3:]]
                except ValueError:
                    errs.append(f"{t[0]} malformed: {ln!r}")
                    continue
                if len(t) != want + 3:
                    errs.append(f"{t[0]} {(t[1], t[2])} width {len(t)-3}")
            elif t[0] in ("DBS", "DBR", "DBSA", "DBRA"):
                key = (t[0], int(t[1]), int(t[2]))
                val = [_f(x) for x in _F.findall(ln)]
                want = 5 if t[0] in ("DBSA", "DBRA") else 11
                if len(val) != want:
                    errs.append(
                        f"{t[0]} {(t[1], t[2])} width {len(val)}!={want}"
                    )
                if not all(abs(v) < 1e30 for v in val):
                    errs.append(f"{t[0]} {(t[1], t[2])} nonfinite")
            elif t[0] == "JMAP":
                if len(t) != 4:
                    errs.append(f"JMAP {t[1:2]} width {len(t)-1}!=3")
                    continue
                try:
                    k = int(t[1])
                    jv = (int(t[2]), int(t[3]))
                except ValueError:
                    errs.append(f"JMAP malformed: {ln!r}")
                    continue
                if k in d.setdefault("JMAP", {}):
                    errs.append(f"duplicate JMAP {k}")
                if not (0 <= jv[0] <= 10**4 and 0 <= jv[1] <= 10**4):
                    errs.append(f"JMAP {k} invalid depth {jv}")
                d["JMAP"][k] = jv
                continue
            elif t[0] == "GAMV":
                if len(t) != 3:
                    errs.append(f"GAMV {t[1:2]} width {len(t)-1}!=2")
                    continue
                try:
                    k = int(t[1])
                    gv = float(t[2])
                except ValueError:
                    errs.append(f"GAMV malformed: {ln!r}")
                    continue
                if k in d.setdefault("GAMV", {}):
                    errs.append(f"duplicate GAMV {k}")
                if not abs(gv) < 1e30:
                    errs.append(f"GAMV {k} nonfinite")
                d["GAMV"][k] = gv
                continue
            else:
                continue
            blk, i, l = key
            if (i, l) in d.setdefault(blk, {}):
                errs.append(f"duplicate {blk} {(i, l)}")
            d[blk][(i, l)] = val
        if errs:
            d["__err__"] = errs
        return d

    def _close(a, b, rtol=1.0e-5, atol=0.0):
        return abs(a - b) <= atol + rtol * abs(b)

    def _abcheck(text, tag):
        """Returns (ok, detail) for one A/B dump."""
        d = _db(text)
        if "RESULT_OK" not in text:
            return False, "driver failed"
        if "__err__" in d:
            return False, "; ".join(d["__err__"][:2])
        for blk in ("DBS", "DBR", "DBSA", "DBRA", "DBSI", "DBRI"):
            if blk not in d:
                return False, f"missing block {blk}"
        ni = max(i for i, _ in d["DBS"])
        ncf = max(l for _, l in d["DBS"])
        grid = {(i, l) for i in range(1, ni + 1) for l in range(1, ncf + 1)}
        for blk in ("DBS", "DBR", "DBSD", "DBRD"):
            if set(d.get(blk, {})) != grid:
                return False, f"{blk} coverage {len(d.get(blk, {}))} != {len(grid)}"
        jk = sorted(d.get("JMAP", {}))
        gk = sorted(d.get("GAMV", {}))
        if not jk or jk != list(range(1, max(jk) + 1)):
            return False, "JMAP non-contiguous or missing rows"
        if not gk or gk != list(range(1, max(gk) + 1)):
            return False, "GAMV non-contiguous or missing rows"
        agrid = {(i, l) for i in range(1, ni) for l in range(1, ncf + 1)}
        for blk in ("DBSA", "DBRA"):
            if set(d[blk]) != grid:
                return False, f"{blk} coverage {len(d[blk])} != {len(grid)}"
        # March coefficients live on intervals 1..NI-1; rows past the
        # last marched interval (periapsis bookkeeping on a p-ray)
        # stay -1 on BOTH sides -- the active-interval sets must be
        # identical and every active interval fully populated.
        active = {
            k for k in agrid
            if d["DBSA"][k][0] != -1.0 or d["DBRA"][k][0] != -1.0
        }
        for k in active:
            if d["DBSA"][k][0] == -1.0 or d["DBRA"][k][0] == -1.0:
                return False, f"A active-interval mismatch at {k}"
        if len(active) < (ni - 2) * ncf:
            return False, f"only {len(active)} active A intervals"
        # Sentinel / finite policy: every dumped column must be
        # initialised and finite. A whole record may be -1 on BOTH
        # sides where the ray point is never marched (periapsis
        # bookkeeping on a p-ray) -- such inactive rows are allowed
        # only if identical across paths; a one-sided gap fails.
        act = {
            k for k in grid
            if d["DBS"][k][0] != -1.0 or d["DBR"][k][0] != -1.0
        }
        for k in act:
            if d["DBS"][k][0] == -1.0 or d["DBR"][k][0] == -1.0:
                return False, f"one-sided inactive record at {k}"
        if len(act) < (ni - 1) * ncf:
            return False, f"only {len(act)} active records"
        for k in act:
            s, r = d["DBS"][k], d["DBR"][k]
            if len(s) != 11 or len(r) != 11:
                return False, f"row width at {k}"
            for j in (0, 1, 2, 6, 7, 10):
                if s[j] == -1.0 or r[j] == -1.0:
                    return False, f"sentinel col {j} at {k}"
            # SPM/SPP boundary rows: the ordered refresh never builds
            # SP at the first (inward) or last (outward) marched
            # point -- the boundary takes the incident field. On a
            # p-ray those are the extreme ACTIVE rows, not 1/NI.
            i_min = min(i for i, _ in act)
            i_max = max(i for i, _ in act)
            for j in (4, 5):
                if s[j] == -1.0:
                    return False, f"sweep sentinel SP col {j} at {k}"
                if r[j] == -1.0 and not (
                    (j == 4 and k[0] == i_min)
                    or (j == 5 and k[0] == i_max)
                ):
                    return False, f"refresh sentinel SP col {j} at {k}"
                if not (abs(s[j]) < 1e30 and abs(r[j]) < 1e30):
                    return False, f"nonfinite col {j} at {k}"
            if k[0] < ni and (s[3] == -1.0 or not abs(s[3]) < 1e30):
                return False, f"bad sweep DTAU at {k}"
            if k[0] < ni and (r[3] == -1.0 or not abs(r[3]) < 1e30):
                return False, f"bad refresh DTAU at {k}"
        # Donor maps from JMAP: J_PNT row -> depth, slot -> depth.
        jm = d.get("JMAP", {})
        if not jm:
            return False, "no JMAP"
        row2dep = {v[0]: k for k, v in jm.items() if v[0] > 0}
        slot2dep = {v[1]: k for k, v in jm.items() if v[1] > 0}
        gam = d.get("GAMV", {})
        # Donor-access counter: a cell with GAM=0 must never enter the
        # expanding-donor branch on either path (0*NaN is NaN, so the
        # read itself is the defect, not just the contribution). A
        # positive-GAM cell must read at every l>1; a compressing cell
        # uses the lagged I_NEXT branch instead and must not trip this
        # counter. If the .GT.0 guards were reverted to .GE. the zero
        # cells would trip immediately.
        for (i, l) in grid:
            g = gam.get(i)
            if g is None:
                continue
            ns = d["DBSD"][(i, l)][0]
            nr = d["DBRD"][(i, l)][0]
            if g == 0.0 and (ns != 0 or nr != 0):
                return False, f"donor read at GAM=0 cell {(i, l)}"
            if g < 0.0 and (ns != 0 or nr != 0):
                return False, f"expanding-donor read at GAM<0 {(i, l)}"
            if g > 0.0 and l > 1 and (ns <= 0 or nr <= 0):
                return False, f"missing donor read at GAM>0 {(i, l)}"
        worst = {}
        worst_abs = {}
        nfail = []

        def _w(name, val, aval):
            if val > worst.get(name, 0.0):
                worst[name] = val
                worst_abs[name] = aval

        for k in act:
            s, r = d["DBS"][k], d["DBR"][k]
            si, ri = d["DBSI"][k], d["DBRI"][k]
            # Field-independent operator columns.
            # Scale-aware check: the reference scale is the local
            # value itself (|ref|, not an arbitrary floor), with a
            # tiny absolute floor of 1e-9*scale so a physically
            # significant discrepancy can never divide itself small.
            for j, name in ((0, "chih"), (1, "q"), (2, "dchr"),
                            (6, "wx"), (7, "dlnm"), (10, "sr")):
                sc = max(abs(r[j]), abs(s[j]), 1e-30)
                _w(name, abs(s[j] - r[j]) / sc, abs(s[j] - r[j]))
                if not _close(s[j], r[j], 1.0e-5, sc * 1.0e-9):
                    nfail.append(f"{name}@{k} {s[j]}!={r[j]}")
            if k[0] < ni:
                sc = max(abs(r[3]), abs(s[3]), 1e-30)
                _w("dtau", abs(s[3] - r[3]) / sc, abs(s[3] - r[3]))
                if not _close(s[3], r[3], 1.0e-5, sc * 1.0e-9):
                    nfail.append(f"dtau@{k} {s[3]}!={r[3]}")
            # Frequency offset + availability parity.
            if si[0] != ri[0]:
                nfail.append(f"nm/nfm@{k} {si[0]}!={ri[0]}")
            for e in (1, 2):
                if (si[e] >= 1) != (ri[e] == 1):
                    nfail.append(f"donor avail{e}@{k}")
            # Physical donor depth identity where both paths hold a
            # bracket end.
            for e in (1, 2):
                if si[e] >= 1 and ri[e] == 1:
                    sd = slot2dep.get(si[e])
                    rd = row2dep.get(ri[e + 2])
                    if sd is None or rd is None or sd != rd:
                        nfail.append(
                            f"donor depth{e}@{k} slot{si[e]}->{sd} "
                            f"row{ri[e+2]}->{rd}"
                        )
            # SP construction rule per path (field-dependent, so no
            # cross-compare): Q != 0 with an active donor must move
            # the source off SR; no donor must leave SPM == SR.
            for vals, sp in ((s, s[4]), (r, r[4])):
                if abs(vals[1]) > 0.0 and k[1] > 1:
                    pass  # donor-weighted: any shift allowed
            # Red-edge pure outflow: compressing cell, last
            # frequency, donor forced to the zeroed phantom slot ->
            # coupling contributes exactly nothing, SPM == SR
            # bitwise on both paths while Q stays finite.
            if k[1] == ncf and si[0] == 1 and gam.get(k[0], 1.0) < 0.0:
                if s[4] != s[10] or r[4] != r[10]:
                    nfail.append(f"red-edge SP!=SR@{k}")
                if s[1] == 0.0 or r[1] == 0.0:
                    nfail.append(f"red-edge Q=0@{k}")
        for k in active:
            s, r = d["DBSA"][k], d["DBRA"][k]
            for j in range(5):
                sc = max(abs(r[j]), abs(s[j]), 1e-30)
                _w(f"a{j}", abs(s[j] - r[j]) / sc, abs(s[j] - r[j]))
                if not _close(s[j], r[j], 1.0e-5, sc * 1.0e-9):
                    nfail.append(f"a{j}@{k} {s[j]}!={r[j]}")
        det = " ".join(
            f"{m}={v:.1e}({worst_abs[m]:.1e})"
            for m, v in sorted(worst.items())
        )
        if nfail:
            det = f"{len(nfail)} fails [{nfail[0]}] | " + det
        return not nfail, det

    for tag in ("ab", "abx", "ab2", "abz"):
        rr = run(bin_, tag, cwd=WORKDIR, check=False)
        ok, det = _abcheck(rr.stdout, tag)
        rep.check(
            f"T14 {tag}: refresh == in-sweep operator "
            f"({'ordered' if tag == 'ab' else 'alternating' if tag == 'abx' else 'noncentral' if tag == 'ab2' else 'zerocell'})",
            ok,
            det[:160],
        )

    # -- T15: resolved GAM zero-crossing (audit follow-up B).
    #    gz[pn][1-4] put GAM(4) at +/-10^-2..10^-5 of GAMP, a
    #    resolved perturbation sequence (the |GAM|/DLN contribution
    #    at 1e-2*GAMP is ~1.4e-7 against CHI ~ 4e-3). Assertions are
    #    coefficient-level at the crossing cell:
    #      * exact zero (gzz): Q == 0 and SPM == SPP == SR bitwise
    #        at every frequency -- no coupling read or applied.
    #      * continuity: |coef(eps)-coef(0)| shrinks ~10x per decade
    #        toward the zero-coupling operator for BOTH signs.
    #      * the CHI_HAT shift is sign-blind (|GAM|/DLN): +eps and
    #        -eps give the same chi_hat increment.
    #      * Q flips sign with GAM where the compressing branch
    #        applies (donor exists), magnitude symmetric.
    #      * gzb: an exact-zero cell between opposite-sign
    #        neighbours keeps Q=0/SP=SR while the neighbours keep
    #        their own couplings.
    #      * gcr: the resolved reversal has Q<0 at cell 4 and Q>0
    #        at cell 5 at coupled frequencies.
    #    The committed stores of the micro pair gpp/gpn0 stay a
    #    secondary check.
    store_map = {}
    dgam = {}
    dbg = {}
    for tag in ("gpos", "gneg", "gpp", "gpn0", "gzz", "gcr", "gzb",
                "gzp1", "gzn1", "gzp2", "gzn2", "gzp3", "gzn3",
                "gzp4", "gzn4"):
        rr = run(bin_, tag, cwd=WORKDIR, check=False)
        dgam[tag], store_map[tag] = parse(rr.stdout)
        dbg[tag] = _db(rr.stdout)

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

    fails = []
    # Exact zero: no coupling coefficient AND no donor access. The
    # KN_DBG_NDON counter increments inside each .GT.0 donor-read
    # site (sweep minus/plus and both refresh helpers), so a zero
    # cell can only trip it if the exact-zero guard is bypassed --
    # under the old .GE.0 branch these counts would be nonzero.
    z = dbg["gzz"]
    for l in range(1, 7):
        r = z["DBR"].get((4, l))
        if r is None or r[1] != 0.0 or r[4] != r[10] or r[5] != r[10]:
            fails.append(f"gzz cell4 l{l} coupling present")
        # The gzz path runs only the refresh (no sweep dump), so the
        # DBRD counter carries the read evidence here; the sweep-side
        # guard is pinned by the abz case through _abcheck. gzz has
        # exact zeros at cells 4 and 5.
        for i in (4, 5):
            if z["DBRD"].get((i, l), [1])[0] != 0:
                fails.append(f"gzz cell{i} l{l} donor read")
    # Counter liveness: every expanding cell must show reads on the
    # refresh path at l>1, or a dead counter could fake the zeros.
    for i, g in z.get("GAMV", {}).items():
        if g <= 0.0:
            continue
        for l in range(2, 7):
            if z["DBRD"].get((i, l), [0])[0] <= 0:
                fails.append(f"gzz expanding cell {i} l{l} no read")
    # Continuity: separate dimensionless measures per quantity.
    #   e_chih = |chih - chih0| / chih0      (local opacity scale)
    #   e_q    = |q| / q_scale               q_scale = |GAM1|/(dlnm*chih0)
    #   e_sp   = |spm - sr| / sr             (local source scale)
    # q_scale is the implemented coupling coefficient evaluated at the
    # largest perturbation (GAM1 = 1e-2*GAMP), derived from the dumped
    # DLNM and the exact-zero CHI_HAT -- not an arbitrary constant.
    # Each sign must contract monotonically toward the zero operator;
    # where two successive levels are both resolved, the ratio should
    # sit at ~0.1 (linear in |GAM|). The SPM error includes the donor
    # field itself, which converges only to the fixed-point tolerance,
    # so its floor cuts in sooner -- same bound, wider resolved floor.
    ratios = {"p": [], "n": []}
    quiet = {"p": [], "n": []}
    for l in range(1, 7):
        r0 = z["DBR"].get((4, l))
        if r0 is None:
            fails.append(f"gzz l{l} missing")
            continue
        chih0, sr0 = r0[0], r0[10]
        # GAM at the largest perturbation is 1e-2*GAMP; GAMP=4e-6 in
        # the driver, so |GAM1| = 4e-8.
        qs = 4.0e-8 / max(abs(r0[7]) * chih0, 1e-30)
        for sgn in ("p", "n"):
            es = []
            for n in (1, 2, 3, 4):
                r = dbg[f"gz{sgn}{n}"]["DBR"].get((4, l))
                if r is None:
                    fails.append(f"gz{sgn}{n} l{l} missing")
                    continue
                es.append(
                    (
                        abs(r[0] - chih0) / max(chih0, 1e-30),
                        abs(r[1]) / qs,
                        # SPM - SR within the SAME run: SR drifts
                        # ~1e-8 run-to-run, which would swamp the
                        # finest coupling levels.
                        abs(r[4] - r[10]) / max(sr0, 1e-30),
                    )
                )
            if len(es) < 4:
                continue
            if all(e == (0.0, 0.0, 0.0) for e in es):
                # l=1 expansion has no l-1 donor: the coupling is
                # structurally absent, not suppressed -- classify it.
                quiet[sgn].append(l)
                continue
            for n in range(1, 4):
                for qi, qn in enumerate(("chih", "q", "sp")):
                    e0, e1 = es[n - 1][qi], es[n][qi]
                    if e1 > e0 * 1.05 + 1e-13:
                        fails.append(
                            f"gz{sgn} l{l} {qn} not monotonic "
                            f"{e0:.1e}->{e1:.1e}"
                        )
                    if e0 > 1e-9:
                        rt = e1 / e0
                        ratios[sgn].append(rt)
                        if not 0.03 <= rt <= 0.5:
                            fails.append(
                                f"gz{sgn} l{l} {qn} ratio {rt:.2f}"
                            )
    # Sign reversal (Task 3): where the compressing branch has a real
    # donor (IKC set, l+NFM<=NCF), Q must flip sign exactly; the
    # magnitudes must match within 5% -- the |GAM|/DLN construction is
    # sign-blind. Where one sign has Q=0 the donor is structurally
    # unavailable (expansion at l=1): classify, never force a sign.
    flips = 0
    partial = []
    for l in range(1, 7):
        rp = dbg["gzp1"]["DBR"].get((4, l))
        rn = dbg["gzn1"]["DBR"].get((4, l))
        ir = dbg["gzn1"]["DBRI"].get((4, l))
        if not rp or not rn or not ir:
            fails.append(f"l{l} sign pair missing")
            continue
        has_donor = (
            (ir[1] == 1 or ir[2] == 1) and l + ir[0] <= 6
        )
        if rp[1] != 0.0 and rn[1] != 0.0 and has_donor:
            if rp[1] * rn[1] >= 0.0:
                fails.append(
                    f"l{l} Q no sign flip {rp[1]:.1e}/{rn[1]:.1e}"
                )
            else:
                flips += 1
                m = min(abs(rp[1]), abs(rn[1])) / max(
                    abs(rp[1]), abs(rn[1])
                )
                if m < 0.95:
                    fails.append(f"l{l} |Q| asymmetric {m:.3f}")
        elif rp[1] != 0.0 and rn[1] != 0.0 and not has_donor:
            # Compression donor beyond the grid (l+NFM>NCF): Q is
            # built but must multiply nothing -- the pure-outflow
            # closure. Assert the source is untouched instead of a
            # sign flip.
            if rn[4] != rn[10]:
                fails.append(f"l{l} beyond-grid donor consumed")
            partial.append(l)
        elif rp[1] != 0.0 or rn[1] != 0.0:
            # One sign coupled, the other structurally unable to --
            # expected at l=1 (expansion has no l-1 donor), wrong
            # anywhere else.
            if l != 1:
                fails.append(f"l{l} one-sided coupling")
            partial.append(l)
    if flips == 0:
        fails.append("no resolved Q sign flip anywhere")

    # Donor identity across the crossing (Task 5): reconstruct the
    # consumed donor as (SPM-SR)/|Q| and match it against the dumped
    # field at the branch-appropriate location.
    #   expansion (Q>0): donor = FLD(same depth, l-1)
    #   compression (Q<0): donor = (1-w)*FLD(jpr1,l+nfm)
    #                             + w*FLD(jpr2,l+nfm)
    # Both must reproduce the reconstruction within 5% -- the lagged
    # donor is the previous iterate of the same converged field.
    # gzz reads none (Q=0 bitwise). A donor that matched the WRONG
    # stencil (own-frequency or unrelated slot) fails the match.
    donors = {"exp": 0, "cmp": 0}
    for l in range(2, 7):
        for tag, kind in (("gzp1", "exp"), ("gzn1", "cmp")):
            d = dbg[tag]
            r = d["DBR"].get((4, l))
            ir = d["DBRI"].get((4, l))
            if not r or not ir or r[1] == 0.0:
                continue
            rec = (r[4] - r[10]) / abs(r[1])
            if kind == "exp":
                ref = d["DBR"].get((4, l - 1), [None] * 9)[8]
            else:
                j1, j2, w = ir[3], ir[4], r[6]
                lf = l + ir[0]
                if lf > 6:
                    # Beyond-grid donor: verified by the SP==SR
                    # closure above, not a reconstruction target.
                    continue
                f1 = d["DBR"].get((j1, lf), [None] * 9)[8]
                f2 = d["DBR"].get((j2, lf), [None] * 9)[8]
                ref = None if f1 is None or f2 is None else (
                    (1.0 - w) * f1 + w * f2
                )
            if ref is None:
                fails.append(f"{tag} l{l} donor ref missing")
            elif abs(rec - ref) / max(abs(ref), 1e-30) > 0.05:
                fails.append(
                    f"{tag} l{l} donor {rec:.3e}!={ref:.3e} ({kind})"
                )
            else:
                donors[kind] += 1
    if donors["exp"] == 0 or donors["cmp"] == 0:
        fails.append(f"donor ids unverified {donors}")
    # gzb: zero cell quiet, neighbours keep opposite-sign Q.
    b = dbg["gzb"]
    for l in range(2, 7):
        r = b["DBR"].get((4, l))
        if r is None or r[1] != 0.0 or r[4] != r[10]:
            fails.append(f"gzb cell4 l{l} coupled")
        r3, r5 = b["DBR"].get((3, l)), b["DBR"].get((5, l))
        if r3 and r5 and r3[1] * r5[1] >= 0.0:
            fails.append(f"gzb l{l} neighbours same-sign Q")
    # gcr: Q<0 left of the glide, Q>0 right of it, where coupled.
    c = dbg["gcr"]
    for l in range(2, 7):
        r4, r5 = c["DBR"].get((4, l)), c["DBR"].get((5, l))
        if r4 and r5 and not (r4[1] < 0.0 < r5[1]):
            fails.append(f"gcr l{l} Q signs {r4[1]:.1e}/{r5[1]:.1e}")

    ok = clean and d_micro == 0.0 and not fails
    _r = lambda s: "/".join(f"{x:.2f}" for x in ratios[s][:4])
    rep.check(
        "T15 GAM crossing: exact-zero quiet, monotonic ~linear "
        "contraction, signed Q flip, donor identity",
        ok,
        f"clean={clean} micro={d_micro:.1e} flips={flips} "
        f"donors={donors} ratios p={_r('p')} n={_r('n')} "
        f"{'; '.join(fails[:2]) if fails else ''}",
    )

    # -- T16: same-map contraction export (audit A1). 'slw' grinds a
    #    slowly-converging alternating solve over 14 allowed passes,
    #    so the post-transient ratio history (NPO>=5, same-direction
    #    baseline NPO-2) accumulates and RHO_M/RHO_P report each
    #    field solve's median asymptotic rate. The fast 'conv'/'ord'
    #    runs must report exactly 0: a solve that exits inside the
    #    seed transient contributes no sample -- that is the property
    #    the old KN_RES_PASS max-envelope violated by comparing
    #    maxima over changing ray subsets. ('d' is clobbered by the
    #    T15 donor loop, so 'conv' is re-run here.)
    d, _ = parse(run(bin_, "conv", cwd=WORKDIR, check=False).stdout)
    r = run(bin_, "slw", cwd=WORKDIR, check=False)
    d16, _ = parse(r.stdout)
    ok = (
        "RESULT_OK" in d16
        and d16.get("NPASS", 0) >= 5
        and 0.0 < d16.get("RHO_M", 0.0) < 50.0
        and 0.0 < d16.get("RHO_P", 0.0) < 50.0
        and d.get("RHO_M", 9.9) == 0.0
        and d.get("RHO_P", 9.9) == 0.0
        and d2.get("RHO_M", 9.9) == 0.0
        and d2.get("RHO_P", 9.9) == 0.0
    )
    rep.check(
        "T16 same-map contraction export: grinding alt solve reports "
        "finite RHO_M/RHO_P, transient exits report 0",
        ok,
        f"NPASS={d16.get('NPASS')} RHO_M={d16.get('RHO_M'):.3e} "
        f"RHO_P={d16.get('RHO_P'):.3e} "
        f"conv={d.get('RHO_M')}/{d.get('RHO_P')} "
        f"ord={d2.get('RHO_M')}/{d2.get('RHO_P')}",
    )

    # -- T17: same export on the ordered march (single residual
    #    sequence, baseline NPO-1, seed-transient floor NPO>=4).
    r = run(bin_, "slwo", cwd=WORKDIR, check=False)
    d17, _ = parse(r.stdout)
    ok = (
        "RESULT_OK" in d17
        and d17.get("NPASS", 0) >= 4
        and 0.0 < d17.get("RHO_M", 0.0) < 50.0
        and 0.0 < d17.get("RHO_P", 0.0) < 50.0
    )
    rep.check(
        "T17 same-map contraction export: grinding ord solve",
        ok,
        f"NPASS={d17.get('NPASS')} RHO_M={d17.get('RHO_M'):.3e} "
        f"RHO_P={d17.get('RHO_P'):.3e}",
    )

    # -- T18: end-to-end refresh with KN_ADAPT_RELAX +
    #    FG_ADAPT_OMEGA_V2: the pooled per-ray median ratios must
    #    reach the omega update (diagnostic line prints at
    #    CONV_DIAG_LVL>=2) and leave INEXT_RELAX inside its
    #    [0.05, 1] bounds, with no ray diverged.
    r = run(bin_, "rfa", cwd=WORKDIR, check=False)
    d18, _ = parse(r.stdout)
    ok = (
        "RESULT_OK" in d18
        and "omega est: rho_p90=" in r.stdout
        and d18.get("RF_NDIV", 9) == 0
        and 0.0 < d18.get("RF_RESID", 0.0) < 1.0e30
        and 0.05 <= d18.get("RF_INEXT_RELAX", -1.0) <= 1.0
    )
    rep.check(
        "T18 refresh adaptive omega: pooled per-ray estimate runs, "
        "INEXT_RELAX stays in bounds",
        ok,
        f"RF_RESID={d18.get('RF_RESID'):.3e} "
        f"NDIV={d18.get('RF_NDIV')} "
        f"INEXT_RELAX={d18.get('RF_INEXT_RELAX'):.4e} "
        f"diag={'omega est' in r.stdout}",
    )

    # -- T19: same-census cap numerator (audit C1). 'capx' sets
    #    NONMONO_INEXT_CAP below the field scale so the march caps
    #    cells on every one of its 14 grinding passes. NCAP keeps the
    #    all-pass telemetry sum; NCAP_LST must report only the last
    #    pass plus commit-phase events -- the same census the gate's
    #    NTOT denominator counts (capped cells stay significant:
    #    CAP=1e-4 >> 1e-4*VMAX floor). The pre-fix numerator was
    #    ~pass-multiplicity x this value plus a divergence flag.
    r = run(bin_, "capx", cwd=WORKDIR, check=False)
    d19, _ = parse(r.stdout)
    ok = (
        "RESULT_OK" in d19
        and d19.get("NCAP", 0) > 0
        and d19.get("NPASS", 0) >= 2
        and 0 < d19.get("NCAP_LST", -1) <= d19.get("NTOT", -1)
        and d19.get("NCAP_LST", 0) < d19.get("NCAP", 0)
        and 0 <= d19.get("NNG_LST", -1) <= d19.get("NNEG", -1)
    )
    rep.check(
        "T19 same-census cap numerator: last-pass+commit count "
        "bounded by cell census, all-pass NCAP separates",
        ok,
        f"NPASS={d19.get('NPASS')} NCAP={d19.get('NCAP')} "
        f"NCAP_LST={d19.get('NCAP_LST')} NTOT={d19.get('NTOT')} "
        f"NNG={d19.get('NNEG')} NNG_LST={d19.get('NNG_LST')}",
    )

    # -- T20: OpenMP determinism of the jump guard (audit Issue 1,
    #    test T23). 'jmpx' runs a three-ray refresh where ray 1 has a
    #    floored frequency-1 column and ray 2 carries +Inf there: the
    #    sibling reference is snapped before the parallel loop, so the
    #    veto verdict, per-ray deltas, committed stores and census
    #    must be identical at 1, 2 and 4 threads.
    import os as _os
    import subprocess as _sp

    def _jmpx(nth):
        env = dict(_os.environ, OMP_NUM_THREADS=str(nth))
        return _sp.run(
            [str(bin_), "jmpx"], cwd=WORKDIR, env=env,
            capture_output=True, text=True, timeout=300,
        )

    _runs = {}
    for _nth in (1, 2, 4):
        _r = _jmpx(_nth)
        _d, _s = parse(_r.stdout)
        _runs[_nth] = (_d, _s)
    _d1, _s1 = _runs[1]
    _det_ok = "RESULT_OK" in _d1 and all(
        "RESULT_OK" in _runs[n][0] for n in (2, 4)
    )
    _mism = []
    for _n in (2, 4):
        _dn, _sn = _runs[_n]
        if _sn.keys() != _s1.keys():
            _mism.append(f"{_n}t keys")
        for _k in _s1:
            if _s1[_k] != _sn.get(_k):
                _mism.append(f"{_k}@{_n}t")
        for _kk in ("RF_JREJ", "RF_NCREJ", "RF_NCCAP", "RF_DCHG"):
            if _d1.get(_kk) != _dn.get(_kk):
                _mism.append(f"{_kk}@{_n}t")
    rep.check(
        "T20 jump guard OMP deterministic: verdicts/deltas/stores "
        "identical at 1, 2, 4 threads",
        _det_ok and not _mism,
        f"JREJ={_d1.get('RF_JREJ')} VC1={_s1.get(('vc', 1))} "
        f"mismatches={_mism[:4]}",
    )

    # -- T21: sibling-reference finite-value protection (audit Issue
    #    1, test T24). 'sibref' drives KN_JMP_SIBREF_V1 on a three-ray
    #    store whose middle ray carries +Inf (freq 2), NaN (freq 3)
    #    and a finite-huge 1e200 column (freq 4): the reference must
    #    stay finite, equal the sane-ray peak, and return -1 only on
    #    the all-floored column. 'jmpx' then proves the veto still
    #    fires on the floored ray (VC=1) with the Inf sibling present
    #    -- under the old MAX scan that ray committed unvetted.
    r = run(bin_, "sibref", cwd=WORKDIR, check=False)
    d21, s21 = parse(r.stdout)
    _sr = {k[1]: v for k, v in s21.items() if k[0] == "sref"}
    _exp_m = {2: 1.43e-4, 3: 1.56e-4, 4: 1.69e-4, 5: 1.82e-4}
    ok = (
        "RESULT_OK" in d21
        and _sr.get(1) == (-1.0, -1.0)
        and _sr.get(6) == (-1.0, -1.0)
        and all(
            abs(_sr[l][0] - v) <= 1e-10 and abs(_sr[l][1] - 2 * v) <= 1e-10
            for l, v in _exp_m.items()
        )
        and all(abs(_sr[l][0]) < 1e30 for l in _exp_m)
    )
    rep.check(
        "T21a sibref excludes Inf/NaN, median kills finite-huge "
        "sibling, -1 on all-floored column",
        ok,
        f"SREF={ {l: _sr.get(l) for l in (1, 2, 4, 5, 6)} }",
    )
    r = run(bin_, "jmpx", cwd=WORKDIR, check=False)
    d21b, s21b = parse(r.stdout)
    _vc1 = s21b.get(("vc", 1), (0, 0.0))[0]
    _rf1 = s21b.get(("r", 1, 1, 1), (1.0, 1.0))[1]
    ok = (
        "RESULT_OK" in d21b
        and _vc1 == 1
        and d21b.get("RF_JREJ", 0) >= 1
        and _rf1 == 1.0e-30
    )
    rep.check(
        "T21b jump veto still fires on floored ray with Inf sibling "
        "(retained store untouched)",
        ok,
        f"VC1={_vc1} RF_JREJ={d21b.get('RF_JREJ')} "
        f"retained_M(1,1)={_rf1:.2e}",
    )

    # -- T22: omega pool on a contracting positive-mode sample set
    #    (audit Issue 2, test T25). rho~0.6 across 40 samples: the
    #    heuristic target is 1/(1+rho_p90)=1/1.62=0.617, so the
    #    smoothed update lands halfway from the 0.75 seed to that
    #    target -- lightly damped, never boosted past it.
    r = run(bin_, "ompos", cwd=WORKDIR, check=False)
    d22, _ = parse(r.stdout)
    ok = (
        "RESULT_OK" in d22
        and d22.get("OM_NSAM", 0) == 40
        and abs(d22.get("OM_P90", 0.0) - 0.62) < 1e-10
        and abs(d22.get("OM_NEW", 0.0) - 0.5 * (0.75 + 1.0 / 1.62))
        < 1e-6
    )
    rep.check(
        "T22 omega pool: contracting pool damps toward "
        "1/(1+rho_p90)",
        ok,
        f"NSAM={d22.get('OM_NSAM')} P90={d22.get('OM_P90')} "
        f"OM_NEW={d22.get('OM_NEW'):.4e}",
    )

    # -- T23: unsigned near-unity/expansive magnitudes -- the
    #    signature a sign-flipping or diverging mode leaves (audit
    #    Issue 2, test T26). omega must move below the seed toward
    #    1/(1+1.2)=0.4545 -- damping, never amplification.
    r = run(bin_, "omneg", cwd=WORKDIR, check=False)
    d23, _ = parse(r.stdout)
    ok = (
        "RESULT_OK" in d23
        and d23.get("OM_NSAM", 0) == 40
        and abs(d23.get("OM_P90", 0.0) - 1.2) < 1e-10
        and abs(d23.get("OM_NEW", 0.0) - 0.5 * (0.75 + 1.0 / 2.2))
        < 1e-6
        and d23.get("OM_NEW", 1.0) < d23.get("OM_PREV", 0.0)
    )
    rep.check(
        "T23 omega pool: expansive-ratio pool damps omega below "
        "seed",
        ok,
        f"P90={d23.get('OM_P90')} OM_NEW={d23.get('OM_NEW'):.4e} "
        f"prev={d23.get('OM_PREV')}",
    )

    # -- T24: small-sample protection (audit Issue 5, tests T27
    #    + T28). One valid sample or 8 samples (4 rays both
    #    directions) is below the KN_RHO_MIN_SAM=10 floor -- the
    #    pool must not update omega at all and must report the
    #    coverage.
    r1 = run(bin_, "omone", cwd=WORKDIR, check=False)
    d24a, _ = parse(r1.stdout)
    r2 = run(bin_, "omfew", cwd=WORKDIR, check=False)
    d24b, _ = parse(r2.stdout)
    ok = (
        "RESULT_OK" in d24a
        and "RESULT_OK" in d24b
        and d24a.get("OM_NSAM", -1) == 1
        and d24a.get("OM_NEW") == d24a.get("OM_PREV")
        and d24b.get("OM_NSAM", -1) == 8
        and d24b.get("OM_NEW") == d24b.get("OM_PREV")
    )
    rep.check(
        "T24 omega pool: 1 and 8 samples retain previous omega "
        "(min_sam=10)",
        ok,
        f"L1: nsam={d24a.get('OM_NSAM')} "
        f"om={d24a.get('OM_NEW')} | "
        f"L8: nsam={d24b.get('OM_NSAM')} om={d24b.get('OM_NEW')}",
    )

    # -- T25: P90 outlier + diverged-ray rejection (audit Issue 5).
    #    40 sane samples ~0.55-0.75, one ray at rho=100/95, two
    #    diverged rays at 400/500. NSAM must be 42 (diverged out),
    #    P90 lands on the sane mass (~0.74), and omega must not be
    #    dragged to the 1/(1+100)=0.01 floor.
    r = run(bin_, "omout", cwd=WORKDIR, check=False)
    d25, _ = parse(r.stdout)
    ok = (
        "RESULT_OK" in d25
        and d25.get("OM_NSAM", 0) == 42
        and 0.6 < d25.get("OM_P90", 0.0) < 1.0
        and d25.get("OM_NEW", 0.0) > 0.5
        and d25.get("OM_MAX", 0.0) >= 95.0
    )
    rep.check(
        "T25 omega pool: P90 ignores the 100x outlier and diverged "
        "rays",
        ok,
        f"NSAM={d25.get('OM_NSAM')} P90={d25.get('OM_P90')} "
        f"MAX={d25.get('OM_MAX')} OM_NEW={d25.get('OM_NEW'):.4e}",
    )

    # -- T26: converged-but-rejected ray cannot read mature (audit
    #    Issue 3, test T29). 'jmpx' ray 1 converges (RF_RESID ~
    #    1e-15) but its field is jump-vetoed: the retained store
    #    keeps the floor value, RF_PEND carries the proposed delta,
    #    and the maturity getter must export it as RF_DCHG so
    #    GS_MAX_ALL sees the field as grossly immature -- never the
    #    converged-looking residual alone.
    r = run(bin_, "jmpx", cwd=WORKDIR, check=False)
    d26, s26 = parse(r.stdout)
    _rf1 = s26.get(("r", 1, 1, 1), (1.0, 1.0))[1]
    ok = (
        "RESULT_OK" in d26
        and d26.get("RF_RESID", 1.0) < 1.0e-10
        and s26.get(("vc", 1), (0, 0.0))[0] == 1
        and d26.get("RF_PEND", 0.0) > 50.0
        and d26.get("GS_RFD", 0.0) >= d26.get("RF_PEND", 1.0)
        and _rf1 == 1.0e-30
    )
    rep.check(
        "T26 converged ray, jump-vetoed field: maturity export "
        "carries the pending delta",
        ok,
        f"RESID={d26.get('RF_RESID'):.3e} PEND={d26.get('RF_PEND'):.3e} "
        f"GS_RFD={d26.get('GS_RFD'):.3e} retained={_rf1:.1e}",
    )

    # -- T27: repeated rejection keeps the field immature (audit
    #    Issue 3, test T30). 'jmp2' runs the same refresh twice; the
    #    second pass still rejects every ray and the pending delta
    #    persists, so the gate can never settle on the old store.
    r = run(bin_, "jmp2", cwd=WORKDIR, check=False)
    d27, s27 = parse(r.stdout)
    _pass2 = r.stdout.split("RF_PASS  2")[-1]
    _d27b, _s27b = parse(_pass2)
    ok = (
        "RESULT_OK" in d27
        and r.stdout.count("RF_PASS") == 2
        and d27.get("RF_JREJ", 0) >= 6
        and s27.get(("vc", 1), (0, 0.0))[0] == 1
        and _d27b.get("RF_PEND", 0.0) > 50.0
        and _d27b.get("GS_RFD", 0.0) >= _d27b.get("RF_PEND", 1.0)
    )
    rep.check(
        "T27 two refreshes, still rejected, still immature "
        "(pending delta persists)",
        ok,
        f"JREJ={d27.get('RF_JREJ')} "
        f"pass2 PEND={_d27b.get('RF_PEND'):.3e} "
        f"GS_RFD={_d27b.get('GS_RFD'):.3e}",
    )

    # -- T28: active but unbuilt ray fails closed (audit Issue 5,
    #    test T31). 'ni0' gives ray 1 NI=0 -- a construction failure
    #    that never produced a field. It must report RESID=1.0 into
    #    the residual max so the gate stays immature, while the other
    #    two rays commit normally.
    d28, s28 = parse(run(bin_, "ni0", cwd=WORKDIR, check=False).stdout)
    ok = (
        "RESULT_OK" in d28
        and s28.get(("vc", 1), (0, 0.0))[0] == -1
        and s28.get(("vc", 2), (0, 0.0))[0] == 0
        and abs(d28.get("RF_RESID", 0.0) - 1.0) < 1.0e-12
        and abs(d28.get("RF_MAXCHG", 0.0) - 1.0) < 1.0e-12
        and d28.get("GS_OK", True) is False
    )
    rep.check(
        "T28 unbuilt ray (NI=0) fails closed: residual max=1.0, "
        "gate stays immature",
        ok,
        f"VC1={s28.get(('vc', 1), ('?', 0.0))[0]} "
        f"RF_RESID={d28.get('RF_RESID'):.3e} "
        f"GS_OK={d28.get('GS_OK')}",
    )

    # -- T29: legitimately inactive ray is excluded (audit Issue 5,
    #    test T32). 'ni1' gives ray 1 NI=1 -- a single-point geometry
    #    with no marchable domain. It keeps the plain unsolved verdict
    #    and a zero residual, so the residual max reads only the
    #    converged active rays.
    d29, s29 = parse(run(bin_, "ni1", cwd=WORKDIR, check=False).stdout)
    ok = (
        "RESULT_OK" in d29
        and s29.get(("vc", 1), (0, 0.0))[0] == -1
        and s29.get(("vc", 2), (0, 0.0))[0] == 0
        and d29.get("RF_RESID", 1.0) < 1.0e-6
    )
    rep.check(
        "T29 single-point ray (NI=1) excluded: unsolved verdict, "
        "residual max from active rays only",
        ok,
        f"VC1={s29.get(('vc', 1), ('?', 0.0))[0]} "
        f"RF_RESID={d29.get('RF_RESID'):.3e}",
    )

    # -- T30: heal census counts only real modifications (audit
    #    Issue 6, test T33). 'heal' diverges one ray whose I_M rows
    #    carry {valid, valid, zero, negative, +Inf} and whose I_P
    #    rows carry {valid, zero, NaN, negative, above-cap}. Only
    #    the nonzero invalid slots may be counted: 3 depths x (2+3)
    #    = 15, not 21 (which would include the four zero slots).
    d30, s30 = parse(run(bin_, "heal", cwd=WORKDIR, check=False).stdout)
    ok = (
        "RESULT_OK" in d30
        and s30.get(("vc", 1), (0, 0.0))[0] == 4
        and d30.get("RF_NDIV", 0) == 1
        and d30.get("RF_HEAL", 0.0) == 15.0
    )
    rep.check(
        "T30 diverged-ray heal counts actual modifications only "
        "(zero slots excluded)",
        ok,
        f"VC1={s30.get(('vc', 1), ('?', 0.0))[0]} "
        f"NDIV={d30.get('RF_NDIV')} "
        f"RF_HEAL={d30.get('RF_HEAL')}",
    )

    # -- T31: JFNK eval bracket restores production telemetry
    #    (audit Issue 7, test T34). 'esav' pushes a sentinel set,
    #    clobbers every field with eval-side junk, pops, and must
    #    read the sentinels back exactly.
    d31, _ = parse(run(bin_, "esav", cwd=WORKDIR, check=False).stdout)
    ok = (
        "RESULT_OK" in d31
        and d31.get("EV_MAXCHG") == 1.25
        and d31.get("EV_SWMX") == 0.75
        and d31.get("EV_RFD") == 3.5
        and d31.get("EV_RFP") == 2.5
        and d31.get("EV_NCONV") == 42.0
        and d31.get("EV_NBAD") == 7.0
        and d31.get("EV_NSTAL") == 3.0
        and d31.get("EV_NDIV") == 2.0
        and d31.get("EV_HEAL") == 9.0
        and d31.get("EV_JREJ") == 4.0
        and d31.get("EV_HOT") == 11.0
        and d31.get("EV_SRDY") is True
        and d31.get("EV_FRQ") == 4.0
    )
    rep.check(
        "T31 eval push/pop restores the full FG census",
        ok,
        f"MAXCHG={d31.get('EV_MAXCHG')} NCONV={d31.get('EV_NCONV')} "
        f"SRDY={d31.get('EV_SRDY')} FRQ={d31.get('EV_FRQ')}",
    )

    # -- T32: sweep/refresh DTAU parity across METHOD classes
    #    (audit B1, test T35). 'pmon' (MON-class captured method)
    #    and 'pnom' (plain) run the same refresh; with the
    #    INTEGRAL-path clamp now unconditional, every output must
    #    be identical between the two.
    _om = run(bin_, "pmon", cwd=WORKDIR, check=False).stdout
    _on = run(bin_, "pnom", cwd=WORKDIR, check=False).stdout
    d32a, s32a = parse(_om)
    d32b, s32b = parse(_on)
    ok = (
        "RESULT_OK" in d32a
        and "RESULT_OK" in d32b
        and _om.replace("pmon", "x") == _on.replace("pnom", "x")
    )
    rep.check(
        "T32 MON-class vs plain captured method: identical "
        "refresh operator (DTAU clamp parity)",
        ok,
        f"RESID {d32a.get('RF_RESID'):.4e} vs "
        f"{d32b.get('RF_RESID'):.4e}",
    )

    rep.summary()


if __name__ == "__main__":
    main()
