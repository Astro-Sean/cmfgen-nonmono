#!/usr/bin/env python3
"""cmfgen_jfnk unit tests: surrogate statistical-equilibrium solve.

Drives the production newton_hook_v1 / gmres_m_v1 through the
se_driver CLI on small surrogate systems (systems/{mild,std,hard})
that reproduce the failure regimes seen in m15_2p82d_SOB runs:

  * stiff SE rows (|F| ~ 1e13) mixed with weak trace-level rows
  * near-proportional Rydberg rows -> near-singular J directions,
    demanded |s|_inf >> step cap
  * floor-level populations (n0 = 1e-30) -> near-null BA rows/cols
  * fixed T/Ne rows, conservation row, structurally-constant row

Tests:
  T0  Fortran vs Python residual identity (F(0) test analogue)
  T1  residual norm anatomy: which rows hold ||F||
  T2  FD JVP error vs eps (production step rule), global + per-class
  T3  demanded Newton step: exact J0 s = -F0 via scipy, |s|_inf, SVD
  T4  preconditioner accuracy: MZDIFF / ZCORE vs expneg, ptcfix, x!=0
  T5  production GMRES on J s = -F(0): pc variants
  T6  production newton_hook_v1 end-to-end: mild should converge,
      hard reproduces the production stall; variant table
  T7  reachable residual floor on the const rows (unsolvable mass)
  ...
  T15 stale-baseline JVP: offset between the stored and evaluated
      residual must appear as delta/h on the offset row only (the
      F_outer-vs-F_eval(0) bug signature)
  T16 structure-flip freeze (systems/flip): a live equation branch
      corrupts the FD column by O(gap/h); replaying the baseline
      branch restores FD-vs-exact agreement (USE_PASSED_REP=T)
  T17 repeat determinism: identical JVP probes are byte-identical
  T18 lagged factor: freeze the synthetic aggregate factor and recover
      the baseline BA derivative without its live low-rank term
  T19 fallback dispatch: preserve a defined caller cap and mode on
      finite stagnation
  T20 evaluation telemetry: report residual callbacks made inside JVPs
      separately from Newton-level residual calls

Usage:  python3 test_se_solver.py [--quick]
Exit status: number of failed tests.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

TD = Path(__file__).resolve().parent
BIN = TD / "bin" / "se_driver"
SYSDIR = TD / "systems"
PLOTDIR = TD / "plots"


def load_spec(name):
    d = SYSDIR / name
    nt, nd, n = np.loadtxt(d / "params.txt", dtype=int)
    a = np.loadtxt(d / "acoll.txt").T
    c = np.loadtxt(d / "cmat.txt").T
    r = np.loadtxt(d / "rvec.txt")
    n0 = np.loadtxt(d / "n0.txt")
    xs = np.loadtxt(d / "xstar.txt")
    mask = np.loadtxt(d / "mask.txt", dtype=int)
    esc = np.loadtxt(d / "escale.txt")
    meta = json.loads((d / "meta.json").read_text())
    return dict(dir=d, nt=int(nt), nd=int(nd), n=int(n), a=a, c=c, r=r,
                n0=n0, xstar=xs, mask=mask, escale=esc, meta=meta)


def se_f(spec, x):
    xx = np.clip(x, -80.0, 80.0)
    f = spec["a"] @ (spec["n0"] * np.exp(xx)) + spec["c"] @ x + spec["r"]
    f = f.copy()
    f[spec["mask"] == 1] = 0.0
    return np.clip(f, -1e20, 1e20)


def se_j(spec, x):
    xx = np.clip(x, -80.0, 80.0)
    j = spec["a"] * (spec["n0"] * np.exp(xx))[None, :] + spec["c"]
    j = j.copy()
    j[spec["mask"] == 1, :] = 0.0
    return j


def run_cli(*args):
    r = subprocess.run([str(BIN)] + [str(a) for a in args],
                       capture_output=True, text=True, check=True,
                       timeout=600)
    return r.stdout


def parse_keyvals(out):
    kv = {}
    for ln in out.splitlines():
        t = ln.split()
        if len(t) == 2 and t[0][0].isalpha() and \
                all(ch.isalnum() or ch == '_' for ch in t[0]):
            try:
                kv[t[0]] = float(t[1])
            except ValueError:
                pass
    return kv


def parse_vec(out, tag):
    rows = {}
    for ln in out.splitlines():
        t = ln.split()
        if len(t) >= 2 and t[0] == tag:
            rows[int(t[1]) - 1] = [float(v) for v in t[2:]]
    return rows


def class_of(spec, i):
    return spec["meta"]["row_class"][i]


def t0_identity(spec):
    """Fortran sef vs python F at x=0 and at a perturbed x."""
    out = run_cli("sef", spec["dir"], 0)
    kv = parse_keyvals(out)
    frows = parse_vec(out, "F")
    f_f = np.array([frows[i][0] for i in range(spec["n"])])
    f_p = se_f(spec, np.zeros(spec["n"]))
    rel = np.linalg.norm(f_f - f_p) / max(np.linalg.norm(f_p), 1e-300)
    ok = rel < 1e-12 and abs(kv["FNORM"] - np.linalg.norm(f_p)) < \
        1e-6 * np.linalg.norm(f_p)
    return ok, f"rel diff {rel:.2e}", dict(f=f_p)


def t1_anatomy(spec, f):
    """Where does ||F|| live?  Prints the top-k table."""
    a = np.abs(f)
    order = np.argsort(-a)
    tot = np.linalg.norm(f)
    top = []
    print("    rank   i    |F_i|        class   share%")
    for kk in range(10):
        i = order[kk]
        share = a[i] ** 2 / tot ** 2 * 100
        top.append((i, a[i], class_of(spec, i), share))
        print(f"    {kk+1:4d} {i:4d}  {a[i]:.3e}  {class_of(spec, i):>6s} "
              f"{share:6.2f}")
    share5 = sum(t[3] for t in top[:5])
    # mass carried by regularized-candidate rows (rowmax < 1e-30 blk)
    j0 = se_j(spec, np.zeros(spec["n"]))
    info = {"top5_share": share5}
    print(f"    top-5 rows carry {share5:.1f}% of ||F||^2; "
          f"max |F_i|/|F| = {a.max()/tot:.2f}")
    return True, f"top5 share {share5:.0f}%", info


def t2_jvp(spec):
    """FD JVP error vs eps, at x=0 and x=x*/2, global and weak-row."""
    rng = np.random.default_rng(5)
    v = rng.uniform(-1, 1, spec["n"])
    v = v / np.linalg.norm(v)
    vfile = TD / "systems" / "_v.txt"
    np.savetxt(vfile, v, fmt="%.16e")
    rows = []
    for tag, x in (("x=0", np.zeros(spec["n"])),
                   ("x=x*/2", 0.5 * spec["xstar"])):
        xfile = TD / "systems" / "_x.txt"
        np.savetxt(xfile, np.clip(x, -80, 80), fmt="%.16e")
        jv_ex = se_j(spec, x) @ v
        for eps in (np.sqrt(np.finfo(float).eps), 1e-6, 1e-4, 1e-2):
            out = run_cli("sejvp", spec["dir"], vfile, eps, xfile)
            kv = parse_keyvals(out)
            rows.append((tag, eps, kv["RELDIFF"],
                         kv["JVFDNORM"], kv["JVEXNORM"]))
    print(f"    {'point':>7s} {'eps':>9s} {'rel|FD-ex|':>11s} "
          f"{'|FD|':>11s} {'|ex|':>11s}")
    for tag, eps, rel, nfd, nex in rows:
        print(f"    {tag:>7s} {eps:9.1e} {rel:11.3e} {nfd:11.3e} "
              f"{nex:11.3e}")
    best = min(r[2] for r in rows[:4])
    ok = best < 0.5  # at x=0 some eps should give usable action
    return ok, f"best rel err at x=0 {best:.2e}", {"rows": rows}


def t3_demand(spec, f):
    """Demanded Newton correction: exact solve of J0 s = -F0."""
    from scipy import linalg
    j0 = se_j(spec, np.zeros(spec["n"]))
    sv = np.linalg.svd(j0, compute_uv=False)
    s, *_ = linalg.lstsq(j0, -f, cond=1e-14)
    s_inf = np.abs(s).max()
    i_top = int(np.abs(s).argmax())
    print(f"    sigma range {sv[-1]:.2e} .. {sv[0]:.2e}  "
          f"cond {sv[0]/max(sv[-1],1e-300):.2e}")
    print(f"    demanded |s|_inf = {s_inf:.3e} at row {i_top} "
          f"({class_of(spec, i_top)}); |s|_2 = {np.linalg.norm(s):.3e}")
    nbig = int((np.abs(s) > 10).sum())
    print(f"    {nbig} of {spec['n']} components exceed step cap 10")
    ok = s_inf > 1e3  # hard problem: demands far beyond cap
    return ok, f"|s|_inf={s_inf:.1e}, {nbig} comps > cap", \
        {"s": s, "sv": sv}


def t4_precond(spec):
    """Preconditioner accuracy: core solve + post-stage effect."""
    rows = []
    # x=0: production chain vs switches
    for expneg, ptc, ptcfix, tag in (
            (1, 0, 0, "prod: expneg on"),
            (0, 0, 0, "expneg off"),
            (0, 1, 1, "expneg off + ptc fix"),
            (1, 1, 1, "prod + ptc")):
        out = run_cli("sepc", spec["dir"], 2, expneg, ptc, ptcfix, 0)
        kv = parse_keyvals(out)
        rows.append(("x=0", tag, kv))
    # x = x*/2 (Newton moved): expneg should hurt
    xfile = TD / "systems" / "_xmid.txt"
    np.savetxt(xfile, np.clip(0.5 * spec["xstar"], -80, 80), fmt="%.16e")
    for expneg, tag in ((1, "expneg on"), (0, "expneg off")):
        out = run_cli("sepc", spec["dir"], 2, expneg, 0, 0, xfile)
        kv = parse_keyvals(out)
        rows.append(("x=x*/2", tag, kv))
    print(f"    {'point':>7s} {'variant':>22s} {'MZDIFF':>9s} "
          f"{'ZINF':>10s} {'ZCORE(post vs core)':>18s}")
    for pt, tag, kv in rows:
        print(f"    {pt:>7s} {tag:>22s} {kv['MZDIFF']:9.2e} "
              f"{kv['ZINF']:10.2e} {kv['ZCORE']:18.2e}")
    # algebraic B^T w apply (production LSQ form) vs direct solve;
    # rowtiny nonzero exercises the reg-row correction branch
    lsq_ok = True
    for rt in (0.0, 1e-30, 1e-6):
        out = run_cli("lsqchk", spec["dir"], 1e-4, rt)
        kv = parse_keyvals(out)
        print(f"    lsqchk rowtiny={rt:8.0e}: RELDIFF {kv['RELDIFF']:.2e}")
        lsq_ok = lsq_ok and kv["RELDIFF"] < 1e-10
    ok = all(r[2]["MZDIFF"] < 0.1 for r in rows) and lsq_ok
    return ok, "core LU solve residual <10%", {"rows": rows}


def t5_gmres(spec):
    """Production GMRES on J s = -F(0): pc variants."""
    curves = {}
    for pc, expneg, ptc, ridge, tag in ((0, 0, 0, 0, "ident"),
                               (1, 0, 0, 0, "jacobi"),
                               (2, 1, 0, 0, "LU prod (expneg)"),
                               (2, 0, 0, 0, "LU noexpneg"),
                               (2, 0, 1, 0, "LU noexpneg + ptc"),
                               (3, 0, 0, 1e-8, "TSVD rcond=1e-8"),
                               (3, 0, 0, 1e-4, "TSVD rcond=1e-4"),
                               (4, 0, 1, 0, "PTC-LU (J+D/dtau)"),
                               (5, 0, 0, 1e-4, "LSQ ridge=1e-4"),
                               (5, 0, 0, 1e-2, "LSQ ridge=1e-2"),
                               (6, 0, 0, 1e-2, "LSQ-alg ridge=1e-2")):
        out = run_cli("segmres", spec["dir"], 10, 30, pc,
                      1e-2, expneg, ptc, 1, ridge)
        rows = [ln.split() for ln in out.splitlines()
                if ln and ln[0].isdigit()]
        ks = np.array([int(r[0]) for r in rows])
        rr = np.array([float(r[2]) for r in rows])
        curves[tag] = (ks, rr / max(rr[0], 1e-300))
        print(f"    {tag:>22s}: res {rr[0]:.2e} -> {rr[-1]:.2e} "
              f"in {ks[-1]} its")
    return True, "curves printed", {"curves": curves}


def t6_newton(spec, xstar):
    """Production newton_hook_v1 end-to-end; variant grid."""
    res = {}
    base = [spec["dir"], 20, 30, 2, 1e-2, 1e-30, 0, 1, 1, 1, 0, 10, 1e-4]
    variants = [
        ("prod config", {}),
        ("expneg off", {7: 0}),
        ("rowtiny 1e-10", {5: 1e-10}),
        ("eqscale", {10: 1}),
        ("cap 100", {11: 100}),
        ("ptc off", {8: 0, 9: 0}),
        ("TSVD rcond=1e-8", {3: 3, 6: 1e-8, 7: 0}),
        ("TSVD rcond=1e-4", {3: 3, 6: 1e-4, 7: 0}),
        ("PTC-LU", {3: 4, 7: 0}),
        ("LSQ ridge=1e-4", {3: 5, 6: 1e-4, 7: 0}),
        ("LSQ ridge=1e-2", {3: 5, 6: 1e-2, 7: 0}),
        ("LSQ-alg ridge=1e-2", {3: 6, 6: 1e-2, 7: 0}),
        ("eps=sqrt(epsm)", {4: np.sqrt(np.finfo(float).eps)}),
    ]
    for tag, patch in variants:
        args = list(base)
        for k, v in patch.items():
            args[k] = v
        try:
            out = run_cli("senewton", *args)
        except subprocess.CalledProcessError as e:
            print(f"    {tag:>24s}: driver failed rc={e.returncode}")
            res[tag] = None
            continue
        kv = parse_keyvals(out)
        hist = [ln.split() for ln in out.splitlines()
                if ln.startswith("SUB")]
        hist = [float(h[1]) for h in hist]
        x_f = np.array([float(ln) for ln in
                        out.split("X\n")[1].splitlines() if ln.strip()])
        xerr = np.linalg.norm(x_f - xstar)
        red = kv.get("FNORM1", np.nan) / kv.get("FNORM0", np.nan)
        res[tag] = dict(kv=kv, hist=hist, red=red, xerr=xerr)
        print(f"    {tag:>24s}: info={int(kv.get('INFO', -1))} "
              f"|F| {kv.get('FNORM0', np.nan):.2e} -> "
              f"{kv.get('FNORM1', np.nan):.2e} (x{red:.3f}) "
              f"nits={int(kv.get('NITS', -1))} "
              f"f/Jv evals={int(kv.get('FEVALS', 0))}/"
              f"{int(kv.get('JVEVALS', 0))} |x-x*|={xerr:.2e}")
    ok = res.get("prod config") is not None
    return ok, "ran variants", res


def t7_floor(spec, f):
    """Residual floor: mass on structurally-constant (unsolvable) rows.
    Even a perfect solver cannot reduce F_i on rows with zero Jacobian.
    """
    j0 = se_j(spec, np.zeros(spec["n"]))
    rowmax = np.abs(j0).max(axis=1)
    blk = spec["nt"]
    const_idx = [i for i in range(spec["n"])
                 if rowmax[i] == 0.0 and abs(f[i]) > 0]
    fconst = np.linalg.norm(f[const_idx]) if const_idx else 0.0
    tot = np.linalg.norm(f)
    print(f"    {len(const_idx)} constant rows, floor ||F|| = {fconst:.3e} "
          f"of {tot:.3e} ({fconst/tot*100:.2f}% of norm)")
    # achievable floor: also rows whose Jacobian is below 1e-30*blkmax
    # (reg candidates at production rowtiny) -- can they be solved?
    reg = []
    for k in range(spec["nd"]):
        lo, hi = k * blk, (k + 1) * blk
        bmax = rowmax[lo:hi].max()
        for i in range(lo, hi):
            if 0 < rowmax[i] < 1e-30 * bmax and abs(f[i]) > 0:
                reg.append(i)
    freg = np.linalg.norm(f[reg]) if reg else 0.0
    print(f"    {len(reg)} rowtiny-reg rows hold {freg:.3e} "
          f"({freg/tot*100:.2f}%)")
    ok = True
    return ok, f"floor={fconst:.2e} ({fconst/tot*100:.1f}%)", \
        {"fconst": fconst, "freg": freg}


def _write_vec(path, arr):
    np.savetxt(path, np.asarray(arr, dtype=float), fmt="%.16e")


def make_partial_mask_spec(spec, n_masked_depths):
    """Spec variant where only the first n_masked_depths depths keep
    their T/Ne mask -- the partial-range FIX_T analog."""
    out = TD / "systems" / ("_pm_" + spec["dir"].name)
    out.mkdir(exist_ok=True)
    for f in spec["dir"].iterdir():
        if f.name != "mask.txt":
            (out / f.name).write_bytes(f.read_bytes())
    mask = spec["mask"].copy()
    nt = spec["nt"]
    for i in range(spec["n"]):
        if mask[i] == 1 and i // nt >= n_masked_depths:
            mask[i] = 0
    np.savetxt(out / "mask.txt", mask, fmt="%d")
    return out


def t8_masked(spec):
    """Masked-T projection semantics: fixed rows are exactly zero in
    F, JVP and the preconditioned step; corrections never move the
    masked components; a partial depth range masks only its rows."""
    n, nt, nd = spec["n"], spec["nt"], spec["nd"]
    midx = np.where(spec["mask"] == 1)[0]
    fidx = np.where(spec["mask"] == 0)[0]
    notes = []
    ok = True
    rng = np.random.default_rng(11)

    # residual projection: perturbed x including masked comps -> the
    # masked rows are still exactly zero (they never carry residual)
    x = rng.uniform(-2, 2, n)
    x[midx] = 5.0
    xf = TD / "systems" / "_xm.txt"
    _write_vec(xf, np.clip(x, -80, 80))
    out = run_cli("sef", spec["dir"], xf)
    frows = parse_vec(out, "F")
    fv = np.array([frows[i][0] for i in range(n)])
    a = np.all(fv[midx] == 0.0) and np.any(fv[fidx] != 0.0)
    ok &= a
    notes.append(f"F masked rows zero: {a}")

    # JVP: production feeds the FD a pre-masked v, so masked columns
    # are never probed; check masked OUTPUT rows are zero both in the
    # FD result and the analytic product
    v = rng.uniform(-1, 1, n)
    v[midx] = 0.0
    vf = TD / "systems" / "_vm.txt"
    _write_vec(vf, v)
    out = run_cli("sejvp", spec["dir"], vf, 1e-4, xf)
    jvrows = parse_vec(out, "JV")
    jvfd = np.array([jvrows[i][0] for i in range(n)])
    jvex = np.array([jvrows[i][1] for i in range(n)])
    a = np.all(jvfd[midx] == 0.0) and np.all(jvex[midx] == 0.0)
    ok &= a
    notes.append(f"JVP masked rows zero: {a}")

    # preconditioned step: masked comps zero, free comps nonzero
    out = run_cli("sepc", spec["dir"], 2, 1, 0, 0, 0)
    zrows = parse_vec(out, "Z")
    zv = np.array([zrows[i][0] for i in range(n)])
    a = np.all(zv[midx] == 0.0) and np.any(zv[fidx] != 0.0)
    ok &= a
    notes.append(f"PC masked comps zero: {a}")

    # Newton never moves masked components: parse final x
    out = run_cli("senewton", spec["dir"], 20, 30, 2, 1e-2, 1e-30,
                  0, 1, 1, 1, 0, 10, 1e-4)
    x_f = np.array([float(ln) for ln in
                    out.split("X\n")[1].splitlines() if ln.strip()])
    a = np.all(x_f[midx] == 0.0)
    ok &= a
    notes.append(f"Newton x masked comps stay 0: {a}")

    # partial range: mask only the first nd//2 depths
    pm_dir = make_partial_mask_spec(spec, nd // 2)
    out = run_cli("sef", pm_dir, xf)
    frows = parse_vec(out, "F")
    fv2 = np.array([frows[i][0] for i in range(n)])
    mlo = np.array([i for i in midx if i // nt < nd // 2])
    mhi = np.array([i for i in midx if i // nt >= nd // 2])
    a = (np.all(fv2[mlo] == 0.0) and mhi.size > 0
         and np.all(fv2[mhi] != 0.0))
    ok &= a
    notes.append(f"partial range masks only its rows: {a}")
    return ok, "; ".join(notes), {}


def t9_jvp_stiff(spec):
    """FD JVP accuracy across conditioning: directions along the
    stiffest, mid and near-null right-singular vectors of J0.  The
    stiff direction must be accurate; the near-null direction shows
    the roundoff floor eps_mach*|F|/(h|Jv|) -- a large error there is
    conditioning, not an implementation bug, so it is checked against
    the predicted bound rather than a fixed tolerance."""
    n = spec["n"]
    rng = np.random.default_rng(7)
    j0 = se_j(spec, np.zeros(n))
    _, sv, vt = np.linalg.svd(j0)
    f0 = se_f(spec, np.zeros(n))
    fnrm = np.linalg.norm(f0)
    vrand = rng.uniform(-1, 1, n)
    vrand /= np.linalg.norm(vrand)
    dirs = {"stiff": vt[0], "mid": vt[n // 2], "null": vt[-1],
            "rand": vrand}
    epss = (1e-2, 1e-6, np.sqrt(np.finfo(float).eps))
    ok = True
    print(f"    {'dir':>6s} {'|Jv|':>11s} " +
          " ".join(f"rel@{e:8.0e}" for e in epss))
    for tag, v in dirs.items():
        v = v / np.linalg.norm(v)
        v[spec["mask"] == 1] = 0.0
        nv = np.linalg.norm(v)
        if nv < 1e-14:
            continue
        v = v / nv
        jv_ex = j0 @ v
        jvn = np.linalg.norm(jv_ex)
        vf = TD / "systems" / "_vj.txt"
        _write_vec(vf, v)
        rels = []
        for eps in epss:
            out = run_cli("sejvp", spec["dir"], vf, eps, 0)
            kv = parse_keyvals(out)
            rels.append(kv["RELDIFF"])
            if not np.isfinite(kv["RELDIFF"]):
                ok = False
        # Roundoff-dominated floor: |FD err| ~ eps_mach*|F|/h with
        # h = eps (x=0, |v|=1).  Threshold per (dir,eps): the fixed
        # truncation tolerance where the direction carries signal,
        # or 100x the predicted roundoff bound where it does not --
        # a bug shows as error far above the bound, conditioning does
        # not.
        bound = [np.finfo(float).eps * fnrm
                 / max(eps * jvn, 1e-300) for eps in epss]
        thresh = [max(t, 100.0 * b)
                  for t, b in zip((1e-2, 1e-4, 1.0), bound)]
        print(f"    {tag:>6s} {jvn:11.3e} " +
              " ".join(f"{r:11.3e}" for r in rels) +
              f"   (floor {bound[-1]:.1e})")
        ok &= all(r < t for r, t in zip(rels, thresh))
    return ok, "dir x eps grid printed", {}


def t10_lsq_grid(spec):
    """LSQ vs LSQ-alg consistency across ridge and regularization
    states; reports componentwise agreement of the two applies."""
    worst = 0.0
    nbit = 0
    ntot = 0
    ok = True
    print(f"    {'ridge':>8s} {'rowtiny':>8s} {'RELDIFF':>11s} "
          f"{'max|dz|':>11s} {'bitid%':>7s}")
    for ridge in (1e-6, 1e-4, 1e-2, 1e-1):
        for rt in (0.0, 1e-30, 1e-6):
            out = run_cli("lsqchk", spec["dir"], ridge, rt)
            kv = parse_keyvals(out)
            zrows = parse_vec(out, "Z")
            z1 = np.array([zrows[i][0] for i in range(spec["n"])])
            z2 = np.array([zrows[i][1] for i in range(spec["n"])])
            nb = int(np.sum(z1 == z2))
            nbit += nb
            ntot += spec["n"]
            mdz = float(np.abs(z1 - z2).max())
            worst = max(worst, kv["RELDIFF"])
            ok &= kv["RELDIFF"] < 1e-10
            print(f"    {ridge:8.0e} {rt:8.0e} {kv['RELDIFF']:11.3e} "
                  f"{mdz:11.3e} {nb/spec['n']*100:6.1f}%")
    return ok, f"worst RELDIFF {worst:.2e}, bit-identical " \
        f"{nbit/ntot*100:.1f}% of comps", {}


def t11_strike(spec):
    """Budget-ladder replica (sestrike): a counted failure raises the
    rung (capped at the argument), a commit resets to 0, a reach-skip
    is exempt -- and Newton runs on EVERY admission (no pause state).
    Requires a spec where the LU admission fails (info=3) and the
    TSVD admission reaches a committed best iterate (info=2) -- only
    'hard' does both."""
    ok = True
    notes = []

    def run(seq, rung_max=4):
        out = run_cli("sestrike", spec["dir"], rung_max, seq)
        adm = {}
        for ln in out.splitlines():
            t = ln.split()
            if t and t[0] == "ADM":
                adm[int(t[1])] = dict(info=int(t[3]), rung=int(t[5]))
        return adm

    # failures climb the ladder one rung each, capped at rung_max;
    # every admission still ran Newton (all ADM rows exist)
    adm = run("FFFFFFFFF")
    a = all(i in adm for i in range(1, 10))
    a &= [adm[i]["rung"] for i in range(1, 7)] == [1, 2, 3, 4, 4, 4]
    ok &= a
    notes.append(f"always-try ladder climbs to cap: {a}")

    # committed update resets the rung to 0
    adm = run("FFSFFFF")
    a = (adm[3]["info"] in (0, 2) and adm[3]["rung"] == 0
         and adm[4]["rung"] == 1 and adm[7]["rung"] == 4)
    ok &= a
    notes.append(f"commit resets rung: {a}")

    # reach-gate skip (info=-2) leaves the rung unchanged
    adm = run("RFFRFF")
    a = (adm[1]["info"] == -2 and adm[1]["rung"] == 0
         and adm[2]["rung"] == 1 and adm[3]["rung"] == 2
         and adm[4]["info"] == -2 and adm[4]["rung"] == 2
         and adm[6]["rung"] == 4)
    ok &= a
    notes.append(f"reach-skip exempt, rung unchanged: {a}")
    return ok, "; ".join(notes), {}


def _auto_admissions(out):
    """Parse setauto output -> {admission: dict}, plus transition lists."""
    adm, trans = {}, []
    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "AUTO":
            adm[int(t[1])] = dict(
                info=int(t[3]), re=float(t[5]), re0=float(t[7]),
                eps=float(t[9]), rung=int(t[11]), alpha=float(t[13]),
                nc=int(t[15]), xt=float(t[17]))
        elif t[0] in ("RELEASE", "REVERT", "PERT"):
            trans.append(ln)
    return adm, trans


def _t_rows(spec):
    # T row is the second-to-last block row (row 35, 0-based)
    return np.arange(spec["nd"]) * spec["nt"] + spec["nt"] - 2


def t12_shock(spec):
    """STELLA-shock mock (FIX_T_AUTO=2 semantics): with adaptive
    damping the T march must stay inside [0, x*_T] (no overshoot,
    never negative), the energy residual RE must fall, and T must
    reach equilibrium; the undamped control jumps the whole deficit
    in one admission -- the excursion the guard exists to bound."""
    tr = _t_rows(spec)
    xs_t = float(np.abs(spec["xstar"][tr]).max())
    if xs_t == 0.0:
        return True, "skipped: no T shock (x*_T=0)", {}
    notes = []
    ok = True

    def run(seq, a0):
        out = run_cli("setauto", spec["dir"], seq, a0, 0.3, 0.1, 0.5, 2)
        return _auto_admissions(out)[0]

    damped = run("D" * 44, 0.05)
    free = run("D", 1.0)

    # The guarded Newton path rejects many early shocks, so the bounded
    # march is longer than the old least-squares-only fixture allowed.
    xts = np.array([damped[i]["xt"] for i in sorted(damped)])
    a = (len(damped) >= 40 and np.all(xts <= xs_t + 0.05)
         and np.all(xts >= -1e-12))
    ok &= a
    notes.append(f"x_T bounded in [0,x*] every admission: {a}")

    # Near the floor the rejected iterates can add noise to RE.  Keep
    # the early march monotone and require a large net reduction by the
    # end rather than allowing local floor noise to hide a runaway.
    res = np.array([damped[i]["re"] for i in sorted(damped)])
    a = (np.all(np.diff(res[:8]) <= 0.05 * res[:7])
         and np.all(res[8:] <= res[0])
         and res[-1] <= 0.01 * res[0])
    ok &= a
    notes.append(f"RE decreases and reaches the floor: {a} "
                 f"({res[0]:.3e} -> {res[-1]:.3e})")

    # T reaches the equilibrium correction
    a = xts[-1] >= 0.97 * xs_t
    ok &= a
    notes.append(f"x_T reached x*={xs_t:.3f}: {a} (xt={xts[-1]:.3f})")

    # the undamped control takes the full deficit in one admission.
    # Adaptive relax (5c4aefe3) raises the first damped step, so the
    # separation is now ~3.5x rather than the pre-adaptive >=4x.
    eps0_d = damped[min(damped)]["eps"]
    eps0_f = free[min(free)]["eps"]
    a = eps0_f >= 3.0 * max(eps0_d, 1e-12)
    ok &= a
    notes.append(f"undamped first step {eps0_f:.3f} >> damped "
                 f"{eps0_d:.3f}: {a}")
    return ok, "; ".join(notes), {}


def t13_tauto_fsm(spec):
    """FIX_T_AUTO ladder FSM: rung 1 holds while committed |dx_T|inf
    exceeds eps; RELEASE fires once RE < rel*RE0 with the commit
    streak met; an outer-iteration T re-shock reverts a rung and
    clears the streak (the production JFNK_SOLVE_CORR_V1 replica)."""
    tr = _t_rows(spec)
    xs_t = float(np.abs(spec["xstar"][tr]).max())
    if xs_t == 0.0:
        return True, "skipped: no T shock (x*_T=0)", {}
    notes = []
    ok = True
    out = run_cli("setauto", spec["dir"], "D" * 44 + "P" + "D" * 4,
                  0.05, 0.3, 0.1, 0.5, 2)
    adm, trans = _auto_admissions(out)

    # The guarded solve rejects several shock iterates before it can
    # commit, so locate the ladder transitions instead of pinning them
    # to a fixed admission index.
    rel = [ln for ln in trans if ln.startswith("RELEASE")]
    rev = [ln for ln in trans if ln.startswith("REVERT")]
    per = [ln for ln in trans if ln.startswith("PERT")]
    rel_i = int(rel[0].split()[1]) if rel else 0
    per_i = int(per[-1].split()[1]) if per else 0
    post_i = per_i + 1

    # rung 1 holds until the first release.
    a = bool(rel) and all(
        adm[i]["rung"] == 1 for i in range(1, rel_i)
    )
    ok &= a
    notes.append(f"rung 1 locked until release: {a}")

    # release to rung 2 fires once RE < rel*RE0 and eps settles.
    a = bool(rel) and adm[rel_i]["rung"] == 2 and \
        adm[rel_i]["alpha"] == 0.2
    ok &= a
    notes.append(f"release 1->2 on RE drop: {a} ({rel[0] if rel else '-'})")

    # The outer-iteration re-shock reverts the highest rung and clears
    # the streak on the next admission.
    a = (bool(rev) and bool(per) and post_i in adm
         and adm[post_i]["rung"] == adm[per_i - 1]["rung"] - 1
         and adm[post_i]["nc"] == 0)
    ok &= a
    notes.append(f"re-shock reverts one rung: {a} "
                 f"({rev[-1] if rev else '-'})")

    # no rung-1 admission may have all three release conditions met
    # (commit streak + RE below rel*RE0 + eps below the threshold)
    viol = [i for i in range(1, rel_i)
            if adm[i]["rung"] == 1 and adm[i]["nc"] >= 2
            and adm[i]["re"] <= 0.3 * adm[i]["re0"]
            and adm[i]["eps"] <= 0.1]
    a = rel_i > 0 and len(viol) == 0
    ok &= a
    notes.append(f"no release without all conditions: {a}")
    return ok, "; ".join(notes), {}


def t14_jvp_tauto(spec):
    """JVP independence under damping: with tauto on, a pure-T
    direction must produce nonzero, FD-vs-exact-consistent T rows in
    J*v (the T columns are never zeroed -- only step acceptance is
    damped); with tauto off the same rows are masked to zero."""
    tr = _t_rows(spec)
    notes = []
    ok = True
    n, nt = spec["n"], spec["nt"]

    v = np.zeros(spec["n"])
    v[tr] = 1.0
    vf = TD / "systems" / "_vT.txt"
    _write_vec(vf, v)

    # damped-but-live mode: T column of J is present and correct
    out = run_cli("sejvp", spec["dir"], vf, 1e-6, 0, 1)
    jv = parse_vec(out, "JV")
    jvfd = np.array([jv[i][0] for i in range(spec["n"])])
    jvex = np.array([jv[i][1] for i in range(spec["n"])])
    kv = parse_keyvals(out)
    a = (np.any(jvfd[tr] != 0.0) and np.any(jvex[tr] != 0.0)
         and kv.get("RELDIFF", 1.0) < 1e-3)
    ok &= a
    notes.append(f"tauto=1: T columns live, reldiff={kv.get('RELDIFF'):.2e}: {a}")

    # hard-mask mode unchanged: same direction gives zero T rows
    out = run_cli("sejvp", spec["dir"], vf, 1e-6, 0)
    jv = parse_vec(out, "JV")
    jvfd0 = np.array([jv[i][0] for i in range(spec["n"])])
    jvex0 = np.array([jv[i][1] for i in range(spec["n"])])
    a = np.all(jvfd0[tr] == 0.0) and np.all(jvex0[tr] == 0.0)
    ok &= a
    notes.append(f"tauto=0: T rows masked zero: {a}")
    return ok, "; ".join(notes), {}


def t15_stale_baseline(spec):
    """Stale subtraction baseline: sejvpst shifts jfnk_fx by delta on
    one row, reproducing the production F_outer-vs-F_eval(0) bug.
    The artifact must equal delta/h, land only on the offset row, and
    dominate the stale JVP error -- the pre-fix colprobe signature."""
    n = spec["n"]
    v = np.zeros(n)
    v[5] = 1.0
    vf = TD / "systems" / "_vstale.txt"
    _write_vec(vf, v)
    eps, delta, krow = 1e-6, 1e-2, 10
    out = run_cli("sejvpst", spec["dir"], vf, eps, 0, delta, krow)
    kv = parse_keyvals(out)
    jv = parse_vec(out, "JV")
    h = kv["H"]
    exp_art = delta / h
    art = np.array([jv[i][0] - jv[i][1] for i in range(n)])
    hot = np.where(np.abs(art) > 1e-3 * exp_art)[0]
    ok_mag = abs(kv["STALEDIFF"] - exp_art) < 2e-2 * exp_art
    ok_loc = list(hot) == [krow - 1]
    ok_dom = kv["STALEEXD"] > 0.9 * exp_art
    ok_sane = kv["CLEANEXD"] < 1e-3 * kv["JVEX"]
    ok = ok_mag and ok_loc and ok_dom and ok_sane
    msg = (f"artifact {kv['STALEDIFF']:.2e} vs d/h {exp_art:.2e}; "
           f"hot rows {[i + 1 for i in hot]}; "
           f"stale-exact {kv['STALEEXD']:.1e}, "
           f"clean-exact {kv['CLEANEXD']:.1e}")
    return ok, msg, {}


def t16_structure_freeze(spec):
    """REPLACE-freeze analogue on systems/flip: a row whose equation
    branch flips when x[icol] crosses a threshold.  A live branch
    corrupts the FD column by O(gap/h); freezing the baseline branch
    (USE_PASSED_REP=T) restores FD-vs-exact agreement."""
    specf = load_spec("flip")
    eps = 1e-6
    live = run_cli("seflip", specf["dir"], 0, 0, eps)
    frz = run_cli("seflip", specf["dir"], 0, 1, eps)
    kvl, kvf = parse_keyvals(live), parse_keyvals(frz)
    rows = parse_vec(live, "JV")
    n = specf["n"]
    diff = np.array([abs(rows[i][0] - rows[i][1]) for i in range(n)])
    irow = int(kvl["FLIPROW"])
    hot = np.where(diff > 1e-3 * diff.max())[0]
    ok_live = kvl["RELDIFF"] > 1.0
    ok_frz = kvf["RELDIFF"] < 1e-3
    ok_loc = list(hot) == [irow - 1]
    ok = ok_live and ok_frz and ok_loc
    msg = (f"live reldiff {kvl['RELDIFF']:.1e} (row {irow}), "
           f"frozen {kvf['RELDIFF']:.1e}; "
           f"hot rows {[i + 1 for i in hot]}")
    return ok, msg, {}


def t17_repeat_determinism(spec):
    """Pristine-vs-repeat identity: the same JVP probe twice must be
    bit-identical (no eval-state contamination between calls)."""
    n = spec["n"]
    v = np.zeros(n)
    v[5] = 1.0
    v[7] = -0.5
    vf = TD / "systems" / "_vrep.txt"
    _write_vec(vf, v)
    a = run_cli("sejvp", spec["dir"], vf, 1e-6, 0)
    b = run_cli("sejvp", spec["dir"], vf, 1e-6, 0)
    ok = (a == b)
    return ok, "byte-identical" if ok else "output drifted", {}


def t18_lagged_factor(spec):
    """A synthetic density factor adds the low-rank derivative omitted
    by the stored BA.  Pinning the factor at baseline must make its FD
    JVP agree with the BA, while the live-factor Jacobian differs."""
    out = run_cli("sefactor", spec["dir"], 1e-6)
    kv = parse_keyvals(out)
    ok_live = kv["LIVE_RELFD"] < 1e-4
    ok_lag = kv["LAG_RELFD"] < 1e-4
    ok_gap = kv["J_GAP"] > 0.1
    ok = ok_live and ok_lag and ok_gap
    msg = (f"live FD err {kv['LIVE_RELFD']:.2e}; "
           f"lagged FD-vs-BA {kv['LAG_RELFD']:.2e}; "
           f"live/BA gap {kv['J_GAP']:.2e}")
    return ok, msg, {}


def t19_fallback_dispatch():
    """Guard the dispatch invariant that finite-stagnation fallback
    starts from a defined caller cap and preserves the FULL/Lambda flags
    unless the eval is catastrophic."""
    source = (TD.parent / "new_main/mod_subs/jfnk_solve_corr_v1.f").read_text()
    init_at = source.index("T1=COR_LIM_IN")
    failure_at = source.index("! info=-2 is the reachability precheck skip")
    failure_end = source.index("FIX_T_AUTO ladder", failure_at)
    failure_block = source[failure_at:failure_end]
    ok = (
        init_at < failure_at
        and "JFNK_BUDGET_RUNG=MIN(JFNK_BUDGET_RUNG+1" in failure_block
        and "LAMBDA_A=LAMBDA_SAVE" in failure_block
        and "COMPUTE_BA_A=COMPUTE_BA_SAVE" in failure_block
        and "T1=JFNK_FALLBACK_CAP" in failure_block
    )
    return ok, "caller cap/mode initialized and restored", {}


def t20_eval_count_source():
    """The OUTGEN counters must distinguish direct Newton F calls from
    all residual callbacks, including those issued inside FD JVPs."""
    root = TD.parent / "new_main/mod_subs"
    files = (root / "jfnk_solve_corr_v1.f", root / "jfnk_cmfgen_interface_v1.f")
    out = [path.read_text() for path in files]
    ok = all(
        "Newton-level F evaluations:" in text
        and "Residual evals incl JVPs:" in text
        and "JFNK_N_RESIDUAL_EVALS" in text
        for text in out
    )
    return ok, "both JFNK reports expose callback and Newton F counts", {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--spec", default=None)
    args = ap.parse_args()
    if not BIN.exists():
        sys.exit(f"build first: bash build.sh")
    names = [args.spec] if args.spec else ["mild", "std", "hard",
                                          "shock"]
    nfail = 0
    for name in names:
        spec = load_spec(name)
        print(f"\n================= spec '{name}' "
              f"(n={spec['n']}, hard={spec['meta']['hard']}) "
              f"=================")
        for tag, fn in (("T0", t0_identity),
                        ("T1", None), ("T2", t2_jvp), ("T3", None),
                        ("T4", t4_precond), ("T5", t5_gmres),
                        ("T6", None), ("T7", None),
                        ("T8", t8_masked), ("T9", t9_jvp_stiff),
                        ("T10", t10_lsq_grid),
                        ("T11", t11_strike if spec["meta"]["hard"]
                         else None),
                        ("T12", t12_shock), ("T13", t13_tauto_fsm),
                        ("T14", t14_jvp_tauto),
                        ("T15", t15_stale_baseline),
                        ("T16", t16_structure_freeze),
                        ("T17", t17_repeat_determinism)):
            label = {"T0": "identity", "T1": "anatomy", "T2": "FD JVP",
                     "T3": "demanded step", "T4": "preconditioner",
                     "T5": "GMRES", "T6": "Newton end-to-end",
                     "T7": "residual floor", "T8": "masked-T",
                     "T9": "JVP stiff", "T10": "LSQ grid",
                     "T11": "budget ladder", "T12": "STELLA shock",
                     "T13": "tauto ladder", "T14": "JVP tauto",
                     "T15": "stale baseline", "T16": "struct freeze",
                     "T17": "repeat determinism"}[tag]
            print(f"  [{tag} {label}]")
            try:
                if tag == "T0":
                    ok, msg, aux = fn(spec)
                    F0 = aux["f"]
                elif tag == "T1":
                    ok, msg, aux = t1_anatomy(spec, F0)
                elif tag == "T3":
                    ok, msg, aux = t3_demand(spec, F0)
                elif tag == "T6":
                    ok, msg, aux = t6_newton(spec, spec["xstar"])
                elif tag == "T7":
                    ok, msg, aux = t7_floor(spec, F0)
                elif tag == "T11" and not spec["meta"]["hard"]:
                    ok, msg = True, "skipped: needs 'hard' (LU fails, " \
                        "TSVD commits)"
                else:
                    ok, msg, aux = fn(spec)
            except subprocess.CalledProcessError as e:
                print("    driver crashed:", e.stdout[-500:],
                      e.stderr[-500:])
                ok, msg = False, "driver crashed"
            print(f"    {'PASS' if ok else 'FAIL'}  ({msg})")
            if not ok:
                nfail += 1
            if args.quick and tag in ("T4", "T5"):
                pass
    if args.spec is None or args.spec == "lag":
        spec = load_spec("lag")
        print("\n================= spec 'lag' (factor derivative) "
              "=================")
        try:
            ok, msg, _ = t18_lagged_factor(spec)
        except subprocess.CalledProcessError as e:
            print("    driver crashed:", e.stdout[-500:], e.stderr[-500:])
            ok, msg = False, "driver crashed"
        print(f"  [T18 lagged factor] {'PASS' if ok else 'FAIL'}  ({msg})")
        if not ok:
            nfail += 1
    if args.spec is None or args.spec == "lag":
        ok, msg, _ = t19_fallback_dispatch()
        print(f"  [T19 fallback dispatch] {'PASS' if ok else 'FAIL'}  ({msg})")
        if not ok:
            nfail += 1
    ok, msg, _ = t20_eval_count_source()
    print(f"  [T20 eval telemetry] {'PASS' if ok else 'FAIL'}  ({msg})")
    if not ok:
        nfail += 1
    print(f"\n{nfail} test(s) failed")
    return nfail


if __name__ == "__main__":
    sys.exit(main())
