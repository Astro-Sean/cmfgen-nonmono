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
    """Strike/cooldown FSM (sestrike replica): 3 strikes -> pause,
    adaptive wait doubling, commit resets, reach-skip exempt.
    Requires a spec where the LU admission fails (info=3) and the
    LSQ-alg admission commits (info=0) -- only 'hard' does both."""
    ok = True
    notes = []

    def run(seq, rearm=2):
        out = run_cli("sestrike", spec["dir"], rearm, seq)
        adm, pauses, rearms = {}, [], []
        for ln in out.splitlines():
            t = ln.split()
            if not t:
                continue
            if t[0] == "ADM":
                adm[int(t[1])] = dict(info=int(t[3]), stag=int(t[5]),
                                      dis=t[7] == "T", cyc=int(t[9]))
            elif t[0] == "PAUSE":
                pauses.append((int(t[1]), int(t[3]), int(t[5])))
            elif t[0] == "REARM":
                rearms.append((int(t[1]), int(t[3])))
        return adm, pauses, rearms

    # strikes accumulate, third disables with cycle=1
    adm, pauses, rearms = run("FFFFFFFFF")
    a = (adm[1]["stag"] == 1 and adm[2]["stag"] == 2
         and adm[3]["dis"] and adm[3]["cyc"] == 1)
    # wait = rearm*2^cyc = 4 -> 3 pauses then rearm at admission 7
    a &= ([p[0] for p in pauses] == [4, 5, 6]
          and pauses[-1][1:] == (3, 4) and rearms == [(7, 4)])
    # second disable doubles the cycle counter
    a &= adm[9]["dis"] and adm[9]["cyc"] == 2
    ok &= a
    notes.append(f"3-strike pause/rearm/doubling: {a}")

    # committed update resets strikes and cycle
    adm, pauses, rearms = run("FFSFFFF")
    a = (adm[3]["info"] == 0 and adm[3]["stag"] == 0
         and adm[3]["cyc"] == 0 and adm[4]["stag"] == 1
         and adm[6]["dis"] and adm[6]["cyc"] == 1
         and pauses == [(7, 1, 4)])
    ok &= a
    notes.append(f"commit resets stag+cycle: {a}")

    # reach-gate skip (info=-2) is exempt from strikes
    adm, pauses, rearms = run("RFFRFF")
    a = (adm[1]["info"] == -2 and adm[1]["stag"] == 0
         and adm[4]["stag"] == 2 and adm[5]["dis"])
    ok &= a
    notes.append(f"reach-skip exempt: {a}")
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

    damped = run("DDDDDDDD", 0.05)
    free = run("D", 1.0)

    # survive the first 5+ admissions with bounded T motion and no
    # overshoot of the equilibrium correction (surrogate T=e^x so
    # positivity is inherent; the excursion bound is the check)
    xts = np.array([damped[i]["xt"] for i in sorted(damped)])
    a = (len(damped) >= 5 and np.all(xts <= xs_t + 0.05)
         and np.all(xts >= -1e-12))
    ok &= a
    notes.append(f"x_T bounded in [0,x*] every admission: {a}")

    # energy residual falls monotonically (small wiggle allowed)
    res = np.array([damped[i]["re"] for i in sorted(damped)])
    a = np.all(np.diff(res) <= 0.05 * res[:-1])
    ok &= a
    notes.append(f"RE monotone decrease: {a} "
                 f"({res[0]:.3e} -> {res[-1]:.3e})")

    # T reaches the equilibrium correction
    a = xts[-1] >= 0.97 * xs_t
    ok &= a
    notes.append(f"x_T reached x*={xs_t:.3f}: {a} (xt={xts[-1]:.3f})")

    # the undamped control takes the full deficit in one admission --
    # >= 4x the first damped step (the STELLA shock excursion)
    eps0_d = damped[min(damped)]["eps"]
    eps0_f = free[min(free)]["eps"]
    a = eps0_f >= 4.0 * max(eps0_d, 1e-12)
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
    out = run_cli("setauto", spec["dir"], "DDDDDDDDPDDDD",
                  0.05, 0.3, 0.1, 0.5, 2)
    adm, trans = _auto_admissions(out)

    # rung 1 held throughout the march (eps > 0.1 blocks release)
    march = [adm[i] for i in range(1, 8)]
    a = all(a_["rung"] == 1 for a_ in march)
    ok &= a
    notes.append(f"rung 1 locked while eps>eps_thr: {a}")

    # release to rung 2 fires once RE < rel*RE0 and eps settles
    rel = [ln for ln in trans if ln.startswith("RELEASE")]
    a = len(rel) >= 1 and adm[8]["rung"] == 2 and \
        adm[8]["alpha"] == 0.2
    ok &= a
    notes.append(f"release 1->2 on RE drop: {a} ({rel[0] if rel else '-'})")

    # the re-shock admission reverts a rung and clears the streak
    rev = [ln for ln in trans if ln.startswith("REVERT")]
    a = len(rev) >= 1 and adm[10]["rung"] == 1 and \
        adm[10]["alpha"] == 0.05 and adm[10]["nc"] == 0
    ok &= a
    notes.append(f"re-shock reverts 2->1: {a} ({rev[0] if rev else '-'})")

    # no rung-1 admission may have all three release conditions met
    # (commit streak + RE below rel*RE0 + eps below the threshold)
    viol = [i for i in range(1, 8)
            if adm[i]["rung"] == 1 and adm[i]["nc"] >= 2
            and adm[i]["re"] <= 0.3 * adm[i]["re0"]
            and adm[i]["eps"] <= 0.1]
    a = len(viol) == 0
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
                        ("T14", t14_jvp_tauto)):
            label = {"T0": "identity", "T1": "anatomy", "T2": "FD JVP",
                     "T3": "demanded step", "T4": "preconditioner",
                     "T5": "GMRES", "T6": "Newton end-to-end",
                     "T7": "residual floor", "T8": "masked-T",
                     "T9": "JVP stiff", "T10": "LSQ grid",
                     "T11": "strike/cooldown", "T12": "STELLA shock",
                     "T13": "tauto ladder", "T14": "JVP tauto"}[tag]
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
                        "LSQ commits)"
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
    print(f"\n{nfail} test(s) failed")
    return nfail


if __name__ == "__main__":
    sys.exit(main())
