#!/usr/bin/env python3
"""Stage-E1 manufactured wavelength-depth block-system battery.

drv_mom_blk writes manufactured row coefficients.  This battery
assembles the dense coupled operator independently, verifies a small
left-preconditioned restarted GMRES prototype against the dense
reference, and measures the lagged frequency iteration on the same
graphs.  No production solver code is exercised at this stage.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, TREE, WORKDIR, Report, compile_driver, run

ND = 5
NCF = 12
CASES = ["exp", "cmp", "ring2", "ring4", "marg", "mix", "depth", "fresh"]


def build():
    return compile_driver(
        "drv_mom_blk",
        [TREE / "tools" / "set_kind_module.f", F90DIR / "drv_mom_blk.f"],
    )


def load_system(tag):
    rows = []
    meta = (None, None)
    with open(WORKDIR / f"{tag}_sys.txt") as f:
        for ln in f:
            if ln.startswith("# ND"):
                p = ln.replace("#", " ").replace(":", " ").split()
                meta = (int(p[-2]), int(p[-1]))
            elif not ln.startswith("#"):
                rows.append(ln.split())
    nd, ncf = meta
    assert (nd, ncf) == (ND, NCF)
    ta = np.zeros((nd, ncf))
    tb = np.zeros_like(ta)
    tc = np.zeros_like(ta)
    x0 = np.zeros_like(ta)
    cup = np.zeros_like(ta)
    xold = np.zeros_like(ta)
    kup = np.zeros((nd, ncf), dtype=int)
    fresh = np.zeros(ncf, dtype=bool)
    for r in rows:
        i = int(r[0]) - 1
        k = int(r[1]) - 1
        ta[i, k], tb[i, k], tc[i, k] = map(float, r[2:5])
        x0[i, k], cup[i, k] = map(float, r[5:7])
        kup[i, k] = int(r[7])
        fresh[k] = r[8].strip() in ("T", ".TRUE.", "True", "1")
        xold[i, k] = float(r[9])
    return ta, tb, tc, x0, cup, kup, fresh, xold


def dense_system(sysdat):
    """Assemble A x=b.  CUP is the positive RHS coefficient, so the
    coupled LHS carries -CUP at the upwind target.  Inactive columns
    use identity rows pinning x to the incoming field."""
    ta, tb, tc, x0, cup, kup, fresh, xold = sysdat
    nd, ncf = tb.shape
    n = nd * ncf
    a = np.zeros((n, n))
    b = np.zeros(n)

    def ix(i, k):
        return k * nd + i

    for k in range(ncf):
        for i in range(nd):
            r = ix(i, k)
            if not fresh[k]:
                a[r, r] = 1.0
                b[r] = xold[i, k]
                continue
            a[r, r] = tb[i, k]
            if i > 0:
                a[r, ix(i - 1, k)] = ta[i, k]
            if i < nd - 1:
                a[r, ix(i + 1, k)] = tc[i, k]
            ku = kup[i, k]
            if ku > 0:
                a[r, ix(i, ku - 1)] -= cup[i, k]
            b[r] = x0[i, k]
    return a, b


def block_preconditioner_inverse(sysdat):
    """Dense inverse of the frequency block-diagonal depth system."""
    ta, tb, tc, x0, cup, kup, fresh, xold = sysdat
    nd, ncf = tb.shape
    n = nd * ncf
    bd = np.zeros((n, n))
    for k in range(ncf):
        for i in range(nd):
            r = k * nd + i
            if not fresh[k]:
                bd[r, r] = 1.0
                continue
            bd[r, r] = tb[i, k]
            if i > 0:
                bd[r, r - 1] = ta[i, k]
            if i < nd - 1:
                bd[r, r + 1] = tc[i, k]
    inv = np.zeros_like(bd)
    for k in range(ncf):
        sl = slice(k * nd, (k + 1) * nd)
        inv[sl, sl] = np.linalg.inv(bd[sl, sl])
    return inv, bd


def gmres_restart(a, b, minv, x0=None, restart=40, tol=1e-14,
                  maxrestart=6):
    """Small dense restarted GMRES prototype, left preconditioned.

    Full orthogonalization with a second Gram-Schmidt pass.  The
    least-squares state is solved directly each iteration; these toy
    systems are at most 60 unknowns, so the cost is immaterial.
    """
    n = len(b)
    x = np.zeros(n) if x0 is None else x0.copy()

    def op(v):
        return minv @ (a @ v)

    bnrm = np.linalg.norm(minv @ b)
    if bnrm == 0.0:
        bnrm = 1.0
    total = 0
    for rst in range(maxrestart):
        z0 = minv @ (b - a @ x)
        beta = np.linalg.norm(z0)
        if beta <= tol * bnrm:
            return x, total, rst, beta / bnrm, True
        vlist = [z0 / beta]
        h = np.zeros((restart + 1, restart))
        y = None
        used = 0
        for j in range(restart):
            w = op(vlist[j])
            # Two full orthogonalization passes keep the small Arnoldi
            # basis orthogonal close to machine precision.
            hj = np.zeros(j + 1)
            for _ in range(2):
                for ii, q in enumerate(vlist):
                    corr = np.dot(q, w)
                    hj[ii] += corr
                    w -= corr * q
            h[: j + 1, j] = hj
            hn = np.linalg.norm(w)
            h[j + 1, j] = hn
            if hn > 0.0:
                vlist.append(w / hn)
            jj = j
            e = np.zeros(jj + 2)
            e[0] = beta
            hh = h[: jj + 2, : jj + 1]
            y, *_ = np.linalg.lstsq(hh, e, rcond=None)
            res = np.linalg.norm(e - hh @ y)
            total += 1
            used = jj + 1
            if res <= tol * bnrm or hn == 0.0:
                break
        vv = np.column_stack(vlist[:used])
        x = x + vv @ y[:used]
        z0 = minv @ (b - a @ x)
        beta = np.linalg.norm(z0)
        if beta <= tol * bnrm:
            return x, total, rst, beta / bnrm, True
    z0 = minv @ (b - a @ x)
    return x, total, maxrestart, np.linalg.norm(z0) / bnrm, False


def lagged_iteration_matrix(a, nd, ncf):
    """Block-Jacobi-in-frequency iteration matrix of the lagged map."""
    n = a.shape[0]
    bd = np.zeros_like(a)
    for k in range(ncf):
        sl = slice(k * nd, (k + 1) * nd)
        bd[sl, sl] = a[sl, sl]
    return np.linalg.solve(bd, bd - a)


def main():
    rep = Report("test_mom_blk")
    exe = build()
    prodexe = compile_driver(
        "drv_mom_blk_solve",
        [
            TREE / "tools" / "set_kind_module.f",
            TREE / "plane" / "mom_blk_gmres_mod.f",
            F90DIR / "drv_mom_blk_solve.f",
        ],
    )

    for mode in CASES:
        tag = f"blk_{mode}"
        r = run(exe, mode, tag, cwd=WORKDIR)
        rep.check(f"{mode}: driver wrote the manufactured system",
                  "RESULT_OK" in r.stdout)
        sysdat = load_system(tag)
        a, b = dense_system(sysdat)
        xd = np.linalg.solve(a, b)
        rep.check(f"{mode}: dense reference finite",
                  np.all(np.isfinite(xd)))
        minv, bd = block_preconditioner_inverse(sysdat)
        xg, nit, nrst, pres, conv = gmres_restart(a, b, minv)
        scale = max(np.max(np.abs(xd)), 1.0)
        err = np.max(np.abs(xg - xd)) / scale
        resid = np.max(np.abs(b - a @ xg)) / max(np.max(np.abs(b)), 1.0)
        rep.check(
            f"{mode}: preconditioned GMRES matches dense reference",
            conv and err < 2.0e-12,
            f"err {err:.2e}, it {nit}, restart {nrst}, pres {pres:.2e}",
        )
        rep.check(
            f"{mode}: final dense residual at machine level",
            resid < 2.0e-12,
            f"max relative residual {resid:.2e}",
        )
        pr = run(prodexe, tag, cwd=WORKDIR)
        rep.check(
            f"{mode}: production GMRES module converged",
            "RESULT_OK" in pr.stdout,
            pr.stdout.strip().splitlines()[-2:],
        )
        xprod = np.zeros(ND * NCF)
        with open(WORKDIR / f"{tag}_xgm.txt") as f:
            for ln in f:
                q = ln.split()
                i = int(q[0]) - 1
                k = int(q[1]) - 1
                xprod[k * ND + i] = float(q[2])
        perr = np.max(np.abs(xprod - xd)) / scale
        prtrue = np.max(np.abs(b - a @ xprod)) / max(np.max(np.abs(b)), 1.0)
        rep.check(
            f"{mode}: production GMRES matches dense reference",
            perr < 5.0e-12,
            f"err {perr:.2e}, dense-resid {prtrue:.2e}",
        )
        if mode != "fresh":
            rep.check(
                f"{mode}: inactive column logic unused as expected",
                np.all(sysdat[6]),
            )
        else:
            frozen = ~sysdat[6]
            nd = sysdat[1].shape[0]
            xf = xg.reshape(nd, len(sysdat[6]), order="F")
            xo = sysdat[7]
            rep.check(
                "fresh: inactive columns remain pinned to old field",
                np.allclose(xf[:, frozen], xo[:, frozen], rtol=0.0,
                            atol=1e-14),
            )

    # Marginal-ring regression anchor: the lagged frequency iteration
    # is the block-Jacobi frequency map.  It is non-contractive on the
    # marginal graph, while the coupled preconditioned GMRES solve
    # still reaches the dense answer.
    sysdat = load_system("blk_marg")
    a, b = dense_system(sysdat)
    miter = lagged_iteration_matrix(a, ND, NCF)
    rho = max(abs(np.linalg.eigvals(miter)))
    minv, _ = block_preconditioner_inverse(sysdat)
    xg, nit, nrst, pres, conv = gmres_restart(a, b, minv)
    xd = np.linalg.solve(a, b)
    err = np.max(np.abs(xg - xd)) / max(np.max(np.abs(xd)), 1.0)
    rep.check(
        "marginal ring: lagged frequency map is non-contractive",
        rho >= 1.0,
        f"rho {rho:.4f}",
    )
    rep.check(
        "marginal ring: coupled GMRES error contraction is exact",
        conv and err < 2.0e-12,
        f"err {err:.2e}, it {nit}, rho {rho:.4f}",
    )

    # Iterated-probe analogue. The 1000/iter amplification measured by
    # drv_mom_nu's iter mode is outer source-feedback amplification
    # through a nearly singular ring; that feedback remains outer-
    # iterated by design. The solver-level question is different: from
    # the same zero initial error state, the lagged frequency iteration
    # amplifies the error while the coupled solve removes it. Pin that
    # contraction directly rather than claiming a linear solve can
    # remove nonlinear H/EPS feedback it does not contain.
    xlag = np.zeros_like(xd)
    lag_err = []
    off = a - bd
    for _ in range(8):
        xlag = np.linalg.solve(bd, b - off @ xlag)
        lag_err.append(
            np.linalg.norm(xlag - xd) / max(np.linalg.norm(xd), 1.0e-30)
        )
    lag_gain = lag_err[-1] / lag_err[-2]
    xblk, bit, brst, bres, bconv = gmres_restart(
        a, b, minv, x0=np.zeros_like(xd)
    )
    blk_gain = (
        np.linalg.norm(xblk - xd) / max(np.linalg.norm(xd), 1.0e-30)
    )
    rep.check(
        "iterated probe analogue: coupled solve contracts solver error",
        lag_gain > 1.0 and bconv and blk_gain < 1.0e-10,
        f"lagged error growth {lag_gain:.3f}/iter, "
        f"block error ratio {blk_gain:.2e} after {bit} iterations",
    )

    rep.summary()


if __name__ == "__main__":
    main()
