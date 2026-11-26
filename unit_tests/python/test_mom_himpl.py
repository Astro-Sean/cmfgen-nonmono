#!/usr/bin/env python3
"""Implicit-H contraction prototype (FIX 4, phase 2 stage 1).

Builds two linear operators in numpy from the drv_mom_augres captured
array dumps:

  * the production-emulated lagged map - solve the captured x block
    (tridiagonal + CUP upwind edges, optionally + |VB|+|VC| diagonal
    mass, i.e. MOM_H_IMPL) with the H/EPS lagged bracket read from the
    previous iterate, then reconstruct H exactly through its
    recurrence. The composite map z -> M z + c has spectral radius
    rho(M). The mass is uncompensated on the RHS (matching
    production), so the massed map keeps a shifted fixed point.
  * the augmented operator - the implemented linear skeleton with the
    lagged bracket and H reconstruction promoted to implicit rows on
    the interleaved unknown z(k) = [x(1),h(1),...,x(ND)]. No mass.

The "hcon" case is the epsmix capture with NU_VB/NU_VC scaled so the
lagged map diverges without mass - the toy analogue of the regime
where production currently reaches for MOM_H_IMPL. b = A_aug z_man
keeps the manufactured field the exact solution, so a converged solve
is checked against the dumped field, not just a residual.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import Report, WORKDIR, run
from test_mom_augres import ND, NCF, build, eval_aug, load

BOOST = 30.0
NX = ND * NCF
NH = (ND - 1) * NCF
NW = 2 * ND - 1
NZ = NCF * NW


def xidx(i, k):
    return k * ND + i


def zxi(i, k):
    return k * NW + 2 * i


def zhi(m, k):
    return k * NW + 2 * m + 1


def eps_donor(c, i, k, nmj, ncf):
    dd = i if c >= 0.0 else i + 1
    if c >= 0.0:
        u = k - 1
    else:
        u = k + nmj[dd] if k + nmj[dd] < ncf else -1
    return dd, u


def assemble(d):
    """Lagged split (a0, lag, bx) and augmented operator (a_aug, b_aug).

    a0: captured tridiagonal + CUP upwind edges, no mass. lag maps the
    block vector [x; h] to the lagged bracket (already LAGSCL-scaled).
    a_aug rows are the implicit coefficients minus the lagged ones
    (moved to the LHS) on the interleaved layout, plus H recurrence
    rows. bx carries NU_SRC; H rows are homogeneous.
    """
    nd, ncf = d["nd"], d["ncf"]
    nmj = d["nmj"]
    a0 = np.zeros((NX, NX))
    lag = np.zeros((NX, NZ))
    bx = np.zeros(NX)
    a_aug = np.zeros((NZ, NZ))
    b_aug = np.zeros(NZ)

    for k in range(ncf):
        outbc = int(d["outbc"][k])
        t1b = d["t1b"][k]
        for i in range(nd):
            r = xidx(i, k)
            zr = zxi(i, k)
            kup = int(d["kup"][i, k])
            bx[r] = d["src"][i, k]
            row = {}
            if i == 0:
                row[zxi(0, k)] = d["tb"][0, k]
                row[zxi(1, k)] = d["tc"][0, k]
                if outbc <= 1:
                    bx[r] = 0.0
                    if 1 <= kup <= ncf:
                        row[zxi(0, kup - 1)] = -d["cup"][0, k]
                else:
                    if 1 <= kup <= ncf:
                        row[zxi(0, kup - 1)] = -d["cup"][0, k]
                    scl = d["lagscl"][0, k]
                    u = int(d["kuh"][0, k])
                    if 1 <= u <= ncf:
                        lag[r, NX + (u - 1) * (nd - 1) + 0] += (
                            -scl * d["hs"][0, k] / t1b
                        )
                    c = d["epsp"][0, k]
                    if c != 0.0:
                        dd, u = eps_donor(c, 0, k, nmj, ncf)
                        if 0 <= u < ncf:
                            # Production subtracts LAG_SCL*T5 from
                            # the boundary XM (sweep :3404), so the
                            # lagged bracket carries -scl*2c/t1b.
                            lag[r, xidx(dd, u)] -= scl * 2.0 * c / t1b
            elif i == nd - 1:
                row[zxi(nd - 2, k)] = d["ta"][nd - 1, k]
                row[zxi(nd - 1, k)] = d["tb"][nd - 1, k]
            else:
                row[zxi(i - 1, k)] = d["ta"][i, k]
                row[zxi(i, k)] = d["tb"][i, k]
                row[zxi(i + 1, k)] = d["tc"][i, k]
                if 1 <= kup <= ncf:
                    row[zxi(i, kup - 1)] = -d["cup"][i, k]
                scl = d["lagscl"][i, k]
                if d["vb"][i, k] != 0.0:
                    u = int(d["kuh"][i - 1, k])
                    if 1 <= u <= ncf:
                        lag[r, NX + (u - 1) * (nd - 1) + i - 1] += (
                            scl * d["vb"][i, k]
                        )
                if d["vc"][i, k] != 0.0:
                    u = int(d["kuh"][i, k])
                    if 1 <= u <= ncf:
                        lag[r, NX + (u - 1) * (nd - 1) + i] += (
                            scl * d["vc"][i, k]
                        )
                for m, sgn in ((i, 1.0), (i - 1, -1.0)):
                    c = d["epsp"][m, k]
                    if c != 0.0:
                        dd, u = eps_donor(c, m, k, nmj, ncf)
                        if 0 <= u < ncf:
                            lag[r, xidx(dd, u)] += sgn * scl * 2.0 * c
            for c_i, v in row.items():
                a_aug[zr, c_i] += v
                a0[r, xidx((c_i % NW) // 2, c_i // NW)] += v
            b_aug[zr] = bx[r]
        for m in range(nd - 1):
            zr = zhi(m, k)
            a_aug[zr, zr] = 1.0
            a_aug[zr, zxi(m + 1, k)] -= d["hu"][m, k]
            a_aug[zr, zxi(m, k)] += d["hl"][m, k]
            u = int(d["kuh"][m, k])
            if 1 <= u <= ncf:
                a_aug[zr, zhi(m, u - 1)] -= d["hs"][m, k]
            c = d["epsp"][m, k]
            if c != 0.0:
                dd, u = eps_donor(c, m, k, nmj, ncf)
                if 0 <= u < ncf:
                    a_aug[zr, zxi(dd, u)] -= 2.0 * c
            e = d["eps"][m, k]
            if e != 0.0:
                dd = m if e >= 0.0 else m + 1
                a_aug[zr, zxi(dd, k)] += 2.0 * e

    # lag is stored on the block layout [x; h]; move its entries onto
    # the interleaved moment rows of a_aug.
    for r in range(NX):
        i, k = r % nd, r // nd
        zr = zxi(i, k)
        for c in range(NZ):
            v = lag[r, c]
            if v == 0.0:
                continue
            if c < NX:
                a_aug[zr, zxi(c % nd, c // nd)] -= v
            else:
                m = (c - NX) % (nd - 1)
                kk = (c - NX) // (nd - 1)
                a_aug[zr, zhi(m, kk)] -= v
    return a0, lag, bx, a_aug, b_aug


def h_recurrence(d):
    """Rx, Rh of h = Rx x + Rh h (implemented recurrence)."""
    nd, ncf = d["nd"], d["ncf"]
    nmj = d["nmj"]
    rx = np.zeros((NH, NX))
    rh = np.zeros((NH, NH))
    for k in range(ncf):
        for m in range(nd - 1):
            r = k * (nd - 1) + m
            rx[r, xidx(m + 1, k)] = d["hu"][m, k]
            rx[r, xidx(m, k)] = -d["hl"][m, k]
            u = int(d["kuh"][m, k])
            if 1 <= u <= ncf:
                rh[r, m + (u - 1) * (nd - 1)] = d["hs"][m, k]
            c = d["epsp"][m, k]
            if c != 0.0:
                dd, u = eps_donor(c, m, k, nmj, ncf)
                if 0 <= u < ncf:
                    rx[r, xidx(dd, u)] += 2.0 * c
            e = d["eps"][m, k]
            if e != 0.0:
                dd = m if e >= 0.0 else m + 1
                rx[r, xidx(dd, k)] -= 2.0 * e
    return rx, rh


def lagged_system(d, mass):
    """(a_impl, g, lx, lh): x_next = a_impl^-1 (b + lx x + lh h),
    h_next = g x_next with g = (I - Rh)^-1 Rx."""
    nd, ncf = d["nd"], d["ncf"]
    a0, lag, _, _, _ = assemble(d)
    rx, rh = h_recurrence(d)
    a_impl = a0.copy()
    if mass:
        for k in range(ncf):
            for i in range(1, nd - 1):
                a_impl[xidx(i, k), xidx(i, k)] += (
                    abs(d["vb"][i, k]) + abs(d["vc"][i, k])
                )
    g = np.linalg.solve(np.eye(NH) - rh, rx)
    return a_impl, g, lag[:, :NX], lag[:, NX:]


def lagged_map_rho(a_impl, g, lx, lh):
    ai = np.linalg.inv(a_impl)
    m = np.zeros((NZ, NZ))
    m[:NX, :NX] = ai @ lx
    m[:NX, NX:] = ai @ lh
    m[NX:, :NX] = g @ (ai @ lx)
    m[NX:, NX:] = g @ (ai @ lh)
    return max(abs(np.linalg.eigvals(m)))


def lagged_iterate(a_impl, g, lx, lh, bx, iters=400):
    ai = np.linalg.inv(a_impl)
    z = np.zeros(NX + NH)
    for _ in range(iters):
        x_new = ai @ (bx + lx @ z[:NX] + lh @ z[NX:])
        z_new = np.concatenate([x_new, g @ x_new])
        step = np.max(np.abs(z_new - z))
        if not np.isfinite(step) or np.max(np.abs(z_new)) > 1.0e30:
            return False, z_new
        if step < 1.0e-12 * max(1.0, np.max(np.abs(z_new))):
            return True, z_new
        z = z_new
    return False, z


def gmres(a, b, m_inv=None, tol=1.0e-11, maxit=400):
    """Restartless left-preconditioned GMRES; returns (x, iters, res)."""
    n = len(b)
    x = np.zeros(n)
    q = []

    def apply(v):
        return v if m_inv is None else m_inv @ v

    r = apply(b - a @ x)
    beta = np.linalg.norm(r)
    if beta < tol:
        return x, 0, beta
    q.append(r / beta)
    h = np.zeros((maxit + 1, maxit))
    for k in range(maxit):
        v = apply(a @ q[k])
        for j in range(k + 1):
            h[j, k] = np.dot(q[j], v)
            v = v - h[j, k] * q[j]
        h[k + 1, k] = np.linalg.norm(v)
        if h[k + 1, k] > 1.0e-14:
            q.append(v / h[k + 1, k])
        e = np.zeros(k + 2)
        e[0] = beta
        hh = h[: k + 2, : k + 1]
        y, *_ = np.linalg.lstsq(hh, e, rcond=None)
        res = np.linalg.norm(e - hh @ y)
        if res < tol:
            sol = np.column_stack(q[: k + 1]) @ y
            return x + sol, k + 1, res
    return x, maxit, res


def block_precond(d, a_aug):
    """Per-frequency dense block inverse (all z rows of one column)."""
    m = np.zeros_like(a_aug)
    for k in range(d["ncf"]):
        sl = slice(k * NW, (k + 1) * NW)
        m[sl, sl] = np.linalg.inv(a_aug[sl, sl])
    return m


def zpack(d):
    z = np.zeros(NZ)
    for k in range(NCF):
        for i in range(ND):
            z[zxi(i, k)] = d["jnext"][i, k]
        for m in range(ND - 1):
            z[zhi(m, k)] = d["hnext"][m, k]
    return z


def main():
    rep = Report("test_mom_himpl")
    exe = build()
    for case in ("exact", "epsmix"):
        r = run(exe, case, f"aug_{case}", cwd=WORKDIR)
        rep.check(f"{case}: driver ran", "RESULT_OK" in r.stdout)

    dex = load("aug_exact")
    deps = load("aug_epsmix")

    # Assembly validation on the unmodified exact capture: the
    # augmented operator annihilates the manufactured field.
    r_ex = eval_aug(dex)
    rep.check(
        "exact: evaluator clean (no clips/caps)",
        r_ex["nclip"] == 0 and r_ex["nepsl"] == 0 and r_ex["nnan"] == 0,
        f"nclip {r_ex['nclip']} nepsl {r_ex['nepsl']} nnan {r_ex['nnan']}",
    )
    _, _, _, a_aug, b_aug = assemble(dex)
    zman = zpack(dex)
    res = np.linalg.norm(a_aug @ zman - b_aug) / max(
        np.linalg.norm(b_aug), 1.0e-30
    )
    rep.check("exact: A_aug z_man = b", res < 1.0e-10, f"res {res:.3e}")

    # hcon: boosted H coupling; b re-synthesised so z_man stays exact.
    d = dict(deps)
    d["vb"] = deps["vb"] * BOOST
    d["vc"] = deps["vc"] * BOOST
    _, _, _, a_aug, b_aug = assemble(d)
    zman = zpack(d)
    b_aug = a_aug @ zman
    bx = np.zeros(NX)
    for k in range(NCF):
        for i in range(ND):
            bx[xidx(i, k)] = b_aug[zxi(i, k)]

    a_impl, g, lx, lh = lagged_system(d, False)
    rho_nm = lagged_map_rho(a_impl, g, lx, lh)
    ok_nm, _ = lagged_iterate(a_impl, g, lx, lh, bx)
    a_impl_m, g, lx, lh = lagged_system(d, True)
    rho_m = lagged_map_rho(a_impl_m, g, lx, lh)
    ok_m, z_m = lagged_iterate(a_impl_m, g, lx, lh, bx)
    rep.check(
        "hcon: lagged diverges w/o mass",
        rho_nm > 1.02 and not ok_nm,
        f"rho {rho_nm:.4f} converged {ok_nm}",
    )
    rep.check(
        "hcon: mass contracts lagged",
        rho_m < 0.98 and ok_m,
        f"rho {rho_m:.4f} converged {ok_m}",
    )

    z_star = np.linalg.solve(a_aug, b_aug)
    res_star = np.linalg.norm(a_aug @ z_star - b_aug) / np.linalg.norm(b_aug)
    rep.check("hcon: augmented solve exact", res_star < 1.0e-10,
              f"res {res_star:.3e}")
    dz = np.max(np.abs(z_star - zman)) / max(1.0, np.max(np.abs(zman)))
    rep.check("hcon: augmented recovers z_man", dz < 1.0e-8,
              f"max rel dz {dz:.3e}")
    if ok_m:
        x_star = np.array(
            [z_star[zxi(i, k)] for k in range(NCF) for i in range(ND)]
        )
        err = np.linalg.norm(z_m[:NX] - x_star) / max(
            np.linalg.norm(x_star), 1.0e-30
        )
        rep.check(
            "hcon: massed fixed point shifted (model error)",
            err > 1.0e-4,
            f"rel err {err:.3e}",
        )

    sol0, it0, res0 = gmres(a_aug, b_aug, maxit=NZ)
    minv = block_precond(d, a_aug)
    sol1, it1, res1 = gmres(a_aug, b_aug, m_inv=minv, maxit=NZ)
    rep.check("hcon: GMRES on A_aug unpreconditioned",
              res0 < 1.0e-9, f"it {it0} res {res0:.3e}")
    rep.check("hcon: GMRES on A_aug block-preconditioned",
              res1 < 1.0e-9, f"it {it1} res {res1:.3e}")
    dz0 = np.max(np.abs(sol0 - zman)) / max(1.0, np.max(np.abs(zman)))
    rep.check("hcon: unpreconditioned GMRES recovers z_man",
              dz0 < 1.0e-7, f"dz {dz0:.3e}")

    print(f"SUMMARY rho_nomass {rho_nm:.4f} rho_mass {rho_m:.4f} "
          f"gmres_np {it0} gmres_bp {it1}")
    rep.summary()


if __name__ == "__main__":
    main()
