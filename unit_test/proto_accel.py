#!/usr/bin/env python3
"""Prototype harness for JFNK acceleration techniques (SOB surrogate).

Replicates the frozen-J surrogate residual
    F(x) = A (n0 e^x) + C x + r            (x = ln n/n0, log-unknown)
and its exact Jacobian
    J(x) = A diag(n0 e^x) + C            (masked rows zeroed)
in numpy (float64), then prototypes:

  T1  log-unknown vs linear-unknown conditioning (already default in
      production via JFNK_LOGX=T -- this script quantifies it).
  T2  GMRES-DR: harmonic-Ritz deflated restarting, recycled across
      Newton steps (Parks et al. 2006; Morgan 2002 GCRO-DR).
  T3  Eisenstat-Walker Choice-2 forcing -- already in production
      (jfnk_eta, jfnk_solver_v1.f90:706); measured here vs fixed tol.
  T4  Schur-complement (n,T) preconditioner vs per-depth block PC.

The default preconditioner replicates the production LSQ path
(JFNK_PRECOND=LSQ): per-depth block of J, near-null rows projected
to identity rows, then ridged normal equations
    (B^T B + lam^2 tr(B^T B)/NT I) z = B^T v
which is the variant that survives the near-singular Rydberg rows.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import lu_factor, lu_solve

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gen_se_system import build_spec, NT, N_NORMAL, RYD, FLOOR, \
    ION_ROW, CONST_ROW, T_ROW, NE_ROW


SPECS = {
    "mild": dict(nd=4, seed=11, hard=False, demand_scale=0.3,
                 ryd_eps=1e-4),
    "std": dict(nd=4, seed=7, hard=False, demand_scale=1.0,
                ryd_eps=1e-6),
    "shock": dict(nd=4, seed=13, hard=False, demand_scale=0.3,
                  ryd_eps=1e-1, tstar=3.0, lam_max=2.0,
                  n0_hi=2.0, ryd_hi=3.0),
    "radpin": dict(nd=4, seed=5, hard=True, demand_scale=1.0,
                   ryd_eps=1e-9, radpin_frac=0.66,
                   radpin_jfac=1e-6, radpin_fmag=1e2),
}


class Surrogate:
    """se_prob_mod replica in numpy (float64)."""

    def __init__(self, name):
        sp = build_spec(**SPECS[name])
        self.a, self.c, self.r = sp["a"], sp["c"], sp["r"]
        self.n0, self.xstar = sp["n0"], sp["xstar"]
        self.mask = sp["mask"].astype(bool)
        self.n, self.nd, self.nt = sp["n"], sp["nd"], sp["nt"]
        self.tauto = True            # T rows live (production regime)
        self.name = name

    def _live(self, i):
        if not self.mask[i]:
            return True
        return self.tauto and (i % self.nt) == T_ROW

    def free_idx(self):
        return np.array([i for i in range(self.n) if self._live(i)])

    def t_idx(self):
        return np.array([i for i in range(self.n)
                         if (i % self.nt) == T_ROW])

    def F(self, x):
        npop = self.n0 * np.exp(np.clip(x, -80.0, 80.0))
        f = self.a @ npop + self.c @ x + self.r
        for i in range(self.n):
            if not self._live(i):
                f[i] = 0.0
        return np.clip(np.nan_to_num(f, nan=1e20), -1e20, 1e20)

    def J(self, x):
        npop = self.n0 * np.exp(np.clip(x, -80.0, 80.0))
        J = self.a * npop[None, :] + self.c
        for i in range(self.n):
            if not self._live(i):
                J[i, :] = 0.0
        return J


# ----------------------------------------------------------------------
# preconditioners (all operate on the full vector, per-depth blocks)
# ----------------------------------------------------------------------

def _reg_rows(blk, rowtiny=1e-30):
    """Production JFNK_ROW_REG: rowmax < rowtiny*blockmax -> identity."""
    blkmax = np.abs(blk).max()
    out = blk.copy()
    for i in range(blk.shape[0]):
        if blkmax > 0 and np.abs(blk[i]).max() < rowtiny * blkmax:
            out[i, :] = 0.0
            out[i, i] = 1.0
    return out


def _lsq_solve(blk, v, lam=1e-2):
    """Ridged normal equations: (B^T B + lam^2 s I)^-1 B^T v."""
    b = _reg_rows(blk)
    g = b.T @ b
    g += (lam ** 2) * (np.trace(g) / blk.shape[0] + 1e-300) * np.eye(len(b))
    return np.linalg.solve(g, b.T @ v)


def _tik_inv_apply(blk, v, rel=1e-10):
    """Tikhonov-regularized inverse apply via SVD (for Schur Jnn)."""
    u, s, vt = np.linalg.svd(blk)
    smax = s[0] if s[0] > 0 else 1.0
    filt = s / (s * s + (rel * smax) ** 2)
    return vt.T @ (filt * (u.T @ v))


def make_pc(prob, x, kind="lsq"):
    J = prob.J(x)
    nt, nd, n = prob.nt, prob.nd, prob.n
    if kind in ("lsq", "lu", "jacobi"):
        def Minv(v):
            z = v.copy()
            for k in range(nd):
                sl = slice(k * nt, (k + 1) * nt)
                blk = J[sl, sl]
                if kind == "lsq":
                    z[sl] = _lsq_solve(blk, z[sl])
                elif kind == "lu":
                    z[sl] = lu_solve(lu_factor(_reg_rows(blk) +
                                               1e-8 * np.eye(nt)),
                                     z[sl])
                else:
                    d = np.diag(_reg_rows(blk))
                    z[sl] = np.where(np.abs(d) > 1e-30,
                                     z[sl] / np.where(np.abs(d) > 1e-30,
                                                      d, 1.0), 0.0)
            return z
        return Minv

    if kind == "schur":
        n_pop = nt - 2                     # 0..nt-3 population rows
        Jnn_lu, JnT, JTn, JTT = [], [], [], []
        for k in range(nd):
            b = k * nt
            blk = J[b:b + nt, b:b + nt]
            Jnn_lu.append(blk[:n_pop, :n_pop])
            JnT.append(blk[:n_pop, T_ROW])
            JTn.append(blk[T_ROW, :n_pop])
            JTT.append(blk[T_ROW, T_ROW])
        S = np.diag(np.array(JTT))
        for k in range(nd):
            S[k, k] -= JTn[k] @ _tik_inv_apply(Jnn_lu[k], JnT[k])
        Sf = lu_factor(S + 1e-10 * np.trace(np.abs(S)) / nd * np.eye(nd)
                       if np.trace(np.abs(S)) > 0 else S)

        def Minv(v):
            z = np.zeros_like(v)
            zn = np.zeros((nd, n_pop))
            for k in range(nd):
                b = k * nt
                zn[k] = _tik_inv_apply(Jnn_lu[k], v[b:b + n_pop])
            trhs = np.array(
                [v[k * nt + T_ROW] - JTn[k] @ zn[k] for k in range(nd)])
            zT = lu_solve(Sf, trhs)
            for k in range(nd):
                b = k * nt
                corr = _tik_inv_apply(Jnn_lu[k], JnT[k] * zT[k])
                z[b:b + n_pop] = zn[k] - corr
                z[b + T_ROW] = zT[k]
            return z
        return Minv
    raise ValueError(kind)


# ----------------------------------------------------------------------
# GMRES: restarted, right-preconditioned, optional recycled subspace
# ----------------------------------------------------------------------

def _harmonic_ritz(V, H, jend, k):
    """Smallest-magnitude harmonic Ritz vectors of the preconditioned
    operator restricted to span(V[:, :jend+1]).  Parks et al. (2006):
    G = H + h_{m+1,m}^2 H^{-T} e_m e_m^T  (upper-Hessenberg H square
    part), then Ritz vectors = V @ eigvecs(G) for smallest |eig(G)|.
    Returns n x k basis in preconditioned-operator space."""
    Hsq = H[:jend + 1, :jend + 1].copy()
    beta2 = H[jend + 1, jend] ** 2
    try:
        f = np.linalg.solve(Hsq.T, np.eye(jend + 1)[:, -1])
    except np.linalg.LinAlgError:
        f = np.linalg.lstsq(Hsq.T, np.eye(jend + 1)[:, -1],
                            rcond=None)[0]
    G = Hsq + beta2 * np.outer(f, np.eye(jend + 1)[:, -1])
    w, Z = np.linalg.eig(G)
    order = np.argsort(np.abs(w))[:k]
    return np.real(V[:, :jend + 1] @ Z[:, order])


def gmres(Jop, Minv, b, m, tol, x0=None, U=None, max_cycles=12,
          keep_k=0):
    """Right-PC restarted GMRES for Jop y = b, x = Minv y.

    If U (n x k, orthonormal, x-space) is given, each cycle first
    solves the Galerkin projection of the residual onto U
    (deflation), then runs an m-step Arnoldi on the deflated
    residual.  If keep_k > 0, returns updated recycled basis from
    harmonic Ritz vectors of the last cycle.
    """
    n = len(b)
    x = np.zeros(n) if x0 is None else x0.copy()
    r = b - Jop(x)
    total_mv = 1 if x0 is not None else 0
    res_hist = [np.linalg.norm(r)]

    def deflate(x, r, U):
        if U is None or U.size == 0:
            return x, r, 0
        AU = np.column_stack([Jop(U[:, j]) for j in range(U.shape[1])])
        c, *_ = np.linalg.lstsq(AU, r, rcond=None)
        return x + U @ c, r - AU @ c, U.shape[1]

    Uout = U
    for _ in range(max_cycles):
        if res_hist[-1] <= tol:
            break
        x, r, mv = deflate(x, r, U)
        total_mv += mv
        res_hist.append(np.linalg.norm(r))
        beta = np.linalg.norm(r)
        if beta <= tol:
            break
        V = np.zeros((n, m + 1))
        H = np.zeros((m + 1, m))
        V[:, 0] = r / beta
        jend = 0
        for j in range(m):
            w = Jop(Minv(V[:, j]))
            total_mv += 1
            for i in range(j + 1):
                H[i, j] = V[:, i] @ w
                w -= H[i, j] * V[:, i]
            H[j + 1, j] = np.linalg.norm(w)
            if H[j + 1, j] > 1e-290:
                V[:, j + 1] = w / H[j + 1, j]
            else:
                jend = j
                break
            e1 = np.zeros(j + 2)
            e1[0] = beta
            y, *_ = np.linalg.lstsq(H[:j + 2, :j + 1], e1, rcond=None)
            jend = j
            if np.linalg.norm(H[:j + 2, :j + 1] @ y - e1) <= tol:
                break
        e1 = np.zeros(jend + 2)
        e1[0] = beta
        y, *_ = np.linalg.lstsq(H[:jend + 2, :jend + 1], e1, rcond=None)
        x = x + Minv(V[:, :jend + 1] @ y)
        r = b - Jop(x)
        total_mv += 1
        res_hist.append(np.linalg.norm(r))
        if keep_k > 0 and jend >= keep_k:
            Uout, _ = np.linalg.qr(_harmonic_ritz(V, H, jend, keep_k))
    return x, total_mv, res_hist, Uout


# ----------------------------------------------------------------------
# Newton drivers and forcing terms
# ----------------------------------------------------------------------

def eta_fixed(k, fn, fn_prev, eta_prev):
    return 1e-6


def eta_ew2(k, fn, fn_prev, eta_prev, gamma=0.9, alpha=1.5,
            emin=1e-4, emax=0.9):
    """Eisenstat-Walker Choice 2 + production safeguard/clamps
    (jfnk_solver_v1.f90:706-715)."""
    if k == 1 or fn_prev <= 0:
        return emax
    e = gamma * (fn / fn_prev) ** alpha
    e = max(e, 0.5 * fn / fn_prev)
    return min(max(e, emin), emax)


def newton(prob, solve, eta_fn, nmax=12, cap=3.0):
    x = np.zeros(prob.n)
    F = prob.F(x)
    hist, tot_mv = [(0, np.linalg.norm(F), 0, 0.0)], 0
    eta_prev, tol_prev = 0.9, np.linalg.norm(F)
    for k in range(1, nmax + 1):
        fn = np.linalg.norm(F)
        if fn < 1e-8 * hist[0][1]:
            break
        eta = eta_fn(k, fn, tol_prev, eta_prev)
        s, mv = solve(prob, x, -F, eta)
        x = x + np.clip(s, -cap, cap)
        F = prob.F(x)
        eta_prev, tol_prev = eta, fn
        tot_mv += mv
        hist.append((k, np.linalg.norm(F), mv, eta))
    return x, hist, tot_mv


# ----------------------------------------------------------------------
# experiments
# ----------------------------------------------------------------------

def exp_t1(prob):
    print(f"\n===== T1 log-vs-linear conditioning [{prob.name}] =====")
    free = prob.free_idx()
    for tag, x in [("x=0 ", np.zeros(prob.n)),
                   ("x*/2", 0.5 * prob.xstar),
                   ("x*  ", prob.xstar)]:
        n = prob.n0 * np.exp(np.clip(x, -80, 80))
        Jx = prob.J(x)
        Jn = Jx * (1.0 / n)[None, :]
        out = []
        for lab, M in [("log", Jx), ("lin", Jn)]:
            sub = M[np.ix_(free, free)]
            s = np.linalg.svd(sub, compute_uv=False)
            s = s[s > 1e-15 * s[0]]
            out.append(f"{lab}: cond={s[0]/s[-1]:.3e}")
        print(f"  {tag}  " + "   ".join(out))


def exp_t3(prob, m=30):
    print(f"\n===== T3 forcing terms [{prob.name}] m={m} =====")
    for lab, efn in [("fixed eta=1e-6", eta_fixed),
                     ("EW choice-2   ", eta_ew2)]:
        def solve(p, x, rhs, eta):
            Minv = make_pc(p, x, "lsq")
            return gmres(lambda v: p.J(x) @ v, Minv, rhs, m,
                         eta * np.linalg.norm(rhs),
                         max_cycles=8)[:2]
        x, hist, mv = newton(prob, solve, efn)
        print(f"  {lab} tot_mv={mv:5d}  |F|f={hist[-1][1]:.3e}")
        print("      (it, |F|, mv, eta):",
              [(h[0], f"{h[1]:.2e}", h[2], f"{h[3]:.1e}")
               for h in hist[1:7]])


def exp_t2(prob, m=30, k=10):
    print(f"\n===== T2 GMRES-DR recycling [{prob.name}] m={m} k={k} =====")
    for lab, use_rec in [("GMRES(m)  ", False), ("GMRES-DR  ", True)]:
        state = {}

        def solve(p, x, rhs, eta):
            Minv = make_pc(p, x, "lsq")
            J = p.J(x)
            s, mv, rh, U = gmres(lambda v: J @ v, Minv, rhs, m,
                                 eta * np.linalg.norm(rhs),
                                 U=state.get("U") if use_rec else None,
                                 keep_k=k if use_rec else 0,
                                 max_cycles=8)
            if use_rec:
                state["U"] = U
            return s, mv
        x, hist, mv = newton(prob, solve, eta_ew2)
        print(f"  {lab} tot_mv={mv:5d}  |F|f={hist[-1][1]:.3e}")
        print("      per-it mv:", [h[2] for h in hist[1:]])


def exp_t4(prob, m=30):
    print(f"\n===== T4 Schur (n,T) PC [{prob.name}] m={m} =====")
    for lab, kind in [("block-LSQ", "lsq"), ("Schur(n,T)", "schur")]:
        def solve(p, x, rhs, eta, kind=kind):
            Minv = make_pc(p, x, kind)
            return gmres(lambda v: p.J(x) @ v, Minv, rhs, m,
                         eta * np.linalg.norm(rhs),
                         max_cycles=8)[:2]
        x, hist, mv = newton(prob, solve, eta_ew2)
        tidx = prob.t_idx()
        print(f"  {lab:10s} tot_mv={mv:5d}  |F|f={hist[-1][1]:.3e}  "
              f"x_T={np.round(x[tidx], 3)} vs {np.round(prob.xstar[tidx], 3)}")
        print("      per-it mv:", [h[2] for h in hist[1:]])


if __name__ == "__main__":
    np.set_printoptions(suppress=True, linewidth=140)
    for name in ["shock", "std", "radpin"]:
        p = Surrogate(name)
        exp_t1(p)
        exp_t3(p)
        exp_t2(p)
        exp_t4(p)
