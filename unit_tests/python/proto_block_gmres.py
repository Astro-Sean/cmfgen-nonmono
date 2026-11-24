"""Phase-C prototype: coupled wavelength-depth block solve vs lagged map.

Manufactured frequency-depth system mimicking the CMF moment structure:
per-frequency depth-tridiagonal blocks (diffusion-like, diagonally
dominant) plus upwind frequency coupling whose direction keys on the
local sign of the velocity-gradient coefficient. A compression band in
the frequency-depth plane reproduces the regime where the lagged
iteration diverges (measured growth c/(1-a) in the drv_mom_nu iter
probe). Questions: (1) spectral radius of the lagged (block-Jacobi in
frequency) iteration; (2) GMRES iterations to convergence on the full
coupled system, unpreconditioned and with the per-frequency tridiagonal
inverse as block-Jacobi preconditioner - the matrix-free analogue of
the existing per-frequency DGTSV solves.
"""
import numpy as np


def build_system(nfreq=24, ndepth=16, coupling=8.0, comp_lo=8, comp_hi=16):
    """Assemble dense coupled operator for the toy system.

    Per frequency l: depth tridiagonal, TB=2.6 on diagonal,
    TA=TC=-0.4 (interior), M-matrix. Frequency coupling of magnitude
    `coupling` reads l-1 where expanding (sign +1) and l+1 where
    compressing (sign -1); the diagonal keeps the implicit upwind
    weight |coupling| (the Baron-Hauschildt chi-hat absorption),
    exactly as PSI does in the production assembly. Compressing rows
    occupy a frequency band [comp_lo, comp_hi) at all depths.
    """
    n = nfreq * ndepth
    a = np.zeros((n, n))

    def idx(l, i):
        return l * ndepth + i

    for l in range(nfreq):
        compressing = comp_lo <= l < comp_hi
        for i in range(ndepth):
            r = idx(l, i)
            a[r, r] = 2.6 + coupling
            if i > 0:
                a[r, idx(l, i - 1)] = -0.4
            if i < ndepth - 1:
                a[r, idx(l, i + 1)] = -0.4
            if compressing:
                if l + 1 < nfreq:
                    a[r, idx(l + 1, i)] = -coupling
                else:
                    a[r, r] -= coupling  # frequency boundary: no neighbour
            else:
                if l - 1 >= 0:
                    a[r, idx(l - 1, i)] = -coupling
                else:
                    a[r, r] -= coupling
    return a, nfreq, ndepth


def block_diag_inverse(a, nfreq, ndepth):
    """Inverse of the frequency block-diagonal part (exact, dense)."""
    n = nfreq * ndepth
    bd = np.zeros((n, n))
    for l in range(nfreq):
        sl = slice(l * ndepth, (l + 1) * ndepth)
        bd[sl, sl] = a[sl, sl]
    return np.linalg.inv(bd), bd


def lagged_spectral_radius(a, nfreq, ndepth):
    """Spectral radius of block-Jacobi in frequency: the lagged map.

    The lagged iteration solves each frequency's depth system exactly
    while reading neighbour frequencies from the previous iterate -
    precisely block-Jacobi in frequency. Iteration matrix
    M = D_b^{-1} (D_b - A); convergence requires rho(M) < 1.
    """
    bdi, bd = block_diag_inverse(a, nfreq, ndepth)
    m = bdi @ (bd - a)
    return max(abs(np.linalg.eigvals(m)))


def gmres(a, b, x0=None, m_inv=None, tol=1e-10, maxit=400):
    """Restartless GMRES, modified Gram-Schmidt, left preconditioning."""
    n = len(b)
    x = np.zeros(n) if x0 is None else x0.copy()

    def apply_m_inv(v):
        return v if m_inv is None else m_inv @ v

    r = apply_m_inv(b - a @ x)
    beta = np.linalg.norm(r)
    if beta < tol:
        return x, 0
    q = [r / beta]
    h = np.zeros((maxit + 1, maxit))
    for k in range(maxit):
        v = apply_m_inv(a @ q[k])
        for j in range(k + 1):
            h[j, k] = np.dot(q[j], v)
            v = v - h[j, k] * q[j]
        h[k + 1, k] = np.linalg.norm(v)
        if h[k + 1, k] > 1e-14:
            q.append(v / h[k + 1, k])
        # residual estimate via the Arnoldi recurrence
        e = np.zeros(k + 2)
        e[0] = beta
        hh = h[: k + 2, : k + 1]
        y, *_ = np.linalg.lstsq(hh, e, rcond=None)
        res = np.linalg.norm(e - hh @ y)
        if res < tol:
            return x, k + 1
    return x, maxit


def run_case(coupling, nfreq=24, ndepth=16):
    a, nf, nd = build_system(nfreq, ndepth, coupling)
    rng = np.random.RandomState(42)
    b = rng.rand(nf * nd)
    rho = lagged_spectral_radius(a, nf, nd)
    _, it_noprec = gmres(a, b)
    bdi, _ = block_diag_inverse(a, nf, nd)
    _, it_prec = gmres(a, b, m_inv=bdi)
    return rho, it_noprec, it_prec


if __name__ == "__main__":
    print("coupling  rho_lagged  gmres_noprec  gmres_blockjacobi")
    for c in [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]:
        rho, i0, i1 = run_case(c)
        verdict = "DIVERGES" if rho >= 1.0 else "contracts"
        print(f"{c:8.1f}  {rho:10.4f}  {i0:12d}  {i1:16d}  lagged {verdict}")
