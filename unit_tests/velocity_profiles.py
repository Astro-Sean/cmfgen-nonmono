#!/usr/bin/env python3
"""Manufactured velocity-profile library for the non-monotonic test suite.

Profiles V001-V012 give every non-monotonic solver path one consistent
set of input velocity fields.  Each profile is a smooth analytic
function of R, so R, V, the exact dV/dR and sigma are all known in
closed form and can be checked against the production discretization.

Conventions (CMFGEN indexing):

- R(1) is the outer radius, R(ND) the inner: R(1) > R(2) > ... > R(ND).
- V is in km/s; R is a normalized radial coordinate (sigma is
  dimensionless, so the absolute R units do not matter).
- sigma = d ln V / d ln R - 1 = (R/V) dV/dR - 1, matching SIGMA_TRUE
  in new_main/subs/recompute_sigma_true.f.  Expansion (dV/dR > 0)
  gives sigma > -1; compression (dV/dR < 0) gives sigma < -1.
- Profiles stay smooth and complex-step safe: only exp, log, tanh and
  powers act on R, so ``complex_step_dvdr`` can rederive dV/dR
  independently of the hand-coded analytic derivative.
"""

import inspect
from dataclasses import dataclass, field

import numpy as np

# ----------------------------------------------------------------------
# Profile functions.  Each takes the radius array (any dtype, complex
# allowed) plus keyword parameters and returns V only; the analytic
# derivative lives in the paired _dV function.  Splitting V and dV/dR
# into separate callables lets the complex-step check stay honest: it
# differentiates V alone and never reads the analytic derivative.
# ----------------------------------------------------------------------


def _v001(r, a=8.0e3):
    """Homologous expansion, V = a R (sigma = 0 exactly)."""
    return a * r


def _d001(r, a=8.0e3):
    return np.full(np.asarray(r).shape, a, dtype=np.asarray(r).dtype)


def _v002(r, vinf=1.2e4, rc=0.02, beta=1.5):
    """Monotonic beta law, V = vinf (1 - rc/R)^beta."""
    return vinf * (1.0 - rc / r) ** beta


def _d002(r, vinf=1.2e4, rc=0.02, beta=1.5):
    return vinf * beta * (rc / r**2) * (1.0 - rc / r) ** (beta - 1.0)


def _v003(r, v0=4.0e3, b=-0.15, r0=1.0):
    """Mild deceleration, V = v0 (R/r0)^b: sigma = b - 1 everywhere."""
    return v0 * (r / r0) ** b


def _d003(r, v0=4.0e3, b=-0.15, r0=1.0):
    return b * _v003(r, v0, b, r0) / r


def _dip(r, amp, r_shock, width):
    """Gaussian depression of V in ln R centred on r_shock.

    u = ln(R/Rs)/w measures the distance in units of the width; the
    dip is a multiplicative factor (1 - amp exp(-u^2/2)) on the
    homologous baseline, so amp < 1 keeps V > 0.
    """
    u = np.log(r / r_shock) / width
    return 1.0 - amp * np.exp(-0.5 * u**2)


def _ddip(r, amp, r_shock, width):
    """d/dR of the dip factor (1 - amp exp(-u^2/2)), u = ln(R/Rs)/w."""
    u = np.log(r / r_shock) / width
    return amp * np.exp(-0.5 * u**2) * u / (width * r)


def _v_dips(r, a, dips):
    """Homologous V = a R carrying one or more Gaussian dips."""
    out = a * r
    for amp, r_shock, width in dips:
        out = out * _dip(r, amp, r_shock, width)
    return out


def _dv_dips(r, a, dips):
    out = np.zeros_like(np.asarray(r), dtype=np.asarray(r).dtype)
    for amp, r_shock, width in dips:
        out = out + _ddip(r, amp, r_shock, width) / _dip(
            r, amp, r_shock, width
        )
    out = (1.0 / r + out) * _v_dips(r, a, dips)
    return out


def _v004(r, a=8.0e3, amp=0.55, r_shock=0.35, width=0.10):
    """Single velocity reversal: one Gaussian dip in ln R."""
    return _v_dips(r, a, [(amp, r_shock, width)])


def _d004(r, a=8.0e3, amp=0.55, r_shock=0.35, width=0.10):
    return _dv_dips(r, a, [(amp, r_shock, width)])


def _v005(r, a=8.0e3, amp=0.9, r_shock=0.35, width=0.06):
    """Strong compression, sigma << -1 inside the dip."""
    return _v_dips(r, a, [(amp, r_shock, width)])


def _d005(r, a=8.0e3, amp=0.9, r_shock=0.35, width=0.06):
    return _dv_dips(r, a, [(amp, r_shock, width)])


def _v006(
    r,
    a=8.0e3,
    amps=(0.5, 0.5, 0.5),
    shocks=(0.60, 0.35, 0.15),
    widths=(0.08, 0.08, 0.08),
):
    """Multiple reversals: three separated Gaussian dips."""
    return _v_dips(r, a, list(zip(amps, shocks, widths)))


def _d006(
    r,
    a=8.0e3,
    amps=(0.5, 0.5, 0.5),
    shocks=(0.60, 0.35, 0.15),
    widths=(0.08, 0.08, 0.08),
):
    return _dv_dips(r, a, list(zip(amps, shocks, widths)))


def _sech2(u):
    """sech^2(u) without overflow for large |u|."""
    e = np.exp(-2.0 * np.abs(u))
    return 4.0 * e / (1.0 + e) ** 2


def _jump(r, a, dj, r_shock, width):
    """Homologous V carrying a tanh deceleration step in ln R.

    u = ln(R/Rs)/w; the multiplicative factor falls from 1 (inner)
    to 1 - dj (outer) across the transition, so V drops by the
    fraction dj near Rs while staying positive for dj < 1.  Away
    from Rs the profile rejoins homologous expansion, which keeps
    the compression region localized and gives two gradient sign
    changes bracketing the shock.
    """
    u = np.log(r / r_shock) / width
    return a * r * (1.0 - dj * 0.5 * (1.0 + np.tanh(u)))


def _djump(r, a, dj, r_shock, width):
    u = np.log(r / r_shock) / width
    f = 1.0 - dj * 0.5 * (1.0 + np.tanh(u))
    df = -dj * 0.5 * _sech2(u) / (width * r)
    return a * f + a * r * df


def _v007(r, a=8.0e3, dj=0.8, r_shock=0.35, width=0.03):
    """Narrow shock transition (tanh in ln R)."""
    return _jump(r, a, dj, r_shock, width)


def _d007(r, a=8.0e3, dj=0.8, r_shock=0.35, width=0.03):
    return _djump(r, a, dj, r_shock, width)


def _v008(r, a=8.0e3, amp=0.8, r_shock=0.35, width=0.40):
    """Broad compression shell: wide Gaussian dip."""
    return _v_dips(r, a, [(amp, r_shock, width)])


def _d008(r, a=8.0e3, amp=0.8, r_shock=0.35, width=0.40):
    return _dv_dips(r, a, [(amp, r_shock, width)])


def _v009(r, v0=5.0e3, eps=3.0e-3, r_shock=0.35, width=0.10):
    """Near-zero velocity gradient: V ~ v0 (1 + eps tanh(u))."""
    u = np.log(r / r_shock) / width
    return v0 * (1.0 + eps * np.tanh(u))


def _d009(r, v0=5.0e3, eps=3.0e-3, r_shock=0.35, width=0.10):
    u = np.log(r / r_shock) / width
    return v0 * eps * (1.0 / np.cosh(u)) ** 2 / (width * r)


def _v010(r, a=8.0e3, amp=0.985, r_shock=0.35, width=0.30):
    """Nearly stagnant region: broad deep dip to a few percent of V."""
    return _v_dips(r, a, [(amp, r_shock, width)])


def _d010(r, a=8.0e3, amp=0.985, r_shock=0.35, width=0.30):
    return _dv_dips(r, a, [(amp, r_shock, width)])


def _v011(
    r,
    a=8.0e3,
    amps=(0.6, 0.6),
    shocks=(0.42, 0.30),
    widths=(0.08, 0.08),
):
    """Two closely spaced reversals (separation ~2 widths)."""
    return _v_dips(r, a, list(zip(amps, shocks, widths)))


def _d011(
    r,
    a=8.0e3,
    amps=(0.6, 0.6),
    shocks=(0.42, 0.30),
    widths=(0.08, 0.08),
):
    return _dv_dips(r, a, list(zip(amps, shocks, widths)))


def _v012(r, a=8.0e3, dj=0.9, r_shock=0.35, width=0.06):
    """Strong smooth velocity jump (10x drop, still resolved)."""
    return _jump(r, a, dj, r_shock, width)


def _d012(r, a=8.0e3, dj=0.9, r_shock=0.35, width=0.06):
    return _djump(r, a, dj, r_shock, width)


# ----------------------------------------------------------------------
# Registry
# ----------------------------------------------------------------------

_PROFILES = {
    "V001": ("monotonic homologous expansion", _v001, _d001),
    "V002": ("monotonic beta-law expansion", _v002, _d002),
    "V003": ("mild deceleration", _v003, _d003),
    "V004": ("single velocity reversal", _v004, _d004),
    "V005": ("strong compression sigma<-1", _v005, _d005),
    "V006": ("multiple velocity reversals", _v006, _d006),
    "V007": ("narrow shock transition", _v007, _d007),
    "V008": ("broad compression shell", _v008, _d008),
    "V009": ("near-zero velocity gradient", _v009, _d009),
    "V010": ("nearly stagnant velocity region", _v010, _d010),
    "V011": ("closely spaced multiple reversals", _v011, _d011),
    "V012": ("strong velocity jump on smooth grid", _v012, _d012),
}


def profile_ids():
    return list(_PROFILES)


@dataclass
class VelocityProfile:
    """One manufactured velocity field on a CMFGEN-ordered grid."""

    pid: str
    title: str
    r: np.ndarray
    v: np.ndarray
    dvdr: np.ndarray
    params: dict = field(default_factory=dict)
    _vfun: object = None
    _dfun: object = None

    @property
    def nd(self):
        return len(self.r)

    @property
    def dlnv_dlnr(self):
        """d ln V / d ln R evaluated with the analytic derivative."""
        return (self.r / self.v) * self.dvdr

    @property
    def sigma(self):
        """CMFGEN SIGMA_TRUE convention: d ln V / d ln R - 1."""
        return self.dlnv_dlnr - 1.0

    @property
    def sigma_min(self):
        return float(np.min(self.sigma))

    @property
    def sigma_max(self):
        return float(np.max(self.sigma))

    def vfun(self, r):
        """Evaluate V on an arbitrary grid (complex-safe)."""
        return self._vfun(np.asarray(r), **self.params)

    def dvfun(self, r):
        """Evaluate the analytic dV/dR on an arbitrary grid."""
        return self._dfun(np.asarray(r), **self.params)

    # -- independent metadata --------------------------------------

    def sign_changes(self):
        """Gradient sign changes: (i, r_cross) with dV/dR(i)*dV/dR(i+1) < 0.

        r_cross is linearly interpolated between the bracketing
        points in ln R.  A run of exactly-zero values collapses to a
        single change (handled by the strict < 0 product test plus
        the near-zero tolerance in callers).
        """
        out = []
        dv = np.real(self.dvdr)
        for i in range(self.nd - 1):
            if dv[i] * dv[i + 1] < 0.0:
                f = abs(dv[i]) / (abs(dv[i]) + abs(dv[i + 1]))
                rc = np.exp(
                    np.log(self.r[i])
                    + f * (np.log(self.r[i + 1]) - np.log(self.r[i]))
                )
                out.append((i, rc))
        return out

    def compression_regions(self):
        """Contiguous index ranges (lo, hi) where sigma < -1 (dV/dR<0).

        Endpoints are interpolated to the sigma = -1 crossing, so the
        returned widths are continuous estimates in ln R and R.
        """
        comp = np.real(self.dvdr) < 0.0
        regions = []
        i = 0
        n = self.nd
        while i < n:
            if not comp[i]:
                i += 1
                continue
            j = i
            while j + 1 < n and comp[j + 1]:
                j += 1
            regions.append((i, j))
            i = j + 1
        return regions

    def compression_widths(self):
        """(dln R, dR) width of each compression region."""
        out = []
        for lo, hi in self.compression_regions():
            out.append(
                (
                    float(np.log(self.r[lo] / self.r[hi])),
                    float(self.r[lo] - self.r[hi]),
                )
            )
        return out

    def perturbed(self, rel_amp, seed, nmodes=3):
        """Small smooth perturbation with a still-exact derivative.

        dV = rel_amp * V * sum_k c_k sin(2 pi k ln(R/r_in)/L + p_k) for
        k = 1..nmodes.  Coefficients come from a seeded generator, so
        the result is deterministic; the derivative stays analytic
        because the perturbation is a closed-form sum of sines.
        """
        rng = np.random.default_rng(seed)
        amps = rng.uniform(-1.0, 1.0, nmodes)
        amps = amps / np.max(np.abs(amps))
        phases = rng.uniform(0.0, 2.0 * np.pi, nmodes)
        ks = np.arange(1, nmodes + 1)
        r_in = float(self.r[-1])
        lr = np.log(self.r[0] / r_in)
        base_vfun, base_dfun = self._vfun, self._dfun
        base_params = dict(self.params)

        def pv(rr, _f=base_vfun, _p=base_params):
            rr = np.asarray(rr)
            x = np.log(rr / r_in) / lr
            s = np.zeros(np.asarray(rr).shape, dtype=np.asarray(rr).dtype)
            for k, c, p in zip(ks, amps, phases):
                s = s + c * np.sin(2.0 * np.pi * k * x + p)
            return _f(rr, **_p) * (1.0 + rel_amp * s)

        def pdv(rr, _f=base_vfun, _d=base_dfun, _p=base_params):
            rr = np.asarray(rr)
            x = np.log(rr / r_in) / lr
            s = np.zeros(np.asarray(rr).shape, dtype=np.asarray(rr).dtype)
            ds = np.zeros(np.asarray(rr).shape, dtype=np.asarray(rr).dtype)
            for k, c, p in zip(ks, amps, phases):
                s = s + c * np.sin(2.0 * np.pi * k * x + p)
                ds = ds + c * np.cos(2.0 * np.pi * k * x + p) * (
                    2.0 * np.pi * k / (rr * lr)
                )
            v = _f(rr, **_p)
            return _d(rr, **_p) * (1.0 + rel_amp * s) + v * rel_amp * ds

        r = self.r
        return VelocityProfile(
            pid=self.pid + "*",
            title=self.title + f" (+{rel_amp:g} perturb)",
            r=r,
            v=pv(r),
            dvdr=pdv(r),
            params=base_params,
            _vfun=pv,
            _dfun=pdv,
        )

    def reversed_view(self):
        """Copy with the array direction reversed (R increasing)."""
        return VelocityProfile(
            pid=self.pid + "R",
            title=self.title + " (reversed grid)",
            r=self.r[::-1].copy(),
            v=self.v[::-1].copy(),
            dvdr=self.dvdr[::-1].copy(),
            params=dict(self.params),
            _vfun=self._vfun,
            _dfun=self._dfun,
        )

    def write_rv(self, path, sigma_col=None):
        """drv_sigma_true fixture format: ND then rows "R V SIGMA"."""
        sig = np.zeros(self.nd) if sigma_col is None else sigma_col
        with open(path, "w") as f:
            f.write(f"{self.nd}\n")
            for a, b, c in zip(self.r, self.v, sig):
                f.write(f"{a:.16e} {b:.16e} {c:.16e}\n")
        return path


def make_profile(pid, nd=80, r_out=1.0, r_in=None, grid="log", **params):
    """Build one registered profile on a fresh grid.

    grid='log' is log-spaced (CMFGEN-style); grid='lin' is uniform.
    Keyword arguments forward to the profile function (e.g. r_shock,
    width, amp/dj/eps per family).  r_in defaults to 0.02, or 0.1 for
    V002 whose beta-law core radius rc=0.02 must stay inside the grid.
    """
    title, vfun, dfun = _PROFILES[pid]
    if r_in is None:
        r_in = 0.1 if pid == "V002" else 0.02
    if grid == "log":
        r = np.geomspace(r_out, r_in, nd)
    elif grid == "lin":
        r = np.linspace(r_out, r_in, nd)
    else:
        raise ValueError(f"unknown grid kind {grid!r}")
    # Filter only the kwargs the profile actually accepts so a sweep
    # can pass a superset (e.g. r_shock) to every family.
    keep = {
        k: v
        for k, v in params.items()
        if k in inspect.signature(vfun).parameters
    }
    return VelocityProfile(
        pid=pid,
        title=title,
        r=r,
        v=vfun(r, **keep),
        dvdr=dfun(r, **keep),
        params=keep,
        _vfun=vfun,
        _dfun=dfun,
    )


def complex_step_dvdr(prof, h=1.0e-30):
    """Independent dV/dR by the complex-step trick.

    Im(V(r + i h))/h is the exact derivative for analytic V to O(h^2);
    with h ~ 1e-30 it is machine-precision and shares nothing with the
    hand-coded analytic derivative, so a formula error in _dXXX shows
    up as a real discrepancy.
    """
    rr = prof.r.astype(complex)
    return np.imag(prof.vfun(rr + 1j * h)) / h


def sigma_true_core_mirror(r, v, rtol=1.0e-12):
    """Independent python mirror of SIGMA_TRUE_CORE's hybrid stencil.

    Re-coded from the documented algorithm (quadratic in ln R on
    same-sign sub-decade stencils, quadratic in linear R at sign
    changes and boundary stencils, |sigma| clipped at 1E6, V == 0 or
    collapsed spacing reported undefined).  Kept in the library so
    every solver test shares one reference stencil; it exists to
    catch coding errors in the Fortran routine, not to restate the
    analytic derivative.
    """
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    nd = len(r)
    s = np.zeros(nd)
    n_undef = 0
    if nd < 3:
        return s, nd

    def clip(x):
        return max(-1.0e6, min(1.0e6, x))

    def log_quad(v0, v1, v2, hp, hm, one_sided=False):
        if one_sided:
            return (
                np.log(abs(v0)) * (2.0 * hp + hm) / (hp * (hp + hm))
                - np.log(abs(v1)) * (hp + hm) / (hp * hm)
                + np.log(abs(v2)) * hp / (hm * (hp + hm))
            )
        return (
            np.log(abs(v0)) * hm / (hp * (hp + hm))
            + np.log(abs(v1)) * (hp - hm) / (hp * hm)
            - np.log(abs(v2)) * hp / (hm * (hp + hm))
        )

    def lin_quad(v0, v1, v2, hp, hm, one_sided=False):
        if one_sided:
            return (
                v0 * (2.0 * hp + hm) / (hp * (hp + hm))
                - v1 * (hp + hm) / (hp * hm)
                + v2 * hp / (hm * (hp + hm))
            )
        return (
            v0 * hm / (hp * (hp + hm))
            + v1 * (hp - hm) / (hp * hm)
            - v2 * hp / (hm * (hp + hm))
        )

    def decade_ok(a, b, c):
        return min(abs(a), abs(b), abs(c)) > 0.1 * max(abs(a), abs(b), abs(c))

    for i in range(1, nd - 1):
        hp = r[i - 1] - r[i]
        hm = r[i] - r[i + 1]
        if min(hp, hm) > rtol * r[i]:
            if (
                v[i - 1] * v[i + 1] > 0.0
                and v[i] * v[i - 1] > 0.0
                and decade_ok(v[i - 1], v[i], v[i + 1])
            ):
                hp = np.log(r[i - 1] / r[i])
                hm = np.log(r[i] / r[i + 1])
                s[i] = log_quad(v[i - 1], v[i], v[i + 1], hp, hm) - 1.0
            elif v[i] != 0.0:
                s[i] = (r[i] / v[i]) * lin_quad(
                    v[i - 1], v[i], v[i + 1], hp, hm
                ) - 1.0
            else:
                n_undef += 1
                continue
            s[i] = clip(s[i])
        else:
            n_undef += 1

    hp = r[0] - r[1]
    hm = r[1] - r[2]
    if min(hp, hm) > rtol * r[0]:
        if v[0] * v[1] > 0.0 and v[1] * v[2] > 0.0 and decade_ok(*v[:3]):
            hp = np.log(r[0] / r[1])
            hm = np.log(r[1] / r[2])
            s[0] = log_quad(v[0], v[1], v[2], hp, hm, True) - 1.0
            s[0] = clip(s[0])
        elif v[0] != 0.0:
            s[0] = (r[0] / v[0]) * lin_quad(v[0], v[1], v[2], hp, hm, True) - 1.0
            s[0] = clip(s[0])
        else:
            n_undef += 1
    else:
        n_undef += 1

    # Inner edge: same one-sided quadratic but mirrored - the boundary
    # point is V(ND) here, so the (2*hp+hm) weight sits on the last
    # stencil point and the overall sign is negated (R runs inward).
    hp = r[nd - 2] - r[nd - 1]
    hm = r[nd - 3] - r[nd - 2]
    if min(hp, hm) > rtol * r[nd - 1]:
        if (
            v[nd - 1] * v[nd - 2] > 0.0
            and v[nd - 2] * v[nd - 3] > 0.0
            and decade_ok(v[nd - 1], v[nd - 2], v[nd - 3])
        ):
            hp = np.log(r[nd - 2] / r[nd - 1])
            hm = np.log(r[nd - 3] / r[nd - 2])
            dv = (
                np.log(abs(v[nd - 3])) * hp / (hm * (hp + hm))
                - np.log(abs(v[nd - 2])) * (hp + hm) / (hp * hm)
                + np.log(abs(v[nd - 1]))
                * (2.0 * hp + hm)
                / (hp * (hp + hm))
            )
            s[nd - 1] = -dv - 1.0
            s[nd - 1] = clip(s[nd - 1])
        elif v[nd - 1] != 0.0:
            dv = (
                v[nd - 3] * hp / (hm * (hp + hm))
                - v[nd - 2] * (hp + hm) / (hp * hm)
                + v[nd - 1] * (2.0 * hp + hm) / (hp * (hp + hm))
            )
            s[nd - 1] = (r[nd - 1] / v[nd - 1]) * (-dv) - 1.0
            s[nd - 1] = clip(s[nd - 1])
        else:
            n_undef += 1
    else:
        n_undef += 1

    return s, n_undef


def parse_driver_sigma(out):
    """Parse drv_sigma_true 'core' output -> (r, v, sigma_true, n_undef)."""
    rows = {}
    n_undef = None
    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "N_UNDEF":
            n_undef = int(t[1])
        elif t[0].isdigit():
            rows[int(t[0]) - 1] = [float(x) for x in t[1:]]
    n = max(rows) + 1
    arr = np.array([rows[i] for i in range(n)])
    return arr[:, 0], arr[:, 1], arr[:, 2], n_undef


def sign_change_crossings(x, s_ref):
    """Indices where a reference curve (e.g. dvdr or sigma+1) crosses
    zero, for comparing discrete sign-change locations."""
    out = []
    for i in range(len(x) - 1):
        if s_ref[i] * s_ref[i + 1] < 0.0:
            f = abs(s_ref[i]) / (abs(s_ref[i]) + abs(s_ref[i + 1]))
            out.append((i, x[i] + f * (x[i + 1] - x[i])))
    return out
