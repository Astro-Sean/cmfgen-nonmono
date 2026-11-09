#!/usr/bin/env python3
"""Regression tests for COMP_OBS_V2 (CMF -> observer-frame spectrum).

Covers the flux-audit items:
  A: zero-velocity limit (VINF=0 -> observer flux equals CMF result)
  E: relativistic transform vs independently evaluated expectation
  F: coverage/non-finite diagnostics (STOP, never silent garbage)
  G: Jy normalization scaling with RMAX

The expectation values are computed in Python directly from the
documented equations:
    FCF_i  = gamma * (1 - mu_i * beta)
    ICF_i  = FCF_i**-3   (DO_FULL_REL=T) else 1
    F_obs  = FLUX_CONST * sum_i w_i * I_CMF(nu_obs*FCF_i) * ICF_i
    FLUX_CONST = 20*pi*(1e18*RMAX/PARSEC)^2   ('IPLUS', RMAX in 1e10 cm)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import F90DIR, TREE, Report, compile_and_link_prod, run

HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "subs/phys_con.f",
    "subs/comp_obs_v2.f",
]

# Driver geometry and grids (must match drv_comp_obs.f90).
MU = [0.0, 0.3, 0.65, 0.95]
WTS = [0.20, 0.30, 0.30, 0.20]
OBS_FREQ = [57.0, 55.0, 53.0, 51.0]
I0, A, B = 2.5, 0.03, 0.1

C_KMS = 2.99792458e5
PARSEC_CM = 3.0856e18


def flux_const(rmax):
    return 20.0 * 3.141592653589793 * (1.0e18 * rmax / PARSEC_CM) ** 2


def expected(obs_freqs, vinf, rmax, do_rel, a, b, i0):
    beta = vinf / C_KMS
    gamma = 1.0 / (1.0 - beta * beta) ** 0.5
    out = []
    for nu_obs in obs_freqs:
        s = 0.0
        for mu, w in zip(MU, WTS):
            fcf = gamma * (1.0 - mu * beta)
            icf = fcf ** -3.0 if do_rel else 1.0
            i_cmf = i0 + a * nu_obs * fcf + b
            s += w * i_cmf * icf
        out.append(flux_const(rmax) * s)
    return out


def parse_obsf(out):
    d = {}
    for ln in out.splitlines():
        t = ln.split()
        if len(t) == 3 and t[0] == "OBSF":
            d[int(t[1])] = float(t[2])
    return d


def build():
    return compile_and_link_prod(
        "cobs", [TREE / h for h in HELPERS] + [F90DIR / "drv_comp_obs.f90"]
    )


def main():
    rep = Report("COMP_OBS_V2")
    binary = build()

    # Case A: zero-velocity limit, constant field.
    # VINF=0 -> gamma=1, FCF=1: I_obs(nu)=I0 -> F=FLUX_CONST*I0*sum(w).
    got = parse_obsf(run(binary, "const").stdout)
    exp = expected(OBS_FREQ, 0.0, 1.0e4, True, 0.0, 0.0, I0)
    dmax = max(abs(got[k + 1] - exp[k]) for k in range(4))
    rep.check(
        "zerov limit: observer flux equals CMF flux",
        dmax <= 1e-9 * max(exp),
        detail=f"max abs diff {dmax:.3e} (flux~{exp[0]:.3e})",
    )

    # Case E1: relativistic transform, linear field, VINF=1e4 km/s.
    got = parse_obsf(run(binary, "lin").stdout)
    exp = expected(OBS_FREQ, 1.0e4, 1.0e4, True, A, B, I0)
    dmax = max(
        abs(got[k + 1] - exp[k]) / exp[k] for k in range(4)
    )
    rep.check(
        "relativistic FCF/ICF transform matches independent eval",
        dmax <= 1e-10,
        detail=f"max rel diff {dmax:.3e}",
    )

    # Case E2: frequency shift alone (DO_FULL_REL=F -> ICF=1).
    got = parse_obsf(run(binary, "norel").stdout)
    exp = expected(OBS_FREQ, 1.0e4, 1.0e4, False, A, B, I0)
    dmax = max(
        abs(got[k + 1] - exp[k]) / exp[k] for k in range(4)
    )
    rep.check(
        "non-rel intensity: frequency transform only (ICF=1)",
        dmax <= 1e-10,
        detail=f"max rel diff {dmax:.3e}",
    )

    # Case G: flux normalization scales as RMAX^2 (rmax2 vs const).
    g1 = parse_obsf(run(binary, "const").stdout)
    g2 = parse_obsf(run(binary, "rmax2").stdout)
    ratio = [g2[k] / g1[k] for k in sorted(g1)]
    dmax = max(abs(r - 4.0) for r in ratio)
    rep.check(
        "Jy normalization scales as RMAX^2",
        dmax <= 1e-9,
        detail=f"flux ratios {['%.6f' % r for r in ratio]}",
    )

    # Case F1: non-finite IPLUS must produce the diagnostic, not a
    # silently stored NaN.
    r = run(binary, "nan", check=False)
    blob = r.stdout + r.stderr
    rep.check(
        "non-finite IPLUS rejected with diagnostic",
        "non-finite CMF intensity" in blob,
        detail=blob.strip().splitlines()[-1] if blob.strip() else "no output",
    )

    # Case F2: OBS_FREQ*FCF(mu=0) above NU_STORE(1) must STOP, not
    # divide by zero into a silent NaN flux.
    r = run(binary, "range", check=False)
    blob = r.stdout + r.stderr
    rep.check(
        "insufficient CMF coverage rejected (no silent extrapolation)",
        "CMF_FREQ outside range" in blob,
        detail=blob.strip().splitlines()[-1] if blob.strip() else "no output",
    )

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
