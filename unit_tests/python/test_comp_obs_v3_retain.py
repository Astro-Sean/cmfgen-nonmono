#!/usr/bin/env python3
"""Regression test for COMP_OBS_V3 NU_STORE retention (09-Nov-2026).

COMP_OBS_V3 carries the same NU_STORE shuffle fix as V2 (retain up to
OBS_FREQ*MAX_FREQ_CONV_FAC) and is the variant cmf_flux_sub_v5 calls,
but no test exercised it.  Same driver geometry as
test_comp_obs_retain.py: more CMF frequencies than NST with a mu<0 ray
so the shuffle runs; the fluxes must match an independent evaluation.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import F90DIR, TREE, Report, compile_and_link_prod, run

HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "subs/phys_con.f",
    "tools/get_lu.f",
    "unix/dir_acc_pars_gen.f",
    "subs/open_dir_acc_v1.f",
    "subs/write_direct_info_v3.f",
    "tools/chk_little_endian.f",
    "subs/comp_obs_v3.f",
]

# Driver geometry and grids (must match drv_comp_obs_v3_retain.f90).
MU = [0.0, 0.3, 0.65, 0.95, -0.5]
WTS = [0.20, 0.25, 0.25, 0.20, 0.10]
VINF = 3.0e4
RMAX = 1.0e4
I0, A, B = 2.5, 0.03, 0.1
N_OBS = 40

C_KMS = 2.99792458e5
PARSEC_CM = 3.0856e18


def obs_freqs():
    return [110.0 - i * (50.0 / (N_OBS - 1)) for i in range(N_OBS)]


def expected():
    beta = VINF / C_KMS
    gamma = 1.0 / (1.0 - beta * beta) ** 0.5
    fc = 20.0 * 3.141592653589793 * (1.0e18 * RMAX / PARSEC_CM) ** 2
    out = []
    for nu_obs in obs_freqs():
        s = 0.0
        for mu, w in zip(MU, WTS):
            fcf = gamma * (1.0 - mu * beta)
            icf = fcf ** -3.0
            i_cmf = I0 + A * nu_obs * fcf + B
            s += w * i_cmf * icf
        out.append(fc * s)
    return out


def parse_obsf(out):
    d = {}
    for ln in out.splitlines():
        t = ln.split()
        if len(t) == 3 and t[0] == "OBSF":
            d[int(t[1])] = float(t[2])
    return d


def main():
    rep = Report("COMP_OBS_V3 NU_STORE retention")
    binary = compile_and_link_prod(
        "cobsr3",
        [TREE / h for h in HELPERS] + [F90DIR / "drv_comp_obs_v3_retain.f90"],
    )

    r = run(binary, check=False)
    blob = r.stdout + r.stderr
    rep.check(
        "shuffle retention completes (no 'CMF_FREQ outside range' STOP)",
        r.returncode == 0 and "CMF_FREQ outside range" not in blob,
        detail=blob.strip().splitlines()[-1] if blob.strip() else "no output",
    )
    got = parse_obsf(r.stdout)
    rep.check("all observer bins produced", len(got) == N_OBS,
              detail=f"got {len(got)}")
    if len(got) == N_OBS:
        exp = expected()
        dmax = max(abs(got[k + 1] - exp[k]) / exp[k] for k in range(N_OBS))
        rep.check(
            "retained-window fluxes match independent evaluation",
            dmax <= 1e-9,
            detail=f"max rel diff {dmax:.3e}",
        )
    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
