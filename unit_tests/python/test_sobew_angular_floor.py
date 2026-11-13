"""Angular T1 census for the non-monotonic Sobolev EW floor.

Regression coverage for the SIGMA_TRUE EW audit
(docs/audits/SOB_EW_SIGMA_TRUE_AUDIT.md):

A. sigma=-1.037 at the zero-gradient direction mu^2=1/|sigma|:
   |1+sigma*mu^2| vanishes, so T1 must equal the -0.2*(1+sigma)=0.0074
   floor arm, not the radial value 0.037.
B. sigma=-1.037 on a radial ray (mu=1): T1=|1+sigma|=0.037.
C. sigma=-155.6 (deepest compression of the m15 2.82d model): the
   -0.2*(1+sigma)=30.92 arm dominates the whole angular quadrature.
D. SIGMA_TRUE_CORE must flag degenerate profiles N_UNDEF>0, and the
   OBS driver must fail closed on them by default (ungated abort,
   ALLOW_SIGMA_UNDEF=T the only override).
E. The production m15 profile keeps N_UNDEF=0 and reproduces the
   recorded sigma_true extremes.
F. A monotonic power-law field is defined everywhere and the
   SIGMA branch is untouched upstream behaviour.
"""

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, TREE, WORKDIR, Report, compile_and_link_prod, run
from test_sob_tau_floor import (
    ND,
    NP,
    beta_predict,
    build,
    exponx,
    parse_rows,
    write_case,
)

SIGMA_TRUE_SRC = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "new_main" / "mod_subs" / "control_variable_mod.f",
    TREE / "unix" / "term_io.f",
    TREE / "new_main" / "subs" / "recompute_sigma_true.f",
]


def beta_at(rows):
    # sobjsim rows are [beta, jbar, znet, vb2] after the index column.
    return rows[:, 0]


def ray_t1(sig, r, p):
    """Production floor formula evaluated per (depth, ray)."""
    mu2 = 1.0 - (p[None, :] / r[:, None]) ** 2
    t1 = np.abs(1.0 + sig * mu2)
    if sig < 0.0:
        t1 = np.maximum(t1, max(-0.2 * (1.0 + sig), 1.0e-4))
    return t1, mu2


def sigma_true_bin():
    return compile_and_link_prod(
        "drv_sigma_true", SIGMA_TRUE_SRC + [F90DIR / "drv_sigma_true.f90"]
    )


def write_rv(name, r, v):
    p = WORKDIR / name
    with open(p, "w") as f:
        f.write(f"{len(r)}\n")
        for a, b in zip(r, v):
            f.write(f"{a:.16e} {b:.16e} {-0.999:.16e}\n")
    return p


def parse_sigma_true(out):
    rows = {}
    n_undef = None
    for ln in out.splitlines():
        t = ln.split()
        if t and t[0] == "N_UNDEF":
            n_undef = int(t[1])
        elif t and t[0].isdigit():
            rows[int(t[0]) - 1] = float(t[3])
    n = max(rows) + 1
    return n_undef, np.array([rows[i] for i in range(n)])


def maxrel(a, b):
    return float(np.max(np.abs(1.0 - a / b)))


def main():
    rep = Report("EW angular T1 floor + undefined-cell policy")
    bins = build()
    sobjsim = bins[("sobjsim", "new")]
    sobewg2 = bins[("sobewg2", "new")]

    # -- A: sigma=-1.037 with a ray on the zero-gradient direction at
    # the innermost depth plus three filler rays.  The zero-dir ray
    # must sit on the -0.2*(1+sigma)=0.0074 floor arm.
    sig = -1.037
    pfrac = [np.sqrt(1.0 - 1.0 / abs(sig)), 0.5, 0.8, 0.95]
    path, r, p, w, gamh = write_case("af_a1", sig, p_frac=pfrac)
    rows, _ = parse_rows(run(sobjsim, "sobjsim", str(path)).stdout)
    t1, mu2 = ray_t1(sig, r, p)
    rep.check(
        "A sobjsim sigma=-1.037 zero-dir: BETA==floored formula",
        np.allclose(beta_at(rows), beta_predict(sig, r, p, w, gamh), rtol=1e-9),
        f"maxrel={maxrel(beta_at(rows), beta_predict(sig, r, p, w, gamh)):.2e}",
    )
    rep.check(
        "A zero-dir ray T1 = 0.0074 floor arm (not 0.037)",
        abs(t1[ND - 1, 0] - 0.0074) < 1.0e-10,
        f"t1={t1[ND-1,0]:.6f} raw={1 + sig * mu2[ND-1,0]:.3e}",
    )
    rows2, head2 = parse_rows(run(sobewg2, "sobewg2", str(path)).stdout)
    rep.check(
        "A sobewg2 sigma=-1.037 zero-dir: EW finite",
        np.isfinite(head2.get("EW", np.nan)) and np.isfinite(rows2).all(),
        f"EW={head2.get('EW'):.4e}",
    )

    # -- B: radial ray p=0 -> mu=1 at every depth -> T1=|1+sigma|=0.037.
    pfrac = [0.0, 0.5, 0.8, 0.95]
    path, r, p, w, gamh = write_case("af_b", sig, p_frac=pfrac)
    rows, _ = parse_rows(run(sobjsim, "sobjsim", str(path)).stdout)
    t1, mu2 = ray_t1(sig, r, p)
    rep.check(
        "B sobjsim sigma=-1.037 radial: T1=0.037",
        np.allclose(beta_at(rows), beta_predict(sig, r, p, w, gamh), rtol=1e-9)
        and np.allclose(t1[:, 0], 0.037, atol=1.0e-12),
        f"t1(radial)={t1[0,0]:.4f} maxrel={maxrel(beta_at(rows), beta_predict(sig, r, p, w, gamh)):.2e}",
    )

    # -- C: sigma=-155.6 (model depth-62 value): -0.2*(1+sigma)=30.92
    # dominates the zero-dir ray; |1+sigma|=154.6 on the radial ray.
    sig = -155.6
    pfrac = [np.sqrt(1.0 - 1.0 / abs(sig)), 0.5, 0.8, 0.95]
    path, r, p, w, gamh = write_case("af_c_zero", sig, p_frac=pfrac)
    rows, _ = parse_rows(run(sobjsim, "sobjsim", str(path)).stdout)
    t1, mu2 = ray_t1(sig, r, p)
    rep.check(
        "C sobjsim sigma=-155.6 zero-dir: floor 30.92 arm",
        np.allclose(beta_at(rows), beta_predict(sig, r, p, w, gamh), rtol=1e-9)
        and abs(t1[ND - 1, 0] - 30.92) < 1.0e-8,
        f"t1={t1[ND-1,0]:.3f} raw={1 + sig * mu2[ND-1,0]:.3e}",
    )
    pfrac = [0.0, 0.2, 0.6, 0.99]
    path, r, p, w, gamh = write_case("af_c_sweep", sig, p_frac=pfrac)
    rows, _ = parse_rows(run(sobjsim, "sobjsim", str(path)).stdout)
    t1, mu2 = ray_t1(sig, r, p)
    rep.check(
        "C sobjsim sigma=-155.6 mu sweep: angular pattern matches",
        np.allclose(beta_at(rows), beta_predict(sig, r, p, w, gamh), rtol=1e-9)
        and np.allclose(t1[:, 0], 154.6, atol=1.0e-10),
        f"t1(radial)={t1[0,0]:.2f} maxrel={maxrel(beta_at(rows), beta_predict(sig, r, p, w, gamh)):.2e}",
    )

    # -- D: degenerate profile flags N_UNDEF>0, and the OBS abort is
    # ungated (fail closed) with ALLOW_SIGMA_UNDEF the only override.
    stbin = sigma_true_bin()
    r_d = np.geomspace(1.0e5, 1.0e3, 25)
    v_d = np.geomspace(50.0, 6.0e3, 25)
    v_d[-1] = 0.0  # stagnant boundary cell -> undefined
    rv = write_rv("af_d_undef.txt", r_d, v_d)
    n_undef, st = parse_sigma_true(run(stbin, "core", str(rv)).stdout)
    rep.check(
        "D sigma_true_core flags stagnant boundary cell",
        n_undef == 1 and st[-1] == 0.0,
        f"N_UNDEF={n_undef} st[-1]={st[-1]}",
    )
    src = (TREE / "obs" / "cmf_flux_sub_v5.f").read_text()
    ungated = re.search(
        r"IF\(\s*\.NOT\.\s*ALLOW_SIGMA_UNDEF\s*\)\s*THEN", src
    )
    gated = re.search(r"ALLOW_SIGMA_UNDEF\s*\.AND\.", src)
    rep.check(
        "D cmf_flux_sub_v5 abort on N_SIG_UNDEF is ungated",
        ungated is not None and gated is None,
        "fail-closed default on every SIGMA_TRUE consumer path",
    )

    # -- E: production m15 2.82d profile -> N_UNDEF=0, recorded range.
    txt = open(
        "/ptmp/sbrennan/SN2024pba/cmfgen_validation/"
        "m15_2p82d_FRESH_sob_dieoff/cmfgen/2.82days/RVTJ"
    ).read().splitlines()

    def block(after_hdr):
        vals = []
        i = after_hdr
        while len(vals) < 70:
            vals += [
                float(x.replace("D", "E")) for x in txt[i].split()
            ]
            i += 1
        return np.array(vals[:70])

    rv = write_rv("af_e_m15.txt", block(13), block(23))
    n_undef, st = parse_sigma_true(run(stbin, "core", str(rv)).stdout)
    rep.check(
        "E m15 profile fully defined, sigma_true range matches",
        n_undef == 0
        and abs(st.min() - (-155.6105)) < 1.0e-3
        and abs(st.max() - (-1.0370)) < 1.0e-3,
        f"N_UNDEF={n_undef} min={st.min():.4f} max={st.max():.4f}",
    )

    # -- F: monotonic power-law V ~ R^-0.7 -> sigma_true = -1.7
    # everywhere except the one-sided boundary stencils.
    r_f = np.geomspace(1.0e5, 1.0e3, 25)
    v_f = 3.0e3 * (r_f / r_f[-1]) ** (-0.7)
    rv = write_rv("af_f_mono.txt", r_f, v_f)
    n_undef, st = parse_sigma_true(run(stbin, "core", str(rv)).stdout)
    rep.check(
        "F monotonic power law: defined, sigma_true=-1.7",
        n_undef == 0
        and np.allclose(st, -1.7, atol=5e-3, rtol=5e-3),
        f"N_UNDEF={n_undef} st range {st.min():.4f}..{st.max():.4f}",
    )

    rep.summary()


if __name__ == "__main__":
    main()
