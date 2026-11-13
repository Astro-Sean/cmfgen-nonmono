#!/usr/bin/env python3
"""Frequency-coverage checks for OBS_FRAME_SUB_V9 (P5).

Runs the production observer-frame routine on manufactured models:

- monotonic outflow, several velocity reversals, negative velocities;
- relativistic mapping and its non-relativistic limit;
- NUM_BANDS partitioning (NOS=2000 -> 2 bands vs NOS=1999 -> 1 band);
- global-grid undercoverage, blue and red edges, which must abort with
  GLOBAL_CMF_GRID_UNDERCOVERAGE rather than emit a truncated spectrum;
- convergence: extending FREQ_CMF outside the required reach must
  leave the spectrum bitwise unchanged;
- a CHI feature anchoring the interpolation brackets must change the
  flux, proving the mapping samples the feature rather than an edge
  value.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import F90DIR, PLOTDIR, TREE, WORKDIR, Report, compile_driver, run

SOURCES = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "subs" / "phys_con.f",
    TREE / "subs" / "mon_interp_fast.f",
    TREE / "subs" / "open_dir_acc_v1.f",
    TREE / "subs" / "write_direct_info_v3.f",
    TREE / "tools" / "chk_little_endian.f",
    TREE / "tools" / "get_lu.f",
    TREE / "unix" / "term_io.f",
    TREE / "unix" / "dir_acc_pars_gen.f",
    TREE / "unix" / "tune.f",
    TREE / "obs" / "obs_frame_sub_v9.f",
    F90DIR / "drv_obs_frame_cov.f",
]


def run_expect_abort(binary, out_dir, mode):
    """Run a mode expected to STOP inside the routine; return output."""
    proc = run(binary, out_dir, mode, cwd=out_dir, check=False)
    return proc.returncode, proc.stdout + proc.stderr


def main():
    rep = Report("OBS frame frequency coverage (P5)")
    out_dir = WORKDIR / "obs_frame_cov"
    out_dir.mkdir(parents=True, exist_ok=True)
    binary = compile_driver("drv_obs_frame_cov", SOURCES)
    run(binary, out_dir, "run", cwd=out_dir)

    mono = np.loadtxt(out_dir / "flux_mono.txt")
    revers = np.loadtxt(out_dir / "flux_revers.txt")
    negv = np.loadtxt(out_dir / "flux_negv.txt")
    rel = np.loadtxt(out_dir / "flux_rel.txt")
    rel_lo = np.loadtxt(out_dir / "flux_rel_lo.txt")
    norel_lo = np.loadtxt(out_dir / "flux_norel_lo.txt")
    b1 = np.loadtxt(out_dir / "flux_b1.txt")
    b2 = np.loadtxt(out_dir / "flux_b2.txt")
    base = np.loadtxt(out_dir / "flux_base.txt")
    ext = np.loadtxt(out_dir / "flux_ext.txt")
    feat = np.loadtxt(out_dir / "flux_feat.txt")

    rep.check(
        "TEST A monotonic outflow: finite positive flux",
        np.all(np.isfinite(mono)) and np.all(mono > 0.0),
        f"flux={mono}",
    )
    rep.check(
        "TEST B reversing velocity: finite positive flux",
        np.all(np.isfinite(revers)) and np.all(revers > 0.0),
        f"flux={revers}",
    )
    rep.check(
        "TEST C negative velocity: finite positive flux",
        np.all(np.isfinite(negv)) and np.all(negv > 0.0),
        f"flux={negv}",
    )
    rep.check(
        "TEST B vs A: reversal changes the mapping result",
        np.any(np.abs(revers - mono) > 0.0),
        f"revers={revers} mono={mono}",
    )
    rep.check(
        "TEST H relativistic run: finite positive flux",
        np.all(np.isfinite(rel)) and np.all(rel > 0.0),
        f"flux={rel}",
    )
    rep.check(
        "TEST H non-relativistic limit: DO_REL at V=1 m/s matches "
        "DO_REL=F (NU_ON_NU0^3 ~ 1e-8)",
        np.allclose(rel_lo, norel_lo, rtol=1.0e-6, atol=0.0),
        f"rel={rel_lo} norel={norel_lo} rel_diff={np.abs(rel_lo - norel_lo)}",
    )
    rep.check(
        "TEST F NUM_BANDS=1 vs 2: identical flux on the shared "
        "observer grid",
        np.array_equal(b1, b2[: b1.size]),
        f"max diff={np.abs(b1 - b2[: b1.size]).max()}",
    )
    rep.check(
        "convergence: extending FREQ_CMF outside the required reach "
        "leaves the spectrum bitwise unchanged",
        np.array_equal(base, ext),
        f"base={base} ext={ext}",
    )
    rep.check(
        "feature inside the Doppler reach changes the flux",
        np.any(np.abs(feat - mono) > 0.0),
        f"feat={feat} mono={mono}",
    )

    rc, out = run_expect_abort(binary, out_dir, "under_blue")
    rep.check(
        "TEST E blue global undercoverage stops the run (bare STOP "
        "returns rc=0, so check the run never completed)",
        "UNEXPECTED" not in out,
        f"rc={rc} tail={out[-300:]}",
    )
    rep.check(
        "TEST E blue abort is named GLOBAL_CMF_GRID_UNDERCOVERAGE",
        "GLOBAL_CMF_GRID_UNDERCOVERAGE" in out,
        out[-600:],
    )
    rep.check(
        "blue diagnostic advises widening MAX_CF",
        "widen MAX_CF" in out,
        out[-600:],
    )

    rc, out = run_expect_abort(binary, out_dir, "under_red")
    rep.check(
        "TEST E red global undercoverage stops the run",
        "UNEXPECTED" not in out,
        f"rc={rc} tail={out[-300:]}",
    )
    rep.check(
        "TEST E red abort is named GLOBAL_CMF_GRID_UNDERCOVERAGE",
        "GLOBAL_CMF_GRID_UNDERCOVERAGE" in out,
        out[-600:],
    )
    rep.check(
        "red diagnostic advises widening MIN_CF",
        "widen MIN_CF" in out,
        out[-600:],
    )

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    x = np.arange(mono.size)
    ax.plot(x, mono, "o-", label="monotonic outflow")
    ax.plot(x, revers, "s-", label="reversing velocity")
    ax.plot(x, negv, "^-", label="negative velocity")
    ax.plot(x, feat, "d--", label="CHI feature in bracket")
    ax.set_xlabel("observer frequency index")
    ax.set_ylabel("formal flux")
    ax.legend()
    fig.savefig(PLOTDIR / "obs_frame_cov.png", dpi=160)
    plt.close(fig)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
