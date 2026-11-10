#!/usr/bin/env python3
"""Regression test for the post-deposition Sobolev ES update.

The production SOB_LINE_ES_UPDATE routine is exercised with controlled
formal-solver and convolution stubs.  This isolates the bookkeeping that
must remove the old scattering source once, solve on the line-bearing
opacity, and install the updated J/J_ES source once.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import F90DIR, PLOTDIR, TREE, WORKDIR, Report, compile_driver, run

ND = 4
NCF = 5

SOURCES = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "subs" / "extend3opac.f",
    TREE / "subs" / "ungrid.f",
    TREE / "obs" / "sob_line_es_update.f",
    F90DIR / "stub_sob_es_fg.f",
    F90DIR / "stub_sob_es_conv.f",
    F90DIR / "drv_sob_line_es.f",
]

CONV_SOURCES = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "subs" / "extend3opac.f",
    TREE / "subs" / "ungrid.f",
    TREE / "obs" / "sob_line_es_update.f",
    TREE / "subs" / "comp_j_conv_v2.f",
    TREE / "subs" / "thomas_rh.f",
    TREE / "subs" / "mon_interp.f",
    TREE / "subs" / "write_direct_info_v3.f",
    TREE / "tools" / "chk_little_endian.f",
    TREE / "unix" / "term_io.f",
    TREE / "unix" / "dir_acc_pars_gen.f",
    F90DIR / "stub_sob_es_fg.f",
    F90DIR / "drv_sob_line_es.f",
]


def base_arrays():
    i = np.arange(1, ND + 1)[:, None]
    k = np.arange(1, NCF + 1)[None, :]
    base = 10.0 + 0.1 * i + 0.01 * k
    chi = 2.0 + 0.01 * i + 0.0 * k
    rj = 3.0 + 0.001 * i + 0.0 * k
    rj_es = 3.5 + 0.001 * i + 0.0001 * k
    ray = 0.02 + 0.001 * i + 0.0 * k
    dust = 0.03 + 0.001 * i + 0.0 * k
    esec = np.full(ND, 0.1)[:, None]
    return base, chi, rj, rj_es, ray, dust, esec


def main():
    rep = Report("SOB post-deposition electron scattering")
    binary = compile_driver("drv_sob_line_es", SOURCES)

    mode0_dir = WORKDIR / "sob_line_es_mode0"
    mode1_dir = WORKDIR / "sob_line_es_mode1"
    mode3_dir = WORKDIR / "sob_line_es_mode3"
    mode0_dir.mkdir(parents=True, exist_ok=True)
    mode1_dir.mkdir(parents=True, exist_ok=True)
    mode3_dir.mkdir(parents=True, exist_ok=True)
    run(binary, mode0_dir, 0)
    run(binary, mode1_dir, 1)
    run(binary, mode3_dir, 3)

    conv_binary = compile_driver("drv_sob_line_es_conv", CONV_SOURCES)
    conv_dir = WORKDIR / "sob_line_es_conv"
    conv_dir.mkdir(parents=True, exist_ok=True)
    run(conv_binary, conv_dir, 2)

    base, _, rj_old, _, ray, dust, esec = base_arrays()
    eta0 = np.loadtxt(mode0_dir / "eta_es.txt")
    rj0 = np.loadtxt(mode0_dir / "rj_es.txt")
    rjes0 = np.loadtxt(mode0_dir / "rj_es_conv.txt")
    expected0 = base + esec * rj_old + (ray + dust) * rj_old

    rep.check(
        "coherent unchanged-J pass restores ETA exactly once",
        np.allclose(eta0, expected0, rtol=1.0e-12, atol=1.0e-14),
        f"max eta err={np.abs(eta0-expected0).max():.3e}",
    )
    rep.check(
        "coherent mode leaves J and J_ES on the formal solution",
        np.allclose(rj0, rj_old, rtol=1.0e-12, atol=0.0)
        and np.allclose(rjes0, rj_old, rtol=1.0e-12, atol=0.0),
    )

    eta1 = np.loadtxt(mode1_dir / "eta_es.txt")
    rj1 = np.loadtxt(mode1_dir / "rj_es.txt")
    rjes1 = np.loadtxt(mode1_dir / "rj_es_conv.txt")
    expected1 = base + esec * 4.25 + (ray + dust) * 4.0

    rep.check(
        "noncoherent pass installs updated J_ES once",
        np.allclose(eta1, expected1, rtol=1.0e-12, atol=1.0e-14),
        f"max eta err={np.abs(eta1-expected1).max():.3e}",
    )
    rep.check(
        "noncoherent pass returns updated J and convolved J_ES",
        np.allclose(rj1, 4.0, rtol=0.0, atol=1.0e-14)
        and np.allclose(rjes1, 4.25, rtol=0.0, atol=1.0e-14),
    )

    diag1 = np.loadtxt(mode1_dir / "ib_out.txt")
    nu = 5.0 - 0.25 * np.arange(NCF)
    emh = np.exp(-nu)
    dbb = 0.5 * nu**4 * emh / (1.0 - emh) ** 2
    rep.check(
        "post-update diffusion boundary uses the new J and opacity",
        np.allclose(diag1[0], 4.0, rtol=0.0, atol=1.0e-14)
        and np.allclose(diag1[1], dbb / (3.0 * 2.04), rtol=1.0e-12),
        f"IB_J={diag1[0]} IB_H={diag1[1]}",
    )

    eta3 = np.loadtxt(mode3_dir / "eta_es.txt")
    rj3 = np.loadtxt(mode3_dir / "rj_es.txt")
    rjes3 = np.loadtxt(mode3_dir / "rj_es_conv.txt")
    rep.check(
        "INC_GRID=F branch gives the same source update",
        np.allclose(eta3, expected1, rtol=1.0e-12, atol=1.0e-14)
        and np.allclose(rj3, 4.0, rtol=0.0, atol=1.0e-14)
        and np.allclose(rjes3, 4.25, rtol=0.0, atol=1.0e-14),
    )

    rj_conv = np.loadtxt(conv_dir / "rj_es.txt")
    rjes_conv = np.loadtxt(conv_dir / "rj_es_conv.txt")
    rep.check(
        "production convolution returns finite positive J_ES",
        np.all(np.isfinite(rjes_conv)) and np.all(rjes_conv > 0.0),
    )
    rep.check(
        "production convolution redistributes the localized J peak",
        np.max(rjes_conv) < np.max(rj_conv)
        and np.any(rjes_conv[:, :-1] > rj_conv[:, :-1]),
        f"J={rj_conv[0]} J_ES={rjes_conv[0]}",
    )
    rep.check(
        "all post-update emissivities remain finite and positive",
        np.all(np.isfinite(eta1)) and np.all(eta1 > 0.0),
    )

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.plot(base[:, 0], label="fixed emissivity")
    ax.plot(expected0[:, 0], label="old scattering source")
    ax.plot(eta1[:, 0], label="updated source")
    ax.set_xlabel("depth index")
    ax.set_ylabel("stored emissivity")
    ax.legend()
    fig.savefig(PLOTDIR / "sob_line_es.png", dpi=160)
    plt.close(fig)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
