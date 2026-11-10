#!/usr/bin/env python3
"""Regression test for resolved Sobolev profile deposition in OBS.

The driver calls ADD_SOB_LINE_PROF from obs/ADD_SOB_LINE_PROF.INC directly.
The assertions cover profile normalization on a nonuniform CMF grid,
per-Hz units, maser/blend-cell drops, and skipping non-SOB or off-grid
lines.
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
NCF = 81
MAX_SIM = 3
N_LINES = 3
C_KMS = 299792.458
F0 = 0.5
VDOP = 50.0

MODULE_SOURCES = [TREE / "tools" / "set_kind_module.f"]


def nu_grid():
    k = np.arange(NCF, dtype=float)
    nu = 0.50035 - 0.0007 * k / (NCF - 1)
    return nu + 2.0e-7 * np.sin(2.0 * np.pi * k / (NCF - 1))


def cell_widths(nu):
    dn = np.empty_like(nu)
    dn[0] = nu[0] - nu[1]
    dn[-1] = nu[-2] - nu[-1]
    dn[1:-1] = 0.5 * (nu[:-2] - nu[2:])
    return dn


def main():
    rep = Report("SOB OBS profile deposition")
    out_dir = WORKDIR / "sob_prof"
    blank_dir = WORKDIR / "sob_prof_blank"
    out_dir.mkdir(parents=True, exist_ok=True)
    blank_dir.mkdir(parents=True, exist_ok=True)

    binary = compile_driver(
        "drv_sob_prof",
        MODULE_SOURCES + [F90DIR / "drv_sob_prof.f"],
        extra_flags=[f"-I{TREE / 'obs'}"],
    )
    result = run(binary, out_dir)
    run(binary, blank_dir, 1)

    eta = np.loadtxt(out_dir / "eta_out.txt", dtype=float)
    chi = np.loadtxt(out_dir / "chi_out.txt", dtype=float)
    diag = np.loadtxt(out_dir / "diag.txt", dtype=float)
    nu = nu_grid()
    dn = cell_widths(nu)
    dnu_d = F0 * VDOP / C_KMS
    support = np.abs(nu - F0) <= 3.0 * dnu_d + 1.0e-14
    weights = np.where(support, np.exp(-(((nu - F0) / dnu_d) ** 2)), 0.0)
    weights /= np.sum(weights * dn)
    prof = 1.0e-15 * weights

    chil = np.array([2.0e-4, 3.0e-4, -1.0e-4, 4.0e-4])
    etal = np.array([0.2, 0.3, 0.4, 0.5])
    clump = np.array([1.0, 2.0, 0.5, 1.5])
    valid = np.array([True, True, False, False])

    eta0 = 0.0
    chi0 = 0.0
    eta_delta = eta - eta0
    chi_delta = chi - chi0

    rep.check(
        "driver completes and all deposited values finite",
        result.returncode == 0
        and np.all(np.isfinite(eta))
        and np.all(np.isfinite(chi)),
    )

    expected_chi = np.zeros((ND, NCF))
    expected_eta = np.zeros((ND, NCF))
    for i in range(ND):
        if valid[i]:
            expected_chi[i, :] = chil[i] * clump[i] * prof
            expected_eta[i, :] = etal[i] * clump[i] * prof

    x_dop = (nu - F0) / dnu_d
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.plot(x_dop, chi_delta[0, :] / (chil[0] * clump[0]), label="depth 1")
    ax.plot(x_dop, chi_delta[1, :] / (chil[1] * clump[1]), label="depth 2")
    ax.plot(x_dop, prof, "k--", lw=1, label="normalized Gaussian")
    ax.set_xlim(-3.2, 3.2)
    ax.set_xlabel(r"$(\nu-\nu_0)/\Delta\nu_D$")
    ax.set_ylabel(r"deposited $\chi_\nu$ / line opacity")
    ax.legend()
    fig.savefig(PLOTDIR / "sob_flux_obs.png", dpi=160)
    plt.close(fig)

    rep.check(
        "deposited opacity and emissivity match normalized profile",
        np.allclose(chi_delta, expected_chi, rtol=2.0e-12, atol=0.0)
        and np.allclose(eta_delta, expected_eta, rtol=2.0e-12, atol=0.0),
        f"max chi err={np.abs(chi_delta-expected_chi).max():.3e}",
    )

    integrals = np.sum(chi_delta * dn[None, :], axis=1)
    expected_integrals = np.where(valid, 1.0e-15 * chil * clump, 0.0)
    rep.check(
        "line opacity integrates to CHIL over frequency",
        np.allclose(integrals, expected_integrals, rtol=2.0e-12, atol=0.0),
        "integrals=" + ",".join(f"{x:.6e}" for x in integrals),
    )

    rep.check(
        "maser and blend-repaired cells inject nothing",
        np.all(chi_delta[2:, :] == 0.0) and np.all(eta_delta[2:, :] == 0.0),
    )

    rep.check(
        "off-grid SOB and BLANK lines leave arrays unchanged",
        np.all(chi_delta[:, ~support] == 0.0)
        and np.all(eta_delta[:, ~support] == 0.0),
        f"support cells={int(support.sum())}",
    )

    n_drop = int(round(diag[0]))
    drop_sum = diag[1]
    inj_sum = diag[2]
    rep.check(
        "drop and injection diagnostics are finite and correct",
        n_drop == 2
        and np.isclose(drop_sum, 0.4 * 0.5 + 0.5 * 1.5)
        and np.isclose(inj_sum, 0.2 * 1.0 + 0.3 * 2.0),
        f"n_drop={n_drop} drop={drop_sum:.6e} inj={inj_sum:.6e}",
    )

    rep.check(
        "deposited line opacity remains positive",
        np.all(chi_delta >= 0.0) and np.all(eta_delta >= 0.0),
    )

    eta_blank = np.loadtxt(blank_dir / "eta_out.txt", dtype=float)
    chi_blank = np.loadtxt(blank_dir / "chi_out.txt", dtype=float)
    diag_blank = np.loadtxt(blank_dir / "diag.txt", dtype=float)
    rep.check(
        "zero-SOB mode is an exact no-op",
        np.all(eta_blank == 0.0)
        and np.all(chi_blank == 0.0)
        and np.all(diag_blank == 0.0),
    )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
