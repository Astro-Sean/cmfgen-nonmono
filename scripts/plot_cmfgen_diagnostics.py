#!/usr/bin/env python3
"""Diagnostic plots for a completed CMFGEN population-solver run.

Reads the files CMFGEN leaves in a run directory after the FINAL /
SINGLE_FINAL job completes and produces a small set of sanity-check
figures in <run_dir>/plots/:

  diag_structure.png    -- V, density, T, Ne: converged model vs hydro input
  diag_regrid.png       -- regridded profile (grid points as markers) with the
                           depth index on the top axis
  diag_convergence.png  -- max population corrections and luminosity vs iteration
  diag_opacity.png      -- mean opacities, radioactive heating, clumping
  diag_radiation.png    -- J/H/K moments and Eddington factor f = K/J
  diag_ionization.png   -- ionization fraction vs radius for each element
  diag_departure.png    -- ground-state departure coefficients (H I, He I)

Data sources (all in the run directory):
  RVTJ          converged atmosphere: R, V, T, ED, opacities, J/H/K, DHYD/DHeI
  SN_HYDRO_DATA hydro input that the model was built from
  POP*          per-element densities and per-ion level populations
  OUTGEN        iteration log (corrections, luminosities)
  WARNINGS      CMFGEN warnings file

Each figure is optional: if its source files are missing the script logs a
note and continues, so a partial run still yields the plots that are
possible.  Exit status is non-zero only if nothing could be plotted.

Usage:
  python3 scripts/plot_cmfgen_diagnostics.py <cmfgen_dir> [--outdir DIR]
"""

import argparse
import re
import sys
from pathlib import Path

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as e:
    print(f"ERROR: Missing required Python package: {e}", file=sys.stderr)
    sys.exit(1)

# Reuse the shared convergence assessment so the classification reported here
# matches check_cmfgen.py / cmfgen_convergence.py.
sys.path.insert(0, str(Path(__file__).parent))
try:
    from cmfgen_convergence import (
        ConvergenceHistory,
        IterationRecord,
        assess_convergence,
    )
    _HAVE_CONV = True
except Exception:
    _HAVE_CONV = False

FLOAT_RE = re.compile(r"^[+-]?[\d.]+[eEdD]?[+-]?\d*$")
LEVELS_RE = re.compile(r"Number of\s+(\w+)\s+levels:\s*(\d+)")

# OUTGEN line patterns (looser than cmfgen_convergence: solver names vary,
# e.g. SOLVEBA_V13 vs SOLVEBA_NONMONO_V1).
ITER_RE = re.compile(r"Current great iteration count is\s+(\d+)")
LUM_RE = re.compile(
    r"Luminosity of star \(d=1,ND\)\(iteration\s+(\d+)\)\s+is:\s+"
    r"([\d.eE+-]+)\s+([\d.eE+-]+)"
)
INC_RE = re.compile(
    r"Maximum % increase at depth\s+(\d+)\s+is\s+([\d.eE+-]+)"
    r".*iteration\s+(\d+)"
)
DEC_RE = re.compile(
    r"Maximum % decrease at depth\s+(\d+)\s+is\s+([\d.eE+-]+)"
    r".*iteration\s+(\d+)"
)
MAXCH_RE = re.compile(r"Maximm changes as returned by\s+\S+\s+is\s+([\d.eE+-]+)")

# Only the last part of OUTGEN is needed once the file gets large; the
# rotation in submit_cmfgen.sh caps it at ~2 GB per cycle anyway.
OUTGEN_MAX_BYTES = 400 * 1024 * 1024

# POP* file -> element label for the ionization figure.
POP_FILES = [
    ("POPHYD", "H"),
    ("POPHE", "He"),
    ("POPOXY", "O"),
    ("POPSIL", "Si"),
    ("POPCAL", "Ca"),
    ("POPIRON", "Fe"),
]

ION_LABELS = {
    "HI": "H I", "HeI": "He I", "OI": "O I", "O2": "O II",
    "SkI": "Si I", "Sk2": "Si II", "CaI": "Ca I", "Ca2": "Ca II",
    "FeI": "Fe I", "Fe2": "Fe II", "FeIII": "Fe III",
}


def _is_float(tok):
    return bool(FLOAT_RE.match(tok))


def _to_float(tok):
    return float(tok.replace("D", "E").replace("d", "e"))


def _section_values(lines, start):
    """Collect whitespace-separated floats from lines[start:] until a
    non-numeric line.  Returns (values, index of first non-numeric line)."""
    vals = []
    i = start
    while i < len(lines):
        toks = lines[i].split()
        if not toks or not all(_is_float(t) for t in toks):
            break
        vals.extend(_to_float(t) for t in toks)
        i += 1
    return vals, i


def read_labeled_sections(path):
    """Parse a CMFGEN labeled-section file (RVTJ, SN_HYDRO_DATA).

    Layout: 'key: value' header lines, then repeated blocks of
    <label line> <N numeric lines>.  Returns (header dict, sections dict).
    """
    lines = Path(path).read_text(errors="replace").splitlines()
    header = {}
    sections = {}
    i = 0
    # Header: 'Key: value' lines up to the first bare label.
    while i < len(lines):
        toks = lines[i].split()
        if ":" in lines[i]:
            key, _, val = lines[i].partition(":")
            header[key.strip()] = val.strip()
            i += 1
            continue
        if toks and not all(_is_float(t) for t in toks):
            break
        i += 1
    # Sections.
    while i < len(lines):
        toks = lines[i].split()
        if toks and not all(_is_float(t) for t in toks):
            label = lines[i].strip()
            vals, j = _section_values(lines, i + 1)
            sections[label] = np.array(vals, dtype=np.float64)
            i = j
        else:
            i += 1
    return header, sections


def read_pop_file(path):
    """Parse a POP* file: element density block followed by per-ion
    level-population blocks stored depth-major (ND x NLEV per ion)."""
    lines = Path(path).read_text(errors="replace").splitlines()
    nd = None
    for line in lines[:10]:
        if line.strip().startswith("ND:"):
            nd = int(line.split(":")[1])
            break
    if nd is None:
        raise ValueError(f"{path}: no ND header")

    # First numeric block after the header = element density (ND values).
    i = 0
    while i < len(lines):
        toks = lines[i].split()
        if toks and all(_is_float(t) for t in toks):
            break
        i += 1
    vals, i = _section_values(lines, i)
    elem = np.array(vals[:nd], dtype=np.float64)

    ions = {}
    while i < len(lines):
        m = LEVELS_RE.search(lines[i])
        if not m:
            i += 1
            continue
        name, nlev = m.group(1), int(m.group(2))
        i += 1
        # Skip the oscillator-date line (a non-numeric line between the
        # levels header and the population data).
        if i < len(lines):
            toks = lines[i].split()
            if toks and not all(_is_float(t) for t in toks):
                i += 1
        vals, i = _section_values(lines, i)
        n = min(len(vals), nd * nlev)
        if n < nd * nlev:
            print(f"  WARNING: {path.name} {name}: expected "
                  f"{nd * nlev} values, got {len(vals)} (truncated block)")
        pop = np.full((nd, nlev), np.nan)
        flat = np.array(vals[:n])
        pop.flat[:n] = flat
        ions[name] = pop
    return nd, elem, ions


def parse_outgen(path):
    """Stream OUTGEN line by line and build per-iteration records.

    A full regex pass over a multi-GB OUTGEN is slow and memory-heavy, so
    for very large files only the tail is scanned.  The iteration history
    of interest (the approach to convergence) is always at the end.
    """
    path = Path(path)
    size = path.stat().st_size
    truncated = size > OUTGEN_MAX_BYTES
    records = {}

    def handle(line):
        m = LUM_RE.search(line)
        if m:
            it = int(m.group(1))
            r = records.setdefault(it, {})
            r["lum_inner"] = float(m.group(2))
            r["lum_outer"] = float(m.group(3))
            return
        m = INC_RE.search(line)
        if m:
            it = int(m.group(3))
            r = records.setdefault(it, {})
            r["max_inc"] = float(m.group(2))
            r["max_inc_depth"] = int(m.group(1))
            handle.last_it = it
            return
        m = DEC_RE.search(line)
        if m:
            it = int(m.group(3))
            r = records.setdefault(it, {})
            r["max_dec"] = float(m.group(2))
            r["max_dec_depth"] = int(m.group(1))
            handle.last_it = it
            return
        m = MAXCH_RE.search(line)
        if m:
            # MAXCH has no iteration number; it follows the inc/dec block.
            it = getattr(handle, "last_it", None)
            if it is not None:
                records.setdefault(it, {})["maxch"] = float(m.group(1))
            return
        m = ITER_RE.search(line)
        if m:
            it = int(m.group(1))
            records.setdefault(it, {})
            handle.last_it = it

    handle.last_it = None

    if truncated:
        with open(path, "rb") as f:
            f.seek(size - OUTGEN_MAX_BYTES)
            data = f.read().decode(errors="replace")
        # Drop the first partial line.
        for line in data.splitlines()[1:]:
            handle(line)
    else:
        with open(path, errors="replace") as f:
            for line in f:
                handle(line)

    iters = []
    for it in sorted(records):
        r = records[it]
        iters.append(
            dict(
                iteration=it,
                max_inc=r.get("max_inc", np.nan),
                max_inc_depth=r.get("max_inc_depth", -1),
                max_dec=r.get("max_dec", np.nan),
                max_dec_depth=r.get("max_dec_depth", -1),
                maxch=r.get("maxch", np.nan),
                lum_inner=r.get("lum_inner", np.nan),
                lum_outer=r.get("lum_outer", np.nan),
            )
        )
    return iters, truncated


def _axis_label_r(ax):
    ax.set_xscale("log")
    ax.set_xlabel(r"Radius ($10^{10}$ cm)  [outer layers at right]")


def _interp_to(r_new, r_old, y_old):
    """Log-log interpolate y(r_old) onto r_new; both grids decrease."""
    xo, yo = r_old[::-1], y_old[::-1]
    xn = r_new[::-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.interp(np.log(xn), np.log(xo), np.log(yo))
    return np.exp(out)[::-1]


def plot_structure(sec_rvtj, sec_hydro, outdir, label):
    need = ["Radius (10^10 cm)", "Velocity (km/s)", "Temperature (10^4K)",
            "Electron density", "Mass Density (gm/cm^3)"]
    if not all(k in sec_rvtj for k in need):
        print("  skip structure: RVTJ sections missing")
        return None
    r = sec_rvtj["Radius (10^10 cm)"]
    v = sec_rvtj["Velocity (km/s)"]
    t = sec_rvtj["Temperature (10^4K)"] * 1e4
    ne = sec_rvtj["Electron density"]
    rho = sec_rvtj["Mass Density (gm/cm^3)"]

    r_h = sec_hydro.get("Radius grid (10^10cm)")
    v_h = sec_hydro.get("Velocity (km/s)")
    t_h = sec_hydro.get("Temperature (10^4 K)")
    rho_h = sec_hydro.get("Density (gm/cm^3)")
    ne_h = sec_hydro.get("Electron density (/cm^3)")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0, 0]
    ax.plot(r, v, "r-", lw=1.5, label="CMFGEN")
    if r_h is not None and v_h is not None:
        ax.plot(r_h, v_h, "b--", lw=1.2, label="Hydro input")
        ax.legend(fontsize=8)
    _axis_label_r(ax)
    ax.set_ylabel("Velocity (km/s)")
    ax.set_title("Velocity")

    ax = axes[0, 1]
    ax.loglog(r, rho, "r-", lw=1.5, label="CMFGEN")
    if r_h is not None and rho_h is not None:
        ax.loglog(r_h, rho_h, "b--", lw=1.2, label="Hydro input")
        ax.legend(fontsize=8)
    _axis_label_r(ax)
    ax.set_ylabel(r"Density (g cm$^{-3}$)")
    ax.set_title("Mass density")

    ax = axes[1, 0]
    ax.loglog(r, t, "r-", lw=1.5, label="CMFGEN T")
    tg = sec_rvtj.get("Grey temperature (10^4K)")
    if tg is not None and len(tg) == len(r):
        ax.loglog(r, tg * 1e4, "g:", lw=1.2, label="Grey T")
    if r_h is not None and t_h is not None:
        ax.loglog(r_h, t_h * 1e4, "b--", lw=1.2, label="Hydro input T")
    ax.legend(fontsize=8)
    _axis_label_r(ax)
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Temperature")

    ax = axes[1, 1]
    ax.loglog(r, ne, "r-", lw=1.5, label="CMFGEN")
    if r_h is not None and ne_h is not None:
        ax.loglog(r_h, ne_h, "b--", lw=1.2, label="Hydro input")
        ax.legend(fontsize=8)
    _axis_label_r(ax)
    ax.set_ylabel(r"$N_e$ (cm$^{-3}$)")
    ax.set_title("Electron density")

    fig.suptitle(f"{label}: converged structure vs hydro input", fontsize=13)
    fig.tight_layout()
    out = outdir / "diag_structure.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def _add_depth_axis(ax, r_cm):
    """Twin x-axis showing the depth index.

    r_cm is the model radius in cm, ordered outer -> inner (decreasing).
    The lower axis is log(R) increasing left to right, so depth 1 (outer
    boundary) sits at the right and ND (inner boundary) at the left.
    Ticks are placed at the actual radius of each labeled depth, so the
    axis is a true 1:1 lookup -- tick crowding reflects grid packing.
    """
    nd = len(r_cm)
    sec = ax.twiny()
    sec.set_xscale("log")
    sec.set_xlim(ax.get_xlim())

    # All depths get a tick.  Labels are placed at evenly spaced
    # target radii snapped to the nearest real depth, so they never
    # overlap regardless of grid packing.  Interior labels prefer
    # multiples of 5 when the rounded index keeps half the target
    # spacing; both endpoints are always shown.
    lgr = np.log10(r_cm)
    n_lab = max(2, min(12, nd))
    targets = np.linspace(lgr[0], lgr[-1], n_lab)
    labeled = []
    for k, t in enumerate(targets):
        d = int(np.argmin(np.abs(lgr - t))) + 1
        if k == 0:
            d = 1
        elif k == n_lab - 1:
            d = nd
        elif labeled:
            d5 = int(round(d / 5.0)) * 5
            if (
                1 < d5 < nd
                and abs(lgr[d5 - 1] - lgr[labeled[-1] - 1])
                >= 0.5 * abs(lgr[-1] - lgr[0]) / n_lab
            ):
                d = d5
        if d not in labeled:
            labeled.append(d)
    labeled = np.array(sorted(labeled))

    minor = np.setdiff1d(np.arange(1, nd + 1), labeled)
    sec.set_xticks(r_cm[labeled - 1])
    sec.set_xticklabels([str(d) for d in labeled])
    sec.set_xticks(r_cm[minor - 1], minor=True)
    sec.xaxis.set_minor_formatter(plt.NullFormatter())
    sec.set_xlabel("Depth index (1 = outer)", fontsize=9)
    sec.tick_params(which="major", labelsize=7)
    sec.tick_params(which="minor", length=2)
    return sec


def _read_new_sn_r_grid(path):
    """Parse NEW_SN_R_GRID (I, R[10^10 cm], V[km/s], SIGMA, ...).

    Written by the structure readers at grid construction, so it carries
    the true (unclamped) SIGMA; RVTJ stores the working sigma, which the
    non-monotonic solver clamps at -0.999.
    """
    if not path.exists():
        return None
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        toks = line.split()
        if len(toks) >= 4 and all(_is_float(t) for t in toks[:4]):
            row = [_to_float(t) for t in toks[:7]]
            row += [np.nan] * (7 - len(row))
            rows.append(row)
    if not rows:
        return None
    a = np.asarray(rows, dtype=np.float64)
    out = {"r": a[:, 1] * 1e10, "v": a[:, 2], "sigma": a[:, 3]}
    if np.isfinite(a[:, 4]).all() and np.isfinite(a[:-1, 6]).all():
        # Diagnostic tau and per-cell log10 increment written by
        # RV_SN_MODEL_SNHYD / RV_SN_MODEL_SNIIN; the increment row is
        # blank on the last depth.
        out["tau"] = a[:, 4]
        out["dlogtau"] = a[:-1, 6]
    return out


def plot_regrid(cdir, sec_rvtj, sec_hydro, outdir, label):
    """The grid CMFGEN actually solved on (RVTJ) vs the hydro input.

    Markers are at every depth so the regridded point distribution is
    visible; the top axis labels depth index for each panel.  Velocity
    and sigma prefer NEW_SN_R_GRID when present, since it holds the
    unclamped sigma the structure reader produced.
    """
    need = ["Radius (10^10 cm)", "Velocity (km/s)", "dlnV/dlnr-1",
            "Mass Density (gm/cm^3)", "Temperature (10^4K)"]
    if not all(k in sec_rvtj for k in need):
        print("  skip regrid: RVTJ sections missing")
        return None

    r = sec_rvtj["Radius (10^10 cm)"] * 1e10
    v = sec_rvtj["Velocity (km/s)"]
    sig = sec_rvtj["dlnV/dlnr-1"]
    rho = sec_rvtj["Mass Density (gm/cm^3)"]
    t = sec_rvtj["Temperature (10^4K)"] * 1e4

    grid = _read_new_sn_r_grid(cdir / "NEW_SN_R_GRID")
    if grid is not None and len(grid["r"]) == len(r):
        r_g, v_g, sig_g = grid["r"], grid["v"], grid["sigma"]
    else:
        r_g, v_g, sig_g = r, v, sig

    r_h = sec_hydro.get("Radius grid (10^10cm)")
    if r_h is not None:
        r_h = r_h * 1e10
    v_h = sec_hydro.get("Velocity (km/s)")
    sig_h = sec_hydro.get("Sigma (dlnV/dlnr-1)")
    rho_h = sec_hydro.get("Density (gm/cm^3)")
    t_h = sec_hydro.get("Temperature (10^4 K)")

    # Rosseland tau on the hydro grid (trapezoid inward from the edge),
    # used for the tau overlay and the tau=2/3 photosphere marker.
    tau_h = None
    r_phot = None
    if (
        r_h is not None
        and rho_h is not None
        and "Kappa (cm^2/gm)" in sec_hydro
    ):
        kap_h = sec_hydro["Kappa (cm^2/gm)"]
        tau_h = np.concatenate(
            ([0.0],
             np.cumsum(0.5 * (kap_h[:-1] * rho_h[:-1]
                              + kap_h[1:] * rho_h[1:])
                       * np.diff(-r_h))))
        if tau_h[-1] > 2.0 / 3.0:
            r_phot = np.exp(
                np.interp(np.log(2.0 / 3.0), np.log(tau_h[1:]),
                          np.log(r_h[1:])))

    fig, axes = plt.subplots(2, 3, figsize=(20, 10))

    def panel(ax, x, y, y_h, ylabel, title, logy=False):
        show = ax.loglog if logy else ax.semilogx
        if r_h is not None and y_h is not None:
            show(r_h, y_h, "b--", lw=1.0, alpha=0.7, label="Hydro input")
        show(x, y, "r-o", ms=3.5, lw=1.0, label="CMFGEN grid")
        if r_phot is not None:
            ax.axvline(r_phot, color="g", ls="-.", lw=1.2,
                       label=r"$\tau_{Ross}=2/3$")
        ax.legend(fontsize=8)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.25, which="both")
        _add_depth_axis(ax, r)

    panel(axes[0, 0], r_g, v_g, v_h, "Velocity (km/s)", "Velocity", logy=True)
    panel(axes[0, 1], r_g, sig_g, sig_h, r"$d\ln V/d\ln R - 1$", "Sigma")
    panel(axes[1, 0], r, rho, rho_h, r"Density (g cm$^{-3}$)",
          "Mass density", logy=True)
    panel(axes[1, 1], r, t, (t_h * 1e4) if t_h is not None else None,
          "Temperature (K)", "Temperature", logy=True)

    # Optical-depth panels: cumulative tau shows where the points sit in
    # the transfer coordinate; the per-cell increment should be flat in
    # the interior (uniform log-tau grid) and is expected to differ in a
    # bounded outer zone, which for VEL_LAW=14 is deliberately log-R.
    if grid is not None and "tau" in grid:
        tau_g = grid["tau"]
        axes[0, 2].loglog(r_g, tau_g, "r-o", ms=3.5, lw=1.0,
                          label="CMFGEN grid")
        if tau_h is not None:
            axes[0, 2].loglog(r_h, tau_h, "b--", lw=1.0, alpha=0.7,
                              label="Hydro input")
        if r_phot is not None:
            axes[0, 2].axvline(r_phot, color="g", ls="-.", lw=1.2,
                               label=r"$\tau_{Ross}=2/3$")
        axes[0, 2].legend(fontsize=8)
        axes[0, 2].set_ylabel(r"$\tau$")
        axes[0, 2].set_title("Optical depth")
        axes[0, 2].grid(alpha=0.25, which="both")
        _add_depth_axis(axes[0, 2], r)

        axes[1, 2].semilogx(r_g[:-1], grid["dlogtau"], "r-o", ms=3.5,
                            lw=1.0)
        axes[1, 2].axhline(np.median(grid["dlogtau"]), color="k",
                           ls=":", lw=1.0, label="median spacing")
        if r_phot is not None:
            axes[1, 2].axvline(r_phot, color="g", ls="-.", lw=1.2,
                               label=r"$\tau_{Ross}=2/3$")
        axes[1, 2].legend(fontsize=8)
        axes[1, 2].set_ylabel(r"$\Delta\log_{10}\tau$ per cell")
        axes[1, 2].set_title("Tau spacing (flat = uniform log-tau grid)")
        axes[1, 2].grid(alpha=0.25, which="both")
        _add_depth_axis(axes[1, 2], r)
    else:
        axes[0, 2].axis("off")
        axes[1, 2].axis("off")

    for ax in axes.flat:
        ax.set_xlabel("Radius (cm)  [outer layers at right]")

    fig.suptitle(
        f"{label}: regridded profile (markers = grid points, ND={len(r)})",
        fontsize=13)
    fig.tight_layout()
    out = outdir / "diag_regrid.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def plot_convergence(iters, truncated, outdir, label):
    if not iters:
        print("  skip convergence: no iteration records in OUTGEN")
        return None
    it = np.array([r["iteration"] for r in iters])
    inc = np.array([r["max_inc"] for r in iters])
    dec = np.array([r["max_dec"] for r in iters])
    maxch = np.array([r["maxch"] for r in iters])
    lin = np.array([r["lum_inner"] for r in iters])
    lout = np.array([r["lum_outer"] for r in iters])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.semilogy(it, np.abs(inc), "r.-", ms=4, lw=0.8, label="max % increase")
    ax.semilogy(it, np.abs(dec), "b.-", ms=4, lw=0.8, label="max % decrease")
    if np.isfinite(maxch).any():
        ax.semilogy(it, np.abs(maxch), "k:", lw=1.2, label="MAXCH (SOLVEBA)")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Max population change (%)")
    ax.set_title("Population corrections")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    if np.isfinite(lin).any():
        ax.plot(it, lin, ".-", ms=4, lw=0.8, label="L inner (d=ND)")
    if np.isfinite(lout).any():
        ax.plot(it, lout, ".-", ms=4, lw=0.8, label="L outer (d=1)")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Luminosity (erg/s)")
    ax.set_title("Luminosity vs iteration")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    title = f"{label}: iteration history"
    if truncated:
        title += f"  (OUTGEN > {OUTGEN_MAX_BYTES // 2**20} MB, tail only)"
    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    out = outdir / "diag_convergence.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def plot_opacity(sec, outdir, label):
    r = sec.get("Radius (10^10 cm)")
    keys = ["Rosseland Mean Opacity", "Flux Mean Opacity",
            "Planck Mean Opacity", "Absorption Mean Opacity"]
    if r is None or not any(k in sec for k in keys):
        print("  skip opacity: RVTJ sections missing")
        return None
    heat = sec.get("Heating: radioactive decay (ergs/cm^3/s)")
    clump = sec.get("Clumping Factor")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    for k in keys:
        arr = sec.get(k)
        if arr is not None and len(arr) == len(r):
            ax.loglog(r, np.abs(arr), lw=1.3, label=k.replace(" Mean Opacity", ""))
    _axis_label_r(ax)
    ax.set_ylabel(r"Mean opacity (cm$^{-1}$)")
    ax.set_title("Mean opacities")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    plotted = False
    if heat is not None and len(heat) == len(r) and np.any(heat > 0):
        ax.loglog(r, heat, "r-", lw=1.3, label="Radioactive heating")
        ax.set_ylabel(r"Heating (erg cm$^{-3}$ s$^{-1}$)")
        plotted = True
    ax2 = None
    if clump is not None and len(clump) == len(r):
        ax2 = ax if not plotted else ax.twinx()
        ax2.semilogy(r, np.maximum(clump, 1e-30), "b--", lw=1.3,
                     label="Clumping factor")
        if plotted:
            ax2.set_ylabel("Clumping factor", color="b")
        else:
            ax2.set_ylabel("Clumping factor")
        plotted = True
    if plotted:
        _axis_label_r(ax)
        ax.set_title("Heating / clumping")
        h1, l1 = ax.get_legend_handles_labels()
        if ax2 is not None and ax2 is not ax:
            h2, l2 = ax2.get_legend_handles_labels()
            h1, l1 = h1 + h2, l1 + l2
        ax.legend(h1, l1, fontsize=8, loc="best")
        ax.grid(alpha=0.3)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "No heating/clumping data", ha="center",
                va="center", transform=ax.transAxes)

    fig.suptitle(f"{label}: opacities and energy input", fontsize=13)
    fig.tight_layout()
    out = outdir / "diag_opacity.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def plot_radiation(sec, outdir, label):
    r = sec.get("Radius (10^10 cm)")
    j = sec.get("J moment of radiation field")
    h = sec.get("H moment of radiation field")
    k = sec.get("K moment of radiation field")
    if r is None or j is None or len(j) != len(r):
        print("  skip radiation: RVTJ J-moment section missing")
        return None

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    for arr, name, c in ((j, "J", "k"), (h, "H", "r"), (k, "K", "b")):
        if arr is not None and len(arr) == len(r):
            ax.semilogy(r, np.abs(arr), c + "-", lw=1.3, label=name)
    _axis_label_r(ax)
    ax.set_ylabel("Moment")
    ax.set_title("Radiation-field moments")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    if k is not None and len(k) == len(r):
        with np.errstate(divide="ignore", invalid="ignore"):
            f = k / np.where(j == 0, np.nan, j)
        ax.plot(r, f, "k-", lw=1.3)
        ax.axhline(1.0 / 3.0, color="gray", ls="--", lw=1,
                   label="1/3 (diffusion limit)")
        ax.set_ylim(0.0, 1.05)
        ax.legend(fontsize=8)
        n_bad = int(np.count_nonzero((f < -0.05) | (f > 1.05)))
        ax.set_title(f"Eddington factor f = K/J ({n_bad} depths outside [0,1])")
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "No K moment", ha="center", va="center",
                transform=ax.transAxes)
    _axis_label_r(ax)
    ax.set_ylabel("f")
    ax.grid(alpha=0.3)

    fig.suptitle(f"{label}: radiation field", fontsize=13)
    fig.tight_layout()
    out = outdir / "diag_radiation.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def plot_ionization(cmfgen_dir, sec, outdir, label):
    r = sec.get("Radius (10^10 cm)")
    if r is None:
        print("  skip ionization: no RVTJ radius")
        return None
    nd = len(r)
    panels = []
    for fname, elem_label in POP_FILES:
        p = Path(cmfgen_dir) / fname
        if not p.exists():
            continue
        try:
            nd_p, elem, ions = read_pop_file(p)
        except Exception as e:
            print(f"  WARNING: could not parse {fname}: {e}")
            continue
        if nd_p != nd:
            print(f"  WARNING: {fname} ND={nd_p} != RVTJ ND={nd}; skipped")
            continue
        if not ions:
            continue
        panels.append((elem_label, elem, ions))
    if not panels:
        print("  skip ionization: no usable POP* files")
        return None

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True)
    for ax, (elem_label, elem, ions) in zip(axes.flat, panels):
        with np.errstate(divide="ignore", invalid="ignore"):
            denom = np.where(elem > 0, elem, np.nan)
            modeled = np.zeros(nd)
            for name, pop in ions.items():
                frac = np.nansum(pop, axis=1) / denom
                modeled += np.nan_to_num(frac)
                ax.plot(r, frac, lw=1.3,
                        label=ION_LABELS.get(name, name))
            rest = 1.0 - modeled
            if np.any(rest > 1e-3):
                ax.plot(r, rest, "k--", lw=1.0,
                        label="higher stages (not in model)")
        ax.set_ylim(-0.02, 1.05)
        ax.set_title(f"{elem_label} ionization fractions")
        ax.legend(fontsize=7, loc="best")
        ax.grid(alpha=0.3)
        ax.set_ylabel("Fraction")
        _axis_label_r(ax)
    for ax in axes.flat[len(panels):]:
        ax.axis("off")

    fig.suptitle(f"{label}: ionization structure from POP files", fontsize=13)
    fig.tight_layout()
    out = outdir / "diag_ionization.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def plot_departure(sec, outdir, label):
    r = sec.get("Radius (10^10 cm)")
    dh = sec.get("DHYD population")
    dhe = sec.get("DHeI population")
    have = [(n, a) for n, a in (("H I (DHYD)", dh), ("He I (DHeI)", dhe))
            if a is not None and r is not None and len(a) == len(r)]
    if not have:
        print("  skip departure: no DHYD/DHeI sections in RVTJ")
        return None
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, arr in have:
        ax.semilogy(r, np.abs(arr), lw=1.3, label=name + " ground-state b")
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    _axis_label_r(ax)
    ax.set_ylabel("Departure coefficient b")
    ax.set_title(f"{label}: ground-state departure coefficients")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = outdir / "diag_departure.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def print_summary(hdr, sec_rvtj, sec_hydro, iters, warnings_path):
    print("\n--- CMFGEN diagnostics summary ---")
    for k in ("ND", "NC", "NP", "NCF", "L(Lsun)", "Mdot(Msun/yr)",
              "H/He abundance", "Was T fixed?", "Species naming convention"):
        if k in hdr:
            print(f"  {k:<28} {hdr[k]}")
    r = sec_rvtj.get("Radius (10^10 cm)")
    rho = sec_rvtj.get("Mass Density (gm/cm^3)")
    if r is not None and rho is not None and len(r) == len(rho) and len(r) > 2:
        rr = r * 1e10
        y = 4.0 * np.pi * rr**2 * rho
        mass = abs(np.sum(0.5 * (y[:-1] + y[1:]) * np.diff(rr))) / 1.989e33
        print(f"  {'Ejecta mass (integ. RVTJ)':<28} {mass:.3f} Msun")
    rho_h = sec_hydro.get("Density (gm/cm^3)")
    r_h = sec_hydro.get("Radius grid (10^10cm)")
    if r_h is not None and rho_h is not None and len(r_h) == len(rho_h):
        rr = r_h * 1e10
        y = 4.0 * np.pi * rr**2 * rho_h
        mass = abs(np.sum(0.5 * (y[:-1] + y[1:]) * np.diff(rr))) / 1.989e33
        print(f"  {'Ejecta mass (hydro input)':<28} {mass:.3f} Msun")
    if iters:
        last = iters[-1]
        print(f"  {'Iterations':<28} {last['iteration']}")
        if np.isfinite(last["lum_inner"]) and np.isfinite(last["lum_outer"]):
            print(f"  {'L inner / L outer (last)':<28} "
                  f"{last['lum_inner']:.4e} / {last['lum_outer']:.4e} erg/s")
        if np.isfinite(last["max_inc"]):
            print(f"  {'Last max % increase':<28} "
                  f"{last['max_inc']:.3e} at depth {last['max_inc_depth']}")
    if _HAVE_CONV and iters:
        hist = ConvergenceHistory()
        for r_i in iters:
            hist.iterations.append(IterationRecord(
                iteration=r_i["iteration"], max_inc=r_i["max_inc"],
                max_inc_depth=r_i["max_inc_depth"], max_dec=r_i["max_dec"],
                max_dec_depth=r_i["max_dec_depth"], maxch=r_i["maxch"],
                lum_inner=r_i.get("lum_inner", np.nan),
                lum_outer=r_i.get("lum_outer", np.nan)))
        try:
            a = assess_convergence(hist, window=10)
            print(f"  {'Convergence assessment':<28} {a.classification}")
        except Exception as e:
            print(f"  Convergence assessment failed: {e}")
    if warnings_path.exists():
        n = sum(1 for _ in open(warnings_path, errors="replace"))
        print(f"  {'WARNINGS lines':<28} {n}")
    print("----------------------------------\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmfgen_dir", help="CMFGEN run directory")
    ap.add_argument("--outdir", default=None,
                    help="Plot output dir (default: <cmfgen_dir>/plots)")
    args = ap.parse_args()

    cdir = Path(args.cmfgen_dir)
    if not cdir.is_dir():
        print(f"ERROR: {cdir} is not a directory", file=sys.stderr)
        sys.exit(1)
    rvtj = cdir / "RVTJ"
    if not rvtj.exists():
        print(f"ERROR: {rvtj} not found (model not finished?)",
              file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir) if args.outdir else cdir / "plots"
    outdir.mkdir(parents=True, exist_ok=True)
    label = f"{cdir.parent.parent.name} {cdir.name}"

    print(f"CMFGEN diagnostics for {cdir}")

    hdr, sec_rvtj = read_labeled_sections(rvtj)
    sec_hydro = {}
    hydro = cdir / "SN_HYDRO_DATA"
    if hydro.exists():
        _, sec_hydro = read_labeled_sections(hydro)
    else:
        print("  WARNING: SN_HYDRO_DATA missing; input overlays skipped")

    iters, truncated = ([], False)
    outgen = cdir / "OUTGEN"
    if outgen.exists():
        iters, truncated = parse_outgen(outgen)
    else:
        print("  WARNING: OUTGEN missing; convergence plot skipped")

    made = []
    for fn in (
        lambda: plot_structure(sec_rvtj, sec_hydro, outdir, label),
        lambda: plot_regrid(cdir, sec_rvtj, sec_hydro, outdir, label),
        lambda: plot_convergence(iters, truncated, outdir, label),
        lambda: plot_opacity(sec_rvtj, outdir, label),
        lambda: plot_radiation(sec_rvtj, outdir, label),
        lambda: plot_ionization(cdir, sec_rvtj, outdir, label),
        lambda: plot_departure(sec_rvtj, outdir, label),
    ):
        try:
            out = fn()
            if out:
                print(f"  wrote {out}")
                made.append(out)
        except Exception as e:
            print(f"  WARNING: plot failed: {e}")

    print_summary(hdr, sec_rvtj, sec_hydro, iters, cdir / "WARNINGS")
    if not made:
        print("ERROR: no diagnostic plots could be produced", file=sys.stderr)
        sys.exit(1)
    print(f"Done: {len(made)} diagnostic plot(s) in {outdir}")


if __name__ == "__main__":
    main()
