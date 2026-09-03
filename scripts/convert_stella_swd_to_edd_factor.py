#!/usr/bin/env python3
"""Convert STELLA .swd_rad radiation data to CMFGEN EDDFACTOR file.

Reads the frequency-dependent radiation moments written by the modified
STELLA .swd_rad output and constructs a binary EDDFACTOR file compatible
with CMFGEN's direct-access read mechanism.

EDDFACTOR file format (from open_rw_eddfactor.f + eddfac_rec_defs_mod.f):
  - Direct-access unformatted binary
  - Record 1 (DIE_CONT_REC): integer access pointer
  - Record 2 (RV_REC): integer R/V coordinate pointer
  - Record 3 (EDD_CONT_REC): access pointer + NCF + NDEXT
  - Record 4 (LINE_CONT_REC): integer access pointer
  - Record 5 (FINISH_REC): nonzero real when all frequencies computed
  - Records 10-12 (INITIAL_RV_REC): R(1:ND), V(1:ND), LANG_COORD(1:ND)
  - Records 15+ (INITIAL_ACCESS_REC): per-frequency RJ(1:ND) + FL

Each per-frequency record contains ND+1 doubles: RJ(1:ND) then FL (Hz).
Frequencies are ordered highest to lowest (CMFGEN convention).
Depth points are ordered outer to inner (CMFGEN convention, depth 1 = outermost).

Usage:
  python convert_stella_swd_to_edd_factor.py SWD_RAD_FILE OUTPUT_DIR [OPTIONS]
    --target-days T    Pick epoch closest to T days (default: last epoch)
    --nd-target N      Target number of CMFGEN depth points (default: from SN_HYDRO_DATA)
    --sn-hydro-data P  Path to SN_HYDRO_DATA for grid mapping (preferred)
    --freq-file P      Path to CMFGEN frequency grid file (optional)
    --plot-dir P       Directory for diagnostic plots

The converter:
  1. Reads .swd_rad for frequency-dependent radiation moments
  2. Maps STELLA depth grid -> CMFGEN depth grid
  3. Converts occupation numbers -> physical J_nu
  4. Interpolates STELLA frequency groups -> CMFGEN frequency grid
  5. Reverses ordering (inner->outer to outer->inner)
  6. Writes EDDFACTOR binary file + EDDFACTOR_INFO
  7. Generates diagnostic plots
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path
from typing import Any

import numpy as np

# Physical constants (CGS)
H_PLANCK = 6.62607015e-27    # erg*s
K_BOLTZ = 1.380649e-16       # erg/K
C_LIGHT = 2.99792458e10      # cm/s
SIGMA_SB = 5.670374419e-5    # erg/cm^2/s/K^4
PI = np.pi

# EDDFACTOR record layout constants (from eddfac_rec_defs_mod.f)
DIE_CONT_REC = 1
RV_REC = 2
EDD_CONT_REC = 3
LINE_CONT_REC = 4
FINISH_REC = 5
INITIAL_RV_REC = 10
INITIAL_ACCESS_REC = 15


# ============================================================================
# .swd_rad reader
# ============================================================================

class SwdRadEpoch:
    """One epoch of frequency-dependent radiation data from .swd_rad."""

    def __init__(self) -> None:
        self.t_day: float = 0.0
        self.nfrus: int = 0
        self.nzon: int = 0
        self.ncnd: int = 0
        self.ufreq: float = 0.0
        self.utp: float = 0.0
        self.ur: float = 0.0
        self.utime: float = 0.0
        self.freqmn: np.ndarray = np.array([])      # (NFRUS,) dimensionless freq
        self.weight: np.ndarray = np.array([])       # (NFRUS,) bandwidth*freq^3
        self.r_cm: np.ndarray = np.array([])         # (NZON,) radius in cm
        self.v_kms: np.ndarray = np.array([])        # (NZON,) velocity in km/s
        self.t_gas: np.ndarray = np.array([])        # (NZON,) gas temperature K
        self.t_rad: np.ndarray = np.array([])        # (NZON,) radiation temperature K
        self.rho: np.ndarray = np.array([])          # (NZON,) density g/cm^3
        # (NZON, NFRUS) — comoving-frame J occupation number from Feau
        self.fradj: np.ndarray = np.array([])
        # (NZON, NFRUS) — comoving-frame H occupation number from Feau
        self.fradh: np.ndarray = np.array([])
        # (NZON, NFRUS) — Eddington factor f = K/J from Feau
        self.eddn: np.ndarray = np.array([])
        # (NZON, NFRUS) — lab-frame J occupation number from Y array
        self.y_j: np.ndarray = np.array([])
        # (NZON, NFRUS) — lab-frame H occupation number from Y array
        self.y_h: np.ndarray = np.array([])

    @property
    def freq_hz(self) -> np.ndarray:
        """Physical frequency in Hz: nu = FREQMN * UFREQ."""
        return self.freqmn * self.ufreq

    @property
    def j_nu_comoving(self) -> np.ndarray:
        """Physical J_nu in erg/cm^2/s/sr/Hz from comoving FradJ.

        J_nu = (2h*nu^3/c^2) * FradJ
        where nu = FREQMN * UFREQ.
        """
        nu = self.freq_hz  # (NFRUS,)
        scale = 2.0 * H_PLANCK * nu**3 / C_LIGHT**2  # (NFRUS,)
        return self.fradj * scale[np.newaxis, :]  # (NZON, NFRUS)

    @property
    def j_nu_lab(self) -> np.ndarray:
        """Physical J_nu from lab-frame Y-array occupation numbers."""
        nu = self.freq_hz
        scale = 2.0 * H_PLANCK * nu**3 / C_LIGHT**2
        return self.y_j * scale[np.newaxis, :]

    @property
    def h_nu_lab(self) -> np.ndarray:
        """Physical H_nu from lab-frame Y-array occupation numbers."""
        nu = self.freq_hz
        scale = 2.0 * H_PLANCK * nu**3 / C_LIGHT**2
        return self.y_h * scale[np.newaxis, :]


def read_swd_rad(path: str) -> dict[float, SwdRadEpoch]:
    """Parse a STELLA .swd_rad file.

    Returns dict keyed by t_day -> SwdRadEpoch.
    """
    epochs: dict[float, SwdRadEpoch] = {}
    current: SwdRadEpoch | None = None
    zone_idx = 0
    freq_idx = 0
    # State: 0=expecting header, 1=reading freq grid, 2=expecting zone header,
    #        3=reading freq data for current zone
    state = 0

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("SWDRAD_V1"):
                # Flush previous epoch
                if current is not None:
                    epochs[current.t_day] = current

                parts = line.split()
                current = SwdRadEpoch()
                current.nfrus = int(parts[1])
                current.t_day = float(parts[2])
                current.nzon = int(parts[3])
                current.ncnd = int(parts[4])
                current.ufreq = float(parts[5])
                current.utp = float(parts[6])
                current.ur = float(parts[7])
                current.utime = float(parts[8])

                # Allocate arrays
                n = current.nfrus
                nz = current.nzon
                current.freqmn = np.zeros(n)
                current.weight = np.zeros(n)
                current.r_cm = np.zeros(nz)
                current.v_kms = np.zeros(nz)
                current.t_gas = np.zeros(nz)
                current.t_rad = np.zeros(nz)
                current.rho = np.zeros(nz)
                current.fradj = np.zeros((nz, n))
                current.fradh = np.zeros((nz, n))
                current.eddn = np.zeros((nz, n))
                current.y_j = np.zeros((nz, n))
                current.y_h = np.zeros((nz, n))

                freq_idx = 0
                zone_idx = 0
                state = 1  # reading freq grid
                continue

            if current is None:
                continue

            parts = line.split()

            if state == 1:
                # Reading frequency grid: NFRUS lines of FREQMN WEIGHT
                if len(parts) >= 2:
                    current.freqmn[freq_idx] = float(parts[0])
                    current.weight[freq_idx] = float(parts[1])
                    freq_idx += 1
                if freq_idx >= current.nfrus:
                    state = 2  # now expect zone headers
                    zone_idx = 0
                    freq_idx = 0
                continue

            if state == 2:
                # Expecting zone header: Km R_cm V_kms T_gas T_rad rho
                if len(parts) >= 6:
                    km = int(parts[0]) - 1  # 0-indexed
                    current.r_cm[km] = float(parts[1])
                    current.v_kms[km] = float(parts[2])
                    current.t_gas[km] = float(parts[3])
                    current.t_rad[km] = float(parts[4])
                    current.rho[km] = float(parts[5])
                    state = 3  # now read freq data for this zone
                    freq_idx = 0
                continue

            if state == 3:
                # Reading frequency data: L FradJ FradH EddN Y_J Y_H
                if len(parts) >= 6:
                    l = int(parts[0]) - 1  # 0-indexed
                    km = zone_idx
                    current.fradj[km, l] = float(parts[1])
                    current.fradh[km, l] = float(parts[2])
                    current.eddn[km, l] = float(parts[3])
                    current.y_j[km, l] = float(parts[4])
                    current.y_h[km, l] = float(parts[5])
                    freq_idx += 1
                    if freq_idx >= current.nfrus:
                        zone_idx += 1
                        freq_idx = 0
                        if zone_idx >= current.nzon:
                            state = 0  # done with this epoch
                        else:
                            state = 2  # expect next zone header
                continue

    # Flush last epoch
    if current is not None:
        epochs[current.t_day] = current

    return epochs


# ============================================================================
# SN_HYDRO_DATA reader (for CMFGEN grid mapping)
# ============================================================================

def read_sn_hydro_data(path: str) -> dict[str, np.ndarray]:
    """Read CMFGEN SN_HYDRO_DATA to get the target depth grid.

    Returns dict with 'r', 'v', 't', 'rho' arrays (outer-to-inner, CGS units).
    """
    data: dict[str, list[float]] = {"r": [], "v": [], "t": [], "rho": []}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("!") or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    data["r"].append(float(parts[0]))
                    data["v"].append(float(parts[1]))
                    data["t"].append(float(parts[2]))
                    data["rho"].append(float(parts[3]))
                except ValueError:
                    continue

    return {k: np.array(v) for k, v in data.items()}


# ============================================================================
# Grid mapping: STELLA -> CMFGEN
# ============================================================================

def map_stella_to_cmfgen_depth(
    stella_r: np.ndarray,
    cmfgen_r: np.ndarray,
) -> np.ndarray:
    """Map CMFGEN depth grid to STELLA zone indices using radius.

    STELLA zones are inner-to-outer (increasing r).
    CMFGEN depths are outer-to-inner (decreasing r).

    Returns array of STELLA zone indices (0-based) for each CMFGEN depth point.
    Uses linear interpolation in radius space.

    Args:
        stella_r: STELLA radii (cm), inner-to-outer (increasing)
        cmfgen_r: CMFGEN radii (cm), outer-to-inner (decreasing)

    Returns:
        indices: array of interpolated STELLA zone indices for each CMFGEN depth
    """
    # STELLA zones: 0, 1, ..., nzon-1 with increasing r
    stella_x = np.arange(len(stella_r))

    # Sort STELLA by radius (should already be increasing)
    sort_idx = np.argsort(stella_r)
    stella_r_sorted = stella_r[sort_idx]
    stella_x_sorted = stella_x[sort_idx]

    # Interpolate: for each CMFGEN radius, find STELLA zone index
    # Clip to valid range
    cmfgen_r_clipped = np.clip(cmfgen_r, stella_r_sorted[0], stella_r_sorted[-1])
    indices = np.interp(cmfgen_r_clipped, stella_r_sorted, stella_x)

    return indices


def interpolate_to_cmfgen_grid(
    stella_data: np.ndarray,
    stella_r: np.ndarray,
    cmfgen_r: np.ndarray,
) -> np.ndarray:
    """Interpolate a per-zone STELLA quantity to CMFGEN depth grid.

    Args:
        stella_data: (NZON,) or (NZON, NFRUS) array
        stella_r: (NZON,) STELLA radii, inner-to-outer
        cmfgen_r: (ND,) CMFGEN radii, outer-to-inner

    Returns:
        (ND,) or (ND, NFRUS) array on CMFGEN grid
    """
    # Sort STELLA by radius
    sort_idx = np.argsort(stella_r)
    r_sorted = stella_r[sort_idx]

    if stella_data.ndim == 1:
        d_sorted = stella_data[sort_idx]
        cmfgen_r_clipped = np.clip(cmfgen_r, r_sorted[0], r_sorted[-1])
        return np.interp(cmfgen_r_clipped, r_sorted, d_sorted)
    else:
        d_sorted = stella_data[sort_idx, :]
        cmfgen_r_clipped = np.clip(cmfgen_r, r_sorted[0], r_sorted[-1])
        result = np.zeros((len(cmfgen_r), stella_data.shape[1]))
        for j in range(stella_data.shape[1]):
            result[:, j] = np.interp(cmfgen_r_clipped, r_sorted, d_sorted[:, j])
        return result


# ============================================================================
# Frequency interpolation
# ============================================================================

def interpolate_frequency_grid(
    stella_j: np.ndarray,
    stella_freq_hz: np.ndarray,
    cmfgen_freq_hz: np.ndarray,
) -> np.ndarray:
    """Interpolate J_nu from STELLA frequency grid to CMFGEN frequency grid.

    STELLA frequencies: ordered low to high (FREQMN increasing)
    CMFGEN frequencies: ordered high to low (by convention)

    Args:
        stella_j: (ND, NFRUS_STELLA) J_nu on STELLA freq grid
        stella_freq_hz: (NFRUS_STELLA,) STELLA frequencies in Hz
        cmfgen_freq_hz: (NCF,) CMFGEN frequencies in Hz

    Returns:
        (ND, NCF) J_nu on CMFGEN freq grid
    """
    nd = stella_j.shape[0]
    ncf = len(cmfgen_freq_hz)
    result = np.zeros((nd, ncf))

    # Sort STELLA frequencies ascending
    sort_idx = np.argsort(stella_freq_hz)
    freq_sorted = stella_freq_hz[sort_idx]

    for i in range(nd):
        j_sorted = stella_j[i, sort_idx]
        # Clip CMFGEN frequencies to STELLA range
        freq_clipped = np.clip(cmfgen_freq_hz, freq_sorted[0], freq_sorted[-1])
        # Interpolate in log space for better behavior
        log_j = np.log10(np.maximum(j_sorted, 1e-300))
        log_result = np.interp(freq_clipped, freq_sorted, log_j)
        result[i, :] = 10.0**log_result

    return result


# ============================================================================
# EDDFACTOR writer
# ============================================================================

def write_edd_factor(
    output_path: str,
    r_cm: np.ndarray,
    v_kms: np.ndarray,
    j_nu: np.ndarray,
    freq_hz: np.ndarray,
) -> None:
    """Write CMFGEN EDDFACTOR binary file + EDDFACTOR_INFO.

    Args:
        output_path: path for EDDFACTOR file
        r_cm: (ND,) radii in cm, outer-to-inner
        v_kms: (ND,) velocities in km/s, outer-to-inner
        j_nu: (ND, NCF) mean intensity in erg/cm^2/s/sr/Hz
              ND depth points (outer-to-inner), NCF frequencies
        freq_hz: (NCF,) frequencies in Hz, ordered highest to lowest
    """
    nd = len(r_cm)
    ncf = len(freq_hz)

    # Convert velocities to cm/s for CMFGEN
    v_cms = v_kms * 1e5

    # LANG_COORD: not used by CMFGEN for spherical models, set to 0
    lang_coord = np.zeros(nd)

    # Determine record size
    # Each record holds ND+1 doubles (RJ(1:ND) + FL)
    # Fortran direct-access RECL is in platform-dependent units
    # We use numpy to write binary and determine sizes
    word_size = 8  # double precision
    record_size_bytes = word_size * (nd + 1)

    # Write EDDFACTOR_INFO companion file
    # Must match CMFGEN's WRITE_DIRECT_INFO_V3 format exactly:
    #   Line 1: ' 12-Apr-2017                                       !INFO format date'
    #   Line 2: ' 23-Jan-2017                                       !File format date'
    #   Line 3: '      {ND}      {RECL}      {WORD_SIZE}      {UNIT_SIZE}      {INT_SIZE}        {LIT_END}'
    # RECL = WORD_SIZE*(ND+1)/UNIT_SIZE (in system units, not bytes)
    # On ifort: WORD_SIZE=8, UNIT_SIZE=4, INT_SIZE=4, LIT_END=T
    info_path = output_path + "_INFO"
    word_size = 8      # double precision bytes
    unit_size = 4      # ifort IOLENGTH units (bytes per RECL unit)
    int_size = 4       # default integer kind bytes
    recl = word_size * (nd + 1) // unit_size  # RECL in system units
    with open(info_path, "w") as f_info:
        f_info.write(" 12-Apr-2017                                       !INFO format date\n")
        f_info.write(" 23-Jan-2017                                       !File format date\n")
        f_info.write(f"{recl:12d}{word_size:12d}{unit_size:12d}{int_size:12d}            T\n")
        f_info.write("        ND        RECL   WORD_SIZE   UNIT_SIZE    INT_SIZE     LIT_END\n")

    # Convert to CMFGEN internal units before writing:
    #   R: cm -> 10^10 cm (divide by 1e10)
    #   V: cm/s -> km/s (divide by 1e5)
    #   Frequency: Hz -> 10^15 Hz (divide by 1e15)
    r_cmfgen = r_cm / 1.0e10
    v_cmfgen = v_cms / 1.0e5
    freq_cmfgen = freq_hz / 1.0e15

    with open(output_path, "wb") as f:
        # Record 1 (DIE_CONT_REC): integer access pointer
        die_cont = INITIAL_ACCESS_REC + ncf
        f.write(struct.pack("i", die_cont))
        _pad_to_record(f, record_size_bytes, 4)

        # Record 2 (RV_REC): integer R/V coordinate pointer
        rv_rec_val = INITIAL_RV_REC
        f.write(struct.pack("i", rv_rec_val))
        _pad_to_record(f, record_size_bytes, 4)

        # Record 3 (EDD_CONT_REC): access pointer + NCF + NDEXT
        f.write(struct.pack("i", INITIAL_ACCESS_REC))
        f.write(struct.pack("i", ncf))
        f.write(struct.pack("i", nd))
        _pad_to_record(f, record_size_bytes, 12)

        # Record 4 (LINE_CONT_REC): integer access pointer
        f.write(struct.pack("i", INITIAL_ACCESS_REC + ncf))
        _pad_to_record(f, record_size_bytes, 4)

        # Record 5 (FINISH_REC): nonzero real when complete
        f.write(struct.pack("d", 1.0))
        _pad_to_record(f, record_size_bytes, 8)

        # Records 6-9: reserved (zeros)
        for _ in range(4):
            f.write(b"\x00" * record_size_bytes)

        # Records 10-12 (INITIAL_RV_REC): R, V, LANG_COORD in CMFGEN units
        f.write(r_cmfgen.astype(np.float64).tobytes())
        _pad_to_record(f, record_size_bytes, nd * 8)

        f.write(v_cmfgen.astype(np.float64).tobytes())
        _pad_to_record(f, record_size_bytes, nd * 8)

        f.write(lang_coord.astype(np.float64).tobytes())
        _pad_to_record(f, record_size_bytes, nd * 8)

        # Records 13-14: padding to reach INITIAL_ACCESS_REC=15
        for _ in range(2):
            f.write(b"\x00" * record_size_bytes)

        # Records 15+ (INITIAL_ACCESS_REC): per-frequency RJ(1:ND) + FL
        for j in range(ncf):
            record_data = np.zeros(nd + 1, dtype=np.float64)
            record_data[:nd] = j_nu[:, j]
            record_data[nd] = freq_cmfgen[j]
            f.write(record_data.tobytes())


def _pad_to_record(f, record_size_bytes: int, written_bytes: int) -> None:
    """Pad current position to next record boundary."""
    remaining = record_size_bytes - (written_bytes % record_size_bytes)
    if remaining < record_size_bytes:
        f.write(b"\x00" * remaining)


# ============================================================================
# Validation
# ============================================================================

def validate_edd_factor(
    j_nu: np.ndarray,
    edd_factor: np.ndarray | None,
    r_cm: np.ndarray,
    freq_hz: np.ndarray,
) -> list[str]:
    """Validate EDDFACTOR contents. Returns list of warnings (empty = OK)."""
    warnings: list[str] = []

    # Check for NaNs/Infs
    if np.any(~np.isfinite(j_nu)):
        warnings.append("J_nu contains NaN or Inf values")
    if edd_factor is not None and np.any(~np.isfinite(edd_factor)):
        warnings.append("Eddington factor contains NaN or Inf values")

    # Check for negative J
    if np.any(j_nu < 0):
        n_neg = np.sum(j_nu < 0)
        warnings.append(f"J_nu has {n_neg} negative values (will be clamped)")

    # Check Eddington factor range [1/3, 1]
    if edd_factor is not None:
        f_min = np.min(edd_factor)
        f_max = np.max(edd_factor)
        if f_min < 0.2:
            warnings.append(f"Eddington factor min={f_min:.4f} < 0.2 (expected >= 1/3)")
        if f_max > 1.5:
            warnings.append(f"Eddington factor max={f_max:.4f} > 1.5 (expected <= 1)")

    # Check ordering
    if not np.all(np.diff(r_cm) <= 0):
        warnings.append("Radii not in decreasing order (expected outer-to-inner)")

    if len(freq_hz) > 1 and not np.all(np.diff(freq_hz) <= 0):
        warnings.append("Frequencies not in decreasing order (expected high-to-low)")

    return warnings


# ============================================================================
# Diagnostic plots
# ============================================================================

def generate_diagnostic_plots(
    epoch: SwdRadEpoch,
    cmfgen_r: np.ndarray,
    j_nu_cmfgen: np.ndarray,
    edd_factor_cmfgen: np.ndarray | None,
    freq_hz_cmfgen: np.ndarray,
    plot_dir: str,
) -> None:
    """Generate diagnostic plots for the EDDFACTOR conversion."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_path = Path(plot_dir)
    plot_path.mkdir(parents=True, exist_ok=True)

    # 1. Eddington factor vs radius
    if edd_factor_cmfgen is not None:
        fig, ax = plt.subplots(figsize=(10, 6))
        for j in [0, len(freq_hz_cmfgen) // 4, len(freq_hz_cmfgen) // 2,
                  3 * len(freq_hz_cmfgen) // 4, len(freq_hz_cmfgen) - 1]:
            if 0 <= j < edd_factor_cmfgen.shape[1]:
                ax.plot(cmfgen_r / 1e15, edd_factor_cmfgen[:, j],
                        label=f"nu={freq_hz_cmfgen[j]:.2e} Hz")
        ax.axhline(1.0 / 3.0, color="k", ls="--", alpha=0.5, label="f=1/3")
        ax.axhline(1.0, color="k", ls=":", alpha=0.5, label="f=1")
        ax.set_xlabel("Radius (1e15 cm)")
        ax.set_ylabel("Eddington factor f")
        ax.set_title("Eddington Factor vs Radius (CMFGEN grid)")
        ax.legend()
        ax.set_ylim(0, 1.5)
        fig.savefig(plot_path / "edd_factor_vs_radius.png", dpi=150)
        plt.close(fig)

    # 2. J_nu profiles at selected frequencies
    fig, ax = plt.subplots(figsize=(10, 6))
    for j in [0, len(freq_hz_cmfgen) // 4, len(freq_hz_cmfgen) // 2,
              3 * len(freq_hz_cmfgen) // 4, len(freq_hz_cmfgen) - 1]:
        if 0 <= j < j_nu_cmfgen.shape[1]:
            ax.plot(cmfgen_r / 1e15, np.log10(np.maximum(j_nu_cmfgen[:, j], 1e-300)),
                    label=f"nu={freq_hz_cmfgen[j]:.2e} Hz")
    ax.set_xlabel("Radius (1e15 cm)")
    ax.set_ylabel("log10 J_nu (erg/cm^2/s/sr/Hz)")
    ax.set_title("J_nu vs Radius (CMFGEN grid)")
    ax.legend()
    fig.savefig(plot_path / "jnu_vs_radius.png", dpi=150)
    plt.close(fig)

    # 3. Radiation temperature vs radius
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(epoch.r_cm / 1e15, epoch.t_gas, "b-", label="T_gas (STELLA)")
    ax.plot(epoch.r_cm / 1e15, epoch.t_rad, "r-", label="T_rad (STELLA)")
    ax.set_xlabel("Radius (1e15 cm)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Temperature Profiles (STELLA)")
    ax.legend()
    ax.set_yscale("log")
    fig.savefig(plot_path / "temperature_vs_radius.png", dpi=150)
    plt.close(fig)

    # 4. STELLA Eddington factor vs radius (raw)
    fig, ax = plt.subplots(figsize=(10, 6))
    for j in [0, epoch.nfrus // 4, epoch.nfrus // 2,
              3 * epoch.nfrus // 4, epoch.nfrus - 1]:
        if 0 <= j < epoch.eddn.shape[1]:
            ax.plot(epoch.r_cm / 1e15, epoch.eddn[:, j],
                    label=f"L={j+1}")
    ax.axhline(1.0 / 3.0, color="k", ls="--", alpha=0.5, label="f=1/3")
    ax.set_xlabel("Radius (1e15 cm)")
    ax.set_ylabel("Eddington factor f")
    ax.set_title("Eddington Factor vs Radius (STELLA raw)")
    ax.legend()
    ax.set_ylim(0, 1.5)
    fig.savefig(plot_path / "stella_edd_factor_vs_radius.png", dpi=150)
    plt.close(fig)


# ============================================================================
# Main conversion function
# ============================================================================

def convert_swd_rad_to_edd_factor(
    swd_rad_path: str,
    output_dir: str,
    target_days: float | None = None,
    sn_hydro_data_path: str | None = None,
    cmfgen_freq_hz: np.ndarray | None = None,
    plot_dir: str | None = None,
) -> dict[str, Any]:
    """Convert STELLA .swd_rad to CMFGEN EDDFACTOR file.

    Args:
        swd_rad_path: Path to .swd_rad file
        output_dir: Directory for EDDFACTOR output
        target_days: Epoch to select (default: last)
        sn_hydro_data_path: Path to SN_HYDRO_DATA for grid mapping
        cmfgen_freq_hz: CMFGEN frequency grid in Hz (optional)
        plot_dir: Directory for diagnostic plots

    Returns:
        dict with conversion metadata and validation results
    """
    # 1. Read .swd_rad
    epochs = read_swd_rad(swd_rad_path)
    if not epochs:
        raise ValueError(f"No epochs found in {swd_rad_path}")

    # 2. Select epoch
    if target_days is not None:
        available = sorted(epochs.keys())
        t_target = target_days
        t_best = min(available, key=lambda t: abs(t - t_target))
        epoch = epochs[t_best]
    else:
        t_best = max(epochs.keys())
        epoch = epochs[t_best]

    print(f"Selected epoch: t={epoch.t_day:.4f} days "
          f"({epoch.nzon} zones, {epoch.nfrus} freq groups)")

    # 3. Get CMFGEN depth grid
    if sn_hydro_data_path:
        hydro = read_sn_hydro_data(sn_hydro_data_path)
        cmfgen_r = hydro["r"]
        cmfgen_v = hydro["v"]
        print(f"CMFGEN grid: {len(cmfgen_r)} depth points from SN_HYDRO_DATA")
    else:
        # Use STELLA grid directly (reversed to outer-to-inner)
        cmfgen_r = epoch.r_cm[::-1].copy()
        cmfgen_v = epoch.v_kms[::-1].copy() * 1e5  # cm/s
        print(f"Using STELLA grid directly: {len(cmfgen_r)} depth points")

    nd = len(cmfgen_r)

    # 4. Get CMFGEN frequency grid
    if cmfgen_freq_hz is None:
        # Use STELLA frequency grid, reversed to high-to-low
        cmfgen_freq_hz = epoch.freq_hz[::-1].copy()
        print(f"Using STELLA frequency grid: {len(cmfgen_freq_hz)} frequencies")
    else:
        print(f"Using provided CMFGEN frequency grid: {len(cmfgen_freq_hz)} frequencies")

    ncf = len(cmfgen_freq_hz)

    # 5. Convert STELLA J to physical units and interpolate to CMFGEN grid
    # Use comoving-frame J (FradJ) as primary source
    j_nu_stella = epoch.j_nu_comoving  # (NZON, NFRUS) in erg/cm^2/s/sr/Hz

    # Check if FradJ has meaningful data (Feau may not have been called)
    if np.all(j_nu_stella == 0):
        print("WARNING: FradJ is all zeros — falling back to lab-frame Y_J")
        j_nu_stella = epoch.j_nu_lab

    # Interpolate to CMFGEN depth grid
    j_nu_on_cmfgen_depth = interpolate_to_cmfgen_grid(
        j_nu_stella, epoch.r_cm, cmfgen_r
    )  # (ND, NFRUS)

    # Interpolate to CMFGEN frequency grid
    j_nu_cmfgen = interpolate_frequency_grid(
        j_nu_on_cmfgen_depth, epoch.freq_hz, cmfgen_freq_hz
    )  # (ND, NCF)

    # 6. Interpolate Eddington factor to CMFGEN grid
    edd_factor_stella = epoch.eddn  # (NZON, NFRUS)
    edd_on_cmfgen_depth = interpolate_to_cmfgen_grid(
        edd_factor_stella, epoch.r_cm, cmfgen_r
    )
    edd_factor_cmfgen = interpolate_frequency_grid(
        edd_on_cmfgen_depth, epoch.freq_hz, cmfgen_freq_hz
    )

    # Clamp Eddington factor to physical range
    edd_factor_cmfgen = np.clip(edd_factor_cmfgen, 1.0 / 3.0, 1.0)

    # 7. Validate
    warnings = validate_edd_factor(j_nu_cmfgen, edd_factor_cmfgen,
                                   cmfgen_r, cmfgen_freq_hz)
    for w in warnings:
        print(f"  WARNING: {w}")

    # 8. Write EDDFACTOR file
    output_path = str(Path(output_dir) / "EDDFACTOR")
    write_edd_factor(output_path, cmfgen_r, cmfgen_v, j_nu_cmfgen, cmfgen_freq_hz)
    print(f"Wrote EDDFACTOR to {output_path}")

    # 9. Generate plots
    if plot_dir:
        generate_diagnostic_plots(
            epoch, cmfgen_r, j_nu_cmfgen, edd_factor_cmfgen,
            cmfgen_freq_hz, plot_dir
        )
        print(f"Diagnostic plots written to {plot_dir}")

    # 10. Return metadata
    return {
        "swd_rad_file": swd_rad_path,
        "selected_epoch_days": epoch.t_day,
        "stella_zones": epoch.nzon,
        "stella_freq_groups": epoch.nfrus,
        "cmfgen_depth_points": nd,
        "cmfgen_freq_points": ncf,
        "j_source": "FradJ (comoving)" if not np.all(epoch.fradj == 0) else "Y_J (lab-frame)",
        "warnings": warnings,
        "edd_factor_min": float(np.min(edd_factor_cmfgen)),
        "edd_factor_max": float(np.max(edd_factor_cmfgen)),
        "j_nu_min": float(np.min(j_nu_cmfgen[j_nu_cmfgen > 0])),
        "j_nu_max": float(np.max(j_nu_cmfgen)),
    }


# ============================================================================
# CLI
# ============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert STELLA .swd_rad to CMFGEN EDDFACTOR file"
    )
    parser.add_argument("swd_rad_file", help="Path to STELLA .swd_rad file")
    parser.add_argument("output_dir", help="Output directory for EDDFACTOR")
    parser.add_argument("--target-days", type=float, default=None,
                        help="Epoch to select (default: last)")
    parser.add_argument("--sn-hydro-data", default=None,
                        help="Path to SN_HYDRO_DATA for grid mapping")
    parser.add_argument("--plot-dir", default=None,
                        help="Directory for diagnostic plots")
    parser.add_argument("--json-output", default=None,
                        help="Write conversion metadata as JSON to this path")

    args = parser.parse_args()

    result = convert_swd_rad_to_edd_factor(
        swd_rad_path=args.swd_rad_file,
        output_dir=args.output_dir,
        target_days=args.target_days,
        sn_hydro_data_path=args.sn_hydro_data,
        plot_dir=args.plot_dir,
    )

    if args.json_output:
        with open(args.json_output, "w") as f:
            json.dump(result, f, indent=2)

    print("\nConversion complete:")
    for k, v in result.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
