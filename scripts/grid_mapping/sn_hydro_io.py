"""Round-trip SN_HYDRO_DATA reader and hard-format validator.

Reads back a CMFGEN ``SN_HYDRO_DATA`` file exactly as CMFGEN's
``rd_sn_data.f`` / ``set_rv_hydro_model_v3.f`` do (header-keyword search,
free-format block reads, outer-to-inner ordering) and enforces the hard
requirements that CMFGEN demands before the file can be used.

Units (confirmed against the local CMFGEN source and prepare_VADAT.py):
    Radius            : 10^10 cm   (set_rv_hydro_model_v3.f L228 comment)
    Velocity          : km/s       (set_rv_hydro_model_v3.f L228 comment)
    Temperature       : 10^4 K     (prepare_VADAT.py L1235)
    Density           : g/cm^3
    Atom/Electron den.: cm^-3
    Rosseland chi     : 10^-10 cm^-1
    Kappa             : cm^2/g
    Mass fractions    : linear (sum to 1 per depth)
    Isotopes          : "<name> <baryon> mass fraction" header + values

The reader is intentionally independent of the writer so it can be used
as a cross-check (round-trip) in tests and in the production pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np


@dataclass
class SnHydroData:
    """Container for a parsed SN_HYDRO_DATA file."""
    n_points: int
    n_mass_fractions: int
    n_isotopes: int
    time_days: float
    r_1e10cm: np.ndarray               # outer->inner, decreasing
    v_kms: np.ndarray
    sigma: np.ndarray
    T_1e4K: np.ndarray
    density_gcm3: np.ndarray
    atom_density: Optional[np.ndarray] = None
    electron_density: Optional[np.ndarray] = None
    rosseland_chi: Optional[np.ndarray] = None
    kappa: Optional[np.ndarray] = None
    clumping_factor: Optional[np.ndarray] = None
    mass_fractions: dict[str, np.ndarray] = field(default_factory=dict)
    isotopes: dict[tuple[str, int], np.ndarray] = field(default_factory=dict)
    headers: dict[str, str] = field(default_factory=dict)  # extra header lines


def _read_numeric_block(lines: list[str], idx: int, n: int) -> tuple[np.ndarray, int]:
    """Read n free-format floats starting after line idx. Returns (values, next_idx)."""
    vals = []
    i = idx
    while i < len(lines):
        s = lines[i].split('!', 1)[0].strip()
        if not s or s.startswith('#'):
            i += 1
            continue
        tokens = s.split()
        try:
            row = [float(t.replace('D', 'E').replace('d', 'e')) for t in tokens]
        except ValueError:
            # Stop if we hit a non-numeric (next section label)
            break
        vals.extend(row)
        i += 1
    if len(vals) != n:
        raise ValueError(f"Expected exactly {n} values, got {len(vals)}")
    return np.asarray(vals, dtype=float), i


def read_sn_hydro_data(path: str | Path) -> SnHydroData:
    """Read a CMFGEN SN_HYDRO_DATA file (mirrors rd_sn_data.f parsing)."""
    path = Path(path)
    with open(path) as fh:
        raw = fh.read()
    lines = raw.splitlines()

    nx = nsp = niso = 0
    old_age = 0.0
    extra_headers: dict[str, str] = {}
    i = 0
    # Header scan (mirrors rd_sn_data.f L106-132)
    while i < len(lines):
        s = lines[i].strip()
        i += 1
        if not s or s.startswith('!'):
            continue
        if 'Number of data points:' in s:
            nx = int(s.split('points:')[1].split()[0])
        elif 'Number of mass fractions:' in s:
            nsp = int(s.split('fractions:')[1].split()[0])
        elif 'Number of isotopes:' in s:
            niso = int(s.split('isotopes:')[1].split()[0])
        elif 'Time(days) since explosion:' in s:
            old_age = float(s.split('explosion:')[1].split()[0]
                            .replace('D', 'E').replace('d', 'e'))
        elif 'Radius grid' in s:
            if nx == 0 or nsp == 0 or old_age == 0.0:
                raise ValueError("Radius grid reached before NX/NSP/time defined")
            i -= 1  # re-process this line in the block scan
            break
        else:
            # Capture optional headers (Energy Lum, Inner boundary ...)
            if ':' in s:
                key, _, val = s.partition(':')
                extra_headers[key.strip()] = val.strip()
        # continue scanning

    if nx < 2:
        raise ValueError(f"Number of data points must be >= 2, got {nx}")
    if nsp < 1:
        raise ValueError(f"Number of mass fractions must be >= 1, got {nsp}")
    if niso < 0:
        raise ValueError(f"Number of isotopes must be >= 0, got {niso}")
    if not np.isfinite(old_age) or old_age <= 0.0:
        raise ValueError(f"Time(days) since explosion must be positive and finite, got {old_age}")

    # Block scan (mirrors rd_sn_data.f L156-218): read labelled blocks until
    # mass fractions begin. Order is not strictly guaranteed by the reader.
    r = v = sigma = T = rho = None
    atom_den = elec_den = clump = chi = kappa = None
    mf: dict[str, np.ndarray] = {}
    iso: dict[tuple[str, int], np.ndarray] = {}

    # We re-scan from the 'Radius grid' line forward.
    seen_sections: set[str] = set()

    def _claim(label: str) -> None:
        if label in seen_sections:
            raise ValueError(f"Duplicate section: {label}")
        seen_sections.add(label)

    while i < len(lines):
        s = lines[i].strip()
        i += 1
        if not s or s.startswith('!'):
            continue
        if s.startswith('Radius grid'):
            _claim('Radius grid')
            r, i = _read_numeric_block(lines, i, nx)
        elif s.startswith('Velocity'):
            _claim('Velocity')
            v, i = _read_numeric_block(lines, i, nx)
        elif s.startswith('Sigma'):
            _claim('Sigma')
            sigma, i = _read_numeric_block(lines, i, nx)
        elif s.startswith('Temperature'):
            _claim('Temperature')
            T, i = _read_numeric_block(lines, i, nx)
        elif s.startswith('Density'):
            _claim('Density')
            rho, i = _read_numeric_block(lines, i, nx)
        elif s.startswith('Electron density'):
            _claim('Electron density')
            elec_den, i = _read_numeric_block(lines, i, nx)
        elif s.startswith('Clumping factor'):
            _claim('Clumping factor')
            clump, i = _read_numeric_block(lines, i, nx)
        elif s.startswith('Atom density'):
            _claim('Atom density')
            atom_den, i = _read_numeric_block(lines, i, nx)
        elif s.startswith('Rosseland'):
            _claim('Rosseland')
            chi, i = _read_numeric_block(lines, i, nx)
        elif s.startswith('Kappa'):
            _claim('Kappa')
            kappa, i = _read_numeric_block(lines, i, nx)
        elif 'mass fraction' in s:
            # Distinguish isotope blocks ("<SPECIES> <baryon> mass fraction")
            # from elemental mass fractions ("<SPECIES> mass fraction").
            parts = s.split()
            is_isotope = False
            if len(parts) >= 4 and parts[-2:] == ['mass', 'fraction']:
                try:
                    baryon = int(parts[-3])
                    spec = parts[-4]
                    is_isotope = True
                except ValueError:
                    pass
            if is_isotope:
                if (spec, baryon) in iso:
                    raise ValueError(f"Duplicate isotope block: {spec} {baryon}")
                vals, i = _read_numeric_block(lines, i, nx)
                iso[(spec, baryon)] = vals
            else:
                spec = s.split('mass fraction')[0].strip()
                if spec in mf:
                    raise ValueError(f"Duplicate mass-fraction block: {spec}")
                vals, i = _read_numeric_block(lines, i, nx)
                mf[spec] = vals
        else:
            tokens = s.split('!', 1)[0].split()
            if tokens:
                numeric = True
                for tok in tokens:
                    try:
                        float(tok.replace('D', 'E').replace('d', 'e'))
                    except ValueError:
                        numeric = False
                        break
                if numeric:
                    raise ValueError(
                        f"Unexpected numeric data outside a section: {s!r}")

    if r is None or v is None or sigma is None or T is None or rho is None:
        raise ValueError("Missing required block (R/V/Sigma/T/Density)")
    if len(mf) != nsp:
        raise ValueError(f"Expected {nsp} mass-fraction blocks, got {len(mf)}")
    if len(iso) != niso:
        raise ValueError(f"Expected {niso} isotope blocks, got {len(iso)}")

    return SnHydroData(
        n_points=nx, n_mass_fractions=nsp, n_isotopes=niso, time_days=old_age,
        r_1e10cm=r, v_kms=v, sigma=sigma, T_1e4K=T, density_gcm3=rho,
        atom_density=atom_den, electron_density=elec_den,
        rosseland_chi=chi, kappa=kappa, clumping_factor=clump,
        mass_fractions=mf, isotopes=iso, headers=extra_headers,
    )


@dataclass
class FormatValidationReport:
    passed: bool = True
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def validate_sn_hydro_data(path: str | Path,
                           require_strict_decreasing_r: bool = True,
                           vel_law: int = 12) -> FormatValidationReport:
    """Hard-format validation of a SN_HYDRO_DATA file.

    Fails (errors) on anything CMFGEN would reject or misinterpret:
      - declared counts vs actual rows
      - non-finite values
      - radius not strictly decreasing (CMFGEN requires R(I) > R(I+1))
      - non-positive density / temperature
      - velocity zero anywhere (sigma divides by V; rv_sn_model_sniin L266)
      - mass fractions not summing to ~1
      - isotope baryon numbers present

    Warnings (non-fatal) for: sigma in the -1.2..-1.0 clamp band, large
    cell ratios, T units sanity.
    """
    rep = FormatValidationReport()
    try:
        d = read_sn_hydro_data(path)
    except Exception as exc:  # noqa: BLE001
        rep.passed = False
        rep.errors.append(f"Unreadable: {exc}")
        return rep

    n = d.n_points

    def _check_finite(name, arr):
        if arr is None:
            return
        if not np.all(np.isfinite(arr)):
            rep.errors.append(f"{name}: {int(np.sum(~np.isfinite(arr)))} non-finite")
            rep.passed = False

    for name, arr in [('R', d.r_1e10cm), ('V', d.v_kms), ('Sigma', d.sigma),
                      ('T', d.T_1e4K), ('Density', d.density_gcm3),
                      ('Atom density', d.atom_density),
                      ('Electron density', d.electron_density),
                      ('Rosseland chi', d.rosseland_chi),
                      ('Kappa', d.kappa),
                      ('Clumping factor', d.clumping_factor)]:
        _check_finite(name, arr)
    for elem, arr in d.mass_fractions.items():
        _check_finite(f'{elem} mass fraction', arr)
    for (iso_name, baryon), arr in d.isotopes.items():
        _check_finite(f'{iso_name} {baryon} isotope', arr)
    if len(d.r_1e10cm) != n:
        rep.errors.append(f"R has {len(d.r_1e10cm)} rows, header says {n}")
        rep.passed = False
    for name, arr in [('V', d.v_kms), ('Sigma', d.sigma), ('T', d.T_1e4K),
                      ('Density', d.density_gcm3)]:
        if arr is not None and len(arr) != n:
            rep.errors.append(f"{name} has {len(arr)} rows, expected {n}")
            rep.passed = False

    # Radius strictly decreasing (outer->inner)
    if require_strict_decreasing_r and len(d.r_1e10cm) >= 2:
        diffs = np.diff(d.r_1e10cm)
        nbad = int(np.sum(diffs >= 0))
        if nbad:
            rep.errors.append(f"R not strictly decreasing ({nbad} violations)")
            rep.passed = False
        ndup = int(np.sum(np.abs(diffs) < 1e-30))
        if ndup:
            rep.errors.append(f"R has {ndup} duplicate adjacent points")
            rep.passed = False

    # Physical signs
    if np.any(d.r_1e10cm <= 0):
        rep.errors.append(f"R has {int(np.sum(d.r_1e10cm<=0))} non-positive")
        rep.passed = False
    if d.atom_density is not None and np.any(d.atom_density <= 0):
        rep.errors.append(f"atom density has {int(np.sum(d.atom_density<=0))} non-positive")
        rep.passed = False
    if d.electron_density is not None and np.any(d.electron_density <= 0):
        rep.errors.append(f"electron density has {int(np.sum(d.electron_density<=0))} non-positive")
        rep.passed = False
    if d.kappa is not None and np.any(d.kappa < 0):
        rep.errors.append(f"kappa has {int(np.sum(d.kappa<0))} negative")
        rep.passed = False
    if d.rosseland_chi is not None and np.any(d.rosseland_chi < 0):
        rep.errors.append(f"Rosseland chi has {int(np.sum(d.rosseland_chi<0))} negative")
        rep.passed = False
    for elem, arr in d.mass_fractions.items():
        if np.any((arr < 0) | (arr > 1)):
            rep.errors.append(f"{elem} mass fraction outside [0,1]")
            rep.passed = False
    for (iso_name, baryon), arr in d.isotopes.items():
        if np.any((arr < 0) | (arr > 1)):
            rep.errors.append(f"{iso_name} {baryon} isotope fraction outside [0,1]")
            rep.passed = False
    if d.density_gcm3 is not None and np.any(d.density_gcm3 <= 0):
        rep.errors.append(f"density has {int(np.sum(d.density_gcm3<=0))} non-positive")
        rep.passed = False
    if d.T_1e4K is not None and np.any(d.T_1e4K <= 0):
        rep.errors.append(f"T has {int(np.sum(d.T_1e4K<=0))} non-positive")
        rep.passed = False

    # Velocity zero -> sigma NaN (VEL_LAW=12 rv_sn_model_sniin L266 divides by V)
    if d.v_kms is not None and np.any(np.abs(d.v_kms) < 1e-30):
        rep.errors.append(f"V has {int(np.sum(np.abs(d.v_kms)<1e-30))} zero points "
                          "(sigma division by zero)")
        rep.passed = False

    # Mass fraction sums
    if d.mass_fractions:
        mf_sum = np.zeros(n)
        for arr in d.mass_fractions.values():
            if len(arr) == n:
                mf_sum += arr
        if np.any(mf_sum <= 0):
            rep.errors.append("mass fractions sum to <= 0 at some depths")
            rep.passed = False
        else:
            max_dev = float(np.max(np.abs(mf_sum - 1.0)))
            if max_dev > 0.05:
                rep.errors.append(f"mass fractions deviate from 1.0 (max={max_dev:.3f})")
                rep.passed = False
            elif max_dev > 0.01:
                rep.warnings.append(f"mass fractions deviate from 1.0 (max={max_dev:.4f})")

    # Sigma clamp band warning
    if d.sigma is not None:
        n_clamp = int(np.sum((d.sigma > -1.2) & (d.sigma < -1.0)))
        if n_clamp:
            rep.warnings.append(f"{n_clamp} sigma points in -1.2..-1.0 (capped to -0.999)")

    # Cell ratio warning
    if len(d.r_1e10cm) >= 3:
        dr = np.abs(np.diff(d.r_1e10cm))
        ratio = np.maximum(dr[:-1] / np.maximum(dr[1:], 1e-30),
                           dr[1:] / np.maximum(dr[:-1], 1e-30))
        if np.any(ratio > 2.0):
            rep.warnings.append(f"max cell ratio {float(np.max(ratio)):.2f} > 2.0")

    return rep
