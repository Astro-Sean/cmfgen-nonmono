"""STELLA .swd snapshot parsing with units and centre/boundary metadata.

Parses STELLA ``*.swd`` files and associated composition files (``.rho``,
``.abn``, ``.xni``) into :class:`SourceShellGrid` objects.

File formats
------------
*.swd* — ASCII, one epoch block per timestep.  Each line is one zone:
    t_day  km  log10(M_ext_Msun)  log10(R_cm)  V_kms
    log10(T_gas_K)  log10(T_rad_K)  log10(rho_gcc)
    log10(P_cgs)  log10(Qv)  log10(eng)  Flum_1e40  WRKX
    [n_bar]  [n_e]  [15 composition values]

*.rho* — final evolved composition (log10 mass fractions)
*.abn* — initial composition (linear mass fractions)
*.xni* — Fortran unformatted, Ni56 mass fraction per zone (inner-to-outer)

All output arrays are **outer-to-inner** (decreasing radius).
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import Optional

import numpy as np

from .shell_model import (
    SourceShellGrid, STELLA_ELEM_NAMES, STELLA_AZ,
    compute_shell_volume, compute_shell_centres, density_from_mass_volume,
    MSUN, AVOGADRO,
)


def normalize_swd_composition(composition: np.ndarray) -> dict[str, np.ndarray]:
    """Normalize STELLA .swd composition from Y_k*A_k to true mass fractions.

    The .swd file writes Yat(ii) = YABUN(ii,Km) * AZ(ii) (number fraction
    times atomic mass).  These sum to <A>, not 1.
    True mass fraction: X_k = Y_k * A_k / <A>.
    """
    comp = np.asarray(composition, dtype=float)
    if comp.ndim != 2 or comp.shape[1] != 15:
        raise ValueError(f"Expected composition shape (nzon, 15), got {comp.shape}")
    total = comp.sum(axis=1, keepdims=True)
    total = np.where(total < 1e-30, 1.0, total)
    mass_fractions = comp / total
    return {name: mass_fractions[:, i] for i, name in enumerate(STELLA_ELEM_NAMES)}


def read_swd(swd_path: str) -> dict[float, dict]:
    """Parse a STELLA *.swd file.

    Returns dict keyed by t_day -> epoch data with arrays sorted outer->inner.
    """
    epochs: dict[float, dict] = {}
    current_t: float | None = None
    buf_km, buf_r, buf_v, buf_T, buf_rho = [], [], [], [], []
    buf_flux, buf_kappa, buf_nbar, buf_ne, buf_comp = [], [], [], [], []
    buf_log_m = []

    def _flush():
        nonlocal current_t, buf_km, buf_r, buf_v, buf_T, buf_rho
        nonlocal buf_flux, buf_kappa, buf_nbar, buf_ne, buf_comp, buf_log_m
        if current_t is not None and buf_r:
            order = np.argsort(buf_r)[::-1]  # outer->inner
            epochs[current_t] = {
                'km':    np.array(buf_km)[order],
                'r':     np.array(buf_r)[order],
                'v':     np.array(buf_v)[order],
                'T':     np.array(buf_T)[order],
                'rho':   np.array(buf_rho)[order],
                'flux':  np.array(buf_flux)[order],
                'kappa': np.array(buf_kappa)[order] if buf_kappa else None,
                'nbar':  np.array(buf_nbar)[order] if buf_nbar else None,
                'ne':    np.array(buf_ne)[order] if buf_ne else None,
                'composition': np.array(buf_comp)[order] if buf_comp else None,
                'log_m': np.array(buf_log_m)[order] if buf_log_m else None,
            }
        buf_km, buf_r, buf_v, buf_T, buf_rho = [], [], [], [], []
        buf_flux, buf_kappa, buf_nbar, buf_ne, buf_comp = [], [], [], [], []
        buf_log_m = []

    with open(swd_path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip()
            if not line.strip():
                continue
            cols = line.split()
            if len(cols) < 13:
                continue
            try:
                t_candidate = float(cols[0])
                km = int(cols[1])
                log_r = float(cols[3])
                v_kms = float(cols[4]) * 1000.0  # .swd unit: 1e8 cm/s -> km/s
                log_T = float(cols[5])
                log_rho = float(cols[7])
                flux_1e40 = float(cols[11])
                kappa_ross = float(cols[12])
            except (ValueError, IndexError):
                # Fixed-width fallback for FORTRAN format concatenation
                try:
                    t_str = line[0:11].strip()
                    t_candidate = float(t_str) if t_str else 0.0
                    km = int(line[11:16])
                    log_r = float(line[28:42])
                    v_kms = float(line[42:51]) * 1000.0
                    log_T = float(line[51:60])
                    log_rho = float(line[69:78])
                    flux_1e40 = float(line[105:116])
                    kappa_ross = float(line[116:127])
                    cols = [line[0:11].strip(), line[11:16].strip(),
                            line[16:28].strip(), line[28:42].strip(),
                            line[42:51].strip(), line[51:60].strip(),
                            line[60:69].strip(), line[69:78].strip(),
                            line[78:87].strip(), line[87:96].strip(),
                            line[96:105].strip(), line[105:116].strip(),
                            line[116:127].strip(), line[127:138].strip(),
                            line[138:149].strip()]
                    for i in range(15):
                        start = 149 + i * 15
                        end = start + 15
                        if end <= len(line):
                            cols.append(line[start:end].strip())
                except (ValueError, IndexError):
                    continue

            if km == 1:
                _flush()
                current_t = t_candidate

            r_cm = 10.0 ** log_r
            T_K = 10.0 ** log_T
            if T_K <= 1.01:
                log_T_gas = float(cols[5])
                T_K = 10.0 ** log_T_gas
            rho = (10.0 ** log_rho) * 1e-6  # .swd rho is log10(rho/1e-6 g/cc)
            flux = flux_1e40 * 1e40
            kappa = kappa_ross

            nbar_val = float(cols[13]) if len(cols) >= 15 else None
            ne_val = float(cols[14]) if len(cols) >= 15 else None

            composition = []
            if len(cols) >= 30:
                for i in range(15):
                    composition.append(float(cols[15 + i]))

            log_m_val = float(cols[2]) if len(cols) >= 3 else None

            buf_km.append(km)
            buf_r.append(r_cm)
            buf_v.append(v_kms)
            buf_T.append(T_K)
            buf_rho.append(rho)
            buf_flux.append(flux)
            buf_kappa.append(kappa)
            if nbar_val is not None:
                buf_nbar.append(nbar_val)
                buf_ne.append(ne_val)
            if composition:
                buf_comp.append(composition)
            if log_m_val is not None:
                buf_log_m.append(log_m_val)

    _flush()
    if not epochs:
        raise ValueError(f"No epochs found in {swd_path}")
    return epochs


def select_epoch(epochs: dict[float, dict], target_days: Optional[float] = None) -> tuple[float, dict]:
    """Select the epoch closest to target_days. Returns (t_day, epoch_data)."""
    available = sorted(epochs.keys())
    if target_days is None:
        t_pick = available[-1]
    elif target_days >= 0:
        positive = [t for t in available if t >= 0]
        if positive:
            t_pick = min(positive, key=lambda t: abs(t - target_days))
            if abs(t_pick - target_days) > max(1.0, 0.1 * (positive[-1] - positive[0])):
                t_pick = positive[-1]
        else:
            t_pick = available[-1]
    else:
        t_pick = min(available, key=lambda t: abs(t - target_days))
    return t_pick, epochs[t_pick]


def read_abn(abn_path: str, nzon: int) -> dict[str, np.ndarray]:
    """Read STELLA *.abn file -> mass fractions per zone (outer->inner)."""
    elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si',
                'S', 'Ar', 'Ca', 'FePeak', 'Ni58', 'Ni56']
    data = {e: [] for e in elements}
    with open(abn_path) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 20 or not cols[0].strip().lstrip('-').isdigit():
                continue
            try:
                vals = [float(cols[c]) for c in range(4, 20)]
            except ValueError:
                continue
            for e, v in zip(elements, vals):
                data[e].append(v)
    for e in elements:
        arr = np.array(data[e])
        if len(arr) < nzon:
            raise ValueError(f"ABN file has {len(arr)} zones, expected {nzon}")
        if len(arr) > nzon:
            arr = arr[:nzon]
        data[e] = arr[::-1]  # inner->outer to outer->inner
    return data


def read_rho(rho_path: str, nzon: int, allow_index_interp: bool = False) -> dict[str, np.ndarray]:
    """Read STELLA .rho file -> final evolved mass fractions (outer->inner).

    Parameters
    ----------
    rho_path : path to the .rho file
    nzon : expected number of zones (from the .swd snapshot)
    allow_index_interp : if True and the .rho zone count differs from nzon,
        fall back to normalized-zone-index interpolation.  This is
        **physically unsafe** (shell positions may be wrong) and is only
        accepted when explicitly requested.  Defaults to False (raise).
    """
    elements = ['Ni', 'H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'S', 'Ar', 'Ca', 'Fe']
    data = {e: [] for e in elements}
    with open(rho_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('zone') or line.startswith('500'):
                continue
            cols = line.split()
            if len(cols) < 22:
                continue
            try:
                vals = [float(cols[c]) for c in range(7, 22)]
            except ValueError:
                continue
            for e, v in zip(elements, vals):
                data[e].append(10**v if v > -10 else 0.0)
    for e in elements:
        arr = np.array(data[e])
        if len(arr) < nzon:
            # Zone-count mismatch: the .rho file has fewer zones than the
            # .swd snapshot.  Interpolating by *normalized zone index*
            # (linspace(0,1)) misplaces composition shells because zone
            # boundaries do not correspond to equal mass/radius fractions.
            # Fail loudly unless the caller explicitly opts in.
            if not allow_index_interp:
                raise ValueError(
                    f".rho file has {len(arr)} zones for element {e} but "
                    f"{nzon} were expected. Zone-count mismatch cannot be "
                    f"resolved safely without a physical coordinate. Pass "
                    f"allow_index_interp=True to accept index-based "
                    f"interpolation (shell positions may be wrong).")
            from scipy.interpolate import interp1d
            r_rho = np.linspace(0, 1, len(arr))
            r_swd = np.linspace(0, 1, nzon)
            arr = interp1d(r_rho, arr, kind='linear', bounds_error=False,
                           fill_value='extrapolate')(r_swd)
        elif len(arr) > nzon:
            arr = arr[:nzon]
        data[e] = arr[::-1]
    return data


def read_xni(xni_path: str, nzon: int) -> np.ndarray:
    """Read STELLA *.xni file -> Ni56 mass fraction per zone (outer->inner).

    Fortran unformatted: [4-byte rec_len] [data] [4-byte rec_len]
    Written inner->outer, reversed to outer->inner.
    """
    with open(xni_path, 'rb') as f:
        data = f.read()
    rec_len = struct.unpack('i', data[:4])[0]
    trailing_len = struct.unpack('i', data[-4:])[0]
    if rec_len != trailing_len:
        raise ValueError(f"XNI record length mismatch: header={rec_len}, trailer={trailing_len}")
    arr = np.frombuffer(data[4:4 + rec_len], dtype=np.float64)
    if len(arr) < nzon:
        raise ValueError(f"XNI file has {len(arr)} zones, expected {nzon}")
    if len(arr) > nzon:
        arr = arr[:nzon]
    return arr[::-1].copy()  # inner->outer to outer->inner, copy to make writable


def build_source_grid(
    epoch_data: dict,
    time_days: float,
    mass_fractions: Optional[dict[str, np.ndarray]] = None,
    ni56_initial: Optional[np.ndarray] = None,
    source_file: str = "",
) -> SourceShellGrid:
    """Build a SourceShellGrid from parsed STELLA epoch data.

    Shell boundaries are constructed from the zone-centre radii by
    placing boundaries at the geometric mean of adjacent zone radii.
    The outermost boundary is r[0] * (1 + dr/(2*r[0])), and the
    innermost boundary is r[-1] * (1 - dr/(2*r[-1])).
    """
    r_centre = np.asarray(epoch_data['r'], dtype=float)
    v_cms = np.asarray(epoch_data['v'], dtype=float) * 1e5  # km/s -> cm/s
    T = np.asarray(epoch_data['T'], dtype=float)
    rho = np.asarray(epoch_data['rho'], dtype=float)
    kappa = np.asarray(epoch_data['kappa'], dtype=float) if epoch_data.get('kappa') is not None else np.full_like(rho, 0.34)
    ne = np.asarray(epoch_data['ne'], dtype=float) if epoch_data.get('ne') is not None else None
    nbar = np.asarray(epoch_data['nbar'], dtype=float) if epoch_data.get('nbar') is not None else None
    log_m = np.asarray(epoch_data['log_m'], dtype=float) if epoch_data.get('log_m') is not None else None

    # Trim trailing zero-density ghost zones (STELLA pads the outer edge with
    # zero-density zones that break conservative remapping and grid selection).
    # These are at the end of the array (innermost in physical space, but
    # outermost in STELLA's outward-to-inward ordering — actually they are
    # the low-velocity outer ghost zones at the array tail).
    n_original = len(r_centre)
    valid_mask = rho > 0
    if not np.all(valid_mask):
        # Find the last contiguous valid zone from the start (outermost)
        # STELLA arrays are outer-to-inner; ghost zones are at the inner end
        # (high index) with rho=0. Trim them.
        last_valid = np.where(valid_mask)[0][-1] if np.any(valid_mask) else -1
        # Also trim leading zeros (outer ghost zones) if any
        first_valid = np.where(valid_mask)[0][0] if np.any(valid_mask) else 0
        if last_valid < n_original - 1 or first_valid > 0:
            n_trimmed = last_valid - first_valid + 1
            print(f"  [build_source_grid] Trimming {n_original - n_trimmed} zero-density ghost zones "
                  f"(indices 0..{first_valid-1} and {last_valid+1}..{n_original-1})")
            r_centre = r_centre[first_valid:last_valid+1]
            v_cms = v_cms[first_valid:last_valid+1]
            T = T[first_valid:last_valid+1]
            rho = rho[first_valid:last_valid+1]
            kappa = kappa[first_valid:last_valid+1]
            if ne is not None:
                ne = ne[first_valid:last_valid+1]
            if nbar is not None:
                nbar = nbar[first_valid:last_valid+1]
            if log_m is not None:
                log_m = log_m[first_valid:last_valid+1]
            if mass_fractions is not None:
                mass_fractions = {k: v[first_valid:last_valid+1] for k, v in mass_fractions.items()}
            if ni56_initial is not None:
                ni56_initial = ni56_initial[first_valid:last_valid+1]
            # Also trim epoch_data composition if present
            if epoch_data.get('composition') is not None:
                comp = epoch_data['composition']
                if isinstance(comp, dict):
                    epoch_data['composition'] = {k: np.asarray(v)[first_valid:last_valid+1] for k, v in comp.items()}

    n = len(r_centre)

    # Construct shell boundaries from zone centres
    # Boundary between zone i and i+1 is at geometric mean
    r_bound = np.empty(n + 1)
    r_bound[1:-1] = np.sqrt(r_centre[:-1] * r_centre[1:])  # between adjacent zones
    # Outer boundary: extrapolate outward
    if n >= 2:
        dr_out = r_centre[0] - r_bound[1]
        r_bound[0] = r_centre[0] + dr_out
        dr_in = r_bound[-2] - r_centre[-1]
        r_bound[-1] = max(r_centre[-1] - dr_in, r_centre[-1] * 0.01)
    else:
        r_bound[0] = r_centre[0] * 1.01
        r_bound[-1] = r_centre[0] * 0.99

    # Ensure strict decrease
    for i in range(len(r_bound) - 1):
        if r_bound[i] <= r_bound[i + 1]:
            r_bound[i] = r_bound[i + 1] * 1.0001

    vol = compute_shell_volume(r_bound)

    # Shell mass from STELLA mass coordinate if available, else from density * volume
    if log_m is not None and len(log_m) == n:
        # log_m = log10(M_ext_Msun) = log10(external mass coordinate in Msun)
        # M_ext is the enclosed mass EXTERIOR to this zone (from outside).
        # Arrays are outer-to-inner, so M_ext increases with index (inward).
        # M_ext[0] = mass exterior to outermost zone (small)
        # M_ext[-1] = mass exterior to innermost zone (total mass)
        # Shell mass[i] = M_ext[i+1] - M_ext[i]  (mass between zones i and i+1)
        # Shell mass[-1] = M_total - M_ext[-1]  (mass interior to innermost zone)
        M_ext = 10.0 ** log_m * MSUN  # grams
        shell_mass = np.empty(n)
        shell_mass[:-1] = M_ext[1:] - M_ext[:-1]  # mass in each shell
        # Innermost shell: use density * volume as approximation
        # (M_ext[-1] is the total enclosed mass, so the innermost shell
        # mass is M_total - M_ext[-1], but we don't know M_total directly.
        # Use the density * volume for the last shell.)
        shell_mass[-1] = rho[-1] * vol[-1]
        shell_mass = np.maximum(shell_mass, 0.0)
        # Fix zones where log_m precision gives zero mass but density is non-zero.
        # This happens in low-mass outer zones where log10 rounding makes adjacent
        # M_ext values identical. Fall back to rho * vol for those zones.
        zero_mass = shell_mass <= 0
        if np.any(zero_mass):
            rho_vol = rho * vol
            n_fixed = int(np.sum(zero_mass & (rho_vol > 0)))
            if n_fixed > 0:
                print(f"  [build_source_grid] Fixing {n_fixed} zero-mass zones with rho*vol fallback "
                      f"(log_m precision loss in low-mass outer zones)")
                shell_mass[zero_mass & (rho_vol > 0)] = rho_vol[zero_mass & (rho_vol > 0)]
    else:
        shell_mass = rho * vol

    # Recompute density from mass/volume for consistency
    density = density_from_mass_volume(shell_mass, vol)

    # Composition
    if mass_fractions is None:
        if epoch_data.get('composition') is not None:
            mass_fractions = normalize_swd_composition(epoch_data['composition'])
        else:
            mass_fractions = {}

    # Baryon density
    baryon_density = nbar if nbar is not None else None

    grid = SourceShellGrid(
        r_bound_cm=r_bound,
        r_centre_cm=r_centre,
        shell_volume_cm3=vol,
        shell_mass_g=shell_mass,
        density_gcm3=density,
        velocity_cms=v_cms,
        temperature_K=T,
        kappa_cgs=kappa,
        mass_fractions=mass_fractions,
        ni56_initial=ni56_initial,
        electron_density=ne,
        baryon_density=baryon_density,
        time_days=time_days,
        source_file=source_file,
        n_zones_original=n,
    )

    # Compute cumulative tau
    from .opacity_tau import compute_tau
    grid.tau_cumulative = compute_tau(r_bound, density, kappa)

    return grid
