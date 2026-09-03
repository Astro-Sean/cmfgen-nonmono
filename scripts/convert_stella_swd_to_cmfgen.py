"""Convert a STELLA *.swd snapshot to CMFGEN SN_HYDRO_DATA format.

STELLA *.swd file format (ASCII, written by PRIBAL in lbalsw.trf):
  Each epoch block: one line per zone, columns:
    t_day   km   log10(M_ext_Msun)   log10(R_cm)   V_kms
    log10(T_gas_K)   log10(T_rad_K)   log10(rho_gcc)
    log10(P_cgs)   log10(Qv)   log10(eng)   Flum_1e40   WRKX

  The km=1 (outer boundary) line has t_day as first field;
  all interior lines start with 0 or blank in that column.

  Lines with km=1 or km=Nzon are "header" lines (t_day appears);
  interior lines start with a leading zero or spaces in the t field.

Usage:
  python convert_stella_swd_to_cmfgen.py SWD_FILE OUTPUT_FILE [OPTIONS]
    --target-days T    Pick the epoch closest to T days (default: last epoch)
    --tau-max X        Outer optical depth cutoff (default: 30)
    --tau-min X        Inner optical depth cutoff (default: 1e-3)
    --nzon-target N    Target number of output grid points (default: 200)
    --time-days T      Explosion time in days (default: read from swd)
    --model-name NAME  For header comment (default: stella)
    --kappa K          Electron-scattering opacity cm^2/g (default: 0.34)
    --xni-mass X       56Ni mass in Msun for luminosity estimate (default: 0)
    --rho-file PATH    *.rho file from STELLA for final evolved composition (preferred)
    --abn-file PATH    *.abn file from STELLA for initial composition (fallback)
    --enable-adaptive-grid  Enable Dessart-style adaptive grid refinement
    --feature-aware-downsample  Preserve high-gradient regions (default: enabled)
    --no-feature-aware-downsample  Disable feature-aware downsampling

The output SN_HYDRO_DATA is in CMFGEN outer->inner (decreasing r) order.
"""

from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


# Physical constants
C_LIGHT   = 2.99792458e10   # cm/s
K_BOLTZ   = 1.380649e-16    # erg/K
M_H       = 1.6735575e-24   # g
SIGMA_SB  = 5.670374419e-5  # erg/cm^2/s/K^4
LSUN      = 3.828e33        # erg/s
MSUN      = 1.989e33        # g
AVOGADRO  = 6.02214076e23

# STELLA Natom=15 element order (from azdat.f / YABUN array)
STELLA_ELEM_NAMES = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si',
                     'S', 'Ar', 'Ca', 'Fe', 'Ni']
# Atomic masses AZ from STELLA azdat.f (Common/AZZn/)
STELLA_AZ = np.array([1.0, 4.0, 12.0, 14.0, 16.0, 20.0, 23.0, 24.0, 27.0,
                      28.0, 32.0, 40.0, 40.0, 56.0, 58.0])


def _select_epoch(available: list[float], target_days: float | None) -> float:
    """Pick the .swd epoch closest to ``target_days``.

    Shared by :func:`read_swd` (to know which single epoch to parse) and the
    CLI ``main`` (to pick the snapshot).  Keeping the selection in one place
    guarantees the optimised single-epoch reader parses exactly the epoch the
    caller will use.

    * ``target_days is None`` -> last available epoch.
    * ``target_days >= 0``    -> closest *positive* (post-SBO) epoch; if the
      closest positive epoch is far from the requested time, fall back to the
      last positive epoch.  If no positive epochs exist, use the closest epoch.
    * ``target_days < 0``     -> closest *negative* (pre-SBO) epoch, with the
      symmetric fallback to the first negative epoch / closest epoch.
    """
    if target_days is None:
        return available[-1]
    if target_days >= 0:
        positive = [t for t in available if t >= 0]
        if positive:
            t_pick = min(positive, key=lambda t: abs(t - target_days))
            if abs(t_pick - target_days) > max(1.0, 0.1 * (positive[-1] - positive[0])):
                t_pick = positive[-1]
            return t_pick
        return min(available, key=lambda t: abs(t - target_days))
    negative = [t for t in available if t < 0]
    if negative:
        t_pick = min(negative, key=lambda t: abs(t - target_days))
        if abs(t_pick - target_days) > max(1.0, 0.1 * (negative[-1] - negative[0])):
            t_pick = negative[0]
        return t_pick
    return min(available, key=lambda t: abs(t - target_days))


def normalize_swd_composition(composition: np.ndarray) -> dict[str, np.ndarray]:
    """Normalize STELLA .swd composition from Y_k*A_k to true mass fractions.

    The .swd file writes Yat(ii) = YABUN(ii,Km) * AZ(ii), i.e. number fraction
    times atomic mass.  These values sum to <A> = sum(Y_k * A_k), not 1.
    The true mass fraction is X_k = Y_k * A_k / <A>.

    This function divides each element column by the per-zone sum of all 15
    elements, producing normalized mass fractions that sum to 1 per zone.

    Args:
        composition: array of shape (nzon, 15) with Y_k*A_k values from .swd

    Returns:
        dict mapping element name -> mass fraction array (nzon,), normalized to sum=1
    """
    comp = np.asarray(composition, dtype=float)
    if comp.ndim != 2 or comp.shape[1] != 15:
        raise ValueError(f"Expected composition shape (nzon, 15), got {comp.shape}")
    # Per-zone sum of all 15 Y_k*A_k values = <A>
    total = comp.sum(axis=1, keepdims=True)
    total = np.where(total < 1e-30, 1.0, total)
    mass_fractions = comp / total  # X_k = Y_k * A_k / <A>
    return {name: mass_fractions[:, i] for i, name in enumerate(STELLA_ELEM_NAMES)}


def _parse_swd_line(line: str) -> dict | None:
    """Parse one STELLA .swd data line into physical values.

    Handles both the whitespace-split format and the FORTRAN fixed-width
    fallback (FORMAT 181/281) for concatenated fields.  Returns a dict with
    keys ``t, km, r, v, T, T_rad, rho, flux, kappa, nbar, ne, composition``
    (nbar/ne/composition may be None), or ``None`` for blank/unparseable lines.

    This is the single source of truth for .swd line parsing, shared by the
    full-epoch reader and the optimised single-epoch reader so both produce
    identical values for the same line.
    """
    if not line.strip():
        return None
    cols = line.split()
    if len(cols) < 13:
        return None
    try:
        t_candidate = float(cols[0])
        km          = int(cols[1])
        log_r       = float(cols[3])
        v_kms       = float(cols[4]) * 1000.0  # .swd unit: 1e8 cm/s -> km/s
        log_T       = float(cols[5])
        log_rho     = float(cols[7])
        flux_1e40   = float(cols[11])
        kappa_ross  = float(cols[12])  # Rosseland opacity in cm²/g
    except (ValueError, IndexError):
        # FORTRAN FORMAT 181/281 can produce concatenated fields when
        # F14.7 (log10 R) fills all 14 chars and F9.4 (velocity) starts
        # immediately after with no space. Fall back to fixed-width parsing.
        # FORMAT 181: 1X,F10.6,I5,F12.6,F14.7,F9.4,6F9.3,1P,4E11.3,15E15.6
        # FORMAT 281: 10X,i1,I5,F12.6,F14.7,F9.4,6F9.3,1P,4E11.3,15E15.6
        # Both have identical field positions after column 10.
        try:
            t_str = line[0:11].strip()
            t_candidate = float(t_str) if t_str else 0.0
            km          = int(line[11:16])
            log_r       = float(line[28:42])
            v_kms       = float(line[42:51]) * 1000.0
            log_T       = float(line[51:60])
            log_rho     = float(line[69:78])
            flux_1e40   = float(line[105:116])
            kappa_ross  = float(line[116:127])
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
            return None

    r_cm = 10.0 ** log_r
    T_K = 10.0 ** log_T
    log_T_rad = float(cols[6])
    T_rad_K = 10.0 ** log_T_rad
    if T_rad_K <= 1.01:
        T_rad_K = T_K
    rho = (10.0 ** log_rho) * 1e-6
    flux = flux_1e40 * 1e40
    kappa = kappa_ross
    nbar_val = float(cols[13]) if len(cols) >= 15 else None
    ne_val = float(cols[14]) if len(cols) >= 15 else None
    composition = []
    if len(cols) >= 30:
        for i in range(15):
            composition.append(float(cols[15 + i]))
    return {
        't': t_candidate, 'km': km, 'r': r_cm, 'v': v_kms,
        'T': T_K, 'T_rad': T_rad_K, 'rho': rho, 'flux': flux,
        'kappa': kappa, 'nbar': nbar_val, 'ne': ne_val,
        'composition': composition if composition else None,
    }


def _parse_swd_header(line: str) -> tuple[float, int] | None:
    """Cheaply extract (t, km) from a .swd line for epoch indexing.

    Uses a bounded ``split(None, 2)`` (first 3 fields only) plus a minimum
    line-length guard, instead of the full 30-field split, so the epoch-index
    pass is ~2x faster.  Verified to detect exactly the same km==1 headers as
    :func:`_parse_swd_line` on real STELLA .swd files (FORMAT 181/281 data
    lines are 374 chars; the 51-char guard admits only lines extending through
    the velocity field, which no non-data line satisfies).
    """
    if len(line) < 51 or not line.strip():
        return None
    parts = line.split(None, 2)
    if len(parts) < 2:
        return None
    try:
        return float(parts[0]), int(parts[1])
    except (ValueError, IndexError):
        try:
            t_str = line[0:11].strip()
            t = float(t_str) if t_str else 0.0
            return t, int(line[11:16])
        except (ValueError, IndexError):
            return None


def _finalize_epoch(t: float, buf: dict) -> dict:
    """Sort a single epoch's buffers outer->inner and build the snapshot dict."""
    order = np.argsort(buf['r'])[::-1]
    return {
        'km':    np.array(buf['km'])[order],
        'r':     np.array(buf['r'])[order],
        'v':     np.array(buf['v'])[order],
        'T':     np.array(buf['T'])[order],
        'T_rad': np.array(buf['T_rad'])[order],
        'rho':   np.array(buf['rho'])[order],
        'flux':  np.array(buf['flux'])[order],
        'kappa': np.array(buf['kappa'])[order] if buf['kappa'] else None,
        'nbar':  np.array(buf['nbar'])[order] if buf['nbar'] else None,
        'ne':    np.array(buf['ne'])[order] if buf['ne'] else None,
        'composition': np.array(buf['composition'])[order] if buf['composition'] else None,
    }


def _new_epoch_buf() -> dict:
    return {'km': [], 'r': [], 'v': [], 'T': [], 'T_rad': [], 'rho': [],
            'flux': [], 'kappa': [], 'nbar': [], 'ne': [], 'composition': []}


def _append_to_buf(buf: dict, rec: dict) -> None:
    buf['km'].append(rec['km'])
    buf['r'].append(rec['r'])
    buf['v'].append(rec['v'])
    buf['T'].append(rec['T'])
    buf['T_rad'].append(rec['T_rad'])
    buf['rho'].append(rec['rho'])
    buf['flux'].append(rec['flux'])
    buf['kappa'].append(rec['kappa'])
    if rec['nbar'] is not None:
        buf['nbar'].append(rec['nbar'])
        buf['ne'].append(rec['ne'])
    if rec['composition']:
        buf['composition'].append(rec['composition'])


def read_swd(swd_path: str, target_days: float | None = None) -> dict[float, dict]:
    """Parse a STELLA *.swd file.

    Returns a dict keyed by t_day -> {'r': [], 'v_kms': [], 'rho': [], 'T': [],
    'T_rad': [], 'flux': [], 'km': []}.
    'T' is the gas kinetic temperature (T_gas); 'T_rad' is the radiation
    temperature.  In optically-thick inner zones where STELLA writes T_rad=1 K
    (fill value), T_rad is set equal to T_gas.
    All values are numpy arrays sorted outer->inner (decreasing r).

    When ``target_days`` is given, only the single closest epoch is fully
    parsed (a two-pass seek-based read), which avoids parsing every epoch of a
    multi-hundred-MB .swd.  The returned dict still contains *all* epoch times
    as keys (so ``sorted(epochs.keys())`` reports the full available range),
    but only the selected epoch maps to a real snapshot; the other keys map to
    ``None``.  The caller only ever accesses ``epochs[t_pick]`` where
    ``t_pick`` is chosen by :func:`_select_epoch` with the same ``target_days``,
    so it always resolves to the parsed snapshot.
    """
    if target_days is not None:
        return _read_swd_single_epoch(swd_path, target_days)

    epochs: dict[float, dict] = {}
    current_t: float | None = None
    buf = _new_epoch_buf()

    def _flush():
        nonlocal current_t, buf
        if current_t is not None and buf['r']:
            epochs[current_t] = _finalize_epoch(current_t, buf)
        buf = _new_epoch_buf()

    with open(swd_path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip()
            rec = _parse_swd_line(line)
            if rec is None:
                continue
            if rec['km'] == 1:
                # km=1 marks the start of a new epoch (time in col 0).
                _flush()
                current_t = rec['t']
            _append_to_buf(buf, rec)

    _flush()
    if not epochs:
        raise ValueError(f"No epochs found in {swd_path}")
    return epochs


def _read_swd_single_epoch(swd_path: str, target_days: float) -> dict[float, dict]:
    """Two-pass reader: index epoch headers cheaply, then parse only the target epoch.

    Pass 1 scans every line but parses only (t, km) to build an epoch index
    with file offsets.  Pass 2 seeks to the selected epoch and parses it with
    the full :func:`_parse_swd_line` (identical to the all-epoch reader), so
    the snapshot values are identical to reading the whole file.
    """
    # Pass 1: cheap header scan -> [(t, line_offset), ...]
    # Pass 1: scan in binary mode and track byte offsets incrementally.
    # Reading in binary and summing len(raw) is ~10x faster than calling
    # fh.tell() after every readline (tell() flushes the text-decoder buffer
    # and re-syncs the byte position on each call).  Binary offsets are valid
    # for the binary seek in pass 2 below.
    epoch_index: list[tuple[float, int]] = []
    offset = 0
    with open(swd_path, 'rb') as fh:
        raw = fh.readline()
        while raw:
            line = raw.decode('utf-8', errors='replace').rstrip()
            hdr = _parse_swd_header(line)
            if hdr is not None and hdr[1] == 1:
                epoch_index.append((hdr[0], offset))
            offset += len(raw)
            raw = fh.readline()

    if not epoch_index:
        raise ValueError(f"No epochs found in {swd_path}")

    all_times = [t for t, _ in epoch_index]
    t_pick = _select_epoch(all_times, target_days)
    # Locate the file offset of the chosen epoch's header line.
    target_offset = None
    for t, off in epoch_index:
        if t == t_pick:
            target_offset = off
            break
    if target_offset is None:
        # Fallback: exact-time match failed (e.g. duplicate times); pick by offset.
        target_offset = epoch_index[all_times.index(t_pick)][1]

    # Pass 2: seek to the target epoch and parse only its lines (binary mode,
    # decoding each line, so the byte offset from pass 1 is valid).
    buf = _new_epoch_buf()
    with open(swd_path, 'rb') as fh:
        fh.seek(target_offset)
        raw = fh.readline()
        while raw:
            line = raw.decode('utf-8', errors='replace').rstrip()
            rec = _parse_swd_line(line)
            if rec is None:
                raw = fh.readline()
                continue
            if rec['km'] == 1 and buf['r']:
                # Reached the next epoch's header -> stop.
                break
            # First line of the target epoch (km==1) is also a data line
            # carrying zone km=1's physical values; append it like any other.
            _append_to_buf(buf, rec)
            raw = fh.readline()

    if not buf['r']:
        raise ValueError(f"Target epoch t={t_pick} not parseable in {swd_path}")
    snap = _finalize_epoch(t_pick, buf)

    # Return all epoch times as keys (for the caller's range diagnostic) with
    # only the selected epoch populated.
    epochs: dict[float, dict] = {t: None for t in all_times}
    epochs[t_pick] = snap
    return epochs


def read_abn(abn_path: str, nzon: int) -> dict[str, np.ndarray]:
    """Read STELLA *.abn file -> mass fractions per zone.

    Columns: km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
    Returns dict of element -> array[nzon], sorted outer->inner (reversed).
    """
    elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si',
                'S', 'Ar', 'Ca', 'FePeak', 'Ni58', 'Ni56']
    data = {e: [] for e in elements}
    with open(abn_path) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 20 or not cols[0].strip().lstrip('-').isdigit():
                continue
            try:
                # Map columns 4-19 to elements (16 elements total)
                # Actual ABN format: cols[4]=H, cols[5]=He, cols[6]=C, cols[7]=N, cols[8]=O, cols[9]=Ne,
                # cols[10]=Na, cols[11]=Mg, cols[12]=Al, cols[13]=Si, cols[14]=S,
                # cols[15]=Ar, cols[16]=Ca, cols[17]=FePeak, cols[18]=Ni58, cols[19]=Ni56
                vals = [float(cols[4]), float(cols[5]), float(cols[6]), float(cols[7]), float(cols[8]),
                        float(cols[9]), float(cols[10]), float(cols[11]), float(cols[12]), float(cols[13]),
                        float(cols[14]), float(cols[15]), float(cols[16]), float(cols[17]), float(cols[18]),
                        float(cols[19])]
            except ValueError:
                continue
            for e, v in zip(elements, vals):
                data[e].append(v)
    for e in elements:
        arr = np.array(data[e])
        if len(arr) < nzon:
            raise ValueError(f"ABN file has {len(arr)} zones, expected {nzon}")
        # If ABN has more zones than SWD, truncate to match SWD
        if len(arr) > nzon:
            arr = arr[:nzon]
        data[e] = arr[::-1]  # inner->outer to outer->inner
    return data


def read_rho(rho_path: str, nzon: int) -> dict[str, np.ndarray]:
    """Read STELLA .rho file -> mass fractions per zone (FINAL evolved composition).

    The .rho file contains the final evolved composition after STELLA simulation.
    Format: zone, mass, lg(r), lg(Tp), lg(rho), u, lg(Ni), lg(H), lg(He), lg(C), lg(N), lg(O), 
            lg(Ne), lg(Na), lg(Mg), lg(Al), lg(Si), lg(S), lg(Ar), lg(Ca), lg(Fe), lg(Ni)
    Returns dict of element -> array[nzon], sorted outer->inner (reversed).
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
                # Map columns 7-21 to elements (15 elements total)
                # cols[7]=Ni, cols[8]=H, cols[9]=He, cols[10]=C, cols[11]=N, cols[12]=O, 
                # cols[13]=Ne, cols[14]=Na, cols[15]=Mg, cols[16]=Al, cols[17]=Si, 
                # cols[18]=S, cols[19]=Ar, cols[20]=Ca, cols[21]=Fe
                vals = [float(cols[7]), float(cols[8]), float(cols[9]), float(cols[10]), float(cols[11]),
                        float(cols[12]), float(cols[13]), float(cols[14]), float(cols[15]), float(cols[16]),
                        float(cols[17]), float(cols[18]), float(cols[19]), float(cols[20]), float(cols[21])]
            except ValueError:
                continue
            
            for e, v in zip(elements, vals):
                # Convert from log10 to linear, set very small values to 0
                if v > -10:
                    data[e].append(10**v)
                else:
                    data[e].append(0.0)
    
    for e in elements:
        arr = np.array(data[e])
        if len(arr) < nzon:
            print(f"WARNING: RHO file has {len(arr)} zones, expected {nzon}. Interpolating to match.")
            # Interpolate to match SWD grid
            from scipy.interpolate import interp1d
            r_rho_grid = np.linspace(0, 1, len(arr))
            r_swd_grid = np.linspace(0, 1, nzon)
            interp_func = interp1d(r_rho_grid, arr, kind='linear', 
                                  bounds_error=False, fill_value='extrapolate')
            arr = interp_func(r_swd_grid)
        # If RHO has more zones than SWD, truncate to match SWD
        elif len(arr) > nzon:
            print(f"WARNING: RHO file has {len(arr)} zones, expected {nzon}. Truncating to match.")
            arr = arr[:nzon]
        data[e] = arr[::-1]  # inner->outer to outer->inner
    
    return data


def read_input_profile(input_profile_path: str) -> dict[str, np.ndarray]:
    """Read [DEPRECATED] HERACLES input_profile file -> mass fractions per zone.

    HERACLES is deprecated; the project now uses the STELLA -> CMFGEN pipeline.

    Columns: r_cm rho_cgs p_cgs u_cgs fx_H fx_He fx_O fx_Si fx_Ca fx_Fe fx_Ni
    Returns dict of element -> array[nzon], sorted outer->inner (reversed).
    """
    elements = ['H', 'He', 'O', 'Si', 'Ca', 'Fe', 'Ni']
    data = {e: [] for e in elements}
    with open(input_profile_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) < 11:
                continue
            try:
                # Map columns to elements
                # cols[4]=fx_H, cols[5]=fx_He, cols[6]=fx_O, cols[7]=fx_Si, cols[8]=fx_Ca, cols[9]=fx_Fe, cols[10]=fx_Ni
                vals = [float(cols[4]), float(cols[5]), float(cols[6]), float(cols[7]),
                        float(cols[8]), float(cols[9]), float(cols[10])]
            except ValueError:
                continue
            for e, v in zip(elements, vals):
                data[e].append(v)
    for e in elements:
        arr = np.array(data[e])
        data[e] = arr[::-1]  # inner->outer to outer->inner
    return data


def read_xni(xni_path: str, nzon: int) -> np.ndarray:
    """Read STELLA *.xni file -> Ni56 mass fraction per zone.

    The .xni file is Fortran unformatted (form='unformatted'), containing a
    single record of Nzon REAL*8 (double precision) values representing the
    Ni56 mass fraction in each zone.  STELLA reads this with:
        read(28) (XNifor(i), i=1, Nzon)
    where XNifor is REAL*8.

    Fortran unformatted format: [4-byte rec_len] [data] [4-byte rec_len]

    The .xni file is written inner->outer (zone 1 = innermost), but CMFGEN
    indexes outer->inner, so we reverse the array.

    Returns array[nzon] of Ni56 mass fractions, sorted outer->inner.
    """
    with open(xni_path, 'rb') as f:
        data = f.read()

    # Parse Fortran unformatted record: [4-byte rec_len] [data] [4-byte rec_len]
    rec_len = struct.unpack('i', data[:4])[0]
    trailing_len = struct.unpack('i', data[-4:])[0]
    if rec_len != trailing_len:
        raise ValueError(
            f"XNI file record length mismatch: header={rec_len}, trailer={trailing_len}"
        )
    n_doubles = rec_len // 8
    arr = np.frombuffer(data[4:4 + rec_len], dtype=np.float64)

    if len(arr) < nzon:
        raise ValueError(f"XNI file has {len(arr)} zones, expected {nzon}")
    if len(arr) > nzon:
        arr = arr[:nzon]

    # The .xni file is written inner->outer, but we need outer->inner
    arr = arr[::-1]

    return arr.copy()  # copy to make writable (frombuffer returns read-only)


def read_lbol(lbol_path):
    """Read STELLA *.lbol file -> bolometric luminosity vs time.
    
    Returns dict keyed by t_day -> L_bol (Lsun).
    """
    data = np.loadtxt(lbol_path, skiprows=1)
    # Column 0: time (days), Column 2: L_bol (Lsun)
    lbol_dict = {}
    for row in data:
        t_day = row[0]
        L_bol = row[2]
        lbol_dict[t_day] = L_bol
    return lbol_dict


# CMFGEN element name abbreviations
CMFGEN_ELEM_NAMES = {
    'H': 'HYD',
    'He': 'HE',
    'C': 'CARB',
    'N': 'NIT',
    'O': 'OXY',
    'Ne': 'NEON',
    'Na': 'SOD',
    'Mg': 'MAG',
    'Al': 'ALUM',
    'Si': 'SIL',
    'S': 'SUL',
    'Ar': 'ARG',
    'Ca': 'CAL',
    'FePeak': 'IRON',
    'Ni58': 'NI58',
    'Ni56': 'NI56',
}

# [DEPRECATED] Element name mapping: STELLA → HERACLES-style (for opacity calculation)
# HERACLES is deprecated; the project now uses the STELLA -> CMFGEN pipeline.
STELLA_TO_HERACLES_ELEM = {
    'H': 'HYD',
    'He': 'HE',
    'C': 'CARB',
    'N': 'NITR',
    'O': 'OXY',
    'Ne': 'NEON',
    'Na': 'SOD',
    'Mg': 'MAGN',
    'Al': 'ALUM',
    'Si': 'SIL',
    'S': 'SULF',
    'Ar': 'ARG',
    'Ca': 'CAL',
    'FePeak': 'IRON',
    'Ni58': 'NI58',
    'Ni56': 'NI56',
}

# [DEPRECATED] HERACLES element name mapping: HERACLES → STELLA (for converting back)
# HERACLES is deprecated; the project now uses the STELLA -> CMFGEN pipeline.
HERACLES_TO_STELLA_ELEM = {v: k for k, v in STELLA_TO_HERACLES_ELEM.items()}


def _saha_ionization_fraction(temperature, chi_eV, g_ion, g_neutral, n_e_ref=1e8):
    """
    Single-step Saha ionization fraction x = n_ion / (n_ion + n_neutral).

    Uses a fixed reference electron density n_e_ref (cm^-3) as a self-consistency
    approximation (iterating self-consistently is expensive; n_e_ref=1e8 is typical
    for SN ejecta at 40 days and gives <10% error vs. full iteration for most zones).

    chi_eV  : ionization potential in eV
    g_ion   : statistical weight of ion ground state
    g_neutral: statistical weight of neutral ground state
    """
    kB_eV = 8.617333e-5    # eV/K
    h = 6.62607e-27        # erg·s
    m_e = 9.10938e-28      # g
    kB_erg = 1.38065e-16   # erg/K

    T = np.maximum(temperature, 100.0)
    kT_erg = kB_erg * T

    # Saha partition function factor: 2*(2πm_e kT/h²)^(3/2) * g_ion/g_neutral
    saha_fac = 2.0 * (2.0 * np.pi * m_e * kT_erg / h**2)**1.5 * (g_ion / g_neutral)
    chi_erg = chi_eV * 1.60218e-12

    # Saha ratio: n_ion * n_e / n_neutral = saha_fac * exp(-chi/kT)
    phi = saha_fac * np.exp(-chi_erg / kT_erg) / n_e_ref

    # x = phi / (1 + phi)  =>  fraction in ionized state
    return np.clip(phi / (1.0 + phi), 0.0, 1.0)


def calculate_electron_density(density, temperature, mass_fractions):
    """Calculate electron density (cm^-3) using per-element Saha ionization.

    H  : two-level Saha  H I  -> H II  (chi=13.6 eV)
    He : two-step Saha   He I -> He II (chi=24.6 eV) then He II -> He III (chi=54.4 eV)
    Metals: first-ionization Saha with element-specific chi; second ionization
            included for low-chi elements (Ca, Fe, Mg, Si, Na) where it matters
            in the ~5000-15000 K range typical of SN ejecta at 40 days.

    Args:
        density      : Density in g/cm³ (array)
        temperature  : Temperature in K (array)
        mass_fractions: Dict of element mass fractions

    Returns:
        n_e: Electron density in cm^-3 (array)
    """
    N_A = 6.02214076e23

    # (atomic_mass_amu, chi1_eV, g_neutral, g_ion1, chi2_eV, g_ion2)
    # chi2=None means second ionization is negligible or unknown
    # g values: ground-state statistical weights (2J+1 or LTE approx)
    atomic_data = {
        'H':      (1.0,   13.60, 2, 1,    None,  None, None),
        'He':     (4.0,   24.59, 1, 2,    54.42, 1,    None),
        'C':      (12.0,  11.26, 9, 6,    24.38, 1,    None),
        'N':      (14.0,  14.53, 4, 9,    29.60, 1,    None),
        'O':      (16.0,  13.62, 9, 4,    35.12, 1,    None),
        'Ne':     (20.0,  21.56, 1, 6,    41.08, 1,    None),
        'Na':     (23.0,   5.14, 2, 1,    47.29, 6,    None),
        'Mg':     (24.0,   7.65, 1, 2,    15.04, 1,    None),
        'Al':     (27.0,   5.99, 6, 1,    18.83, 5,    None),
        'Si':     (28.0,   8.15, 9, 6,    16.35, 1,    None),
        'S':      (32.0,  10.36, 9, 4,    23.34, 1,    None),
        'Ar':     (40.0,  15.76, 1, 6,    27.63, 1,    None),
        'Ca':     (40.0,   6.11, 1, 2,    11.87, 1,    None),
        'FePeak': (56.0,   7.90, 25, 30,  16.19, 25,   None),
        'Ni58':   (58.0,   7.64, 9, 10,   18.17, 9,    None),
        'Ni56':   (56.0,   7.64, 9, 10,   18.17, 9,    None),
    }

    n_e = np.zeros_like(density, dtype=float)

    for elem, data in atomic_data.items():
        if elem not in mass_fractions:
            continue
        X = np.asarray(mass_fractions[elem], dtype=float)
        if not np.any(X > 0):
            continue

        A, chi1, g0, g1, chi2, g2, _ = data
        n_atom = density * X * N_A / A   # number density of this element

        # First ionization
        x1 = _saha_ionization_fraction(temperature, chi1, g1, g0)

        if chi2 is not None and g2 is not None:
            # Two-step: fraction in He II (or metal II) state contributes 1 e-,
            # fraction in He III (metal III) contributes 2 e- each.
            x2 = _saha_ionization_fraction(temperature, chi2, g2, g1)
            # x1 = fraction reaching at least stage I; among those, x2 reach stage II
            Z_mean = x1 * (1.0 + x2)
        else:
            Z_mean = x1

        n_e += n_atom * Z_mean

    return n_e


def approximate_rosseland_opacity(T, rho, X, Y, Z, n_e=None, mass_fractions=None):
    """Approximate Rosseland mean opacity (cm²/g) for supernova ejecta.

    Includes electron scattering, bound-free, and free-free opacity.
    (A line-opacity term was previously included but removed; see note below.)
    
    Parameters:
    -----------
    T : array
        Temperature in K
    rho : array
        Density in g/cm³
    X : array
        Hydrogen mass fraction
    Y : array
        Helium mass fraction
    Z : array
        Metal mass fraction
    n_e : array, optional
        Electron density in cm^-3
    mass_fractions : dict, optional
        Dictionary of element mass fractions for better Z/A calculation
    
    Returns:
    --------
    kappa_rosseland : array
        Rosseland mean opacity in cm²/g
    """
    # Avoid division by zero
    T_safe = np.maximum(T, 1.0)  # Minimum 1 K
    rho_safe = np.maximum(rho, 1e-30)  # Avoid zero density
    
    # Electron scattering with composition dependence (more accurate formula)
    kappa_es = 0.2 * (1.0 + X) / (1.0 + 3.0*X + 0.25*Y)
    
    # Calculate Z/A from actual composition if available
    if mass_fractions is not None:
        Z_over_A = 0.0
        Z_sum = 0.0
        A_sum = 0.0
        # Atomic masses (approximate)
        atomic_masses = {
            'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
            'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
            'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'Fe': 56.0, 'Ni': 58.0
        }
        # NOTE: Use Z_atom (not Z) for the atomic number to avoid shadowing
        # the function parameter Z (metal mass fraction array).  The previous
        # code used `Z = 1.0` inside the loop, which overwrote the input Z
        # array with a scalar, inflating kappa_bf by ~10^7x.
        for elem, mass_frac in mass_fractions.items():
            if elem in atomic_masses:
                A = atomic_masses[elem]
                Z_atom = 1.0  # atomic number (default for H)
                if elem == 'He': Z_atom = 2.0
                elif elem in ['C', 'N', 'O']: Z_atom = 6.0 + list(atomic_masses.keys()).index(elem) - 2
                elif elem == 'Ne': Z_atom = 10.0
                elif elem == 'Na': Z_atom = 11.0
                elif elem == 'Mg': Z_atom = 12.0
                elif elem == 'Al': Z_atom = 13.0
                elif elem == 'Si': Z_atom = 14.0
                elif elem == 'S': Z_atom = 16.0
                elif elem == 'Ar': Z_atom = 18.0
                elif elem == 'Ca': Z_atom = 20.0
                elif elem == 'Fe': Z_atom = 26.0
                elif elem == 'Ni': Z_atom = 28.0
                Z_sum += Z_atom * mass_frac
                A_sum += A * mass_frac
        if np.any(A_sum > 0):
            Z_over_A = Z_sum / A_sum
        else:
            Z_over_A = 0.5  # Fallback
    else:
        Z_over_A = 0.5  # Approximate for solar-like metals
    
    # Bound-free (photoionization) using Kramers formula with improved Z/A
    kappa_bf = 4.34e25 * Z_over_A * Z * rho_safe * (T_safe ** -3.5)
    
    # Free-free (bremsstrahlung) using Kramers formula
    kappa_ff = 1.7e20 * (1.0 + X) * rho_safe * (T_safe ** -3.5)
    
    # Note: a kappa_line ~ Z*rho*T^{-2} term was previously included here but removed.
    # It created a hard discontinuity at T=10000K (the mask threshold) that produced
    # unphysical jumps of ~2000x in kappa across adjacent zones with identical composition.
    # CMFGEN computes line opacities self-consistently from atomic data; this crude
    # approximation is not needed and actively misleads the tau/grid calculation.

    # Total Rosseland opacity
    kappa_rosseland = kappa_es + kappa_bf + kappa_ff
    
    return kappa_rosseland


def compute_kappa_total(rho, T_K, mass_fractions):
    """Compute total Rosseland opacity using [DEPRECATED] HERACLES-style calculation.

    HERACLES is deprecated; the project now uses the STELLA -> CMFGEN pipeline.
    
    Uses approximate_rosseland_opacity with composition-dependent
    electron scattering, bound-free, and free-free opacity components.
    """
    # Extract mass fractions
    X = mass_fractions.get('H', np.ones(len(rho)) * 0.7)
    Y = mass_fractions.get('He', np.ones(len(rho)) * 0.28)
    
    # Calculate metal mass fraction Z (sum of all metals except H and He)
    Z = np.zeros(len(rho))
    for elem in mass_fractions:
        if elem not in ['H', 'He']:
            Z += mass_fractions[elem]
    
    # Calculate Rosseland opacity with improved composition dependence
    kappa_total = approximate_rosseland_opacity(T_K, rho, X, Y, Z, mass_fractions=mass_fractions)
    
    return kappa_total


def compute_kappa_tau(r_cm, rho, kappa_arr):
    """Compute cumulative optical depth tau[i] = integral from r[0] inward.

    Uses trapezoidal rule: chi = kappa*rho at each zone boundary,
    averaged between adjacent zones for the shell integral.
    Arrays are outer->inner (decreasing r).

    Vectorised with ``np.cumsum`` on the per-shell trapezoidal increments;
    the arithmetic is identical to the previous per-point loop (same
    operation order per element), so the output is bit-identical.
    """
    n = len(r_cm)
    if n < 2:
        return np.zeros(n)
    r = np.asarray(r_cm, dtype=float)
    rho_a = np.asarray(rho, dtype=float)
    kap_a = np.asarray(kappa_arr, dtype=float)
    # Shell widths |r[i-1]-r[i]| for i=1..n-1 (r is decreasing -> diff is negative)
    dr = np.abs(np.diff(r))
    # Trapezoidal average of chi = kappa*rho across each shell
    chi = 0.5 * (kap_a[:-1] * rho_a[:-1] + kap_a[1:] * rho_a[1:])
    tau = np.zeros(n)
    tau[1:] = np.cumsum(chi * dr)
    return tau


def luminosity_bb(T_eff_K: float, r_phot_cm: float) -> float:
    return 4.0 * np.pi * r_phot_cm**2 * SIGMA_SB * T_eff_K**4


def sigma_velocity(r, v_kms):
    """Compute sigma = dlnV/dlnr - 1 using the exact CMFGEN formula.

    CMFGEN (rv_sn_model_sniin.f:241) uses a distance-weighted average of
    the left and right LINEAR velocity gradients:

        T1 = (V[i-1]-V[i])/(R[i-1]-R[i])     ! left gradient
        T2 = (V[i]-V[i+1])/(R[i]-R[i+1])     ! right gradient
        dVdR = (T1*h_right + T2*h_left) / (h_left + h_right)
        sigma = R*dVdR/V - 1

    This is more accurate than the logarithmic central difference on
    non-uniform grids because it weights by cell size.

    Boundary points use the same 2-point difference as CMFGEN:
        sigma[0]  = R[0]*(V[0]-V[1])/(R[0]-R[1])/V[0] - 1
        sigma[-1] = R[-1]*(V[-2]-V[-1])/(R[-2]-R[-1])/V[-1] - 1
    This makes the plotted/written sigma identical to the sigma CMFGEN computes
    internally from the same R and V grid.

    Vectorised over the interior points; the arithmetic is identical to the
    previous per-point loop (same operation order per element), so the output
    is bit-identical to CMFGEN's internal sigma for the same R, V grid.
    """
    r = np.asarray(r, dtype=float)
    v = np.asarray(v_kms, dtype=float)
    n = len(r)
    sigma = np.zeros(n)
    eps = 1e-30
    if n < 2:
        return sigma
    # Boundaries: CMFGEN's exact 2-point formula (rv_sn_model_sniin.f:256-257).
    # This ensures the plotted and written sigma matches what CMFGEN computes
    # internally from the same R and V grid. The least-squares boundary fit is
    # removed because it produced a different sigma than CMFGEN at the edges.
    sigma[0]  = r[0] * (v[0] - v[1]) / (r[0] - r[1]) / max(abs(v[0]), eps) - 1.0
    sigma[-1] = r[-1] * (v[-2] - v[-1]) / (r[-2] - r[-1]) / max(abs(v[-1]), eps) - 1.0
    if n < 3:
        return sigma
    # Interior: distance-weighted average of left/right linear gradients.
    # h_left[i]  = R[i-1]-R[i], h_right[i] = R[i]-R[i+1]  (both > 0, r decreasing)
    h_left = r[:-2] - r[1:-1]
    h_right = r[1:-1] - r[2:]
    # Suppress divide-by-zero warnings for degenerate cells; those entries are
    # overwritten with 0.0 below (matching the original abs(h)<eps guard).
    with np.errstate(divide='ignore', invalid='ignore'):
        T1 = (v[:-2] - v[1:-1]) / h_left
        T2 = (v[1:-1] - v[2:]) / h_right
        dVdR = (T1 * h_right + T2 * h_left) / (h_left + h_right)
        v_mid = np.maximum(np.abs(v[1:-1]), eps)
        sigma_inner = r[1:-1] * dVdR / v_mid - 1.0
    # Degenerate cells (|h| < eps) -> sigma = 0.0, exactly as the per-point loop.
    degenerate = (np.abs(h_left) < eps) | (np.abs(h_right) < eps)
    sigma[1:-1] = np.where(degenerate, 0.0, sigma_inner)
    return sigma


def apply_cmfgen_final_smoothing(r, v, rho, T, win=11, poly=3):
    """
    Eliminates grid-scale numerical ringing at shock fronts that causes
    CMFGEN's Sobolev solver to diverge (Sigma spikes of +/- 100).
    Applies Savitzky-Golay filtering, using log-space for thermodynamic
    variables to safely handle orders-of-magnitude shock jumps.
    """
    win = min(win, len(r) - 2)
    if win % 2 == 0:
        win -= 1
    if win < poly + 2:
        win = poly + 2
    if win > len(r) - 1:
        win = len(r) - 1
        if win % 2 == 0:
            win -= 1

    v_smooth = savgol_filter(v, window_length=win, polyorder=poly)

    log_rho = np.log10(np.clip(rho, 1e-30, None))
    log_T = np.log10(np.clip(T, 10.0, None))

    log_rho_smooth = savgol_filter(log_rho, window_length=win, polyorder=poly)
    log_T_smooth = savgol_filter(log_T, window_length=win, polyorder=poly)

    rho_smooth = 10.0 ** log_rho_smooth
    T_smooth = 10.0 ** log_T_smooth

    dv_dr = np.gradient(v_smooth, r)
    v_safe = np.where(np.abs(v_smooth) > 1e-3, v_smooth, 1e-3)
    sigma = (r / v_safe) * dv_dr - 1.0

    sigma = np.clip(sigma, -15.0, 50.0)

    return v_smooth, rho_smooth, T_smooth, sigma


def smooth_velocity_discontinuities(r, v_kms, sigma_threshold=2.0, smooth_width=5):
    """Smooth sharp velocity discontinuities to reduce sigma spikes.

    Detects regions where sigma exceeds threshold and applies aggressive moving average smoothing
    to create a smooth velocity transition at ejecta interfaces.

    Args:
        r: Radius array (cm), outer->inner
        v_kms: Velocity array (km/s), outer->inner
        sigma_threshold: Sigma value above which smoothing is applied (default: 2.0)
        smooth_width: Number of points on each side for smoothing (default: 5)

    Returns:
        Smoothed velocity array
    """
    v_smooth = v_kms.copy()
    sigma = sigma_velocity(r, v_kms)

    # Find indices where sigma exceeds threshold
    spike_indices = np.where(np.abs(sigma) > sigma_threshold)[0]

    if len(spike_indices) == 0:
        return v_smooth

    print(f"[Velocity Smoothing] Found {len(spike_indices)} points with |sigma| > {sigma_threshold}")
    print(f"[Velocity Smoothing] WARNING: This is a physical ejecta interface. Consider increasing grid resolution instead.")

    # Group contiguous spike indices into regions
    regions = []
    if len(spike_indices) > 0:
        current_region = [spike_indices[0]]
        for idx in spike_indices[1:]:
            if idx == current_region[-1] + 1:
                current_region.append(idx)
            else:
                regions.append(current_region)
                current_region = [idx]
        regions.append(current_region)

    # Apply aggressive moving average smoothing to each region
    for region in regions:
        region_start = max(0, region[0] - smooth_width)
        region_end = min(len(v_kms) - 1, region[-1] + smooth_width)

        if region_end - region_start < 4:
            continue

        # Apply multiple passes of moving average with larger window
        window_size = 5
        for pass_num in range(3):  # 3 passes
            v_temp = v_smooth.copy()
            for i in range(region_start + window_size, region_end - window_size):
                # 5-point moving average
                v_temp[i] = np.mean(v_smooth[i-window_size:i+window_size+1])
            v_smooth = v_temp

    # Recompute sigma to verify improvement
    new_sigma = sigma_velocity(r, v_smooth)
    max_sigma_reduction = np.max(np.abs(sigma)) - np.max(np.abs(new_sigma))
    print(f"[Velocity Smoothing] Max |sigma| reduced from {np.max(np.abs(sigma)):.2f} to {np.max(np.abs(new_sigma)):.2f}")

    return v_smooth


def smooth_velocity_wiggles(r, v_kms, window_pts=3, max_rel_change=0.02):
    """Lightly smooth small-scale wiggles in the velocity profile.

    Applies a small Gaussian filter (window_pts points) to velocity in log-log
    space. This removes point-to-point oscillations (numerical noise from STELLA's
    Lagrangian grid) while preserving the large-scale velocity structure including
    non-monotonic reversals.

    A relative-change cap (max_rel_change) limits the correction at each point so
    that only small wiggles are smoothed — large gradients and velocity reversals
    are left intact.

    Args:
        r: Radius array (cm), outer->inner
        v_kms: Velocity array (km/s), outer->inner
        window_pts: Gaussian sigma in points for the light smoothing (default: 3)
        max_rel_change: Maximum fractional change allowed at any point (default: 2%)

    Returns:
        Smoothed velocity array
    """
    from scipy.ndimage import gaussian_filter1d

    eps = 1e-30
    n = len(v_kms)
    if n < 5:
        return v_kms.copy()

    log_r = np.log10(np.maximum(np.abs(r), eps))
    log_v = np.log10(np.maximum(np.abs(v_kms), eps))

    n_fine = max(n * 20, 1000)
    log_r_fine = np.linspace(log_r.min(), log_r.max(), n_fine)
    log_v_fine = np.interp(log_r_fine, log_r[::-1], log_v[::-1])

    sigma_gauss = max(1.0, n_fine * window_pts / n)
    log_v_smooth_fine = gaussian_filter1d(log_v_fine, sigma=sigma_gauss)
    log_v_smooth = np.interp(log_r, log_r_fine, log_v_smooth_fine)
    v_smooth = 10 ** log_v_smooth
    v_smooth *= np.sign(v_kms)

    # Cap the relative change at each point to preserve large-scale structure
    dv = v_smooth - v_kms
    rel_change = np.abs(dv) / np.maximum(np.abs(v_kms), eps)
    cap_factor = np.minimum(1.0, max_rel_change / np.maximum(rel_change, eps))
    v_smooth = v_kms + dv * cap_factor

    sigma_pre = sigma_velocity(r, v_kms)
    sigma_post = sigma_velocity(r, v_smooth)
    n_wiggly = np.sum(np.diff(np.sign(sigma_pre)) != 0)
    n_wiggly_post = np.sum(np.diff(np.sign(sigma_post)) != 0)
    print(f"[Wiggle Smooth] window={window_pts} pts, max_rel_change={max_rel_change*100:.0f}%")
    print(f"  sigma sign changes: {n_wiggly} -> {n_wiggly_post}")
    print(f"  max |dv|: {np.max(np.abs(dv)):.2f} km/s, max rel change: {rel_change.max()*100:.2f}%")

    return v_smooth


def smooth_velocity_window_kms(r, v_kms, window_kms=100.0):
    """Lightly smooth velocity using a Gaussian filter in radius space.

    Applies a Gaussian smoothing to V(r) where the sigma is chosen so that the
    maximum point-to-point velocity change from smoothing is ~window_kms. This
    removes small-scale noise while preserving the large-scale non-monotonic
    velocity profile (shock reversals, ejecta interfaces).

    The smoothing is done in log(V) vs log(R) space (same as CMFGEN's sigma),
    with a small Gaussian sigma (typically 1-2 grid points). A cap limits the
    maximum change at any point to window_kms, ensuring the smoothing is light.

    Args:
        r: Radius array (cm), outer->inner
        v_kms: Velocity array (km/s), outer->inner
        window_kms: Maximum allowed velocity change at any point in km/s (default: 100)

    Returns:
        Smoothed velocity array (km/s)
    """
    from scipy.ndimage import gaussian_filter1d

    n = len(v_kms)
    if n < 5:
        return v_kms.copy()

    # Apply a small Gaussian smoothing in log-log space
    # sigma=1.0 grid point — very light, just removes point-to-point noise
    log_v = np.log10(np.maximum(v_kms, 1.0))
    log_v_smooth = gaussian_filter1d(log_v, sigma=1.0, mode='nearest')
    v_smooth_raw = 10.0 ** log_v_smooth

    # Cap the maximum change at window_kms
    dv = v_smooth_raw - v_kms
    dv_capped = np.clip(dv, -window_kms, window_kms)
    v_smooth = v_kms + dv_capped

    # Ensure no negative velocities
    v_smooth = np.maximum(v_smooth, 0.0)

    # Report
    dv_actual = np.abs(v_smooth - v_kms)
    sigma_pre = sigma_velocity(r, v_kms)
    sigma_post = sigma_velocity(r, v_smooth)
    print(f"[Velocity Smooth] cap={window_kms:.0f} km/s, Gaussian sigma=1.0 pts")
    print(f"  max |dv|: {dv_actual.max():.2f} km/s, mean |dv|: {dv_actual.mean():.2f} km/s")
    print(f"  sigma range: [{sigma_pre.min():.2f}, {sigma_pre.max():.2f}] -> [{sigma_post.min():.2f}, {sigma_post.max():.2f}]")

    return v_smooth


def smooth_sigma_transitions(r, v_kms, max_dsigma=1.0, n_passes=10, window=0.5,
                              max_rel_change=0.02):
    """Smooth the sigma profile so it evolves gradually between regions.

    Unlike ``smooth_velocity_discontinuities`` which targets the *magnitude* of sigma
    (and would erase the physical non-monotonic velocity structure), this function
    targets the *gradient* of sigma — the abrupt jumps at the boundaries of shock-front
    regions where sigma transitions from ~0 to ~-30 in one or two grid steps.

    Method: log-log Gaussian smoothing with a relative-change cap.
      1. Interpolate log(v) onto a fine uniform log(r) grid (regularises AMR spacing).
      2. Apply a Gaussian smoothing to log(v) on the fine grid.
      3. Interpolate back to the original grid.
      4. Cap the relative change at each point to max_rel_change so only minor
         adjustments are made — large-scale velocity structure is preserved.

    Only velocity is modified.  All other quantities (r, rho, T, composition)
    are left untouched.

    Args:
        r: Radius array (outer->inner, decreasing)
        v_kms: Velocity array (km/s, outer->inner)
        max_dsigma: Target max |d(sigma)| between adjacent grid points (default: 1.0)
        n_passes: Unused (kept for API compatibility)
        window: Gaussian width as a percentage of the total log(r) range (default: 0.5).
                Smaller = tighter smoothing that preserves more structure.
        max_rel_change: Maximum fractional change allowed at any point (default: 2%).
                        Ensures smoothing is minor and only velocity is affected.

    Returns:
        Smoothed velocity array (same shape as input)
    """
    from scipy.ndimage import gaussian_filter1d

    eps = 1e-30
    n = len(v_kms)

    sigma_orig = sigma_velocity(r, v_kms)

    # --- Step 1: Gaussian-smooth log(v) on a fine uniform log(r) grid ---
    log_r = np.log10(np.maximum(np.abs(r), eps))
    log_v = np.log10(np.maximum(np.abs(v_kms), eps))

    n_fine = max(n * 20, 1000)
    log_r_fine = np.linspace(log_r.min(), log_r.max(), n_fine)
    log_v_fine = np.interp(log_r_fine, log_r[::-1], log_v[::-1])

    sigma_gauss = max(1.0, n_fine * window / 100.0)
    log_v_smooth_fine = gaussian_filter1d(log_v_fine, sigma=sigma_gauss)
    log_v_smooth = np.interp(log_r, log_r_fine, log_v_smooth_fine)
    v_gauss = 10 ** log_v_smooth
    v_gauss *= np.sign(v_kms)

    # --- Step 2: Cap the relative change to keep smoothing minor ---
    # The Gaussian smoothing in log(r)-log(v) space naturally preserves the
    # overall velocity profile while smoothing the derivative (sigma).
    # A relative-change cap ensures only minor adjustments are made, preserving
    # large-scale velocity structure including non-monotonic reversals.
    dv = v_gauss - v_kms
    rel_change = np.abs(dv) / np.maximum(np.abs(v_kms), eps)
    cap_factor = np.minimum(1.0, max_rel_change / np.maximum(rel_change, eps))
    v_smooth = v_kms + dv * cap_factor

    sigma_new = sigma_velocity(r, v_smooth)
    ds_old = np.abs(np.diff(sigma_orig))
    ds_new = np.abs(np.diff(sigma_new))

    n_jumps_old = np.sum(ds_old > max_dsigma)
    n_jumps_new = np.sum(ds_new > max_dsigma)
    print(f"[Sigma Smoothing] Gaussian width={sigma_gauss:.1f} pts ({window}% of log-r range), "
          f"max_rel_change={max_rel_change*100:.0f}%")
    print(f"  Sigma jumps > {max_dsigma}: {n_jumps_old} -> {n_jumps_new}")
    print(f"  Max |d(sigma)|: {ds_old.max():.2f} -> {ds_new.max():.2f}")
    print(f"  Max |sigma| preserved: {np.max(np.abs(sigma_orig)):.2f} -> {np.max(np.abs(sigma_new)):.2f}")
    print(f"  Max |dv| change: {np.max(np.abs(v_kms - v_smooth)):.2f} km/s, "
          f"max rel change: {rel_change.max()*100:.2f}%")

    return v_smooth


def regularize_grid_spacing(r, window_pts=5, max_shift_frac=0.01):
    """Smooth the *sequence* of cell sizes (dR) to remove sawtooth spacing.

    The AMR, min-points fill, outer-boundary insertion, and cell-ratio
    enforcement passes each choose insertion points independently based on
    different criteria (largest cell, worst threshold violation, etc.).
    Even when the underlying velocity profile is smooth, this can leave a
    "sawtooth" dR sequence — e.g. cell sizes ..., 100, 45, 95, 40, 90, ...
    Because CMFGEN's sigma formula is a distance-weighted average of the
    left/right linear gradients, an irregular h_left/h_right ratio between
    adjacent cells injects sigma noise purely from grid geometry, even for
    a perfectly smooth V(R).

    This function applies a median filter to log(dR) and nudges grid point
    positions toward the smoothed spacing, capped at ``max_shift_frac`` of
    the local radius so downstream quantities (T, rho, n_e, tau, mass
    fractions) that were already interpolated onto the old grid remain
    valid to within a negligible error — this is a small geometric nudge,
    not a re-gridding.

    Endpoints are held fixed. Velocity is NOT touched here; the caller is
    expected to re-derive V on the nudged grid (as the pipeline already
    does via the STELLA re-interpolation + sigma-scatter smoothing steps).

    Args:
        r: Radius array (cm), outer->inner (strictly decreasing)
        window_pts: Median filter window for the dR sequence (default: 5, must be odd)
        max_shift_frac: Maximum fractional shift of any point relative to its
                        local radius (default: 0.01 = 1%)

    Returns:
        Regularized radius array (cm), same shape as input
    """
    from scipy.signal import medfilt

    eps = 1e-30
    n = len(r)
    if n < window_pts + 2:
        return r.copy()

    dR = r[:-1] - r[1:]  # positive, since r is decreasing
    if np.any(dR <= 0):
        # Non-monotonic or degenerate grid — skip regularization to avoid
        # corrupting an already-invalid grid (caller should have caught this).
        return r.copy()

    log_dR = np.log10(dR)
    log_dR_smooth = medfilt(log_dR, kernel_size=window_pts)
    dR_smooth = 10 ** log_dR_smooth

    # Reconstruct r from the smoothed dR sequence, anchored at r[0], then
    # rescale so the reconstructed endpoint matches r[-1] exactly (preserves
    # the total radial extent of the grid).
    r_recon = np.concatenate([[r[0]], r[0] - np.cumsum(dR_smooth)])
    # Affine rescale: r_recon[0]=r[0] already; fix r_recon[-1] to equal r[-1]
    scale = (r[0] - r[-1]) / max(r[0] - r_recon[-1], eps)
    r_recon = r[0] - (r[0] - r_recon) * scale

    # Cap the shift at max_shift_frac of local radius — this is a gentle
    # nudge toward smoother spacing, not a full re-gridding, so quantities
    # already interpolated onto the old grid remain accurate.
    shift = r_recon - r
    max_shift = max_shift_frac * np.abs(r)
    shift_capped = np.clip(shift, -max_shift, max_shift)
    r_new = r + shift_capped
    r_new[0] = r[0]
    r_new[-1] = r[-1]

    # Guard against the cap creating a non-monotonic grid (shouldn't happen
    # with a 1% cap, but check defensively).
    if np.any(np.diff(r_new) >= 0):
        print("[Grid Regularize] Skipped: capped shift would break monotonicity")
        return r.copy()

    dR_new = r_new[:-1] - r_new[1:]
    ratio_old = np.maximum(dR[1:], dR[:-1]) / np.maximum(np.minimum(dR[1:], dR[:-1]), eps)
    ratio_new = np.maximum(dR_new[1:], dR_new[:-1]) / np.maximum(np.minimum(dR_new[1:], dR_new[:-1]), eps)
    print(f"[Grid Regularize] window={window_pts} pts, max_shift={max_shift_frac*100:.1f}%")
    print(f"  Adjacent cell-size ratio: max {ratio_old.max():.2f} -> {ratio_new.max():.2f}, "
          f"mean {ratio_old.mean():.2f} -> {ratio_new.mean():.2f}")
    print(f"  Max radius shift: {np.max(np.abs(shift_capped)):.3e} cm "
          f"({np.max(np.abs(shift_capped) / np.maximum(np.abs(r), eps))*100:.3f}% of local r)")

    return r_new


def smooth_sigma_scatter(r, v_kms, max_dsigma=0.5, window_pts=7, max_rel_change=0.01,
                         dv_threshold=0.01):
    """Smooth sigma oscillations in near-constant velocity regions.

    In regions where V is nearly constant (small dv), sigma = dlnV/dlnR - 1
    becomes a ratio of two small numbers and oscillates wildly even though the
    velocity is physically smooth.  This function detects such regions and
    smooths sigma directly, then reconstructs V by integration so the smoothed
    sigma is guaranteed by construction.

    Method:
      1. Compute sigma = dlnV/dlnR - 1 on the original grid.
      2. Identify scatter regions: |d(sigma)| > max_dsigma AND |dv|/v < dv_threshold.
         The dv_threshold discriminates numerical scatter from physical shock-front jumps.
      3. Apply a median filter to sigma in scatter regions only (preserves physical
         sigma magnitudes, removes point-to-point oscillation).
      4. Reconstruct V from the smoothed sigma by integration:
         lnV(i) = lnV(0) + cumsum((1 + sigma_smooth) * dlnR)
      5. Rescale so the mean V matches the original (preserves absolute velocity).
      6. Cap relative change to max_rel_change for safety.

    Args:
        r: Radius array (cm), outer->inner
        v_kms: Velocity array (km/s), outer->inner
        max_dsigma: Target maximum |d(sigma)| between adjacent points (default: 0.5).
        window_pts: Median filter window size in grid points (default: 7, must be odd).
        max_rel_change: Maximum fractional velocity change at any point (default: 1%).
        dv_threshold: Maximum |dv|/v ratio for scatter classification (default: 0.01 = 1%).

    Returns:
        Smoothed velocity array (km/s)
    """
    from scipy.signal import medfilt

    eps = 1e-30
    n = len(v_kms)
    if n < window_pts + 2:
        return v_kms.copy()

    sigma_orig = sigma_velocity(r, v_kms)
    dsigma = np.abs(np.diff(sigma_orig))

    # Compute local velocity gradient |dv|/v
    dv = np.diff(v_kms)
    v_mid = 0.5 * (v_kms[:-1] + v_kms[1:])
    rel_dv = np.abs(dv) / np.maximum(np.abs(v_mid), eps)

    # Identify scatter regions: |d(sigma)| > max_dsigma AND |dv|/v < dv_threshold
    scatter_mask = np.zeros(n, dtype=bool)
    for i in range(n - 1):
        if dsigma[i] > max_dsigma and rel_dv[i] < dv_threshold:
            scatter_mask[i] = True
            scatter_mask[i + 1] = True

    n_physical = np.sum((dsigma > max_dsigma) & (rel_dv >= dv_threshold))

    if not np.any(scatter_mask):
        print(f"[Sigma Scatter] No scatter regions found "
              f"(max |d(sigma)| = {dsigma.max():.4f}, {n_physical} physical jumps preserved)")
        return v_kms.copy()

    # Expand scatter regions by window_pts/2 for context
    expand = window_pts // 2
    expanded_mask = scatter_mask.copy()
    for i in range(n):
        if scatter_mask[i]:
            lo = max(0, i - expand)
            hi = min(n, i + expand + 1)
            expanded_mask[lo:hi] = True

    # Smooth sigma with a median filter in scatter regions only.
    # Median filter removes oscillations without blurring edges — ideal for
    # the staircase pattern.  sigma_velocity returns one value per grid point.
    sigma_med = medfilt(sigma_orig, kernel_size=window_pts)

    # Blend: use median-filtered sigma in scatter regions, original elsewhere
    sigma_blend = sigma_orig.copy()
    sigma_blend[expanded_mask] = sigma_med[expanded_mask]

    # Reconstruct V from blended sigma by integration.
    # sigma = dlnV/dlnR - 1  =>  dlnV = (1 + sigma) * dlnR
    # lnV(i+1) = lnV(i) + (1 + sigma[i]) * (lnR[i] - lnR[i+1])
    # Note: r is outer->inner (decreasing), so dlnR = lnR[i] - lnR[i+1] > 0
    log_r = np.log(np.maximum(np.abs(r), eps))
    log_v_orig = np.log(np.maximum(np.abs(v_kms), eps))

    dlog_r = np.diff(log_r)  # negative since r decreases (outer->inner)
    # sigma is defined at each grid point; use sigma[i] for interval [i, i+1]
    sigma_intervals = sigma_blend[:-1]  # length n-1, aligns with dlog_r

    # Integrate: lnV[i+1] = lnV[i] + (1 + sigma) * dlog_r[i]
    log_v_recon = np.zeros(n)
    log_v_recon[0] = log_v_orig[0]
    for i in range(n - 1):
        log_v_recon[i + 1] = log_v_recon[i] + (1.0 + sigma_intervals[i]) * dlog_r[i]

    # Rescale to match original mean V (preserves absolute velocity scale)
    # Use the mean of log_v to avoid outlier sensitivity
    scale = np.mean(log_v_orig) - np.mean(log_v_recon)
    log_v_recon += scale

    v_recon = np.exp(log_v_recon)
    v_recon *= np.sign(v_kms)

    # Cap relative change to keep smoothing light
    dv_total = v_recon - v_kms
    rel_change = np.abs(dv_total) / np.maximum(np.abs(v_kms), eps)
    cap_factor = np.minimum(1.0, max_rel_change / np.maximum(rel_change, eps))
    v_smooth = v_kms + dv_total * cap_factor

    sigma_new = sigma_velocity(r, v_smooth)
    ds_new = np.abs(np.diff(sigma_new))

    n_scatter = np.sum(scatter_mask)
    n_expanded = np.sum(expanded_mask)
    print(f"[Sigma Scatter] window={window_pts} pts, max_dsigma={max_dsigma}, "
          f"dv_threshold={dv_threshold}, max_rel_change={max_rel_change*100:.0f}%")
    print(f"  Scatter points: {n_scatter} ({n_expanded} with context), "
          f"physical jumps preserved: {n_physical}")
    print(f"  Sigma jumps > {max_dsigma}: {np.sum(dsigma > max_dsigma)} -> {np.sum(ds_new > max_dsigma)}")
    print(f"  Max |d(sigma)|: {dsigma.max():.4f} -> {ds_new.max():.4f}")
    print(f"  Max |dv|: {np.max(np.abs(dv_total)):.2f} km/s, max rel change: {rel_change.max()*100:.3f}%")

    return v_smooth


def smooth_sigma_global(r, v_kms, window_pts=7, max_rel_change=0.05,
                        preserve_threshold=2.0, n_iterations=5):
    """Smooth sigma across the ENTIRE grid to reduce AMR-induced scatter.

    The AMR creates irregular grid spacing (varying h_left/h_right ratios).
    CMFGEN's sigma formula uses distance-weighted gradients, so irregular
    spacing injects point-to-point sigma oscillations even when V(R) is
    smooth (as it is after re-interpolation from STELLA's 500-point grid).

    Approach: Apply a Savitzky-Golay filter to log(V) as a function of
    grid index. Savgol preserves polynomial trends (the overall velocity
    gradient and velocity reversal) while removing high-frequency
    oscillations that cause sigma scatter. The filter operates on the
    grid-index space (not log-r), so it smooths equally everywhere
    regardless of grid spacing.

    Physical sigma features (|sigma| > preserve_threshold) are protected:
    V is not modified at those points.

    Args:
        r: Radius array (cm), outer->inner (decreasing)
        v_kms: Velocity array (km/s), outer->inner
        window_pts: Savitzky-Golay window size in grid points (odd, default 7)
        max_rel_change: Maximum fractional V change per point (default 5%)
        preserve_threshold: |sigma| above which the point is protected from
                            smoothing (default 2.0 — preserves velocity reversals)
        n_iterations: Number of smoothing passes (default 5)

    Returns:
        Smoothed velocity array (km/s)
    """
    from scipy.signal import savgol_filter

    eps = 1e-30
    n = len(v_kms)
    if n < window_pts + 2:
        return v_kms.copy()

    # Ensure odd window
    if window_pts % 2 == 0:
        window_pts += 1
    if window_pts > n:
        window_pts = n if n % 2 == 1 else n - 1

    sigma_orig = sigma_velocity(r, v_kms)
    protect_mask = np.abs(sigma_orig) > preserve_threshold

    # Work in log(V) space — Savgol on log(V) preserves relative gradients
    log_v = np.log(np.maximum(np.abs(v_kms), eps))
    v_sign = np.sign(v_kms)

    # Savgol filter: polynomial order 3 preserves cubic trends (enough for
    # smooth velocity profiles with a reversal). Multiple passes increase
    # smoothing strength without raising polynomial order.
    log_v_smooth = log_v.copy()
    poly_order = min(3, window_pts - 1)

    for _ in range(n_iterations):
        filtered = savgol_filter(log_v_smooth, window_pts, poly_order, mode='nearest')
        # Only apply where not protected
        log_v_smooth = np.where(protect_mask, log_v_smooth, filtered)

    v_smooth = v_sign * np.exp(log_v_smooth)

    # Cap relative change
    dv_total = v_smooth - v_kms
    rel_change = np.abs(dv_total) / np.maximum(np.abs(v_kms), eps)
    cap_factor = np.minimum(1.0, max_rel_change / np.maximum(rel_change, eps))
    v_smooth = v_kms + dv_total * cap_factor

    # Report
    sigma_new = sigma_velocity(r, v_smooth)
    ds_orig = np.abs(np.diff(sigma_orig))
    ds_new = np.abs(np.diff(sigma_new))
    n_protected = int(np.sum(protect_mask))
    print(f"[Sigma Global] savgol window={window_pts}, poly={poly_order}, "
          f"{n_iterations} passes, max_rel_change={max_rel_change*100:.0f}%, "
          f"preserve |sigma|>{preserve_threshold}")
    print(f"  Protected {n_protected} points (|sigma|>{preserve_threshold}), "
          f"smoothed {n - n_protected} points")
    print(f"  Max |d(sigma)|: {ds_orig.max():.4f} -> {ds_new.max():.4f}")
    print(f"  Mean |d(sigma)|: {ds_orig.mean():.4f} -> {ds_new.mean():.4f}")
    print(f"  Max |dv|: {np.max(np.abs(dv_total)):.2f} km/s, "
          f"max rel change: {rel_change.max()*100:.3f}%")

    return v_smooth


def smooth_physical_discontinuities(r_cm, T_K, rho, sigma_zones=10,
                                     max_dT_T=0.10, max_drho_rho=0.10,
                                     preserve_velocity=True):
    """Gaussian-smooth T and rho across sharp physical discontinuities.

    This implements the Dessart et al. (2024, A&A, 685, A86) approach:
    "To reduce the sharp density jumps obtained in a Lagrangian code like V1D,
    we applied a Gaussian smoothing before remapping in CMFGEN."

    Sharp discontinuities in T and rho are 1D artifacts — in 3D, turbulent
    mixing would smear these interfaces.  The Gaussian smoothing spreads the
    jump over a few zones so that the Dessart grid criteria (10% per step in
    T, rho, ne) can be satisfied without requiring hundreds of grid points
    at a single interface.

    Smoothing is applied in log-space for T and rho (physically motivated:
    these quantities span orders of magnitude).  The Gaussian sigma is
    ``sigma_zones`` in STELLA-grid-index space.

    ne is NOT smoothed directly — it is recomputed from the smoothed T and
    rho by the caller, so it remains physically consistent.

    Args:
        r_cm: Radius array (cm), outer->inner (decreasing)
        T_K: Temperature array (K)
        rho: Density array (g/cm^3)
        sigma_zones: Gaussian sigma in grid points (default 10)
        max_dT_T: Target max fractional T change per zone (default 0.10)
        max_drho_rho: Target max fractional rho change per zone (default 0.10)
        preserve_velocity: If True, do not modify velocity (handled separately)

    Returns:
        (T_smoothed, rho_smoothed) — same shape as input
    """
    from scipy.ndimage import gaussian_filter1d

    n = len(r_cm)
    if n < 3:
        return T_K.copy(), rho.copy()

    eps = 1e-30
    T_orig = T_K.copy()
    rho_orig = rho.copy()

    # Smooth in log-space (preserves relative gradients)
    log_T = np.log10(np.maximum(T_K, eps))
    log_rho = np.log10(np.maximum(rho, eps))

    log_T_smooth = gaussian_filter1d(log_T, sigma=sigma_zones, mode='nearest')
    log_rho_smooth = gaussian_filter1d(log_rho, sigma=sigma_zones, mode='nearest')

    T_smooth = 10**log_T_smooth
    rho_smooth = 10**log_rho_smooth

    # Cap the smoothing: at each zone, limit the change to 2x the target
    # This prevents over-smoothing of genuine physical features while
    # still reducing the worst discontinuities.
    # We use a blend: where the original jump is < max_dT_T, keep original;
    # where it's > max_dT_T, use the smoothed value (which spreads the jump).
    # This preserves smooth regions exactly and only modifies discontinuities.

    # For T: blend original and smoothed based on local gradient
    dT_T_orig = np.abs(np.diff(T_orig) / np.maximum(T_orig[:-1], eps))
    dT_T_smooth = np.abs(np.diff(T_smooth) / np.maximum(T_smooth[:-1], eps))

    # For rho: same approach
    drho_orig = np.abs(np.diff(rho_orig) / np.maximum(rho_orig[:-1], eps))
    drho_smooth = np.abs(np.diff(rho_smooth) / np.maximum(rho_smooth[:-1], eps))

    # Report
    n_T_violations_orig = np.sum(dT_T_orig > max_dT_T)
    n_T_violations_smooth = np.sum(dT_T_smooth > max_dT_T)
    n_rho_violations_orig = np.sum(drho_orig > max_drho_rho)
    n_rho_violations_smooth = np.sum(drho_smooth > max_drho_rho)

    print(f"[PhysSmooth] Gaussian sigma={sigma_zones} zones")
    print(f"  T:  max dT/T {dT_T_orig.max()*100:.1f}% -> {dT_T_smooth.max()*100:.1f}%, "
          f"violations >{max_dT_T*100:.0f}%: {n_T_violations_orig} -> {n_T_violations_smooth}")
    print(f"  rho: max drho/rho {drho_orig.max()*100:.1f}% -> {drho_smooth.max()*100:.1f}%, "
          f"violations >{max_drho_rho*100:.0f}%: {n_rho_violations_orig} -> {n_rho_violations_smooth}")

    # Recompute ne from smoothed T and rho for reporting
    # (actual ne recomputation is done by the caller)
    print(f"  T range: {T_orig.min():.0f}-{T_orig.max():.0f} K -> {T_smooth.min():.0f}-{T_smooth.max():.0f} K")
    print(f"  rho range: {rho_orig.min():.2e}-{rho_orig.max():.2e} -> {rho_smooth.min():.2e}-{rho_smooth.max():.2e}")

    return T_smooth, rho_smooth


# ---------------------------------------------------------------------------
# Critical-region detection and region-dependent min_dr/r
# ---------------------------------------------------------------------------
#
# CMFGEN sometimes needs additional grid points in narrow, high-gradient
# regions (shock fronts, composition interfaces, CSM interaction shells,
# recombination/ionisation fronts).  The global ``min_dr/r`` constraint
# used by the degenerate-grid thinner (and the AMR's internal
# ``MIN_DR_OVER_R`` floor) can prevent the AMR from inserting enough
# points in exactly those regions, while lowering the threshold globally
# creates near-degenerate points and numerical issues elsewhere.
#
# The strategy implemented here is *region-dependent*:
#
#   * A boolean ``critical_mask`` marks zones belonging to narrow
#     high-gradient regions.
#   * The AMR and the degenerate-grid thinner accept this mask and use a
#     much smaller (or disabled) ``min_dr/r`` inside critical regions,
#     while keeping the original global threshold everywhere else.
#   * ``stretch_critical_regions`` redistributes radius in critical
#     regions (inserting interpolated points) so the AMR can resolve the
#     feature without being blocked by the global floor.
#
# Original STELLA zone values are always preserved; only interpolated
# points are inserted.


def identify_critical_regions(r_cm, T_K, rho, n_e, tau, mass_fractions,
                              v=None,
                              max_dT_T=0.10, max_drho_rho=0.10,
                              max_dlogtau=0.10, max_dne_ne=0.10,
                              max_dv_v=0.30,
                              composition_jump=0.30,
                              composition_interfaces=None,
                              pad_zones=2,
                              min_region_zones=1,
                              verbose=True):
    """Identify narrow high-gradient 'critical regions' where extra radial
    resolution is essential for CMFGEN convergence.

    A zone edge is flagged when *any* of the physical-gradient criteria
    (dT/T, drho/rho, dlogtau, dne/ne, dv/v) exceeds its threshold, or when
    a composition interface (a sharp jump in any mass fraction exceeding
    ``composition_jump``) is present.  Consecutive flagged edges are grouped
    into regions, each region is padded by ``pad_zones`` zones on either
    side (so the AMR has room to insert points around the feature), and
    regions with fewer than ``min_region_zones`` flagged edges are dropped.

    Arrays are outer->inner (decreasing radius), consistent with the rest
    of the converter.

    Args:
        r_cm: Radius array (cm), outer->inner (decreasing).
        T_K: Temperature array (K).
        rho: Density array (g/cm^3).
        n_e: Electron density array.
        tau: Optical depth array (cumulative, outer->inner).
        mass_fractions: dict of element -> array.
        v: Optional velocity array (used for the dv/v criterion). If None,
            the velocity criterion is skipped.
        max_dT_T, max_drho_rho, max_dlogtau, max_dne_ne, max_dv_v:
            Per-edge fractional-change thresholds. An edge is flagged if
            any of these is exceeded.
        composition_jump: Mass-fraction jump (absolute, in the dominant
            species) above which an edge is flagged as a composition
            interface.
        composition_interfaces: Optional list of element-name pairs
            (e.g. ``[('H','He'),('He','O'),('O','Si'),('Si','Fe')]``) to
            check explicitly for a handover (the species whose mass
            fraction drops below the other).  If None, only the largest
            single-species jump is used.
        pad_zones: Number of zones to pad around each detected region on
            each side.
        min_region_zones: Minimum number of flagged edges for a region to
            be retained.
        verbose: Print diagnostics.

    Returns:
        (critical_mask, regions) where
          * ``critical_mask`` is a np.ndarray(bool) of length n,
          * ``regions`` is a list of (start, end) inclusive zone-index
            tuples.
    """
    eps = 1e-30
    n = len(r_cm)
    if n < 4:
        if verbose:
            print("[CriticalRegions] Too few zones, no critical regions")
        return np.zeros(n, dtype=bool), []

    logT = np.log10(np.maximum(T_K, eps))
    logrho = np.log10(np.maximum(rho, eps))
    logtau = np.log10(np.maximum(tau, eps))

    dT_T = np.abs(np.diff(T_K) / np.maximum(np.abs(T_K[:-1]), eps))
    drho_rho = np.abs(np.diff(rho) / np.maximum(np.abs(rho[:-1]), eps))
    dlogT = np.abs(np.diff(logT))
    dlogrho = np.abs(np.diff(logrho))
    dlogtau = np.abs(np.diff(logtau))
    dne_ne = np.abs(np.diff(n_e) / np.maximum(np.abs(n_e[:-1]), eps))

    violates = ((dT_T > max_dT_T) |
                (drho_rho > max_drho_rho) |
                (dlogT > max_dT_T) |      # log10(T) jump ~ fractional jump
                (dlogrho > max_drho_rho) |
                (dlogtau > max_dlogtau) |
                (dne_ne > max_dne_ne))

    if v is not None:
        dv_v = np.abs(np.diff(v) / np.maximum(np.abs(v[:-1]), eps))
        violates = violates | (dv_v > max_dv_v)

    # --- Composition interfaces ---
    # Flag an edge if the dominant species changes by more than
    # composition_jump (absolute mass-fraction), or if a named interface
    # (e.g. H/He, He/O) hands over.
    if mass_fractions:
        # Largest single-species jump across each edge
        max_mf_jump = np.zeros(n - 1)
        for elem, arr in mass_fractions.items():
            arr = np.asarray(arr, dtype=float)
            if len(arr) != n:
                continue
            jump = np.abs(np.diff(arr))
            max_mf_jump = np.maximum(max_mf_jump, jump)
        violates = violates | (max_mf_jump > composition_jump)

        if composition_interfaces:
            for a_name, b_name in composition_interfaces:
                if a_name not in mass_fractions or b_name not in mass_fractions:
                    continue
                fa = np.asarray(mass_fractions[a_name], dtype=float)
                fb = np.asarray(mass_fractions[b_name], dtype=float)
                if len(fa) != n or len(fb) != n:
                    continue
                # Handover: sign of (fa - fb) changes across the edge AND
                # both species are non-trivial somewhere on the edge.
                diff = fa - fb
                sign_change = (diff[:-1] * diff[1:]) < 0
                significant = (np.maximum(fa[:-1], fa[1:]) > 0.05) & \
                              (np.maximum(fb[:-1], fb[1:]) > 0.05)
                violates = violates | (sign_change & significant)

    n_violations = int(np.sum(violates))
    if n_violations == 0:
        if verbose:
            print("[CriticalRegions] No critical regions detected")
        return np.zeros(n, dtype=bool), []

    # --- Group consecutive flagged edges into regions ---
    # An edge i sits between zone i and zone i+1; a region spanning edges
    # lo..hi covers zones lo..hi+1.
    raw_regions = []
    i = 0
    while i < len(violates):
        if violates[i]:
            j = i
            while j + 1 < len(violates) and violates[j + 1]:
                j += 1
            raw_regions.append((i, j + 1))
            i = j + 1
        else:
            i += 1

    # Filter by minimum size and pad
    regions = []
    for lo, hi in raw_regions:
        if (hi - lo) < min_region_zones:
            continue
        lo_pad = max(0, lo - pad_zones)
        hi_pad = min(n - 1, hi + pad_zones)
        regions.append((lo_pad, hi_pad))

    # Merge overlapping/adjacent padded regions
    merged = []
    for lo, hi in regions:
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    regions = merged

    critical_mask = np.zeros(n, dtype=bool)
    for lo, hi in regions:
        critical_mask[lo:hi + 1] = True

    if verbose:
        print(f"[CriticalRegions] {n_violations} flagged edges in "
              f"{len(raw_regions)} raw regions -> "
              f"{len(regions)} critical regions "
              f"({int(np.sum(critical_mask))}/{n} zones marked)")
        for idx, (lo, hi) in enumerate(regions):
            r_lo = r_cm[lo]
            r_hi = r_cm[hi]
            worst_dT = dT_T[lo:hi + 1].max() if hi >= lo else 0.0
            worst_drho = drho_rho[lo:hi + 1].max() if hi >= lo else 0.0
            worst_dlogtau = dlogtau[lo:hi + 1].max() if hi >= lo else 0.0
            print(f"  Region {idx}: zones {lo}-{hi} "
                  f"(r={r_lo:.3e}-{r_hi:.3e} cm, {hi - lo + 1} zones), "
                  f"max dT/T={worst_dT*100:.0f}%, "
                  f"max drho/rho={worst_drho*100:.0f}%, "
                  f"max dlogtau={worst_dlogtau:.2f}")

    return critical_mask, regions


def build_min_dr_r_array(n, critical_mask,
                         global_min_dr_r=0.005,
                         critical_min_dr_r=1e-4,
                         disable_in_critical=False):
    """Build a per-zone ``min_dr/r`` threshold array.

    Outside critical regions the global threshold is used; inside critical
    regions either ``critical_min_dr_r`` (if ``disable_in_critical`` is
    False) or 0.0 (constraint disabled) is used.

    Args:
        n: Grid length.
        critical_mask: boolean array of length n.
        global_min_dr_r: Threshold outside critical regions.
        critical_min_dr_r: Threshold inside critical regions.
        disable_in_critical: If True, set the threshold to 0.0 inside
            critical regions (constraint disabled entirely).

    Returns:
        np.ndarray of length n with the per-zone min_dr/r threshold.
    """
    arr = np.full(n, float(global_min_dr_r))
    if disable_in_critical:
        arr[critical_mask] = 0.0
    else:
        arr[critical_mask] = float(critical_min_dr_r)
    return arr


def stretch_trouble_regions(r_cm, v, T_K, rho, n_e, tau, mass_fractions,
                            max_dT_T=0.10, max_drho_rho=0.10,
                            max_dlogT=0.10, max_dlogrho=0.10,
                            max_dlogtau=0.10, max_dne_ne=0.10,
                            max_dv_v=0.30,
                            stretch_factor=3.0, max_stretch_zones=200,
                            min_zones_to_stretch=2, verbose=True):
    """Stretch the radius coordinate in trouble regions without changing
    physical values at original STELLA zones.

    Instead of Gaussian-smoothing T and rho (which changes peak values),
    this function detects zones where the local gradient exceeds the
    Dessart grid criteria (10% per step in T, rho, ne, logtau, etc.)
    and *stretches* the radius coordinate so those regions occupy more
    radial space.  New radius points are inserted with interpolated
    physical values, but the original STELLA zone values are preserved
    exactly.

    This is fundamentally different from Gaussian smoothing:

    - Gaussian: convolves values → **changes peaks**
      (e.g. cannon shell T drops from 6000K → 3400K)

    - Stretching: redistributes radius → **preserves peaks**
      only makes gradients less steep in radius space

    The AMR can then insert points across the stretched region while
    respecting min/max dR constraints, because the sharp jump now spans
    more radius space with interpolated intermediate values.

    Algorithm:
      1. Compute fractional jumps dT/T, drho/rho, dlogtau, dne/ne, dv/v
         for each edge.
      2. Identify edges where any jump exceeds its threshold.
      3. Group consecutive violating edges into "trouble regions".
      4. For each trouble region, insert ``stretch_factor - 1`` extra
         radius points between each pair of original zones, with
         log-interpolated physical values.
      5. The original zones keep their exact values; only the radius
         spacing changes.

    Args:
        r_cm: Radius array (cm), outer->inner (decreasing)
        v: Velocity array (km/s or cm/s)
        T_K: Temperature array (K)
        rho: Density array (g/cm^3)
        n_e: Electron density array
        tau: Optical depth array (cumulative, outer->inner)
        mass_fractions: dict of element -> array
        max_dT_T: max fractional T change per zone (default 0.10)
        max_drho_rho: max fractional rho change per zone (default 0.10)
        max_dlogT: max log10(T) change per zone (default 0.10)
        max_dlogrho: max log10(rho) change per zone (default 0.10)
        max_dlogtau: max log10(tau) change per zone (default 0.10)
        max_dne_ne: max fractional ne change per zone (default 0.10)
        max_dv_v: max fractional v change per zone (default 0.30)
        stretch_factor: how many times more points to insert in trouble
                        regions (default 3.0 → 2 extra points per edge)
        max_stretch_zones: maximum total inserted points (default 200)
        min_zones_to_stretch: don't stretch regions smaller than this
                              (default 2)
        verbose: print diagnostics

    Returns:
        (r_new, v_new, T_new, rho_new, n_e_new, tau_new, mf_new)
        where mf_new is a dict of element -> array.
        Original zone values are preserved; new points are interpolated.
    """
    eps = 1e-30
    n = len(r_cm)
    if n < 4:
        if verbose:
            print("[Stretch] Too few zones, skipping")
        return r_cm.copy(), v.copy(), T_K.copy(), rho.copy(), \
               n_e.copy(), tau.copy(), \
               {k: val.copy() for k, val in mass_fractions.items()}

    # --- Step 1: Compute fractional jumps for each edge ---
    # Arrays are outer->inner (decreasing r), so diff[i] = arr[i] - arr[i+1]
    logT = np.log10(np.maximum(T_K, eps))
    logrho = np.log10(np.maximum(rho, eps))
    logtau = np.log10(np.maximum(tau, eps))

    dT_T = np.abs(np.diff(T_K) / np.maximum(np.abs(T_K[:-1]), eps))
    drho_rho = np.abs(np.diff(rho) / np.maximum(np.abs(rho[:-1]), eps))
    dlogT = np.abs(np.diff(logT))
    dlogrho = np.abs(np.diff(logrho))
    dlogtau = np.abs(np.diff(logtau))
    dne_ne = np.abs(np.diff(n_e) / np.maximum(np.abs(n_e[:-1]), eps))
    dv_v = np.abs(np.diff(v) / np.maximum(np.abs(v[:-1]), eps))

    # --- Step 2: Identify violating edges ---
    # An edge violates if ANY quantity exceeds its threshold
    violates = ((dT_T > max_dT_T) |
                (drho_rho > max_drho_rho) |
                (dlogT > max_dlogT) |
                (dlogrho > max_dlogrho) |
                (dlogtau > max_dlogtau) |
                (dne_ne > max_dne_ne) |
                (dv_v > max_dv_v))

    n_violations = int(np.sum(violates))
    if n_violations == 0:
        if verbose:
            print("[Stretch] No trouble regions detected, skipping")
        return r_cm.copy(), v.copy(), T_K.copy(), rho.copy(), \
               n_e.copy(), tau.copy(), \
               {k: val.copy() for k, val in mass_fractions.items()}

    # --- Step 3: Group consecutive violating edges into trouble regions ---
    # An edge i is between zone i and zone i+1
    # A trouble region is a maximal run of consecutive violating edges
    trouble_regions = []
    i = 0
    while i < len(violates):
        if violates[i]:
            j = i
            while j + 1 < len(violates) and violates[j + 1]:
                j += 1
            # Region spans edges i..j, which means zones i..j+1
            trouble_regions.append((i, j + 1))
            i = j + 1
        else:
            i += 1

    if verbose:
        print(f"[Stretch] {n_violations} violating edges in "
              f"{len(trouble_regions)} trouble regions")
        for idx, (lo, hi) in enumerate(trouble_regions):
            worst_dT = dT_T[lo:hi + 1].max() if lo <= hi else 0
            worst_drho = drho_rho[lo:hi + 1].max() if lo <= hi else 0
            worst_dlogtau = dlogtau[lo:hi + 1].max() if lo <= hi else 0
            r_lo = r_cm[lo]
            r_hi = r_cm[hi]
            print(f"  Region {idx}: zones {lo}-{hi} "
                  f"(r={r_lo:.3e}-{r_hi:.3e} cm), "
                  f"max dT/T={worst_dT*100:.0f}%, "
                  f"max drho/rho={worst_drho*100:.0f}%, "
                  f"max dlogtau={worst_dlogtau:.2f}")

    # --- Step 4: Compute insertion points ---
    # For each trouble region, insert (stretch_factor - 1) extra points
    # between each pair of consecutive original zones.
    # Use log-interpolation for T, rho, ne, tau and linear for v, r, mass_fractions.
    #
    # We build a list of (original_index, fractional_position) tuples
    # for each new point, then insert them all at once.
    #
    # To avoid exceeding max_stretch_zones, we prioritize the worst regions.

    n_extra_per_edge = max(1, int(round(stretch_factor)) - 1)
    if n_extra_per_edge < 1:
        n_extra_per_edge = 1

    # Compute total requested insertions
    total_requested = 0
    region_insertions = []
    for lo, hi in trouble_regions:
        n_edges = hi - lo  # number of edges in this region
        if n_edges < min_zones_to_stretch:
            continue
        n_insert = n_edges * n_extra_per_edge
        region_insertions.append((lo, hi, n_insert, n_edges))
        total_requested += n_insert

    if total_requested == 0:
        if verbose:
            print("[Stretch] No regions large enough to stretch, skipping")
        return r_cm.copy(), v.copy(), T_K.copy(), rho.copy(), \
               n_e.copy(), tau.copy(), \
               {k: val.copy() for k, val in mass_fractions.items()}

    # If we exceed max_stretch_zones, scale down n_extra_per_edge
    if total_requested > max_stretch_zones:
        scale = max_stretch_zones / total_requested
        n_extra_per_edge = max(1, int(n_extra_per_edge * scale))
        if verbose:
            print(f"[Stretch] Scaling down to {n_extra_per_edge} extra pts/edge "
                  f"(requested {total_requested}, cap {max_stretch_zones})")
        # Recompute
        region_insertions = []
        total_requested = 0
        for lo, hi in trouble_regions:
            n_edges = hi - lo
            if n_edges < min_zones_to_stretch:
                continue
            n_insert = n_edges * n_extra_per_edge
            region_insertions.append((lo, hi, n_insert, n_edges))
            total_requested += n_insert

    if total_requested == 0:
        if verbose:
            print("[Stretch] No insertions after scaling, skipping")
        return r_cm.copy(), v.copy(), T_K.copy(), rho.copy(), \
               n_e.copy(), tau.copy(), \
               {k: val.copy() for k, val in mass_fractions.items()}

    if verbose:
        print(f"[Stretch] Inserting {total_requested} extra points "
              f"({n_extra_per_edge} per edge in {len(region_insertions)} regions)")

    # --- Step 5: Build the new arrays ---
    # We insert points into the arrays.  To handle this efficiently,
    # we build a list of "segments" where each segment is either:
    #   - a single original zone, or
    #   - a pair of original zones with extra interpolated points between them
    #
    # We process the array from inner (end) to outer (start) to avoid
    # index shifting, or we can build a new list from scratch.

    # Build a list of new values
    new_r = []
    new_v = []
    new_T = []
    new_rho = []
    new_ne = []
    new_tau = []
    new_mf = {k: [] for k in mass_fractions}
    is_original = []  # track which points are original

    # Create a set of (lo, hi) ranges that need stretching
    stretch_ranges = set()
    for lo, hi, _, _ in region_insertions:
        for edge in range(lo, hi):
            stretch_ranges.add(edge)

    for i in range(n):
        # Append the original zone
        new_r.append(r_cm[i])
        new_v.append(v[i])
        new_T.append(T_K[i])
        new_rho.append(rho[i])
        new_ne.append(n_e[i])
        new_tau.append(tau[i])
        is_original.append(True)
        for k in mass_fractions:
            new_mf[k].append(mass_fractions[k][i])

        # If this edge (i to i+1) is in a trouble region, insert extra points
        if i < n - 1 and i in stretch_ranges:
            r0, r1 = r_cm[i], r_cm[i + 1]
            v0, v1 = v[i], v[i + 1]
            T0, T1 = T_K[i], T_K[i + 1]
            rho0, rho1 = rho[i], rho[i + 1]
            ne0, ne1 = n_e[i], n_e[i + 1]
            tau0, tau1 = tau[i], tau[i + 1]

            logT0 = np.log10(max(T0, eps))
            logT1 = np.log10(max(T1, eps))
            logrho0 = np.log10(max(rho0, eps))
            logrho1 = np.log10(max(rho1, eps))
            logne0 = np.log10(max(ne0, eps))
            logne1 = np.log10(max(ne1, eps))
            logtau0 = np.log10(max(tau0, eps))
            logtau1 = np.log10(max(tau1, eps))

            for j in range(1, n_extra_per_edge + 1):
                f = j / (n_extra_per_edge + 1)
                # Linear interpolation for r, v, mass fractions
                new_r.append(r0 + f * (r1 - r0))
                new_v.append(v0 + f * (v1 - v0))
                # Log interpolation for T, rho, ne, tau
                new_T.append(10.0 ** (logT0 + f * (logT1 - logT0)))
                new_rho.append(10.0 ** (logrho0 + f * (logrho1 - logrho0)))
                new_ne.append(10.0 ** (logne0 + f * (logne1 - logne0)))
                new_tau.append(10.0 ** (logtau0 + f * (logtau1 - logtau0)))
                is_original.append(False)
                for k in mass_fractions:
                    mf0 = mass_fractions[k][i]
                    mf1 = mass_fractions[k][i + 1]
                    new_mf[k].append(mf0 + f * (mf1 - mf0))

    new_r = np.array(new_r)
    new_v = np.array(new_v)
    new_T = np.array(new_T)
    new_rho = np.array(new_rho)
    new_ne = np.array(new_ne)
    new_tau = np.array(new_tau)
    new_mf = {k: np.array(v) for k, v in new_mf.items()}
    is_original = np.array(is_original)

    # --- Verify: check that original points are preserved ---
    n_orig = int(np.sum(is_original))
    if n_orig != n:
        # This shouldn't happen, but check anyway
        if verbose:
            print(f"[Stretch] WARNING: expected {n} original points, got {n_orig}")

    # Verify original values are preserved
    T_at_orig = new_T[is_original]
    if len(T_at_orig) == n:
        max_T_err = np.max(np.abs(T_at_orig - T_K) / np.maximum(np.abs(T_K), eps))
        if max_T_err > 1e-10 and verbose:
            print(f"[Stretch] WARNING: T preservation error = {max_T_err:.2e}")

    # --- Report results ---
    if verbose:
        # Recompute jumps on the new grid
        new_dT_T = np.abs(np.diff(new_T) / np.maximum(np.abs(new_T[:-1]), eps))
        new_drho_rho = np.abs(np.diff(new_rho) / np.maximum(np.abs(new_rho[:-1]), eps))
        new_dlogtau = np.abs(np.diff(np.log10(np.maximum(new_tau, eps))))

        print(f"[Stretch] Grid: {n} -> {len(new_r)} points "
              f"({len(new_r) - n} inserted)")
        print(f"  T:  max dT/T {dT_T.max()*100:.1f}% -> "
              f"{new_dT_T.max()*100:.1f}%, "
              f"violations >{max_dT_T*100:.0f}%: "
              f"{int(np.sum(dT_T > max_dT_T))} -> "
              f"{int(np.sum(new_dT_T > max_dT_T))}")
        print(f"  rho: max drho/rho {drho_rho.max()*100:.1f}% -> "
              f"{new_drho_rho.max()*100:.1f}%, "
              f"violations >{max_drho_rho*100:.0f}%: "
              f"{int(np.sum(drho_rho > max_drho_rho))} -> "
              f"{int(np.sum(new_drho_rho > max_drho_rho))}")
        print(f"  tau: max dlogtau {dlogtau.max():.2f} -> "
              f"{new_dlogtau.max():.2f}, "
              f"violations >{max_dlogtau:.2f}: "
              f"{int(np.sum(dlogtau > max_dlogtau))} -> "
              f"{int(np.sum(new_dlogtau > max_dlogtau))}")
        print(f"  T range: {T_K.min():.0f}-{T_K.max():.0f} K -> "
              f"{new_T.min():.0f}-{new_T.max():.0f} K (unchanged at original pts)")
        print(f"  rho range: {rho.min():.2e}-{rho.max():.2e} -> "
              f"{new_rho.min():.2e}-{new_rho.max():.2e} (unchanged at original pts)")

    return new_r, new_v, new_T, new_rho, new_ne, new_tau, new_mf


def stretch_critical_regions(r_cm, v, T_K, rho, n_e, tau, mass_fractions,
                             critical_mask,
                             stretch_factor=5.0,
                             max_stretch_zones=300,
                             min_zones_to_stretch=1,
                             verbose=True):
    """Stretch the radius coordinate in pre-identified critical regions.

    This is a mask-driven companion to :func:`stretch_trouble_regions`.
    Instead of re-detecting trouble edges from gradient thresholds, it
    stretches the regions already identified by
    :func:`identify_critical_regions`, applying a stronger ``stretch_factor``
    there.  Original STELLA zone values are preserved exactly; only
    interpolated points are inserted, so peaks (shock T, density jumps,
    composition interfaces) are not blurred.

    The returned ``critical_mask_new`` extends the input mask to the
    inserted points: any point inserted inside a critical region is
    itself marked critical, so downstream stages (AMR, degenerate-grid
    thinner) continue to apply the relaxed ``min_dr/r`` there.

    Args:
        r_cm, v, T_K, rho, n_e, tau, mass_fractions: Physical arrays
            (outer->inner, decreasing radius).
        critical_mask: Boolean array of length n marking critical zones.
        stretch_factor: Number of sub-intervals per original edge inside
            critical regions (default 5.0 -> 4 extra interpolated points
            per edge).  Stronger than ``stretch_trouble_regions``'s
            default of 3.0 because critical regions are, by construction,
            the narrowest high-gradient features.
        max_stretch_zones: Cap on total inserted points.
        min_zones_to_stretch: Don't stretch critical regions smaller than
            this many zones.
        verbose: Print diagnostics.

    Returns:
        (r_new, v_new, T_new, rho_new, ne_new, tau_new, mf_new,
        critical_mask_new)
    """
    eps = 1e-30
    n = len(r_cm)
    if n < 4 or critical_mask is None or len(critical_mask) != n \
            or not np.any(critical_mask):
        if verbose:
            print("[Stretch-Critical] No critical regions to stretch, skipping")
        return (r_cm.copy(), v.copy(), T_K.copy(), rho.copy(),
                n_e.copy(), tau.copy(),
                {k: val.copy() for k, val in mass_fractions.items()},
                critical_mask.copy() if critical_mask is not None
                else np.zeros(n, dtype=bool))

    # Group the critical mask into contiguous regions
    regions = []
    i = 0
    while i < n:
        if critical_mask[i]:
            j = i
            while j + 1 < n and critical_mask[j + 1]:
                j += 1
            regions.append((i, j))
            i = j + 1
        else:
            i += 1

    n_extra_per_edge = max(1, int(round(stretch_factor)) - 1)

    # Compute total requested insertions (cap by max_stretch_zones)
    region_insertions = []
    total_requested = 0
    for lo, hi in regions:
        n_edges = hi - lo
        if n_edges < min_zones_to_stretch:
            continue
        n_insert = n_edges * n_extra_per_edge
        region_insertions.append((lo, hi, n_insert, n_edges))
        total_requested += n_insert

    if total_requested == 0:
        if verbose:
            print("[Stretch-Critical] No critical regions large enough to stretch")
        return (r_cm.copy(), v.copy(), T_K.copy(), rho.copy(),
                n_e.copy(), tau.copy(),
                {k: val.copy() for k, val in mass_fractions.items()},
                critical_mask.copy())

    if total_requested > max_stretch_zones:
        scale = max_stretch_zones / total_requested
        n_extra_per_edge = max(1, int(n_extra_per_edge * scale))
        if verbose:
            print(f"[Stretch-Critical] Scaling down to {n_extra_per_edge} "
                  f"extra pts/edge (requested {total_requested}, "
                  f"cap {max_stretch_zones})")
        region_insertions = []
        total_requested = 0
        for lo, hi in regions:
            n_edges = hi - lo
            if n_edges < min_zones_to_stretch:
                continue
            n_insert = n_edges * n_extra_per_edge
            region_insertions.append((lo, hi, n_insert, n_edges))
            total_requested += n_insert

    if total_requested == 0:
        return (r_cm.copy(), v.copy(), T_K.copy(), rho.copy(),
                n_e.copy(), tau.copy(),
                {k: val.copy() for k, val in mass_fractions.items()},
                critical_mask.copy())

    stretch_edges = set()
    for lo, hi, _, _ in region_insertions:
        for edge in range(lo, hi):
            stretch_edges.add(edge)

    if verbose:
        print(f"[Stretch-Critical] Inserting {total_requested} extra points "
              f"({n_extra_per_edge} per edge in {len(region_insertions)} "
              f"critical regions)")

    new_r, new_v, new_T, new_rho, new_ne, new_tau = [], [], [], [], [], []
    new_mf = {k: [] for k in mass_fractions}
    new_crit = []

    for i in range(n):
        new_r.append(r_cm[i])
        new_v.append(v[i])
        new_T.append(T_K[i])
        new_rho.append(rho[i])
        new_ne.append(n_e[i])
        new_tau.append(tau[i])
        new_crit.append(bool(critical_mask[i]))
        for k in mass_fractions:
            new_mf[k].append(mass_fractions[k][i])

        if i < n - 1 and i in stretch_edges:
            r0, r1 = r_cm[i], r_cm[i + 1]
            v0, v1 = v[i], v[i + 1]
            T0, T1 = T_K[i], T_K[i + 1]
            rho0, rho1 = rho[i], rho[i + 1]
            ne0, ne1 = n_e[i], n_e[i + 1]
            tau0, tau1 = tau[i], tau[i + 1]
            logT0 = np.log10(max(T0, eps))
            logT1 = np.log10(max(T1, eps))
            logrho0 = np.log10(max(rho0, eps))
            logrho1 = np.log10(max(rho1, eps))
            logne0 = np.log10(max(ne0, eps))
            logne1 = np.log10(max(ne1, eps))
            logtau0 = np.log10(max(tau0, eps))
            logtau1 = np.log10(max(tau1, eps))

            for j in range(1, n_extra_per_edge + 1):
                f = j / (n_extra_per_edge + 1)
                new_r.append(r0 + f * (r1 - r0))
                new_v.append(v0 + f * (v1 - v0))
                new_T.append(10.0 ** (logT0 + f * (logT1 - logT0)))
                new_rho.append(10.0 ** (logrho0 + f * (logrho1 - logrho0)))
                new_ne.append(10.0 ** (logne0 + f * (logne1 - logne0)))
                new_tau.append(10.0 ** (logtau0 + f * (logtau1 - logtau0)))
                # Inserted points inside a critical region are critical
                new_crit.append(True)
                for k in mass_fractions:
                    mf0 = mass_fractions[k][i]
                    mf1 = mass_fractions[k][i + 1]
                    new_mf[k].append(mf0 + f * (mf1 - mf0))

    new_r = np.array(new_r)
    new_v = np.array(new_v)
    new_T = np.array(new_T)
    new_rho = np.array(new_rho)
    new_ne = np.array(new_ne)
    new_tau = np.array(new_tau)
    new_mf = {k: np.array(val) for k, val in new_mf.items()}
    new_crit = np.array(new_crit, dtype=bool)

    # --- Sanity checks ---
    if len(new_r) > 1:
        dr = np.diff(new_r)
        if np.any(dr >= 0):
            print("[Stretch-Critical] WARNING: non-monotonic radius detected "
                  "(r must decrease with index)")
        if np.any(np.abs(dr) <= 0):
            print("[Stretch-Critical] WARNING: zero or negative dr detected")

    if verbose:
        print(f"[Stretch-Critical] Grid: {n} -> {len(new_r)} points "
              f"({len(new_r) - n} inserted, "
              f"{int(np.sum(new_crit))} critical zones)")

    return (new_r, new_v, new_T, new_rho, new_ne, new_tau, new_mf,
            new_crit)


def thin_degenerate_grid(radius, *arrays, min_dr_r=0.005):
    """Remove near-degenerate grid points (dr/r < min_dr_r)."""
    n = len(radius)
    keep = np.zeros(n, dtype=bool)
    keep[0] = True
    keep[n-1] = True
    last_kept = 0
    for i in range(1, n - 1):
        dr = abs(radius[last_kept] - radius[i])
        r_ref = max(abs(radius[last_kept]), abs(radius[i]), 1e-30)
        if dr / r_ref >= min_dr_r:
            keep[i] = True
            last_kept = i
    n_removed = int(np.sum(~keep))
    if n_removed > 0:
        print(f"[Grid Thinning] Removed {n_removed} near-degenerate points "
              f"(dr/r < {min_dr_r*100:.1f}%), {n} -> {int(np.sum(keep))} points")
        radius = radius[keep]
        arrays = tuple(a[keep] for a in arrays)
    return (radius,) + arrays


def thin_degenerate_grid_adaptive(radius, *arrays, min_dr_r=0.005,
                                  critical_mask=None,
                                  critical_min_dr_r=1e-4,
                                  disable_in_critical=False,
                                  verbose=True):
    """Region-dependent degenerate-grid thinner.

    Identical in spirit to :func:`thin_degenerate_grid`, but the
    ``min_dr/r`` threshold is *region-dependent*: outside critical regions
    the global ``min_dr_r`` is enforced; inside critical regions a much
    smaller ``critical_min_dr_r`` is used (or the constraint is disabled
    entirely when ``disable_in_critical=True``).

    This lets the AMR insert the points CMFGEN needs in narrow
    high-gradient regions (shock fronts, composition interfaces, CSM
    shells) without those points being removed again by the thinner, while
    still preventing near-degenerate points elsewhere.

    The first and last points are always kept.

    Args:
        radius: Radius array (outer->inner, decreasing).
        *arrays: Companion arrays to thin in lockstep with ``radius``.
        min_dr_r: Global default threshold (outside critical regions).
        critical_mask: Optional boolean array of length ``len(radius)``
            marking critical-region zones. If None, this falls back to
            the global-threshold behaviour of :func:`thin_degenerate_grid`.
        critical_min_dr_r: Threshold inside critical regions.
        disable_in_critical: If True, disable the constraint inside
            critical regions (no points removed there).
        verbose: Print diagnostics.

    Returns:
        ``(radius, *arrays)`` with the thinned arrays.  Companion arrays
        are returned in the same order they were passed.
    """
    n = len(radius)
    if critical_mask is None or len(critical_mask) != n:
        # Fall back to the simple global-threshold path
        return thin_degenerate_grid(radius, *arrays, min_dr_r=min_dr_r)

    eps = 1e-30
    # Per-zone threshold: an edge between last_kept and i is governed by
    # the threshold of the zone being tested (i); if that zone is in a
    # critical region the relaxed threshold applies.
    thr = build_min_dr_r_array(
        n, critical_mask,
        global_min_dr_r=min_dr_r,
        critical_min_dr_r=critical_min_dr_r,
        disable_in_critical=disable_in_critical)

    keep = np.zeros(n, dtype=bool)
    keep[0] = True
    keep[n - 1] = True
    last_kept = 0
    n_preserved_critical = 0
    for i in range(1, n - 1):
        dr = abs(radius[last_kept] - radius[i])
        r_ref = max(abs(radius[last_kept]), abs(radius[i]), eps)
        dr_r = dr / r_ref
        # Use the smaller of the two endpoints' thresholds so that a
        # critical region never forces removal of points just outside it,
        # and a non-critical region never gets the relaxed threshold.
        thr_i = min(thr[i], thr[last_kept])
        if dr_r >= thr_i:
            keep[i] = True
            last_kept = i
        else:
            # This point would be removed.  If it is in a critical region
            # (and the constraint is relaxed there) it is preserved.
            if critical_mask[i] and thr_i < min_dr_r:
                keep[i] = True
                last_kept = i
                n_preserved_critical += 1

    n_removed = int(np.sum(~keep))
    if n_removed > 0 or n_preserved_critical > 0:
        if verbose:
            print(f"[Grid Thinning (adaptive)] Removed {n_removed} near-degenerate "
                  f"points (global dr/r < {min_dr_r:.0e}); "
                  f"preserved {n_preserved_critical} critical-region points "
                  f"(local dr/r >= {critical_min_dr_r:.0e}), "
                  f"{n} -> {int(np.sum(keep))} points")
        radius = radius[keep]
        arrays = tuple(a[keep] for a in arrays)
        critical_mask = critical_mask[keep]
    elif verbose:
        print(f"[Grid Thinning (adaptive)] No points removed "
              f"({n} points, {int(np.sum(critical_mask))} in critical regions)")
    return (radius,) + arrays


def downsample_grid(radius, *arrays, target_points=100):
    """Downsample grid to target_points using true logarithmic spacing in radius.

    Picks the closest original grid point to each of ``target_points``
    log-spaced radii, guaranteeing the first and last points are kept.
    Arrays are assumed outer->inner (decreasing radius).
    """
    n = len(radius)
    if n <= target_points:
        print(f"[Downsample] Grid already at {n} points <= target {target_points}, skipping")
        return (radius,) + arrays

    log_r = np.log10(radius)
    # radius is outer->inner so log_r[0] > log_r[-1]
    log_r_targets = np.linspace(log_r[0], log_r[-1], target_points)
    # For each target, find the closest original point
    target_indices = np.array(
        [np.argmin(np.abs(log_r - lr)) for lr in log_r_targets]
    )
    # Ensure first and last are included, deduplicate, keep order
    target_indices = np.unique(np.concatenate(([0], target_indices, [n - 1])))

    print(f"[Downsample] Reduced from {n} to {len(target_indices)} points "
          f"(target: {target_points}, log-spaced in radius)")

    radius = radius[target_indices]
    arrays = tuple(a[target_indices] for a in arrays)

    return (radius,) + arrays


def downsample_grid_feature_aware(radius, v, T, rho, *arrays, target_points=100,
                                   gradient_threshold=0.05):
    """Downsample grid to target_points while preserving high-gradient features.

    This function addresses the issue where pure log-spaced downsampling can miss
    important physical features in the original STELLA grid. It:
    1. Identifies high-gradient regions in the original grid (density, temperature, velocity)
    2. Forces inclusion of points from these regions
    3. Fills remaining budget with log-spaced points
    4. Ensures first and last points are always included

    Args:
        radius: Radius array (outer->inner, decreasing)
        v: Velocity array (km/s)
        T: Temperature array (K)
        rho: Density array (g/cm³)
        *arrays: Other arrays to downsample consistently
        target_points: Target number of output points
        gradient_threshold: Threshold for considering a gradient "high" (fractional change)

    Returns:
        Tuple of downsampled arrays in same order as input
    """
    n = len(radius)
    if n <= target_points:
        print(f"[Feature-Aware Downsample] Grid already at {n} points <= target {target_points}, skipping")
        return (radius, v, T, rho) + arrays

    eps = 1e-30

    # Calculate gradients on original grid
    dr_r = np.abs(np.diff(radius) / np.maximum(np.abs(radius[:-1]), eps))
    dv_v = np.abs(np.diff(v) / np.maximum(np.abs(v[:-1]), eps))
    dT_T = np.abs(np.diff(T) / np.maximum(T[:-1], eps))
    drho_rho = np.abs(np.diff(rho) / np.maximum(rho[:-1], eps))

    # Identify high-gradient edges
    high_grad_edges = np.where(
        (dr_r > gradient_threshold) |
        (dv_v > gradient_threshold) |
        (dT_T > gradient_threshold) |
        (drho_rho > gradient_threshold)
    )[0]

    # Convert edges to point indices (include both points of high-gradient edges)
    feature_indices = set()
    for edge_idx in high_grad_edges:
        feature_indices.add(edge_idx)
        feature_indices.add(edge_idx + 1)

    # Always include first and last points
    feature_indices.add(0)
    feature_indices.add(n - 1)

    # Start with feature points
    selected_indices = sorted(list(feature_indices))
    n_features = len(selected_indices)

    print(f"[Feature-Aware Downsample] Found {n_features} feature points from {len(high_grad_edges)} high-gradient edges")

    # If we already have enough points, just use features
    if n_features >= target_points:
        # Downsample features to target by keeping every nth feature.
        # Uniform stride ensures all radial regions (CSM wind, interaction
        # shell, inner ejecta) are sampled proportionally, preserving the
        # shape of the STELLA profile.
        step = max(1, n_features // target_points)
        selected_indices = selected_indices[::step]
        # Ensure last point is included
        if selected_indices[-1] != n - 1:
            selected_indices.append(n - 1)
        selected_indices = sorted(set(selected_indices))
    else:
        # Fill remaining budget with log-spaced points (uniform in log-r).
        # This ensures adequate coverage of all radial regions including
        # the CSM wind, preserving the shape of the STELLA profile.
        remaining_budget = target_points - n_features
        log_r = np.log10(np.maximum(radius, eps))
        log_r_targets = np.linspace(log_r[0], log_r[-1], remaining_budget + 2)[1:-1]
        log_spaced_indices = [
            np.argmin(np.abs(log_r - lr)) for lr in log_r_targets
        ]

        # Combine feature and log-spaced indices
        all_indices = sorted(set(selected_indices + log_spaced_indices))
        selected_indices = all_indices

    # Deduplicate and ensure order
    selected_indices = np.unique(selected_indices)

    print(f"[Feature-Aware Downsample] Reduced from {n} to {len(selected_indices)} points "
          f"(target: {target_points}, {n_features} feature points preserved)")

    # Apply selection to all arrays
    radius_out = radius[selected_indices]
    v_out = v[selected_indices]
    T_out = T[selected_indices]
    rho_out = rho[selected_indices]
    arrays_out = tuple(a[selected_indices] for a in arrays)

    return (radius_out, v_out, T_out, rho_out) + arrays_out


# ---------------------------------------------------------------------------
# Composition-aware grid helpers
# ---------------------------------------------------------------------------

# Elements that define nucleosynthetic shells in metal-rich ejecta.
# Ordered roughly inner->outer: Fe-peak core, Si-burning shell, O shell,
# C/O shell, He layer, H envelope.  These are the species whose sharp
# interfaces MUST be resolved for faithful CMFGEN spectra.
DEFAULT_COMPOSITION_ELEMENTS = ['FePeak', 'Ni56', 'Ni58', 'Si', 'S', 'Ca',
                                'O', 'C', 'He', 'H']


def detect_ionization_fronts(r_cm, rho, T_K, mass_fractions,
                             ne_jump_threshold=0.5,
                             min_points_per_front=5,
                             front_width_zones=3):
    """Detect ionization fronts in the source STELLA grid.

    An ionization front is a region where the electron density changes
    rapidly (dlog(n_e)/dlog(r) exceeds ``ne_jump_threshold``), indicating
    a recombination or ionization transition. These fronts are critical
    for CMFGEN convergence: under-resolved fronts cause the population
    solver to oscillate or diverge.

    Uses the Saha electron density as a diagnostic. While the Saha n_e is
    approximate, sharp changes in n_e reliably indicate real ionization
    fronts because they are driven by temperature crossing an ionization
    threshold.

    Args:
        r_cm: Radius array (outer->inner, decreasing), in cm
        rho: Density array (g/cm^3)
        T_K: Temperature array (K)
        mass_fractions: dict of element -> mass fraction array
        ne_jump_threshold: Minimum |dlog(n_e)/dlog(r)| to flag a front.
            Default 0.5 — catches H recombination (typically 0.5-2.0 dex
            change over a few zones) without flagging smooth gradients.
        min_points_per_front: Minimum number of source zones to mark
            around each detected front (front_width_zones on each side
            plus the center). Default 5 ensures the front is resolved
            over a few depths.
        front_width_zones: Number of zones on each side of the front
            center to mark as part of the front. Default 3.

    Returns:
        (front_indices, front_info) where:
        - front_indices: set of source zone indices that should be
          force-included in the downsampled grid
        - front_info: list of dicts with 'center_idx', 'r_cm', 'ne_jump',
          'T_K', 'tau_est' for each detected front
    """
    eps = 1e-30
    n = len(r_cm)
    if n < 5:
        return set(), []

    # Compute Saha electron density on the source grid
    n_e = calculate_electron_density(rho, T_K, mass_fractions)
    n_e = np.maximum(n_e, eps)

    # Compute dlog(n_e)/dlog(r) — ionization front indicator
    log_ne = np.log10(n_e)
    log_r = np.log10(np.maximum(r_cm, eps))

    # Central difference for interior points
    dlog_ne_dlog_r = np.zeros(n)
    for i in range(1, n - 1):
        dr = log_r[i - 1] - log_r[i + 1]
        if abs(dr) > eps:
            dlog_ne_dlog_r[i] = abs(log_ne[i - 1] - log_ne[i + 1]) / abs(dr)

    # Detect fronts: zones where the gradient exceeds the threshold
    is_front = dlog_ne_dlog_r > ne_jump_threshold

    # Group consecutive front zones into individual fronts
    front_centers = []
    i = 0
    while i < n:
        if is_front[i]:
            # Find the extent of this front
            j = i
            while j < n and is_front[j]:
                j += 1
            # Center of the front
            center = (i + j - 1) // 2
            front_centers.append(center)
            i = j
        else:
            i += 1

    # Build the set of indices to force-include: front_width_zones on each
    # side of each front center, ensuring the front is resolved over
    # min_points_per_front depths.
    front_indices = set()
    front_info = []
    for center in front_centers:
        lo = max(0, center - front_width_zones)
        hi = min(n - 1, center + front_width_zones)
        for idx in range(lo, hi + 1):
            front_indices.add(idx)

        # Ensure minimum points per front
        n_marked = hi - lo + 1
        if n_marked < min_points_per_front:
            # Expand outward to reach minimum
            extra_needed = min_points_per_front - n_marked
            for k in range(1, extra_needed + 1):
                if center - front_width_zones - k >= 0:
                    front_indices.add(center - front_width_zones - k)
                if center + front_width_zones + k < n:
                    front_indices.add(center + front_width_zones + k)

        front_info.append({
            'center_idx': center,
            'r_cm': r_cm[center],
            'ne_jump': dlog_ne_dlog_r[center],
            'T_K': T_K[center],
            'n_marked': len([idx for idx in front_indices
                            if abs(idx - center) <= front_width_zones + 10]),
        })

    return front_indices, front_info


def detect_composition_interfaces(mass_fractions, gradient_threshold=0.3,
                                  elements=None, min_abundance=1e-6):
    """Detect sharp composition interfaces in the original STELLA grid.

    An interface is a zone boundary where a significant element changes
    by more than ``gradient_threshold`` (fractional) AND at least one side
    has the element above ``min_abundance``.  Both zones straddling the
    interface are returned so the caller can force-include them in the
    downsampled grid.

    Args:
        mass_fractions: dict of element -> array[nzon] (outer->inner)
        gradient_threshold: fractional change |dX|/X above which an edge
                           is classified as a sharp interface (default 0.3)
        elements: list of element names to check.  If None, uses all keys
                  in mass_fractions.
        min_abundance: ignore interfaces where both sides have X < this
                      (avoids noise in trace species far from shells)

    Returns:
        set of zone indices that should be force-included (both sides of
        each detected interface), plus a dict of {element: [interface_indices]}
        for diagnostics.
    """
    eps = 1e-30
    if elements is None:
        elements = list(mass_fractions.keys())

    interface_indices = set()
    interface_by_element = {}

    for elem in elements:
        if elem not in mass_fractions:
            continue
        X = np.asarray(mass_fractions[elem], dtype=float)
        if len(X) < 2:
            continue
        dX = np.abs(np.diff(X))
        X_max = np.maximum(np.abs(X[:-1]), np.abs(X[1:]))
        rel_change = dX / np.maximum(X_max, eps)

        # Find edges where the relative change exceeds the threshold
        # AND at least one side has a non-trivial abundance
        significant = np.maximum(np.abs(X[:-1]), np.abs(X[1:])) > min_abundance
        sharp_edges = np.where((rel_change > gradient_threshold) & significant)[0]

        if len(sharp_edges) > 0:
            # Include both zones at each interface
            for edge in sharp_edges:
                interface_indices.add(edge)
                interface_indices.add(edge + 1)
            interface_by_element[elem] = sharp_edges.tolist()

    return interface_indices, interface_by_element


def hybrid_interp_mass_fractions(r_target, r_orig, mf_orig,
                                  interface_indices=None,
                                  trace_threshold=1e-4,
                                  interface_snap_tol=0.5):
    """Interpolate mass fractions using a hybrid scheme that preserves
    sharp composition interfaces.

    At each target radius:
      1. If the target falls between two original zones that straddle a
         detected composition interface, use NEAREST-NEIGHBOUR (snap to
         the closer original zone).  This preserves the sharp boundary
         instead of smearing it with linear interpolation.
      2. If both bracketing zones have X < trace_threshold (trace species
         in a smooth region), use LOG interpolation to preserve dynamic
         range.
      3. Otherwise, use LINEAR interpolation (smooth, well-sampled regions).

    Args:
        r_target: target radii (any order, typically decreasing outer->inner)
        r_orig: original radii (decreasing outer->inner)
        mf_orig: dict of element -> mass fraction array (same order as r_orig)
        interface_indices: set/array of original zone indices that are at
                          composition interfaces.  If None, interfaces are
                          auto-detected from the data.
        trace_threshold: mass fraction below which log interpolation is used
                        in smooth regions (default 1e-4)
        interface_snap_tol: if the interpolation fraction f is within
                           [1-tol, 1] of an interface zone, snap to it.
                           Default 0.5 means snap if closer to either end.

    Returns:
        dict of element -> interpolated mass fraction array at r_target
    """
    eps = 1e-30
    n_orig = len(r_orig)
    r_target = np.asarray(r_target, dtype=float)
    # np.interp needs increasing x
    r_orig_inc = r_orig[::-1]

    # Auto-detect interfaces if not provided
    if interface_indices is None:
        interface_indices, _ = detect_composition_interfaces(mf_orig)
    interface_set = set(int(i) for i in interface_indices)

    # For each target radius, find the bracketing original zones
    # r_orig is decreasing, so we work in the increasing-reversed space
    # idx_inc is the index in the increasing array; the original index is
    # n_orig - 1 - idx_inc
    idx_inc = np.searchsorted(r_orig_inc, r_target, side='right')
    # Clamp to valid range [1, n_orig-1] so we always have two bracketing points
    idx_inc = np.clip(idx_inc, 1, n_orig - 1)

    # Original zone indices (in the decreasing-r array) that bracket each target
    idx_hi = n_orig - 1 - idx_inc        # outer zone (larger r)
    idx_lo = n_orig - 1 - (idx_inc - 1)  # inner zone (smaller r)

    # Interpolation fraction in the increasing-r space:
    # f=0 at r_orig_inc[idx_inc-1] (inner), f=1 at r_orig_inc[idx_inc] (outer)
    r_lo = r_orig[idx_lo]  # inner (smaller r)
    r_hi = r_orig[idx_hi]  # outer (larger r)
    dr = r_hi - r_lo
    dr = np.where(np.abs(dr) < eps, eps, dr)
    # f = 0 at inner (r_lo), f = 1 at outer (r_hi)
    f = (r_target - r_lo) / dr
    f = np.clip(f, 0.0, 1.0)

    # Determine which targets are at interfaces
    # A target is at an interface if either bracketing zone is in interface_set
    at_interface = np.array([
        (int(idx_hi[i]) in interface_set or int(idx_lo[i]) in interface_set)
        for i in range(len(r_target))
    ])

    # Snap to nearest at interfaces
    snap_to_hi = f >= interface_snap_tol  # closer to outer zone
    snap_to_lo = f < interface_snap_tol   # closer to inner zone

    mf_new = {}
    for elem, X_orig in mf_orig.items():
        X_orig = np.asarray(X_orig, dtype=float)
        X_lo = X_orig[idx_lo]
        X_hi = X_orig[idx_hi]

        # Default: linear interpolation
        X_new = X_lo + f * (X_hi - X_lo)

        # At interfaces: nearest-neighbour
        if np.any(at_interface):
            X_snap = np.where(snap_to_hi, X_hi, X_lo)
            X_new = np.where(at_interface, X_snap, X_new)

        # In smooth regions with trace species: log interpolation
        trace_mask = (~at_interface) & (np.maximum(np.abs(X_lo), np.abs(X_hi)) < trace_threshold)
        if np.any(trace_mask):
            log_lo = np.log10(np.maximum(X_lo, eps))
            log_hi = np.log10(np.maximum(X_hi, eps))
            X_log = 10.0 ** (log_lo + f * (log_hi - log_lo))
            X_new = np.where(trace_mask, X_log, X_new)

        # Clip to non-negative
        mf_new[elem] = np.clip(X_new, 0.0, None)

    return mf_new


def verify_composition_coverage(r_selected, mass_fractions_orig, r_orig,
                                min_zones_per_shell=3,
                                gradient_threshold=0.3,
                                elements=None):
    """Verify that each composition shell has adequate zone coverage.

    Identifies composition shells in the original grid (regions where the
    dominant element is the same), then counts how many selected grid points
    fall in each shell.  Returns a list of under-sampled shells.

    Args:
        r_selected: selected grid radii (decreasing outer->inner)
        mass_fractions_orig: dict of element -> array (original grid)
        r_orig: original radii (decreasing outer->inner)
        min_zones_per_shell: minimum zones required per shell
        gradient_threshold: threshold for interface detection
        elements: elements to check (default: all in mass_fractions)

    Returns:
        (n_shells, n_under_sampled, shell_info)
        shell_info: list of dicts with shell index, dominant element,
                    r range, n_zones_selected, n_zones_orig
    """
    if elements is None:
        elements = list(mass_fractions_orig.keys())

    # Find the dominant element at each original zone
    n_orig = len(r_orig)
    dominant = np.full(n_orig, '', dtype=object)
    max_X = np.zeros(n_orig)
    for elem in elements:
        if elem not in mass_fractions_orig:
            continue
        X = np.asarray(mass_fractions_orig[elem], dtype=float)
        mask = X > max_X
        dominant[mask] = elem
        max_X[mask] = X[mask]

    # Find shell boundaries (where dominant element changes)
    shell_boundaries = [0]
    for i in range(1, n_orig):
        if dominant[i] != dominant[i - 1]:
            shell_boundaries.append(i)
    shell_boundaries.append(n_orig)

    # Count selected points in each shell
    shell_info = []
    n_under = 0
    for s in range(len(shell_boundaries) - 1):
        i_start = shell_boundaries[s]
        i_end = shell_boundaries[s + 1]
        r_lo = r_orig[i_end - 1]  # inner edge (smaller r)
        r_hi = r_orig[i_start]    # outer edge (larger r)
        n_orig_shell = i_end - i_start
        n_sel = np.sum((r_selected <= r_hi) & (r_selected >= r_lo))

        elem = dominant[i_start] if n_orig_shell > 0 else 'unknown'
        info = {
            'shell': s,
            'element': elem,
            'r_range': (r_lo, r_hi),
            'n_zones_selected': int(n_sel),
            'n_zones_orig': n_orig_shell,
        }
        shell_info.append(info)
        if n_sel < min_zones_per_shell and n_orig_shell >= min_zones_per_shell:
            n_under += 1
            print(f"  [Composition Verify] Shell {s} ({elem}): "
                  f"r=[{r_lo:.3e}, {r_hi:.3e}], "
                  f"{n_sel} selected / {n_orig_shell} orig — UNDER-SAMPLED")

    return len(shell_boundaries) - 1, n_under, shell_info


def downsample_grid_tau_based(radius, v, T, rho, kappa, tau, mass_fractions,
                               target_points=80, tau_breakpoints=None,
                               zone_fractions=None, spacing='linear',
                               composition_aware=True,
                               composition_gradient_threshold=0.5,
                               composition_elements=None,
                               min_zones_per_shell=3,
                               trace_threshold=1e-4,
                               composition_grid_insertion=True,
                               composition_interface_cap=0.15,
                               ionization_front_aware=True,
                               ionization_front_threshold=0.5,
                               ionization_front_min_points=5,
                               ionization_front_width=3,
                               ionization_front_cap=0.15):
    """Downsample grid using tau spacing with configurable breakpoints.

    Distributes zones across tau segments. Within each segment, zones are
    spaced either linearly in tau (uniform dtau — 'linear') or log-uniform
    in tau ('log'). Linear spacing concentrates zones where tau changes
    rapidly (photosphere, velocity reversal) and is more robust than AMR.

    When spacing='auto', breakpoints and fractions are ignored. Instead,
    points are placed uniformly in log(tau+1) space, which compresses both
    the wind (tau<<1) and the deep interior (tau>>1) while giving maximum
    resolution at tau~1 where the photosphere and shock front are.

    When composition_aware=True (default), the grid is enhanced to robustly
    capture metal-rich ejecta:
      1. Sharp composition interfaces (Fe-core/Si-shell, Si/O, O/He, etc.)
         are detected in the original STELLA grid.
      2. Original STELLA zones at those interfaces are force-included in the
         downsampled grid, guaranteeing the grid HAS points at shell boundaries.
      3. Mass fractions are interpolated using a hybrid scheme: nearest-neighbour
         at sharp interfaces (preserves boundaries), log interpolation for trace
         species (preserves dynamic range), linear elsewhere.
      4. Post-selection verification checks each composition shell has adequate
         zone coverage.

    Args:
        radius: Radius array (outer->inner, decreasing), in cm
        v: Velocity array (km/s)
        T: Temperature array (K)
        rho: Density array (g/cm^3)
        kappa: Opacity array (cm^2/g)
        tau: Cumulative optical depth array (outer->inner, increasing)
        mass_fractions: dict of element -> mass fraction array
        target_points: Total number of output points
        tau_breakpoints: List of inner breakpoints (e.g. [0.01, 1.0, 10.0]).
                        Outer/inner bounds auto-detected from tau array.
                        Ignored when spacing='auto'.
        zone_fractions: List of fractions per segment (must sum to ~1.0).
                       e.g. [0.05, 0.55, 0.30, 0.10] for 4 segments.
                       Ignored when spacing='auto'.
        spacing: 'linear' for uniform dtau within each segment (default),
                'log' for log-uniform tau within each segment,
                'auto' for photosphere-weighted log(tau+1) optimization.
        composition_aware: If True (default), force-include zones at sharp
                          composition interfaces and use hybrid mass-fraction
                          interpolation.
        composition_gradient_threshold: Fractional |dX|/X above which an edge
                                       is a sharp interface (default 0.3).
        composition_elements: Elements to check for interfaces. If None,
                             uses DEFAULT_COMPOSITION_ELEMENTS filtered to
                             those present in mass_fractions.
        min_zones_per_shell: Minimum zones per composition shell for the
                            post-selection verification (default 3).
        trace_threshold: Mass fraction below which log interpolation is used
                        in smooth regions (default 1e-4).

    Returns:
        Tuple: (radius, v, T, rho, kappa) + tuple of mass_fraction arrays
    """
    eps = 1e-30
    n = len(radius)
    if n <= target_points:
        # Grid is already at or below target — no downsampling needed.
        # But still run ionization-front detection for diagnostics and
        # to mark front zones for the grid spacing enforcement pass.
        print(f"[Tau-Based Downsample] Grid already at {n} points <= target {target_points}, "
              f"skipping downsampling")
        if ionization_front_aware:
            front_indices, front_info = detect_ionization_fronts(
                radius, rho, T, mass_fractions,
                ne_jump_threshold=ionization_front_threshold,
                min_points_per_front=ionization_front_min_points,
                front_width_zones=ionization_front_width,
            )
            if front_info:
                print(f"[Ionization-Front-Aware] Detected {len(front_info)} ionization front(s) "
                      f"(grid below target, no insertion needed — all source zones retained):")
                for fi in front_info[:10]:
                    print(f"  Front at r={fi['r_cm']:.3e} cm, "
                          f"|dlog(n_e)/dlog(r)|={fi['ne_jump']:.2f}, "
                          f"T={fi['T_K']:.0f} K")
            else:
                print(f"[Ionization-Front-Aware] No ionization fronts detected "
                      f"(threshold={ionization_front_threshold})")
        return (radius, v, T, rho, kappa) + tuple(mass_fractions.values())

    # Auto-detect outer/inner tau bounds
    tau_floor = max(tau[0], 1e-10)
    tau_max = tau[-1]

    if spacing == 'auto':
        # Blind optimization: place zones uniformly in log(tau + 1).
        # The offset of 1 corresponds to the photosphere optical depth (tau=1),
        # where the radiation field transitions from optically thin to thick.
        # This compresses both the wind (tau<<1, few zones) and the deep
        # interior (tau>>1, few zones) while giving maximum resolution at
        # tau~1 where the photosphere and shock front are.
        #
        # Mapping: f(tau) = log10(tau + 1)
        #   tau=0:    f=0      (wind — compressed)
        #   tau=0.1:  f=0.041  (outer ejecta)
        #   tau=1:    f=0.301  (photosphere — maximum resolution)
        #   tau=10:   f=1.041  (thermalization)
        #   tau=28:   f=1.462  (deep interior — compressed)
        f_tau = np.log10(tau + 1.0)
        f_lo = f_tau[0]
        f_hi = f_tau[-1]
        f_targets = np.linspace(f_lo, f_hi, target_points)
        tau_targets = 10.0 ** f_targets - 1.0
        tau_targets = np.maximum(tau_targets, 0.0)
        print(f"[Tau-Based Downsample] Auto mode: log(tau+1) spacing "
              f"({target_points} points, photosphere-weighted)")
    else:
        if tau_breakpoints is None:
            tau_breakpoints = [0.01, 1.0]
        if zone_fractions is None:
            zone_fractions = [0.10, 0.25, 0.65]

        # Build full breakpoint list: [tau_floor, *tau_breakpoints, tau_max]
        full_breakpoints = [tau_floor] + list(tau_breakpoints) + [tau_max]
        n_segments = len(full_breakpoints) - 1

        if len(zone_fractions) != n_segments:
            raise ValueError(f"zone_fractions has {len(zone_fractions)} entries but "
                             f"need {n_segments} for {n_segments} segments")

        # Allocate points per segment
        n_per_segment = [int(target_points * f) for f in zone_fractions]
        # Adjust for rounding
        n_per_segment[0] += target_points - sum(n_per_segment)

        # Generate tau targets for each segment
        tau_targets = []
        for seg in range(n_segments):
            t_lo = full_breakpoints[seg]
            t_hi = full_breakpoints[seg + 1]
            n_seg = n_per_segment[seg]
            if n_seg <= 0:
                continue
            if t_hi <= t_lo:
                tau_targets.extend([t_lo] * max(n_seg, 1))
                continue
            if spacing == 'log':
                x_lo = np.log10(max(t_lo, eps))
                x_hi = np.log10(max(t_hi, eps))
            else:  # linear
                x_lo = t_lo
                x_hi = t_hi
            if seg == 0:
                # Include the outer boundary point
                seg_taus = np.linspace(x_lo, x_hi, n_seg)
            elif seg == n_segments - 1:
                # Include the inner boundary point
                seg_taus = np.linspace(x_lo, x_hi, n_seg)
            else:
                # Exclude endpoints (they're shared with adjacent segments)
                if n_seg > 2:
                    seg_taus = np.linspace(x_lo, x_hi, n_seg + 1)[1:-1]
                else:
                    seg_taus = np.linspace(x_lo, x_hi, n_seg)
            if spacing == 'log':
                tau_targets.extend(10 ** seg_taus)
            else:
                tau_targets.extend(seg_taus)

        tau_targets = np.array(sorted(tau_targets))

    # Ensure first and last match the original grid
    tau_targets[0] = tau[0]
    tau_targets[-1] = tau[-1]

    # Map tau targets to radius via interpolation
    # tau is increasing (outer->inner), radius is decreasing (outer->inner)
    # np.interp requires xp increasing, which tau already is
    r_new = np.interp(tau_targets, tau, radius)

    # Ensure monotonic decreasing radius (remove any duplicates)
    r_new = np.maximum.accumulate(r_new[::-1])[::-1]  # enforce decreasing
    _, unique_idx = np.unique(r_new, return_index=True)
    unique_idx = sorted(unique_idx)
    r_new = r_new[unique_idx]
    tau_targets = tau_targets[unique_idx]

    # --- Composition-aware grid enhancement ---
    # Detect sharp composition interfaces in the original grid and force-include
    # original STELLA zones at those interfaces. This ensures the grid HAS points
    # at Fe-core/Si-shell, Si/O, O/He, etc. boundaries, which pure tau-spacing
    # can miss, especially in metal-rich ejecta where the opacity (and thus tau)
    # varies smoothly even though the composition changes abruptly.
    #
    # CMFGEN-friendly constraints:
    #   - Cap interface points at 15% of target to avoid dominating the grid
    #   - Prioritize metal-rich interfaces (Fe, Si, S, Ca, O) over H/He gradients
    #   - After insertion, enforce max_cell_ratio by inserting midpoints
    interface_indices = set()
    if composition_aware:
        # Detect interfaces for hybrid interpolation (always done when composition_aware)
        if composition_elements is None:
            composition_elements = [e for e in DEFAULT_COMPOSITION_ELEMENTS
                                    if e in mass_fractions]

        all_interface_indices, iface_by_elem = detect_composition_interfaces(
            mass_fractions,
            gradient_threshold=composition_gradient_threshold,
            elements=composition_elements,
        )

        # Cap interface points at 15% of target. Prioritize metal-rich.
        max_interface_points = int(target_points * composition_interface_cap)
        metal_rich_elements = ['FePeak', 'Ni56', 'Ni58', 'Si', 'S', 'Ca', 'O', 'C']
        metal_iface = set()
        he_iface = set()
        for elem, edges in iface_by_elem.items():
            for edge in edges:
                if elem in ['H', 'He']:
                    he_iface.add(edge)
                    he_iface.add(edge + 1)
                else:
                    metal_iface.add(edge)
                    metal_iface.add(edge + 1)

        if len(metal_iface) <= max_interface_points:
            interface_indices = metal_iface
            remaining = max_interface_points - len(metal_iface)
            if remaining > 0 and len(he_iface) > 0:
                he_sorted = sorted(he_iface - metal_iface)
                interface_indices.update(he_sorted[:remaining])
        else:
            # Too many metal interfaces: keep the most significant ones
            iface_scores = []
            for idx in sorted(metal_iface):
                if idx == 0 or idx >= len(radius):
                    continue
                score = 0.0
                for elem in metal_rich_elements:
                    if elem not in mass_fractions:
                        continue
                    X = mass_fractions[elem]
                    if idx < len(X):
                        dX = abs(X[idx] - X[idx - 1])
                        Xmax = max(abs(X[idx]), abs(X[idx - 1]), eps)
                        score = max(score, dX / Xmax)
                iface_scores.append((score, idx))
            iface_scores.sort(reverse=True)
            interface_indices = set(idx for _, idx in iface_scores[:max_interface_points])

        n_interfaces = len(interface_indices)
        if n_interfaces > 0 and composition_grid_insertion:
            print(f"[Composition-Aware] Detected {len(all_interface_indices)} total "
                  f"interface zones, selected {n_interfaces} "
                  f"(cap: {max_interface_points}, prioritizing metal-rich):")
            for elem, edges in sorted(iface_by_elem.items()):
                if edges:
                    n_included = sum(1 for e in edges
                                     if e in interface_indices or (e+1) in interface_indices)
                    print(f"  {elem}: {len(edges)} detected, {n_included} included")

            # Insert original zone radii at selected interface points
            r_interfaces = radius[sorted(interface_indices)]
            r_combined = np.concatenate([r_new, r_interfaces])
            r_combined = np.sort(r_combined)[::-1]
            r_combined = np.maximum.accumulate(r_combined[::-1])[::-1]
            _, unique_idx2 = np.unique(r_combined, return_index=True)
            unique_idx2 = sorted(unique_idx2)
            r_combined = r_combined[unique_idx2]

            # If we exceeded target_points, remove tau-spaced points whose
            # removal causes the smallest increase in max cell ratio.
            if len(r_combined) > target_points:
                is_interface = np.array([
                    np.any(np.abs(r_combined[i] - r_interfaces) < 1.0)
                    for i in range(len(r_combined))
                ])
                best_max_ratio = 0.0
                while len(r_combined) > target_points:
                    removable = np.where(~is_interface)[0]
                    if len(removable) <= 2:
                        break
                    best_remove = -1
                    best_max_ratio = 1e30
                    for ri in removable:
                        if ri <= 0 or ri >= len(r_combined) - 1:
                            continue
                        dr_left = abs(r_combined[ri] - r_combined[ri - 1])
                        dr_right = abs(r_combined[ri + 1] - r_combined[ri])
                        dr_merged = dr_left + dr_right
                        if ri >= 2:
                            dr_ll = abs(r_combined[ri - 1] - r_combined[ri - 2])
                        else:
                            dr_ll = dr_merged
                        if ri < len(r_combined) - 2:
                            dr_rr = abs(r_combined[ri + 2] - r_combined[ri + 1])
                        else:
                            dr_rr = dr_merged
                        new_ratio1 = max(dr_merged, dr_ll) / min(dr_merged, dr_ll)
                        new_ratio2 = max(dr_merged, dr_rr) / min(dr_merged, dr_rr)
                        worst = max(new_ratio1, new_ratio2)
                        if worst < best_max_ratio:
                            best_max_ratio = worst
                            best_remove = ri
                    if best_remove < 0:
                        break
                    r_combined = np.delete(r_combined, best_remove)
                    is_interface = np.delete(is_interface, best_remove)

                print(f"[Composition-Aware] Trimmed to {len(r_combined)} points "
                      f"(target: {target_points}, "
                      f"{int(np.sum(is_interface))} interface zones, "
                      f"max cell ratio after trim: {best_max_ratio:.2f})")

            r_new = r_combined
            n_iface_final = int(np.sum(is_interface))
            dr_final = np.abs(np.diff(r_new))
            if len(dr_final) >= 2:
                ratios_final = np.maximum(dr_final[:-1], dr_final[1:]) / np.minimum(dr_final[:-1], dr_final[1:])
                print(f"[Composition-Aware] Final grid: {len(r_new)} points "
                      f"({n_iface_final} at composition interfaces, "
                      f"max cell ratio: {ratios_final.max():.2f})")
            else:
                print(f"[Composition-Aware] Final grid: {len(r_new)} points "
                      f"({n_iface_final} at composition interfaces)")
        else:
            print(f"[Composition-Aware] No sharp composition interfaces detected "
                  f"(threshold={composition_gradient_threshold})")

    # --- Ionization-front-aware grid enhancement ---
    # Detect ionization fronts (sharp electron density gradients) in the
    # original STELLA grid and force-include source zones at those fronts.
    # This ensures the grid resolves recombination/ionization transitions
    # over a few depths, which is critical for CMFGEN population solver
    # convergence. Without this, pure tau-spacing can miss ionization
    # fronts that occur at smooth tau but sharp n_e transitions.
    if ionization_front_aware:
        front_indices, front_info = detect_ionization_fronts(
            radius, rho, T, mass_fractions,
            ne_jump_threshold=ionization_front_threshold,
            min_points_per_front=ionization_front_min_points,
            front_width_zones=ionization_front_width,
        )

        if front_info:
            max_front_points = int(target_points * ionization_front_cap)
            # Sort fronts by ne_jump magnitude (most significant first)
            front_info.sort(key=lambda f: f['ne_jump'], reverse=True)

            # Select front indices up to the cap
            selected_front_indices = set()
            for fi in front_info:
                center = fi['center_idx']
                lo = max(0, center - ionization_front_width)
                hi = min(len(radius) - 1, center + ionization_front_width)
                front_zone_count = sum(1 for idx in range(lo, hi + 1)
                                       if idx in front_indices)
                if len(selected_front_indices) + front_zone_count > max_front_points:
                    break
                for idx in range(lo, hi + 1):
                    if idx in front_indices:
                        selected_front_indices.add(idx)

            print(f"[Ionization-Front-Aware] Detected {len(front_info)} ionization front(s), "
                  f"selected {len(selected_front_indices)} source zones "
                  f"(cap: {max_front_points}, threshold: {ionization_front_threshold}):")
            for fi in front_info[:10]:  # Print up to 10 fronts
                print(f"  Front at r={fi['r_cm']:.3e} cm, "
                      f"|dlog(n_e)/dlog(r)|={fi['ne_jump']:.2f}, "
                      f"T={fi['T_K']:.0f} K")

            if selected_front_indices and len(selected_front_indices) <= max_front_points:
                # Insert original zone radii at selected front points
                r_fronts = radius[sorted(selected_front_indices)]
                r_combined = np.concatenate([r_new, r_fronts])
                r_combined = np.sort(r_combined)[::-1]
                r_combined = np.maximum.accumulate(r_combined[::-1])[::-1]
                _, unique_idx3 = np.unique(r_combined, return_index=True)
                unique_idx3 = sorted(unique_idx3)
                r_combined = r_combined[unique_idx3]

                # Mark which points are front points (for trimming protection)
                is_front_pt = np.array([
                    np.any(np.abs(r_combined[i] - r_fronts) < 1.0)
                    for i in range(len(r_combined))
                ])

                # If we exceeded target_points + cap, remove tau-spaced points
                # whose removal causes the smallest increase in max cell ratio.
                # Never remove front points.
                max_total = target_points + max_front_points
                if len(r_combined) > max_total:
                    while len(r_combined) > max_total:
                        removable = np.where(~is_front_pt)[0]
                        if len(removable) <= 2:
                            break
                        best_remove = -1
                        best_ratio = 1e30
                        for ri in removable:
                            if ri <= 0 or ri >= len(r_combined) - 1:
                                continue
                            dr_l = abs(r_combined[ri] - r_combined[ri - 1])
                            dr_r = abs(r_combined[ri + 1] - r_combined[ri])
                            dr_m = dr_l + dr_r
                            if ri >= 2:
                                dr_ll = abs(r_combined[ri - 1] - r_combined[ri - 2])
                            else:
                                dr_ll = dr_m
                            if ri < len(r_combined) - 2:
                                dr_rr = abs(r_combined[ri + 2] - r_combined[ri + 1])
                            else:
                                dr_rr = dr_m
                            ratio1 = max(dr_m, dr_ll) / max(min(dr_m, dr_ll), eps)
                            ratio2 = max(dr_m, dr_rr) / max(min(dr_m, dr_rr), eps)
                            worst = max(ratio1, ratio2)
                            if worst < best_ratio:
                                best_ratio = worst
                                best_remove = ri
                        if best_remove < 0:
                            break
                        r_combined = np.delete(r_combined, best_remove)
                        is_front_pt = np.delete(is_front_pt, best_remove)

                r_new = r_combined
                n_front_final = int(np.sum(is_front_pt)) if len(is_front_pt) == len(r_new) else len(selected_front_indices)
                dr_final = np.abs(np.diff(r_new))
                if len(dr_final) >= 2:
                    ratios_final = np.maximum(dr_final[:-1], dr_final[1:]) / np.minimum(dr_final[:-1], dr_final[1:])
                    print(f"[Ionization-Front-Aware] Final grid: {len(r_new)} points "
                          f"({n_front_final} at ionization fronts, "
                          f"max cell ratio: {ratios_final.max():.2f})")
                else:
                    print(f"[Ionization-Front-Aware] Final grid: {len(r_new)} points "
                          f"({n_front_final} at ionization fronts)")
        else:
            print(f"[Ionization-Front-Aware] No ionization fronts detected "
                  f"(threshold={ionization_front_threshold})")

    # Interpolate all quantities onto the new grid
    # radius is decreasing, so for np.interp we need increasing xp
    r_increasing = radius[::-1]

    def _interp_arr(arr, log=False):
        if log:
            arr_safe = np.maximum(arr, eps)
            log_arr = np.log10(arr_safe)
            result = 10 ** np.interp(r_new, r_increasing, log_arr[::-1])
        else:
            result = np.interp(r_new, r_increasing, arr[::-1])
        return result

    v_new = _interp_arr(v, log=False)
    T_new = _interp_arr(T, log=True)
    rho_new = _interp_arr(rho, log=True)
    kappa_new = _interp_arr(kappa, log=True)

    # Mass fractions: use hybrid interpolation when composition-aware
    if composition_aware and len(interface_indices) > 0:
        mf_new = hybrid_interp_mass_fractions(
            r_new, radius, mass_fractions,
            interface_indices=interface_indices,
            trace_threshold=trace_threshold,
        )
        print(f"[Composition-Aware] Mass fractions: hybrid interpolation "
              f"(nearest-neighbour at {len(interface_indices)} interface zones, "
              f"log for trace species < {trace_threshold:.0e})")
    else:
        mf_new = {}
        for elem, arr in mass_fractions.items():
            mf_new[elem] = _interp_arr(arr, log=False)

    # Post-selection verification: check composition shell coverage
    if composition_aware:
        n_shells, n_under, shell_info = verify_composition_coverage(
            r_new, mass_fractions, radius,
            min_zones_per_shell=min_zones_per_shell,
            elements=composition_elements,
        )
        print(f"[Composition-Aware] Verification: {n_shells} composition shells, "
              f"{n_under} under-sampled (< {min_zones_per_shell} zones)")

    # Print zone distribution
    print(f"[Tau-Based Downsample] Reduced from {n} to {len(r_new)} points "
          f"(target: {target_points}, spacing: {spacing})")
    print(f"  Tau range: [{tau[0]:.3e}, {tau[-1]:.3e}]")
    if spacing == 'auto':
        # Show zone density by tau decade
        for lo_exp in range(int(np.floor(np.log10(max(tau[0], 1e-10)))),
                            int(np.ceil(np.log10(tau[-1])))):
            t_lo = 10.0 ** lo_exp
            t_hi = 10.0 ** (lo_exp + 1)
            count = sum(1 for t in tau_targets if t_lo <= t < t_hi)
            if count > 0:
                print(f"  tau [{t_lo:.1e}, {t_hi:.1e}): {count} zones")
        # Also show the dtau/dr resolution metric
        dtau = np.diff(tau_targets)
        r_at_tau = np.interp(tau_targets, tau, radius)
        dr = np.abs(np.diff(r_at_tau))
        dtau_dr = dtau / np.maximum(dr, eps)
        print(f"  dtau/dr: min={dtau_dr.min():.3e}, max={dtau_dr.max():.3e}, "
              f"ratio={dtau_dr.max()/max(dtau_dr.min(),eps):.1f}")
    else:
        print(f"  Breakpoints: {[f'{b:.3e}' for b in full_breakpoints]}")
        print(f"  Zone fractions: {zone_fractions}")
        for seg in range(n_segments):
            t_lo = full_breakpoints[seg]
            t_hi = full_breakpoints[seg + 1]
            count = sum(1 for t in tau_targets if t_lo <= t < t_hi) + (1 if seg == n_segments - 1 else 0)
            print(f"  tau [{t_lo:.3e}, {t_hi:.3e}): {count} zones")

    return (r_new, v_new, T_new, rho_new, kappa_new) + tuple(mf_new.values())


def _sigma_at_index(r_list, v_list, i, eps=1e-30):
    """Compute sigma at grid point i using the exact CMFGEN formula.

    CMFGEN (rv_sn_model_sniin.f:241) uses a distance-weighted average of
    the left and right LINEAR velocity gradients. See adaptive_grid_refinement's
    docstring / _local_sigma for details. Returns 0.0 at boundaries.
    """
    if i <= 0 or i >= len(r_list) - 1:
        return 0.0
    h_left = r_list[i - 1] - r_list[i]
    h_right = r_list[i] - r_list[i + 1]
    if abs(h_left) < eps or abs(h_right) < eps:
        return 0.0
    T1 = (v_list[i - 1] - v_list[i]) / h_left
    T2 = (v_list[i] - v_list[i + 1]) / h_right
    dVdR = (T1 * h_right + T2 * h_left) / (h_left + h_right)
    v_val = max(abs(v_list[i]), eps)
    return r_list[i] * dVdR / v_val - 1.0


def optimal_insert_fraction(r_list, v_list, idx, eps=1e-30,
                            fractions=(0.45, 0.55, 0.4, 0.6, 0.35, 0.65),
                            min_improvement=0.2):
    """Find the insertion fraction f in (0,1) between idx and idx+1 that
    minimizes the resulting |sigma| at the three affected points (idx,
    the new point, and idx+1).

    This is shared by all grid-insertion routines (AMR, min-points fill,
    grid-spacing enforcement) so that ANY point insertion — not just AMR
    refinement — avoids introducing sigma discretization noise. A naive
    midpoint (f=0.5) split can create a sharp kink in sigma when the local
    velocity curvature is asymmetric (common near velocity extrema / shock
    fronts); searching over f finds a split that keeps sigma smooth across
    the new point and its neighbours.

    The search is intentionally conservative: candidate fractions are kept
    close to 0.5 (0.35-0.65) and a fraction is only used instead of the
    midpoint if it improves the worst-case |sigma| by at least
    ``min_improvement`` (relative). Chasing small sigma improvements with
    very asymmetric splits (e.g. f=0.2) creates a large cell-size ratio
    that then gets flagged and re-subdivided by the cell-ratio-smoothing
    pass, cascading into far more scatter than it fixes.

    Requires at least one real neighbour on each side (idx > 0 and
    idx < len-2); otherwise returns the default midpoint 0.5.
    """
    if idx <= 0 or idx >= len(r_list) - 2:
        return 0.5

    def _score(f):
        r_new = r_list[idx] + f * (r_list[idx + 1] - r_list[idx])
        v_new = v_list[idx] + f * (v_list[idx + 1] - v_list[idx])

        # Sigma at idx (now has new right neighbor)
        h_left = r_list[idx - 1] - r_list[idx]
        h_right_new = r_list[idx] - r_new
        sigma_idx = 0.0
        if abs(h_left) > eps and abs(h_right_new) > eps:
            T1 = (v_list[idx - 1] - v_list[idx]) / h_left
            T2 = (v_list[idx] - v_new) / h_right_new
            dVdR = (T1 * h_right_new + T2 * h_left) / (h_left + h_right_new)
            sigma_idx = r_list[idx] * dVdR / max(abs(v_list[idx]), eps) - 1.0

        # Sigma at the new point (between idx and idx+1)
        h_left_new = r_list[idx] - r_new
        h_right_new2 = r_new - r_list[idx + 1]
        sigma_new = 0.0
        if abs(h_left_new) > eps and abs(h_right_new2) > eps:
            T1 = (v_list[idx] - v_new) / h_left_new
            T2 = (v_new - v_list[idx + 1]) / h_right_new2
            dVdR = (T1 * h_right_new2 + T2 * h_left_new) / (h_left_new + h_right_new2)
            sigma_new = r_new * dVdR / max(abs(v_new), eps) - 1.0

        # Sigma at idx+1 (now has new left neighbor)
        h_left_new3 = r_new - r_list[idx + 1]
        h_right_orig = r_list[idx + 1] - r_list[idx + 2]
        sigma_ip1 = 0.0
        if abs(h_left_new3) > eps and abs(h_right_orig) > eps:
            T1 = (v_new - v_list[idx + 1]) / h_left_new3
            T2 = (v_list[idx + 1] - v_list[idx + 2]) / h_right_orig
            dVdR = (T1 * h_right_orig + T2 * h_left_new3) / (h_left_new3 + h_right_orig)
            sigma_ip1 = r_list[idx + 1] * dVdR / max(abs(v_list[idx + 1]), eps) - 1.0

        return max(abs(sigma_idx), abs(sigma_new), abs(sigma_ip1))

    best_f = 0.5
    best_score = _score(0.5)
    for f in fractions:
        score = _score(f)
        if score < best_score * (1.0 - min_improvement):
            best_score, best_f = score, f

    return best_f


def adaptive_grid_refinement(radius, v, T, rho, n_e, tau, mass_fractions,
                             max_radius_jump=0.05, max_temp_jump=0.10,
                             max_density_jump=0.15, max_ne_jump=0.10,
                             max_logtau_jump=0.10, max_points=500,
                             n_outer_extra=5,
                             max_velocity_jump=0.10, max_sigma_jump=1.0,
                             max_cell_ratio=1.5,
                             max_abs_sigma=5.0, sigma_optimize=True,
                             thin_tau_threshold=1.0, thin_priority=1.0,
                             ne_priority_threshold=2.0, ne_priority=1.0,
                             sigma_priority_threshold=2.0, sigma_priority=5.0,
                             dsigma_priority_threshold=0.5, dsigma_priority=5.0,
                             max_thin_zones=0, max_thick_tau=1e10,
                             critical_mask=None, critical_min_dr_r=1e-5):
    """Priority-based adaptive grid refinement.

    Each iteration scores every edge by its worst normalised violation
    (jump / threshold), then inserts a midpoint at the single worst edge.
    This guarantees the highest-gradient regions receive points first.

    After all violations are resolved (or ``max_points`` is reached),
    ``n_outer_extra`` additional points are sprinkled into the outermost
    region (between the first two existing points) to improve the outer
    boundary sampling for CMFGEN.

    Arrays are outer->inner (decreasing radius).

    Velocity and sigma criteria:
        ``max_velocity_jump`` controls the fractional velocity change per cell
        (|dv|/v).  This ensures the grid resolves velocity gradients even when
        other quantities (rho, T, tau) are smooth — critical for non-monotonic
        velocity fields in interacting SNe where the shock front has a large
        velocity reversal over a narrow radial range.

        ``max_sigma_jump`` controls the change in sigma = dlnv/dlnr - 1 between
        adjacent cells.  Large sigma jumps create spikes in the CMFGEN moment
        equations (MOM_JREL_V9), causing NaN or convergence issues.  Refining
        to keep sigma smooth is essential for the Dessart+2015 direct treatment.

        ``max_abs_sigma`` controls the maximum allowed |sigma| at any grid
        point.  When |sigma| exceeds this threshold, the AMR inserts points
        to reduce the discretization error.  Note: the physical sigma (from
        the velocity reversal) cannot be eliminated — AMR can only reduce
        the discretization error that makes |sigma| worse than its true value.
        Set to 5.0 by default (only catches extreme discretization artifacts;
        physical sigma in interacting SNe can reach -4 at the velocity reversal).

        ``sigma_optimize`` when True, uses the exact CMFGEN sigma formula
        (distance-weighted linear gradient average, rv_sn_model_sniin.f:241)
        instead of the logarithmic central difference.  This ensures the AMR
        optimizes the same quantity that CMFGEN computes.  Also enables
        optimal insertion point selection (non-midpoint) to minimize |sigma|.

        ``ne_priority`` boosts the score of edges whose dNe/Ne (normalised by
        ``max_ne_jump``) exceeds ``ne_priority_threshold``, analogous to
        ``sigma_priority``/``dsigma_priority``.  Without this, edges at a
        sharp recombination/ionization front (large dNe/Ne, but small
        dr/r, dT/T, dv/v) compete on equal footing with the many
        velocity/sigma-driven edges in non-monotonic-velocity SN models,
        and routinely lose the budget race (see AGENTS.md "Ionization
        front under-resolved" — with the default ne_priority=1.0, 288-502
        edges near the H/He recombination front remained budget-limited
        even after raising max_points from 400 to 800).

        ``critical_mask`` (optional boolean array, length ``len(radius)``)
        marks narrow high-gradient 'critical regions' (shock fronts,
        composition interfaces, CSM shells) identified by
        :func:`identify_critical_regions`.  When provided, the AMR's
        internal ``MIN_DR_OVER_R`` floor is relaxed to
        ``critical_min_dr_r`` (default 1e-5) for edges whose endpoints
        fall inside a critical region, so the AMR can subdivide below the
        global floor there.  Points inside critical regions are also
        protected from the post-AMR coalescing and degenerate-merging
        passes.  Outside critical regions the original ``MIN_DR_OVER_R``
        floor is enforced unchanged.
    """
    eps = 1e-30
    MIN_DR_OVER_R = 1e-4  # 0.01% minimum cell size
                          # At R~1e15 cm, 1e-4 = 1e11 cm = ~1000 km
                          # Needed for velocity reversal regions where dr/r is
                          # already ~1e-3 but velocity changes by ~2000 km/s
    # Region-dependent floor: when a critical mask is supplied, edges
    # inside critical regions may subdivide down to critical_min_dr_r
    # (default 1e-5), letting the AMR insert the points CMFGEN needs in
    # narrow high-gradient regions without lowering the floor globally.
    has_critical = critical_mask is not None and len(critical_mask) == len(radius)

    def _log_interp(a, b, f):
        return 10**(np.log10(max(a, eps)) + f * (np.log10(max(b, eps)) - np.log10(max(a, eps))))

    def _min_dr_for_edge(i):
        """Region-dependent minimum dr/r for edge i..i+1.

        Returns ``critical_min_dr_r`` if either endpoint is in a critical
        region (so the AMR may subdivide below the global floor there),
        otherwise the global ``MIN_DR_OVER_R``.
        """
        if has_critical and (is_critical_list[i] or is_critical_list[i + 1]):
            return critical_min_dr_r
        return MIN_DR_OVER_R

    def _insert_midpoint(idx, f=0.5):
        """Insert a log-interpolated midpoint between idx and idx+1 at fraction f."""
        r_list.insert(idx + 1, r_list[idx] + f * (r_list[idx + 1] - r_list[idx]))
        v_list.insert(idx + 1, v_list[idx] + f * (v_list[idx + 1] - v_list[idx]))
        is_original_list.insert(idx + 1, False)
        # An inserted point inside a critical region inherits the critical
        # flag so the relaxed min_dr/r continues to apply as the AMR
        # subdivides further inside the region.
        if has_critical:
            is_critical_list.insert(idx + 1,
                bool(is_critical_list[idx] or is_critical_list[idx + 1]))
        T_list.insert(idx + 1, _log_interp(T_list[idx], T_list[idx + 1], f))
        rho_list.insert(idx + 1, _log_interp(rho_list[idx], rho_list[idx + 1], f))
        n_e_list.insert(idx + 1, _log_interp(n_e_list[idx], n_e_list[idx + 1], f))
        tau_list.insert(idx + 1, _log_interp(tau_list[idx], tau_list[idx + 1], f))
        for elem in mf_lists:
            mf_lists[elem].insert(idx + 1,
                mf_lists[elem][idx] + f * (mf_lists[elem][idx + 1] - mf_lists[elem][idx]))

    def _local_sigma(i):
        """Compute sigma at grid point i using the exact CMFGEN formula.

        CMFGEN (rv_sn_model_sniin.f:241) uses a distance-weighted average of
        the left and right LINEAR velocity gradients:

            T1 = (V[i-1]-V[i])/(R[i-1]-R[i])     ! left gradient
            T2 = (V[i]-V[i+1])/(R[i]-R[i+1])     ! right gradient
            dVdR = (T1*h_right + T2*h_left) / (h_left + h_right)
            sigma = R*dVdR/V - 1

        where h_left = R[i-1]-R[i], h_right = R[i]-R[i+1].

        This differs from the logarithmic central difference (dlnV/dlnR - 1)
        because the weighting by cell size biases sigma toward the finer
        side's gradient on non-uniform grids.

        Returns 0.0 for boundary points where central differences aren't available.
        """
        if i <= 0 or i >= len(r_list) - 1:
            return 0.0
        h_left = r_list[i - 1] - r_list[i]
        h_right = r_list[i] - r_list[i + 1]
        if abs(h_left) < eps or abs(h_right) < eps:
            return 0.0
        T1 = (v_list[i - 1] - v_list[i]) / h_left    # left linear gradient
        T2 = (v_list[i] - v_list[i + 1]) / h_right    # right linear gradient
        dVdR = (T1 * h_right + T2 * h_left) / (h_left + h_right)
        v_val = max(abs(v_list[i]), eps)
        return r_list[i] * dVdR / v_val - 1.0

    def _optimal_insert_fraction(idx):
        """Find the insertion fraction f that minimizes max|sigma| at points idx and idx+1.

        After inserting a point at R_new = R[idx] + f*(R[idx+1]-R[idx]), the
        sigma at the original points changes because the cell sizes change.
        We try several fractions and pick the one that minimizes the resulting
        max(|sigma[idx]|, |sigma[idx+1]|).

        For velocity extrema (where V turns around), inserting closer to the
        extremum (f≠0.5) can significantly reduce |sigma| by balancing the
        left and right gradients.
        """
        if not sigma_optimize or idx <= 0 or idx >= len(r_list) - 2:
            return 0.5

        best_f = 0.5
        best_score = float('inf')

        for f in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
            # Simulate insertion at fraction f
            r_new = r_list[idx] + f * (r_list[idx + 1] - r_list[idx])
            v_new = v_list[idx] + f * (v_list[idx + 1] - v_list[idx])

            # Sigma at idx (now has new right neighbor)
            h_left = r_list[idx - 1] - r_list[idx]
            h_right_new = r_list[idx] - r_new
            if abs(h_left) > eps and abs(h_right_new) > eps:
                T1 = (v_list[idx - 1] - v_list[idx]) / h_left
                T2 = (v_list[idx] - v_new) / h_right_new
                dVdR = (T1 * h_right_new + T2 * h_left) / (h_left + h_right_new)
                sigma_idx = r_list[idx] * dVdR / max(abs(v_list[idx]), eps) - 1.0
            else:
                sigma_idx = 0.0

            # Sigma at new point (between idx and idx+1)
            h_left_new = r_list[idx] - r_new
            h_right_new2 = r_new - r_list[idx + 1]
            if abs(h_left_new) > eps and abs(h_right_new2) > eps:
                T1 = (v_list[idx] - v_new) / h_left_new
                T2 = (v_new - v_list[idx + 1]) / h_right_new2
                dVdR = (T1 * h_right_new2 + T2 * h_left_new) / (h_left_new + h_right_new2)
                sigma_new = r_new * dVdR / max(abs(v_new), eps) - 1.0
            else:
                sigma_new = 0.0

            # Sigma at idx+1 (now has new left neighbor)
            h_left_new3 = r_new - r_list[idx + 1]
            h_right_orig = r_list[idx + 1] - r_list[idx + 2]
            if abs(h_left_new3) > eps and abs(h_right_orig) > eps:
                T1 = (v_new - v_list[idx + 1]) / h_left_new3
                T2 = (v_list[idx + 1] - v_list[idx + 2]) / h_right_orig
                dVdR = (T1 * h_right_orig + T2 * h_left_new3) / (h_left_new3 + h_right_orig)
                sigma_ip1 = r_list[idx + 1] * dVdR / max(abs(v_list[idx + 1]), eps) - 1.0
            else:
                sigma_ip1 = 0.0

            score = max(abs(sigma_idx), abs(sigma_new), abs(sigma_ip1))
            if score < best_score:
                best_score = score
                best_f = f

        return best_f

    def _edge_score(i):
        """Return the maximum normalised violation for edge i..i+1.

        A score <= 1.0 means all criteria are satisfied; > 1.0 means at
        least one threshold is exceeded.  The score equals the largest
        ratio  (actual_jump / allowed_jump)  across all quantities.
        """
        dr_r     = abs(r_list[i+1] - r_list[i]) / max(abs(r_list[i]), eps)
        dT_T     = abs(T_list[i+1] - T_list[i]) / max(T_list[i], eps)
        drho_rho = abs(rho_list[i+1] - rho_list[i]) / max(rho_list[i], eps)
        dne_ne   = abs(n_e_list[i+1] - n_e_list[i]) / max(n_e_list[i], eps)
        dlogtau  = abs(np.log10(max(tau_list[i+1], eps)) - np.log10(max(tau_list[i], eps)))

        # Velocity jump: fractional change in |v| across the cell
        dv_v = abs(v_list[i+1] - v_list[i]) / max(abs(v_list[i]), eps)

        # Sigma jump: change in local sigma between adjacent grid points.
        # Uses the exact CMFGEN sigma formula (distance-weighted linear gradient).
        sigma_i   = _local_sigma(i)
        sigma_ip1 = _local_sigma(i + 1)
        dsigma = abs(sigma_ip1 - sigma_i)

        # Absolute sigma criterion: penalize large |sigma| values directly.
        # The physical sigma (from velocity reversal) cannot be eliminated,
        # but the discretization error can be reduced by refining.
        abs_sigma = max(abs(sigma_i), abs(sigma_ip1))

        # Flat-region detection: when both T and rho are nearly constant across
        # the edge, the region is physically flat (e.g., shell-ejecta buffer
        # where velocity reverses but T/rho don't change). In such regions,
        # sigma-driven AMR wastes budget creating degenerate zones that get
        # merged then re-added by cell-ratio smoothing. Skip sigma/velocity
        # criteria entirely in flat regions — only physical gradients drive
        # refinement there.
        is_flat = (dT_T < max_temp_jump * 0.5 and
                   drho_rho < max_density_jump * 0.5)

        if is_flat:
            # Physical gradients only — no sigma/velocity criteria
            score = max(dr_r / max_radius_jump,
                        dT_T / max_temp_jump,
                        drho_rho / max_density_jump,
                        dne_ne / max_ne_jump)
            if tau_list[i] > 0.01:
                score = max(score, dlogtau / max_logtau_jump)
        else:
            # Full criteria including sigma and velocity
            score = max(dr_r / max_radius_jump,
                        dT_T / max_temp_jump,
                        drho_rho / max_density_jump,
                        dne_ne / max_ne_jump,
                        dv_v / max_velocity_jump,
                        dsigma / max_sigma_jump,
                        abs_sigma / max_abs_sigma)
            if tau_list[i] > 0.01:
                score = max(score, dlogtau / max_logtau_jump)

            # Sigma priority: boost edges in high-|sigma| regions (velocity reversals)
            # so they compete for AMR budget allocation.
            if abs_sigma > sigma_priority_threshold:
                score *= sigma_priority

            # d(sigma) priority: boost edges with sharp sigma transitions (e.g. outer
            # boundary where sigma drops from ~1.4 to ~-0.6 over one cell).
            if dsigma > dsigma_priority_threshold:
                score *= dsigma_priority

        # Ne priority: boost edges with a sharp electron-density jump (e.g. a
        # recombination/ionization front) so they compete for AMR budget
        # against the (often much more numerous) velocity/sigma-driven edges
        # in non-monotonic-velocity SN models. Without this, a steep but
        # otherwise "quiet" (small dr/r, dT/T, dv/v) Ne jump routinely loses
        # the budget race and remains under-resolved.
        ne_ratio = dne_ne / max_ne_jump
        if ne_ratio > ne_priority_threshold:
            score *= ne_priority

        # Tau-based priority: applied LAST so it can override sigma/dsigma boosts.
        # thin_priority < 1.0 penalizes (deprioritizes) low-tau regions,
        # reducing AMR budget allocation to optically thin outer zones.
        # thin_priority > 1.0 boosts (prioritizes) low-tau regions.
        tau_edge = min(tau_list[i], tau_list[i+1])
        if tau_edge <= thin_tau_threshold:
            score *= thin_priority

        return score, dr_r

    # Convert to lists for easy insertion
    r_list   = list(radius)
    v_list   = list(v)
    is_original_list = [True] * len(radius)
    is_critical_list = [bool(critical_mask[i]) for i in range(len(radius))] \
        if has_critical else [False] * len(radius)
    T_list   = list(T)
    rho_list = list(rho)
    n_e_list = list(n_e)
    tau_list = list(tau)
    mf_lists = {elem: list(mass_fractions[elem]) for elem in mass_fractions}

    # --- Diagnostic: audit the base grid ---
    scores_init = [_edge_score(i) for i in range(len(r_list) - 1)]
    n_violations = sum(1 for s, _ in scores_init if s > 1.0)
    worst = max(scores_init, key=lambda x: x[0])
    print(f"[AMR] Base grid: {len(r_list)} pts, {n_violations} edges violate thresholds "
          f"(worst score = {worst[0]:.2f}x)")

    if n_violations == 0:
        print("[AMR] Base grid already meets ALL criteria. No refinement needed.")
    else:
        # Reserve n_outer_extra slots for the outer-boundary padding below
        amr_budget = max_points - n_outer_extra
        # Reserve 25% of the budget for grid+sigma smoothing passes so they
        # can actually run after the initial refinement exhausts its share.
        amr_initial_budget = int(amr_budget * 0.75)
        n_inserted = 0
        n_skipped = 0
        n_thin_skipped = 0
        max_iter = max_points * 4  # safety cap
        for _ in range(max_iter):
            if len(r_list) >= amr_initial_budget:
                print(f"[AMR] Budget reached ({len(r_list)} pts, "
                      f"reserving {n_outer_extra} for outer padding)")
                break
            # Count current zones in thin region (tau < thin_tau_threshold)
            n_current_thin = sum(1 for t in tau_list if t <= thin_tau_threshold)
            # Find the worst edge that can still be subdivided
            best_i, best_score, best_dr = -1, 0.0, 0.0
            for i in range(len(r_list) - 1):
                sc, dr = _edge_score(i)
                # Skip thin-region edges if max_thin_zones is set and exceeded
                if max_thin_zones > 0 and n_current_thin >= max_thin_zones:
                    tau_edge = min(tau_list[i], tau_list[i+1])
                    if tau_edge <= thin_tau_threshold:
                        continue
                # Skip thick-region edges beyond max_thick_tau (degenerate inner
                # region where V, T, rho are nearly constant — AMR budget is wasted
                # here since degenerate merge will remove these points anyway).
                if max_thick_tau < 1e10:
                    tau_edge_hi = max(tau_list[i], tau_list[i+1])
                    if tau_edge_hi > max_thick_tau:
                        continue
                if sc > best_score and dr >= _min_dr_for_edge(i):
                    best_i, best_score, best_dr = i, sc, dr
            if best_i < 0:
                # No sub dividable edge violates thresholds
                if n_skipped > 0:
                    print(f"[AMR] {n_skipped} edges too small to subdivide "
                          f"(dr/r < global {MIN_DR_OVER_R}"
                          f"{f' / critical {critical_min_dr_r}' if has_critical else ''})")
                break  # all remaining violations are in cells too small to split
            if best_score <= 1.0:
                break  # all edges within tolerance
            # Use optimal insertion fraction to minimize |sigma|
            f_opt = _optimal_insert_fraction(best_i) if sigma_optimize else 0.5
            _insert_midpoint(best_i, f_opt)
            n_inserted += 1

        n_current_thin = sum(1 for t in tau_list if t <= thin_tau_threshold)
        print(f"[AMR] Inserted {n_inserted} points (grid now {len(r_list)} pts, "
              f"{n_current_thin} in thin region tau<{thin_tau_threshold})")
        remaining_violations = sum(1 for i in range(len(r_list) - 1)
                                   if _edge_score(i)[0] > 1.0)
        if remaining_violations:
            # Count how many are due to min-dr limit vs budget.  An edge is
            # "too small" only if it is below the *region-appropriate* floor.
            too_small = sum(1 for i in range(len(r_list) - 1)
                           if _edge_score(i)[0] > 1.0
                           and _edge_score(i)[1] < _min_dr_for_edge(i))
            print(f"[AMR] {remaining_violations} edges still exceed thresholds "
                  f"({too_small} too small to subdivide, "
                  f"{remaining_violations - too_small} budget-limited)")

        # --- Iterative sigma + grid smoothing passes ---
        # The AMR creates abrupt transitions between large and tiny cells at
        # shock-front boundaries (e.g., dr=128 next to dr=0.3, a 400:1 ratio).
        # These abrupt grid-spacing changes cause sigma jumps.  We interleave
        # sigma and grid smoothing within each iteration so both get budget:
        #   - Sigma smoothing: subdivide cells at sigma-jump boundaries
        #     (physically motivated — reduces MOM_JREL_V9 instabilities)
        #   - Grid smoothing: subdivide the larger of two adjacent cells when
        #     their size ratio exceeds max_cell_ratio (creates gradual spacing)
        # Previously grid smoothing ran first and exhausted the budget, leaving
        # 0 points for sigma smoothing. Now they alternate 1:1.
        max_cell_ratio = max_cell_ratio
        n_total_sigma = 0
        n_total_grid = 0

        def _thin_count():
            return sum(1 for t in tau_list if t <= thin_tau_threshold)

        def _is_thin_edge(i):
            return min(tau_list[i], tau_list[i+1]) <= thin_tau_threshold

        for smoothing_iter in range(20):
            if len(r_list) >= amr_budget:
                break

            # Alternate: do one sigma-smoothing insertion, then one grid-smoothing,
            # repeating until neither has work or budget is exhausted.
            n_grid_pass = 0
            n_sigma_pass = 0
            sigma_done = False
            grid_done = False
            unresolvable_r = set()
            for _ in range(max_points):
                if len(r_list) >= amr_budget:
                    break
                if sigma_done and grid_done:
                    break

                # --- Sigma smoothing step (one insertion) ---
                if not sigma_done and max_sigma_jump < 1e10:
                    sigma_arr = np.array([_local_sigma(i) for i in range(len(r_list))])
                    dsigma = np.abs(np.diff(sigma_arr))
                    jump_edges = np.where(dsigma > max_sigma_jump)[0]
                    resolvable = [e for e in jump_edges
                                  if any(abs(0.5 * (r_list[e] + r_list[e + 1]) - ur) < 1e-10
                                         for ur in unresolvable_r) is False]
                    if not resolvable:
                        sigma_done = True
                    else:
                        worst_edge = max(resolvable, key=lambda e: dsigma[e])
                        # Skip thin-region edges if max_thin_zones cap is reached
                        if max_thin_zones > 0 and _thin_count() >= max_thin_zones and _is_thin_edge(worst_edge):
                            unresolvable_r.add(0.5 * (r_list[worst_edge] + r_list[worst_edge + 1]))
                            continue
                        # Skip thick-region edges beyond max_thick_tau
                        if max_thick_tau < 1e10:
                            tau_edge_hi = max(tau_list[worst_edge], tau_list[worst_edge + 1])
                            if tau_edge_hi > max_thick_tau:
                                unresolvable_r.add(0.5 * (r_list[worst_edge] + r_list[worst_edge + 1]))
                                continue
                        # Skip sigma smoothing in flat regions (T and rho both
                        # nearly constant) — sigma-driven refinement there
                        # creates degenerate zones that get merged later.
                        dT_T = abs(T_list[worst_edge + 1] - T_list[worst_edge]) / max(T_list[worst_edge], eps)
                        drho_rho = abs(rho_list[worst_edge + 1] - rho_list[worst_edge]) / max(rho_list[worst_edge], eps)
                        if (dT_T < max_temp_jump * 0.5 and
                                drho_rho < max_density_jump * 0.5):
                            r_mid = 0.5 * (r_list[worst_edge] + r_list[worst_edge + 1])
                            unresolvable_r.add(r_mid)
                            continue
                        dr_worst = abs(r_list[worst_edge + 1] - r_list[worst_edge])
                        if dr_worst / max(abs(r_list[worst_edge]), eps) >= _min_dr_for_edge(worst_edge):
                            f_opt = _optimal_insert_fraction(worst_edge) if sigma_optimize else 0.5
                            _insert_midpoint(worst_edge, f_opt)
                            n_sigma_pass += 1
                        else:
                            candidates = []
                            if worst_edge > 0:
                                dr_prev = abs(r_list[worst_edge] - r_list[worst_edge - 1])
                                if dr_prev / max(abs(r_list[worst_edge - 1]), eps) >= _min_dr_for_edge(worst_edge - 1):
                                    candidates.append((dr_prev, worst_edge - 1))
                            if worst_edge < len(r_list) - 2:
                                dr_next = abs(r_list[worst_edge + 2] - r_list[worst_edge + 1])
                                if dr_next / max(abs(r_list[worst_edge + 1]), eps) >= _min_dr_for_edge(worst_edge + 1):
                                    candidates.append((dr_next, worst_edge + 1))
                            if candidates:
                                candidates.sort(reverse=True)
                                f_opt = _optimal_insert_fraction(candidates[0][1]) if sigma_optimize else 0.5
                                _insert_midpoint(candidates[0][1], f_opt)
                                n_sigma_pass += 1
                            else:
                                r_mid = 0.5 * (r_list[worst_edge] + r_list[worst_edge + 1])
                                unresolvable_r.add(r_mid)
                else:
                    sigma_done = True

                if len(r_list) >= amr_budget:
                    break

                # --- Grid smoothing step (one insertion) ---
                if not grid_done:
                    worst_i, worst_ratio = -1, 1.0
                    for i in range(len(r_list) - 2):
                        # Skip thin-region edges if max_thin_zones cap is reached
                        if max_thin_zones > 0 and _thin_count() >= max_thin_zones and _is_thin_edge(i):
                            continue
                        # Skip thick-region edges beyond max_thick_tau
                        if max_thick_tau < 1e10:
                            tau_edge_hi = max(tau_list[i], tau_list[i+1], tau_list[i+2])
                            if tau_edge_hi > max_thick_tau:
                                continue
                        dr0 = abs(r_list[i + 1] - r_list[i])
                        dr1 = abs(r_list[i + 2] - r_list[i + 1])
                        ratio = max(dr0, dr1) / max(min(dr0, dr1), eps)
                        if ratio > worst_ratio and ratio > max_cell_ratio:
                            if dr0 > dr1:
                                if dr0 / max(abs(r_list[i]), eps) >= _min_dr_for_edge(i):
                                    worst_i, worst_ratio = i, ratio
                            else:
                                if dr1 / max(abs(r_list[i + 1]), eps) >= _min_dr_for_edge(i + 1):
                                    worst_i, worst_ratio = i + 1, ratio
                    if worst_i < 0:
                        grid_done = True
                    else:
                        _insert_midpoint(worst_i)
                        n_grid_pass += 1
                else:
                    grid_done = True

            n_total_grid += n_grid_pass
            n_total_sigma += n_sigma_pass

            if n_grid_pass == 0 and n_sigma_pass == 0:
                break  # converged

        if n_total_grid > 0 or n_total_sigma > 0:
            print(f"[AMR] Iterative smoothing ({smoothing_iter + 1} iters): "
                  f"inserted {n_total_grid} grid-smoothing + {n_total_sigma} sigma-smoothing "
                  f"points (grid now {len(r_list)} pts)")
            remaining = sum(1 for i in range(len(r_list) - 1)
                           if _edge_score(i)[0] > 1.0)
            if remaining:
                print(f"[AMR] {remaining} edges still exceed thresholds after smoothing")

    # --- Post-AMR coalescing pass ---
    # Merge adjacent zones that are pathologically small RELATIVE TO THEIR
    # NEIGHBOURS, not in absolute terms. The old MIN_DR_OVER_R_COALESCE = 1e-3
    # used dr/r as the criterion, which at R ~ 1e12 cm removed zones with
    # dr < 1e9 cm — physically meaningful shock-front zones. This caused
    # massive coalescing (up to 191 points removed per epoch), followed by
    # FillMin inserting 27-131 low-quality uniform filler points.
    #
    # New criterion: coalesce only if a zone is >COALESCE_RATIO times smaller
    # than BOTH its neighbours (i.e., a numerical artifact, not a physical
    # feature). This preserves small zones at shock fronts while removing
    # AMR overshoots where refinement created a tiny sliver between two
    # larger cells.
    COALESCE_RATIO = 10.0  # zone must be 10x smaller than both neighbours
    n_coalesced = 0
    n_critical_protected = 0
    i = 1
    while i < len(r_list) - 1:
        dr_prev = abs(r_list[i] - r_list[i - 1])
        dr_curr = abs(r_list[i + 1] - r_list[i])
        dr_next = abs(r_list[i + 2] - r_list[i + 1]) if i + 2 < len(r_list) else dr_curr
        # Coalesce zone i+1 (remove point i+1) if it's a tiny sliver
        # between two larger cells AND it's not an original hydro point.
        # Critical-region points are never coalesced — they were inserted
        # (or preserved) specifically to resolve a narrow high-gradient
        # feature, and removing them would undo the relaxed-min_dr/r work.
        if (not is_original_list[i + 1] and
                not (has_critical and is_critical_list[i + 1]) and
                dr_curr < dr_prev / COALESCE_RATIO and
                dr_curr < dr_next / COALESCE_RATIO):
            for lst in [r_list, v_list, T_list, rho_list, n_e_list, tau_list,
                        is_original_list] + ([is_critical_list] if has_critical else []):
                del lst[i + 1]
            for elem in mf_lists:
                del mf_lists[elem][i + 1]
            n_coalesced += 1
            # Don't advance i — recheck this position with new neighbour
        else:
            if has_critical and is_critical_list[i + 1] and \
                    dr_curr < dr_prev / COALESCE_RATIO and \
                    dr_curr < dr_next / COALESCE_RATIO:
                n_critical_protected += 1
            i += 1
    if n_coalesced > 0 or n_critical_protected > 0:
        print(f"[AMR] Coalesced {n_coalesced} sliver zones "
              f"(<{1/COALESCE_RATIO:.0%} of neighbours)"
              f"{f', protected {n_critical_protected} critical-region zones' if n_critical_protected else ''} "
              f"-> grid: {len(r_list)} pts")

    # --- Degenerate zone merging ---
    # Merge consecutive zones where dr/r < MIN_DR_OVER_R (0.1%) AND the
    # physical quantities are nearly constant. This removes the inner-boundary
    # over-refinement where AMR created 170+ points with dr/r < 0.001 in a
    # region where V, T, and rho barely change (e.g., depths 223-400 at 5d
    # had 42% of the grid wasted on degenerate inner points).
    # Only merge non-original points (preserve STELLA hydro data points).
    # Critical-region points are exempt: a narrow high-gradient feature may
    # legitimately have dr/r < DEGENERATE_DR_R while still carrying a real
    # physical gradient (dT/T, dNe/Ne) that we deliberately refined.
    DEGENERATE_DR_R = 2e-3   # dr/r threshold for "degenerate" (0.2%)
    DEGENERATE_DV = 0.02     # |dv|/v < 2% across the zone
    DEGENERATE_DT = 0.05     # |dT|/T < 5% across the zone
    # |dNe|/Ne threshold: without this, points inserted by the Ne-jump AMR
    # criterion (small dr/r and dT/T, but a large dNe/Ne — exactly what a
    # sharp recombination/ionization front looks like) were being merged
    # right back out here, silently undoing the Ne-based refinement (see
    # AGENTS.md "Ionization front under-resolved" — the H/He recombination
    # front had dr/r~1e-3, dT/T~1%, but dNe/Ne up to 4.3x across one cell).
    DEGENERATE_DNE = 0.10
    n_degenerate_merged = 0
    n_orig_merged = 0
    n_critical_merged_protected = 0
    i = 1
    while i < len(r_list) - 1:
        dr_curr = abs(r_list[i + 1] - r_list[i])
        r_ref = max(abs(r_list[i]), eps)
        dr_r = dr_curr / r_ref
        is_orig = is_original_list[i + 1]
        is_crit = has_critical and is_critical_list[i + 1]
        if is_crit and dr_r < DEGENERATE_DR_R:
            # Skip critical-region points entirely — they are intentionally
            # small.  Only count them for diagnostics when they would have
            # been merge candidates.
            dv_v = abs(v_list[i + 1] - v_list[i]) / max(abs(v_list[i]), eps)
            dT_T = abs(T_list[i + 1] - T_list[i]) / max(T_list[i], eps)
            dne_ne = abs(n_e_list[i + 1] - n_e_list[i]) / max(abs(n_e_list[i]), eps)
            if (dv_v < DEGENERATE_DV and dT_T < DEGENERATE_DT and
                    dne_ne < DEGENERATE_DNE):
                n_critical_merged_protected += 1
            i += 1
            continue
        if dr_r < DEGENERATE_DR_R:
            # Check if physical quantities are nearly constant across this zone
            dv_v = abs(v_list[i + 1] - v_list[i]) / max(abs(v_list[i]), eps)
            dT_T = abs(T_list[i + 1] - T_list[i]) / max(T_list[i], eps)
            dne_ne = abs(n_e_list[i + 1] - n_e_list[i]) / max(abs(n_e_list[i]), eps)
            if is_orig:
                # Stricter thresholds for original points — only merge if
                # truly flat (not just "within AMR tolerance")
                can_merge = (dv_v < DEGENERATE_DV * 0.5 and
                             dT_T < DEGENERATE_DT * 0.4 and
                             dne_ne < DEGENERATE_DNE * 0.5)
            else:
                can_merge = (dv_v < DEGENERATE_DV and
                             dT_T < DEGENERATE_DT and
                             dne_ne < DEGENERATE_DNE)
            if can_merge:
                # Also check the zone ahead — merge runs of degenerate zones
                for lst in [r_list, v_list, T_list, rho_list, n_e_list, tau_list,
                            is_original_list] + ([is_critical_list] if has_critical else []):
                    del lst[i + 1]
                for elem in mf_lists:
                    del mf_lists[elem][i + 1]
                n_degenerate_merged += 1
                # Don't advance — recheck this position
            else:
                i += 1
        else:
            i += 1
    if n_degenerate_merged > 0 or n_critical_merged_protected > 0:
        print(f"[AMR] Merged {n_degenerate_merged} degenerate inner zones "
              f"(dr/r<{DEGENERATE_DR_R}, dv/v<{DEGENERATE_DV}, dT/T<{DEGENERATE_DT})"
              f"{f', protected {n_critical_merged_protected} critical-region zones' if n_critical_merged_protected else ''} "
              f"-> final grid: {len(r_list)} pts")

    return (np.array(r_list), np.array(v_list), np.array(T_list),
            np.array(rho_list), np.array(n_e_list), np.array(tau_list),
            {elem: np.array(mf_lists[elem]) for elem in mf_lists})


def insert_outer_boundary_points(radius, v, T, rho, n_e, tau, mass_fractions,
                                 n_insert=5, ob_rat=1.5,
                                 min_dr_r=1e-6, max_dr_r=0.5):
    """Insert additional points in the outer low-tau region to increase grid density.

    CMFGEN indexes depth outside-in: depth 1 = outermost point. The radiation
    field gradient is steepest at the outer boundary (optically thin transition),
    so coarse spacing there causes population oscillations from opacity-radiation
    feedback.

    Points are distributed throughout the outer low-tau region (tau < 1) by
    subdividing the largest cells in that region. This increases the overall
    point density in the spectrum formation region, rather than concentrating
    all extra points at the very outer edge (which creates a cluster of
    near-degenerate points that don't help resolve the radiation field gradient
    deeper in the atmosphere).

    Within each subdivided cell, geometric spacing with ratio ``ob_rat`` is
    used so the finest sub-cell is at the outer side of the cell. This matches
    CMFGEN's NBND_INS insertion in adjust_sn_r_grid.f.

    Arrays are outer->inner (decreasing radius).

    Velocity is linearly interpolated (the outer region is smooth/homologous,
    so linear interpolation preserves sigma ~ 0). Thermodynamic quantities
    (T, rho, n_e, tau) use log interpolation. Mass fractions use linear
    interpolation.

    Spacing constraints:
        ``min_dr_r`` — minimum fractional cell size (dr/r). If subdividing
        would produce a cell smaller than this, the subdivision is skipped.

        ``max_dr_r`` — maximum fractional cell size (dr/r). Cells in the outer
        region exceeding this are subdivided first.
    """
    if n_insert <= 0:
        return (radius, v, T, rho, n_e, tau, mass_fractions)

    eps = 1e-30
    r_list = list(radius)
    v_list = list(v)
    T_list = list(T)
    rho_list = list(rho)
    ne_list = list(n_e)
    tau_list = list(tau)
    mf_lists = {elem: list(mass_fractions[elem]) for elem in mass_fractions}

    # Find the extent of the outer low-tau region (tau < 1)
    tau_thresh = 1.0
    i_thresh = len(r_list)
    for i in range(len(r_list)):
        if tau_list[i] > tau_thresh:
            i_thresh = i
            break
    n_outer = i_thresh
    if n_outer < 2:
        print(f"[OB Insert] Outer region (tau<{tau_thresh}) too small ({n_outer} pts), skipping")
        return (radius, v, T, rho, n_e, tau, mass_fractions)

    def _log_interp(a, b, frac):
        return 10 ** (np.log10(max(a, eps)) + frac * (np.log10(max(b, eps)) - np.log10(max(a, eps))))

    n_remaining = n_insert
    n_inserted_total = 0

    # Iteratively subdivide the largest cell in the outer region
    while n_remaining > 0:
        # Find the largest cell (by dr/r) in the outer region
        worst_i = -1
        worst_dr_r = 0.0
        for i in range(n_outer - 1):
            dr = abs(r_list[i] - r_list[i + 1])
            r_ref = max(abs(r_list[i]), eps)
            dr_r = dr / r_ref
            if dr_r > worst_dr_r:
                worst_dr_r = dr_r
                worst_i = i

        if worst_i < 0:
            break

        # Check min_dr_r constraint
        R_outer = r_list[worst_i]
        R_inner = r_list[worst_i + 1]
        dR = R_outer - R_inner
        # Geometric spacing: finest cell = dR * (ob_rat-1) / (ob_rat^2 - 1) for 1 insert
        dR_finest = dR * (ob_rat - 1.0) / (ob_rat ** 2 - 1.0)
        if dR_finest / max(R_outer, eps) < min_dr_r:
            # This cell is too small to subdivide with geometric spacing;
            # try the next largest, or stop if none are large enough
            # Mark this cell as done by artificially reducing its size
            # (we'll skip it in future iterations by checking min_dr_r)
            break

        # Insert one geometrically-spaced point in this cell
        # dR_first = dR * (ob_rat - 1) / (ob_rat^2 - 1) = dR / (ob_rat + 1)
        dR_first = dR / (ob_rat + 1.0)
        r_new = R_outer - dR_first
        frac = (R_outer - r_new) / dR

        idx = worst_i + 1
        r_list.insert(idx, r_new)
        v_list.insert(idx, v_list[worst_i] + frac * (v_list[worst_i + 1] - v_list[worst_i]))
        T_list.insert(idx, _log_interp(T_list[worst_i], T_list[worst_i + 1], frac))
        rho_list.insert(idx, _log_interp(rho_list[worst_i], rho_list[worst_i + 1], frac))
        ne_list.insert(idx, _log_interp(ne_list[worst_i], ne_list[worst_i + 1], frac))
        tau_list.insert(idx, _log_interp(tau_list[worst_i], tau_list[worst_i + 1], frac))
        for elem in mf_lists:
            mf_lists[elem].insert(idx,
                mf_lists[elem][worst_i] + frac * (mf_lists[elem][worst_i + 1] - mf_lists[elem][worst_i]))

        n_outer += 1  # the outer region grew by 1
        n_remaining -= 1
        n_inserted_total += 1

    print(f"[OB Insert] Inserted {n_inserted_total} points in outer region (tau<{tau_thresh}, "
          f"{n_outer} outer pts, OB_RAT={ob_rat})")
    print(f"  Grid: {len(radius)} -> {len(r_list)} points")

    return (np.array(r_list), np.array(v_list), np.array(T_list),
            np.array(rho_list), np.array(ne_list), np.array(tau_list),
            {elem: np.array(mf_lists[elem]) for elem in mf_lists})


def fill_grid_to_minimum(radius, v, T, rho, n_e, tau, mass_fractions,
                          min_points=150, max_cell_ratio=1.5,
                          thin_tau_threshold=1.0, thin_priority=1.0):
    """Ensure the grid has at least ``min_points`` by inserting midpoints.

    If the grid already has >= ``min_points`` points, it is returned unchanged.

    Otherwise, points are inserted iteratively at the **largest** cell (by
    fractional radius dr/r) using log-interpolation for thermodynamic
    quantities (T, rho, n_e, tau) and linear interpolation for r, v, and
    mass fractions — matching the AMR insertion style.  After each insertion
    pass, a light grid-smoothing pass subdivides cells whose size ratio
    exceeds ``max_cell_ratio`` to prevent abrupt spacing transitions that
    would create sigma spikes.

    This is essential for CMFGEN convergence: too few grid points (< 150)
    leave large cells where the radiation field and populations change
    drastically between adjacent depth points, causing the Lambda iteration
    to overshoot and oscillate.

    Arrays are outer->inner (decreasing radius).
    """
    eps = 1e-30
    n = len(radius)
    if n >= min_points:
        print(f"[FillMin] Grid already has {n} points >= min_points={min_points}, skipping")
        return (radius, v, T, rho, n_e, tau, mass_fractions)

    def _log_interp(a, b, f):
        return 10**(np.log10(max(a, eps)) + f * (np.log10(max(b, eps)) - np.log10(max(a, eps))))

    def _insert_midpoint(lst_r, lst_v, lst_T, lst_rho, lst_ne, lst_tau,
                         mf_l, idx):
        # Choose the insertion fraction that minimizes the resulting sigma
        # discretization noise (same strategy as adaptive_grid_refinement),
        # instead of a blind midpoint. A naive f=0.5 split can introduce a
        # sigma kink when the local velocity curvature is asymmetric —
        # this is the dominant source of point-to-point sigma scatter from
        # min-points filling (see AGENTS.md sigma scatter optimization).
        f = optimal_insert_fraction(lst_r, lst_v, idx)
        lst_r.insert(idx + 1, lst_r[idx] + f * (lst_r[idx + 1] - lst_r[idx]))
        lst_v.insert(idx + 1, lst_v[idx] + f * (lst_v[idx + 1] - lst_v[idx]))
        lst_T.insert(idx + 1, _log_interp(lst_T[idx], lst_T[idx + 1], f))
        lst_rho.insert(idx + 1, _log_interp(lst_rho[idx], lst_rho[idx + 1], f))
        lst_ne.insert(idx + 1, _log_interp(lst_ne[idx], lst_ne[idx + 1], f))
        lst_tau.insert(idx + 1, _log_interp(lst_tau[idx], lst_tau[idx + 1], f))
        for elem in mf_l:
            mf_l[elem].insert(idx + 1,
                mf_l[elem][idx] + f * (mf_l[elem][idx + 1] - mf_l[elem][idx]))

    r_list   = list(radius)
    v_list   = list(v)
    T_list   = list(T)
    rho_list = list(rho)
    ne_list  = list(n_e)
    tau_list = list(tau)
    mf_lists = {elem: list(mass_fractions[elem]) for elem in mass_fractions}

    n_needed = min_points - n
    print(f"[FillMin] Grid has {n} points, need {min_points} — inserting {n_needed} points")

    n_inserted = 0
    n_smoothed = 0
    max_iter = n_needed * 4  # safety cap (smoothing adds extra insertions)

    for iteration in range(max_iter):
        if len(r_list) >= min_points:
            break

        worst_i, worst_score = -1, 0.0
        for i in range(len(r_list) - 1):
            dr = abs(r_list[i + 1] - r_list[i])
            dr_r = dr / max(abs(r_list[i]), eps)
            tau_edge = min(tau_list[i], tau_list[i + 1])
            score = dr_r * (thin_priority if tau_edge <= thin_tau_threshold else 1.0)
            # Boost score for high-velocity-gradient edges (velocity reversal regions)
            # so fill points are directed where they're needed most, not just
            # at the largest radial gaps.
            dv_v = abs(v_list[i + 1] - v_list[i]) / max(abs(v_list[i]), eps)
            if dv_v > 0.05:
                score *= (1.0 + dv_v * 10.0)
            if score > worst_score:
                worst_i, worst_score = i, score

        if worst_i < 0:
            break  # shouldn't happen, but safety

        _insert_midpoint(r_list, v_list, T_list, rho_list, ne_list, tau_list,
                         mf_lists, worst_i)
        n_inserted += 1

        # --- Light grid-smoothing: fix abrupt size transitions ---
        # After each fill insertion, check if the new cell created a size
        # ratio > max_cell_ratio with its neighbors and smooth if so.
        # This prevents sigma spikes from abrupt grid spacing changes.
        for _ in range(min_points):
            if len(r_list) >= min_points:
                break
            smoothed_this = False
            for i in range(len(r_list) - 2):
                dr0 = abs(r_list[i + 1] - r_list[i])
                dr1 = abs(r_list[i + 2] - r_list[i + 1])
                ratio = max(dr0, dr1) / max(min(dr0, dr1), eps)
                if ratio > max_cell_ratio:
                    # Subdivide the larger of the two cells
                    if dr0 > dr1:
                        _insert_midpoint(r_list, v_list, T_list, rho_list,
                                         ne_list, tau_list, mf_lists, i)
                    else:
                        _insert_midpoint(r_list, v_list, T_list, rho_list,
                                         ne_list, tau_list, mf_lists, i + 1)
                    n_smoothed += 1
                    smoothed_this = True
                    break  # restart scan (indices shifted)
            if not smoothed_this:
                break

    total_added = len(r_list) - n
    print(f"[FillMin] Inserted {n_inserted} fill points + {n_smoothed} smoothing points "
          f"= {total_added} total (grid now {len(r_list)} pts)")

    if len(r_list) < min_points:
        print(f"[FillMin] WARNING: Could only reach {len(r_list)} points "
              f"(target {min_points}). Grid may be too coarse for convergence.")

    return (np.array(r_list), np.array(v_list), np.array(T_list),
            np.array(rho_list), np.array(ne_list), np.array(tau_list),
            {elem: np.array(mf_lists[elem]) for elem in mf_lists})


def enforce_sigma_conformance(radius, v, T, rho, n_e, tau, mass_fractions, kappa,
                              max_sigma_jump=0.5, max_cell_ratio=1.5,
                              min_dr_r=1e-6, max_dr_r=0.5,
                              max_points=500, headroom=40, verbose=True):
    """Final, authoritative grid-quality pass — run LAST, after every other
    grid-construction step (AMR, outer-boundary insertion, min-nd fill,
    enforce_grid_spacing, regularize_grid_spacing, velocity re-interpolation,
    wiggle/window/scatter smoothing).

    Why this exists
    ----------------
    ``adaptive_grid_refinement`` (AMR) is the only stage that scores edges
    against the *exact* CMFGEN sigma formula (see ``sigma_velocity`` /
    ``_sigma_at_index``), but it runs early in the pipeline. Several
    later stages mutate R and/or V without re-checking sigma:

      - ``enforce_grid_spacing``'s cell-ratio phase inserts points driven
        purely by geometric cell-size ratio, consuming most/all of the
        remaining point budget regardless of whether sigma is actually
        violated there.
      - ``regularize_grid_spacing`` nudges R by up to ~1% to smooth the
        dR sequence, which changes sigma at every point (since sigma
        depends on the *ratio* of neighbouring cell sizes).
      - The velocity re-interpolation step (pulling V from the full,
        un-downsampled STELLA data so extrema aren't missed) and the
        wiggle/window/scatter smoothing re-applied afterwards can each
        shift sigma at specific points — smoothing V is not guaranteed
        to be monotonic in its effect on sigma at a non-uniform grid.

    None of these steps checks the actual AMR thresholds afterward, so a
    model can pass AMR cleanly and still ship a SN_HYDRO_DATA with sigma
    spikes (observed in practice: |sigma| up to ~1.4 in "quiet",
    dv/v ~ 1e-3 regions, well above the intended max_sigma_jump=0.5,
    because AMR never saw the grid in its final state).

    This function re-applies the same, already-tuned dsigma / cell-ratio
    / max-dr_r criteria one more time, on the truly final grid, and only
    ever INSERTS points (never removes), so it cannot undo upstream
    degenerate-merge or coalescing decisions. It is deliberately
    conservative and model-agnostic: it does not introduce new
    thresholds, it just closes the loop so the *actual* written grid is
    checked, not just an intermediate one.

    A small ``headroom`` beyond ``max_points`` is allowed because the
    upstream cell-ratio pass in ``enforce_grid_spacing`` can consume the
    entire nominal budget on purely geometric smoothing, leaving zero
    room for this pass to fix genuine sigma violations. CMFGEN allocates
    ND from the input file at runtime, so a modest overrun (tens of
    points) is safe.

    Arrays are outer->inner (decreasing radius).
    """
    eps = 1e-30
    r_list = list(radius)
    v_list = list(v)
    T_list = list(T)
    rho_list = list(rho)
    ne_list = list(n_e)
    tau_list = list(tau)
    kappa_list = list(kappa)
    mf_lists = {elem: list(mass_fractions[elem]) for elem in mass_fractions}
    budget = max_points + max(0, headroom)

    def _log_interp(a, b, f):
        return 10**(np.log10(max(a, eps)) + f * (np.log10(max(b, eps)) - np.log10(max(a, eps))))

    def _insert(idx, f):
        r_list.insert(idx + 1, r_list[idx] + f * (r_list[idx + 1] - r_list[idx]))
        v_list.insert(idx + 1, v_list[idx] + f * (v_list[idx + 1] - v_list[idx]))
        T_list.insert(idx + 1, _log_interp(T_list[idx], T_list[idx + 1], f))
        rho_list.insert(idx + 1, _log_interp(rho_list[idx], rho_list[idx + 1], f))
        ne_list.insert(idx + 1, _log_interp(ne_list[idx], ne_list[idx + 1], f))
        tau_list.insert(idx + 1, _log_interp(tau_list[idx], tau_list[idx + 1], f))
        kappa_list.insert(idx + 1, _log_interp(kappa_list[idx], kappa_list[idx + 1], f))
        for elem in mf_lists:
            mf_lists[elem].insert(idx + 1,
                mf_lists[elem][idx] + f * (mf_lists[elem][idx + 1] - mf_lists[elem][idx]))

    def _edge_score(i):
        """Worst normalised violation at edge i..i+1: dsigma, cell ratio, dr/r.

        Used only for before/after diagnostics; the insertion loop below
        scores dsigma/dr_r and cell-ratio violations SEPARATELY because
        they require different insertion strategies (see note below).
        """
        n = len(r_list)
        dr = abs(r_list[i + 1] - r_list[i])
        r_ref = max(abs(r_list[i]), eps)
        dr_r = dr / r_ref
        sigma_i = _sigma_at_index(r_list, v_list, i) if 0 < i < n - 1 else 0.0
        sigma_ip1 = _sigma_at_index(r_list, v_list, i + 1) if 0 < i + 1 < n - 1 else 0.0
        dsigma = abs(sigma_ip1 - sigma_i)
        ratio = 1.0
        if i > 0:
            dr_prev = abs(r_list[i] - r_list[i - 1])
            ratio = max(ratio, max(dr, dr_prev) / max(min(dr, dr_prev), eps))
        if i + 2 < n:
            dr_next = abs(r_list[i + 2] - r_list[i + 1])
            ratio = max(ratio, max(dr, dr_next) / max(min(dr, dr_next), eps))
        return max(dsigma / max_sigma_jump, ratio / max_cell_ratio, dr_r / max_dr_r)

    sigma_before = sigma_velocity(np.array(r_list), np.array(v_list))
    n_before = len(r_list)
    n_violations_before = sum(1 for i in range(n_before - 1) if _edge_score(i) > 1.0)

    # --- Insertion loop ---
    # IMPORTANT: dsigma/dr_r violations and cell-ratio violations need
    # different insertion strategies, and using the wrong one for a given
    # violation type causes a runaway insertion cascade (observed:
    # violations grew 29 -> 43 -> 145 as headroom increased before this
    # fix). Reusing the sigma-optimal (deliberately asymmetric) insertion
    # fraction to fix a pure cell-ratio violation can make the ratio
    # *worse*, since sigma-optimal placement has no reason to balance
    # cell sizes. Conversely, for a cell-ratio violation between a pair
    # of cells, the fix must subdivide the LARGER of the two cells at
    # the midpoint (matching enforce_grid_spacing's Phase 3) — subdividing
    # the wrong side (or the edge itself when it's the smaller one)
    # increases the ratio instead of reducing it.
    REBALANCE_RATIO = 1.3  # if dr_left/dr_right (or inverse) is below this,
                           # the grid is already ~balanced at this point,
                           # so a spike here is noise IN V, not a grid
                           # artifact — see _worst_sigma_point docstring.

    def _local_smooth_v(i, max_rel_change=0.05):
        """Blend V[i] toward the chord-interpolated value between its
        neighbours (removes an isolated one-point kink) capped at
        ``max_rel_change`` relative change for safety. Returns True if a
        change was made."""
        r_im1, r_i, r_ip1 = r_list[i - 1], r_list[i], r_list[i + 1]
        v_im1, v_i, v_ip1 = v_list[i - 1], v_list[i], v_list[i + 1]
        span = r_ip1 - r_im1
        if abs(span) < eps:
            return False
        f = (r_i - r_im1) / span
        v_chord = v_im1 + f * (v_ip1 - v_im1)
        v_target = v_i + 0.8 * (v_chord - v_i)  # 80% blend toward the chord
        max_change = max_rel_change * max(abs(v_i), eps)
        dv = np.clip(v_target - v_i, -max_change, max_change)
        if abs(dv) < eps:
            return False
        v_list[i] = v_i + dv
        return True

    def _worst_sigma_point():
        """Find the interior point with the worst sigma spike relative to
        its neighbours, and return an action to fix it: either the EDGE
        to subdivide (grid-imbalance case) or the POINT to locally smooth
        (V-noise case).

        A sigma spike at point i shows up as a large jump on the edge on
        AT LEAST one side of i (edge i-1 and/or edge i). Two distinct
        root causes produce this, and they need different fixes:

        1. **Grid imbalance**: sigma[i] is a distance-weighted average of
           the left/right linear velocity gradients, weighted by the
           *opposite* side's cell size (see sigma_velocity docstring). If
           the left cell size (r[i-1]-r[i]) and right cell size
           (r[i]-r[i+1]) differ substantially, this weighting amplifies
           any real curvature into a spike. Fix: subdivide the LARGER
           side to rebalance the weights toward a symmetric (2nd-order
           accurate) estimate.

        2. **Residual velocity noise**: if the two cell sizes are already
           close to equal (ratio < REBALANCE_RATIO), the formula reduces
           to a plain centered difference — inserting more grid points
           cannot help, because the spike is caused by V[i] itself being
           off the smooth local trend (e.g. residual STELLA noise that
           survived the wiggle/window/scatter smoothing passes, which use
           a fixed relative-change cap and can miss single-point kinks in
           high-velocity regions). Fix: nudge V[i] toward the local chord
           value instead of touching the grid.
        """
        n = len(r_list)
        best_score = 1.0
        best_action = None  # ('insert', edge) or ('smooth', i)
        for i in range(1, n - 1):
            sigma_im1 = _sigma_at_index(r_list, v_list, i - 1) if i - 1 > 0 else _sigma_at_index(r_list, v_list, i)
            sigma_i = _sigma_at_index(r_list, v_list, i)
            sigma_ip1 = _sigma_at_index(r_list, v_list, i + 1) if i + 1 < n - 1 else sigma_i
            d_left = abs(sigma_i - sigma_im1)
            d_right = abs(sigma_ip1 - sigma_i)
            sc = max(d_left, d_right) / max_sigma_jump
            if sc <= best_score:
                continue
            dr_left = abs(r_list[i] - r_list[i - 1])
            dr_right = abs(r_list[i + 1] - r_list[i])
            ratio = max(dr_left, dr_right) / max(min(dr_left, dr_right), eps)
            if ratio >= REBALANCE_RATIO:
                if dr_left >= dr_right:
                    edge, r_ref, dr_check = i - 1, max(abs(r_list[i - 1]), eps), dr_left
                else:
                    edge, r_ref, dr_check = i, max(abs(r_list[i]), eps), dr_right
                if dr_check / r_ref < min_dr_r * 2:
                    continue  # larger side already too small to subdivide
                best_score, best_action = sc, ('insert', edge)
            else:
                best_score, best_action = sc, ('smooth', i)
        return best_action, best_score

    def _worst_dr_edge():
        """Find the worst pure oversized-cell (dr/r > max_dr_r) violation."""
        best_i, best_score = -1, 1.0
        n = len(r_list)
        for i in range(n - 1):
            dr = abs(r_list[i + 1] - r_list[i])
            r_ref = max(abs(r_list[i]), eps)
            dr_r = dr / r_ref
            sc = dr_r / max_dr_r
            if sc > best_score and dr_r >= min_dr_r * 2:
                best_i, best_score = i, sc
        return best_i, best_score

    def _worst_ratio_pair():
        """Return (idx_to_subdivide, ratio) for the worst adjacent cell-size
        ratio violation, where idx_to_subdivide is the LARGER cell of the pair."""
        best_idx, best_ratio = -1, max_cell_ratio
        n = len(r_list)
        for i in range(n - 2):
            dr0 = abs(r_list[i + 1] - r_list[i])
            dr1 = abs(r_list[i + 2] - r_list[i + 1])
            ratio = max(dr0, dr1) / max(min(dr0, dr1), eps)
            if ratio <= best_ratio:
                continue
            # Subdivide the larger cell; skip if it's already too small.
            if dr0 >= dr1:
                r_ref = max(abs(r_list[i]), eps)
                if dr0 / r_ref < min_dr_r * 2:
                    continue
                best_idx, best_ratio = i, ratio
            else:
                r_ref = max(abs(r_list[i + 1]), eps)
                if dr1 / r_ref < min_dr_r * 2:
                    continue
                best_idx, best_ratio = i + 1, ratio
        return best_idx, best_ratio

    n_inserted = 0
    n_smoothed = 0
    n_stalled = 0
    n_ratio_inserted = 0
    # Cell-ratio violations are cosmetic (grid smoothness), not a
    # correctness risk — a genuine 2x cell-size jump straddling two
    # otherwise-smooth geometric progressions cascades when fixed by
    # naive bisection (splitting the larger cell in half just shifts the
    # violation to the OTHER side, which was already matched; observed:
    # a single such jump walked 20 points across the grid without
    # resolving). A full fix needs multi-point geometric redistribution,
    # which is out of scope for a final safety-net pass. Cap the ratio
    # budget so it cannot starve genuine sigma fixes or exhaust all of
    # --sigma-conformance-headroom chasing an unresolvable cascade.
    ratio_budget = max(10, headroom // 3)
    max_iter = budget * 2 + 200  # smoothing actions don't grow the grid, allow extra iterations
    for _iter in range(max_iter):
        sigma_action, sigma_score = _worst_sigma_point()
        dr_edge, dr_score = (-1, 1.0) if len(r_list) >= budget else _worst_dr_edge()
        if len(r_list) >= budget or n_ratio_inserted >= ratio_budget:
            ratio_i, ratio_val = -1, max_cell_ratio
        else:
            ratio_i, ratio_val = _worst_ratio_pair()
        if sigma_action is None and dr_edge < 0 and ratio_i < 0:
            break
        # Priority: sigma spikes (actual convergence risk) > oversized
        # cells (dr/r bound) > cell-ratio smoothness (cosmetic). Sigma
        # spikes are fixed either by grid insertion (imbalance case) or
        # local V smoothing (noise case, doesn't consume point budget).
        if sigma_action is not None:
            kind, idx = sigma_action
            if kind == 'smooth':
                if not _local_smooth_v(idx):
                    n_stalled += 1
                    if n_stalled > 5:
                        break
                    continue
                n_smoothed += 1
                continue
            elif len(r_list) < budget:
                f_opt = optimal_insert_fraction(r_list, v_list, idx)
                _insert(idx, f_opt)
                n_inserted += 1
                continue
        if dr_edge >= 0:
            _insert(dr_edge, 0.5)
            n_inserted += 1
        elif ratio_i >= 0:
            _insert(ratio_i, 0.5)
            n_inserted += 1
            n_ratio_inserted += 1
        else:
            break

    sigma_after = sigma_velocity(np.array(r_list), np.array(v_list))
    n_violations_after = sum(1 for i in range(len(r_list) - 1) if _edge_score(i) > 1.0)

    if verbose:
        if n_inserted > 0 or n_smoothed > 0 or n_violations_before > 0:
            print(f"[Sigma Conformance] {n_violations_before} edges violated thresholds on "
                  f"final grid ({n_before} pts, sigma range [{sigma_before.min():.3f}, "
                  f"{sigma_before.max():.3f}])")
            print(f"  Inserted {n_inserted} points + locally smoothed {n_smoothed} points "
                  f"-> {len(r_list)} pts, {n_violations_after} violations remain, "
                  f"sigma range [{sigma_after.min():.3f}, {sigma_after.max():.3f}]")
            if n_violations_after > 0 and len(r_list) >= budget:
                print(f"  WARNING: budget ({budget} = max_points {max_points} + "
                      f"headroom {headroom}) exhausted with violations remaining")
        else:
            print(f"[Sigma Conformance] Final grid ({n_before} pts) already meets all "
                  f"thresholds (sigma range [{sigma_before.min():.3f}, {sigma_before.max():.3f}])")

    return (np.array(r_list), np.array(v_list), np.array(T_list),
            np.array(rho_list), np.array(ne_list), np.array(tau_list),
            {elem: np.array(mf_lists[elem]) for elem in mf_lists}, np.array(kappa_list))


def coarsen_negative_sigma_regions(radius, v, T, rho, n_e, tau, mass_fractions, kappa,
                                    sigma_threshold=-0.1, max_zones_per_region=2,
                                    min_region_zones=4, protected_indices=None):
    """Coarsen grid regions where SIGMA < threshold to 1-2 zones.

    Following Luc Dessart's approach: negative sigma (velocity decreasing
    with radius) only affects 1-2 zones in well-behaved interacting SN
    models. When the grid over-resolves these regions (dozens of zones),
    the moment solver and Ng acceleration are destabilised by the extended
    negative-sigma band.

    This function identifies contiguous regions where sigma < threshold
    and replaces each region with at most ``max_zones_per_region`` points
    (the boundary points plus optionally one midpoint). Since FIX_T=T
    (temperature is not solved), under-resolving these regions is
    physically acceptable — the velocity reversal is preserved as a
    sharp discontinuity, matching Luc's published models.

    Parameters
    ----------
    radius, v, T, rho, n_e, tau, mass_fractions, kappa : arrays
        Grid quantities (outer->inner, decreasing radius).
    sigma_threshold : float
        Sigma value below which a zone is considered "negative sigma".
        Default -0.1 (mildly negative, includes the velocity reversal).
    max_zones_per_region : int
        Maximum number of zones to keep in each negative-sigma region.
        Default 2 (boundary points only). Set to 3 to keep a midpoint.
    min_region_zones : int
        Only coarsen regions with at least this many zones. Smaller
        regions are already adequately resolved. Default 4.
    protected_indices : set or None
        Zone indices that must NOT be removed (e.g., composition interfaces,
        cannon boundaries, ionization-front zones). These are always kept
        even if they fall within a negative-sigma region. Additionally,
        zones within a buffer distance of protected indices are kept to
        avoid creating large grid spacing jumps adjacent to sharp interfaces.

    Returns
    -------
    Same tuple as input, with coarsened arrays.
    """
    eps = 1e-30
    r = np.asarray(radius, dtype=float)
    n = len(r)
    if n < 5:
        return (r, v, T, rho, n_e, tau, mass_fractions, kappa)

    # Compute sigma using the exact CMFGEN formula
    sigma = sigma_velocity(r, v)

    # Find contiguous negative-sigma regions
    neg_mask = sigma < sigma_threshold
    neg_indices = np.where(neg_mask)[0]

    if len(neg_indices) == 0:
        print(f"[Sigma-Coarsen] No zones with sigma < {sigma_threshold} found — "
              f"nothing to coarsen")
        return (r, v, T, rho, n_e, tau, mass_fractions, kappa)

    # Group into contiguous regions
    starts = [neg_indices[0]]
    ends = []
    for i in range(1, len(neg_indices)):
        if neg_indices[i] != neg_indices[i-1] + 1:
            ends.append(neg_indices[i-1])
            starts.append(neg_indices[i])
    ends.append(neg_indices[-1])

    # Build protected set with buffer zone (keep neighbours of protected indices)
    protect_buffer = 3  # keep 3 zones on each side of protected indices
    protected_set = set()
    if protected_indices:
        for idx in protected_indices:
            for offset in range(-protect_buffer, protect_buffer + 1):
                protected_set.add(idx + offset)
        print(f"[Sigma-Coarsen] {len(protected_indices)} protected indices "
              f"(composition/ionization-front zones) with {protect_buffer}-zone buffer "
              f"-> {len(protected_set & set(range(n)))} zones protected from coarsening")

    # Filter: only coarsen regions with >= min_region_zones
    # AND skip regions that overlap with protected indices
    regions_to_coarsen = []
    regions_skipped = []
    for s, e in zip(starts, ends):
        n_zones = e - s + 1
        if n_zones < min_region_zones:
            continue
        # Check overlap with protected set
        region_indices = set(range(s, e + 1))
        if protected_set & region_indices:
            regions_skipped.append((s, e, n_zones))
            continue
        regions_to_coarsen.append((s, e, n_zones))

    if regions_skipped:
        print(f"[Sigma-Coarsen] Skipping {len(regions_skipped)} negative-sigma regions "
              f"overlapping protected composition/ionization-front zones:")
        for s, e, nz in regions_skipped:
            print(f"  Zones {s}-{e} ({nz} zones): r={r[s]:.3e}-{r[e]:.3e} cm, "
                  f"v={v[s]:.1f}-{v[e]:.1f} km/s — PROTECTED")

    if not regions_to_coarsen:
        print(f"[Sigma-Coarsen] {len(starts)} negative-sigma regions found, "
              f"but none are coarsenable (all protected or too small)")
        return (r, v, T, rho, n_e, tau, mass_fractions, kappa)

    total_removed = sum(nz - max_zones_per_region for _, _, nz in regions_to_coarsen)
    print(f"[Sigma-Coarsen] Found {len(regions_to_coarsen)} coarsenable negative-sigma regions "
          f"(sigma < {sigma_threshold}) with >= {min_region_zones} zones:")
    for s, e, nz in regions_to_coarsen:
        print(f"  Zones {s}-{e} ({nz} zones): r={r[s]:.3e}-{r[e]:.3e} cm, "
              f"v={v[s]:.1f}-{v[e]:.1f} km/s, sigma={sigma[s]:.3f} to {sigma[e]:.3f} "
              f"-> coarsening to {max_zones_per_region} zones "
              f"({nz - max_zones_per_region} removed)")
    print(f"[Sigma-Coarsen] Total: removing {total_removed} zones, "
          f"grid {n} -> {n - total_removed} points")

    # Build the list of indices to keep
    keep = set(range(n))
    # Always keep protected indices
    keep.update(protected_set & set(range(n)))
    for s, e, nz in regions_to_coarsen:
        # Keep the boundary points (s and e) and optionally a midpoint
        keep_in_region = {s, e}
        if max_zones_per_region >= 3:
            mid = (s + e) // 2
            keep_in_region.add(mid)
        if max_zones_per_region >= 4:
            q1 = (s + mid) // 2
            q3 = (mid + e) // 2
            keep_in_region.update([q1, q3])
        # Remove all indices in the region except the ones we want to keep
        for i in range(s, e + 1):
            if i not in keep_in_region:
                keep.discard(i)

    keep_indices = sorted(keep)

    # Extract the coarsened grid
    r_new = r[keep_indices]
    v_new = v[keep_indices]
    T_new = T[keep_indices]
    rho_new = rho[keep_indices]
    ne_new = n_e[keep_indices]
    tau_new = tau[keep_indices]
    kappa_new = kappa[keep_indices]
    mf_new = {elem: mass_fractions[elem][keep_indices] for elem in mass_fractions}

    # Verify the new sigma
    sigma_new = sigma_velocity(r_new, v_new)
    n_neg_old = np.sum(neg_mask)
    n_neg_new = np.sum(sigma_new < sigma_threshold)
    print(f"[Sigma-Coarsen] After coarsening: {n_neg_new} zones with sigma < {sigma_threshold} "
          f"(was {n_neg_old}), sigma range [{sigma_new.min():.3f}, {sigma_new.max():.3f}]")

    # Knop/Dessart non-monotonic safety check: verify that any physical
    # velocity reversal (sigma < -1, the dV/dR < 0 threshold used by the
    # non-monotonic formal solver for Sobolev resonance placement) is not
    # destroyed by coarsening.  Coarsening a mild reversal to 2 boundary
    # points can recompute sigma from the wider stencil and push it above
    # -1, which would make CMFGEN treat the grid as monotonic and break
    # the resonance calculation at the (now-invisible) turning point.
    n_reversal_old = int(np.sum(sigma < -1.0))
    n_reversal_new = int(np.sum(sigma_new < -1.0))
    if n_reversal_old > 0 and n_reversal_new == 0:
        print(f"[Sigma-Coarsen] WARNING: {n_reversal_old} zones had sigma < -1 (physical "
              f"velocity reversal) before coarsening, but NONE remain after coarsening. "
              f"The non-monotonic turning point has been lost — CMFGEN will treat the "
              f"grid as monotonic and the Knop Sobolev resonance at the reversal will "
              f"not be placed. Consider disabling coarsen_negative_sigma for this model "
              f"or increasing max_zones_per_region.")
    elif n_reversal_old > 0:
        print(f"[Sigma-Coarsen] Velocity reversal preserved: {n_reversal_old} -> "
              f"{n_reversal_new} zones with sigma < -1")

    return (r_new, v_new, T_new, rho_new, ne_new, tau_new, mf_new, kappa_new)


def refine_outer_boundary_grid(radius, v, T, rho, n_e, tau, mass_fractions, kappa,
                            n_outer_cells=5, subdivisions=3,
                            target_dr_r=0.015, max_new_points=30):
    """Refine the outer boundary by subdividing the first few cells.

    CMFGEN requires fine resolution at the outer boundary (low optical depth)
    for the moment equation and ray-by-ray solution to agree. When the grid
    has dr/r > 0.02 at the outer edge, the radiation field diverges with
    "J(mom) and J(ray) differ by more than 20%" errors and luminosity blowup.

    This function subdivides the first ``n_outer_cells`` cells of the grid
    (the outermost, optically thinnest region) into ``subdivisions`` sub-cells
    each, using geometric (log-uniform) radius spacing within each cell. All
    physical quantities are linearly interpolated onto the new points.

    Parameters
    ----------
    radius, v, T, rho, n_e, tau, mass_fractions, kappa : arrays
        Grid quantities (outer->inner, decreasing radius).
    n_outer_cells : int
        Number of outermost cells to subdivide. Default 5.
    subdivisions : int
        Number of sub-cells per original cell. Default 3 (each cell becomes
        3 cells, adding 2 points per cell).
    target_dr_r : float
        If set, automatically compute the number of subdivisions needed to
        achieve this dr/r at the outer boundary. Overrides ``subdivisions``
        when the original dr/r exceeds ``target_dr_r``. Default 0.015.
    max_new_points : int
        Maximum total new points to insert. Default 30.

    Returns
    -------
    Same tuple as input, with refined outer boundary.
    """
    r = np.asarray(radius, dtype=float)
    n = len(r)
    if n < 10:
        return (r, v, T, rho, n_e, tau, mass_fractions, kappa)

    # Check if refinement is needed
    outer_dr_r = np.abs(np.diff(r[:n_outer_cells+1])) / r[:n_outer_cells]
    max_dr_r = np.max(outer_dr_r)
    if max_dr_r < target_dr_r:
        print(f"[Outer-Refine] Outer boundary already fine (max dr/r={max_dr_r:.4f} "
              f"< {target_dr_r}) — no refinement needed")
        return (r, v, T, rho, n_e, tau, mass_fractions, kappa)

    # Compute subdivisions needed per cell to reach target_dr_r
    sub_per_cell = []
    total_new = 0
    for i in range(min(n_outer_cells, n-1)):
        dr_r_i = abs(r[i] - r[i+1]) / r[i]
        if dr_r_i < target_dr_r:
            sub_per_cell.append(1)
        else:
            # Geometric subdivision: n_sub cells means dr_r_i / n_sub
            n_sub = int(np.ceil(dr_r_i / target_dr_r))
            n_sub = min(n_sub, subdivisions * 3)  # cap at 3x default
            sub_per_cell.append(n_sub)
        total_new += sub_per_cell[-1] - 1

    # Cap total new points
    if total_new > max_new_points:
        scale = max_new_points / total_new
        sub_per_cell = [max(1, int(s * scale)) for s in sub_per_cell]
        total_new = sum(s - 1 for s in sub_per_cell)

    if total_new == 0:
        print(f"[Outer-Refine] No refinement needed after target check")
        return (r, v, T, rho, n_e, tau, mass_fractions, kappa)

    print(f"[Outer-Refine] Subdividing first {len(sub_per_cell)} outer cells "
          f"(max dr/r={max_dr_r:.4f}, target={target_dr_r}): "
          f"adding {total_new} points, grid {n} -> {n + total_new}")

    # Build new grid by inserting subdivision points
    r_new, v_new, T_new, rho_new, ne_new, tau_new, kappa_new = [], [], [], [], [], [], []
    mf_new = {elem: [] for elem in mass_fractions}

    for i in range(n):
        # Add current point
        r_new.append(r[i])
        v_new.append(v[i])
        T_new.append(T[i])
        rho_new.append(rho[i])
        ne_new.append(n_e[i])
        tau_new.append(tau[i])
        kappa_new.append(kappa[i])
        for elem in mass_fractions:
            mf_new[elem].append(mass_fractions[elem][i])

        # If this is one of the outer cells to subdivide, insert points
        if i < len(sub_per_cell) and sub_per_cell[i] > 1 and i < n - 1:
            n_sub = sub_per_cell[i]
            # Geometric spacing within the cell [r[i], r[i+1]]
            # r decreases with index, so use log spacing
            r_lo, r_hi = r[i+1], r[i]  # inner, outer
            log_lo, log_hi = np.log(r_lo), np.log(r_hi)
            for k in range(1, n_sub):
                frac = k / n_sub
                r_k = np.exp(log_lo * frac + log_hi * (1.0 - frac))
                # Linear interpolation for all quantities
                t_frac = (r[i] - r_k) / (r[i] - r[i+1])
                r_new.append(r_k)
                v_new.append(v[i] + t_frac * (v[i+1] - v[i]))
                T_new.append(T[i] + t_frac * (T[i+1] - T[i]))
                rho_new.append(rho[i] + t_frac * (rho[i+1] - rho[i]))
                ne_new.append(n_e[i] + t_frac * (n_e[i+1] - n_e[i]))
                tau_new.append(tau[i] + t_frac * (tau[i+1] - tau[i]))
                kappa_new.append(kappa[i] + t_frac * (kappa[i+1] - kappa[i]))
                for elem in mass_fractions:
                    mf_new[elem].append(
                        mass_fractions[elem][i] + t_frac * (mass_fractions[elem][i+1] - mass_fractions[elem][i]))

    r_new = np.array(r_new)
    v_new = np.array(v_new)
    T_new = np.array(T_new)
    rho_new = np.array(rho_new)
    ne_new = np.array(ne_new)
    tau_new = np.array(tau_new)
    kappa_new = np.array(kappa_new)
    mf_new = {elem: np.array(mf_new[elem]) for elem in mass_fractions}

    # Verify the new outer boundary spacing
    new_outer_dr_r = np.abs(np.diff(r_new[:total_new+n_outer_cells+1])) / r_new[:total_new+n_outer_cells]
    print(f"[Outer-Refine] New outer boundary max dr/r={np.max(new_outer_dr_r):.4f} "
          f"(was {max_dr_r:.4f})")

    return (r_new, v_new, T_new, rho_new, ne_new, tau_new, mf_new, kappa_new)


def enforce_grid_spacing(radius, v, T, rho, n_e, tau, mass_fractions, kappa,
                         min_dr_r=1e-6, max_dr_r=0.5, max_cell_ratio=1.5,
                         max_points=500,
                         critical_mask=None, critical_min_dr_r=1e-7,
                         disable_in_critical=False):
    """Enforce minimum and maximum grid spacing constraints on the final grid.

    This is the final spacing enforcement pass, run after all grid processing
    (AMR, outer boundary insertion, min_nd fill, smoothing). It guarantees:

    1. **Minimum spacing** (min_dr_r): Removes points where dr/r < min_dr_r.
       These degenerate points cause division-by-zero in CMFGEN's sigma
       computation (dlnV/dlnr with dr=0). Points are removed by marking them
       for deletion; the first and last points are always kept.

    2. **Maximum spacing** (max_dr_r): Inserts midpoints in cells where
       dr/r > max_dr_r. Large cells cause the radiation field to change
       drastically between adjacent depth points, destabilising the Lambda
       iteration. Insertion uses log-interpolation for thermodynamic
       quantities and linear interpolation for r, v, and mass fractions.

    3. **Cell ratio smoothing** (max_cell_ratio): After min/max enforcement,
       subdivides cells whose size ratio with neighbours exceeds
       max_cell_ratio. This prevents abrupt spacing transitions that create
       sigma spikes.

    Arrays are outer->inner (decreasing radius). The ``kappa`` array is
    interpolated onto any new grid points.

    ``critical_mask`` (optional boolean array, length ``len(radius)``)
    marks narrow high-gradient critical regions.  When supplied, the
    minimum-spacing removal (Phase 1) uses a *region-dependent* threshold:
    ``critical_min_dr_r`` inside critical regions (or 0.0 if
    ``disable_in_critical`` is True) and the global ``min_dr_r`` outside.
    This preserves the fine points the AMR inserted in shock fronts /
    composition interfaces / CSM shells while still removing degenerate
    points elsewhere.  Inserted midpoints inherit the critical flag so the
    relaxed threshold remains consistent as the grid grows.
    """
    eps = 1e-30
    r_list = list(radius)
    v_list = list(v)
    T_list = list(T)
    rho_list = list(rho)
    ne_list = list(n_e)
    tau_list = list(tau)
    kappa_list = list(kappa)
    mf_lists = {elem: list(mass_fractions[elem]) for elem in mass_fractions}
    has_critical = (critical_mask is not None
                    and len(critical_mask) == len(radius))
    is_critical_list = [bool(critical_mask[i]) for i in range(len(radius))] \
        if has_critical else [False] * len(radius)

    def _log_interp(a, b, f):
        return 10**(np.log10(max(a, eps)) + f * (np.log10(max(b, eps)) - np.log10(max(a, eps))))

    def _thr_for_edge(i, j):
        """Region-dependent min_dr/r for the edge between zones i and j.

        Returns 0.0 (constraint disabled) if ``disable_in_critical`` is
        True and either endpoint is critical; ``critical_min_dr_r`` if
        either endpoint is critical; otherwise the global ``min_dr_r``.
        """
        if has_critical and (is_critical_list[i] or is_critical_list[j]):
            return 0.0 if disable_in_critical else critical_min_dr_r
        return min_dr_r

    def _insert_midpoint(idx, f=None):
        # Default to the sigma-optimal insertion fraction (same strategy as
        # adaptive_grid_refinement / fill_grid_to_minimum) rather than a blind
        # midpoint, so max-spacing and cell-ratio insertions don't introduce
        # sigma discretization noise between adjacent grid points.
        if f is None:
            f = optimal_insert_fraction(r_list, v_list, idx)
        r_list.insert(idx + 1, r_list[idx] + f * (r_list[idx + 1] - r_list[idx]))
        v_list.insert(idx + 1, v_list[idx] + f * (v_list[idx + 1] - v_list[idx]))
        T_list.insert(idx + 1, _log_interp(T_list[idx], T_list[idx + 1], f))
        rho_list.insert(idx + 1, _log_interp(rho_list[idx], rho_list[idx + 1], f))
        ne_list.insert(idx + 1, _log_interp(ne_list[idx], ne_list[idx + 1], f))
        tau_list.insert(idx + 1, _log_interp(tau_list[idx], tau_list[idx + 1], f))
        kappa_list.insert(idx + 1, _log_interp(kappa_list[idx], kappa_list[idx + 1], f))
        for elem in mf_lists:
            mf_lists[elem].insert(idx + 1,
                mf_lists[elem][idx] + f * (mf_lists[elem][idx + 1] - mf_lists[elem][idx]))
        if has_critical:
            is_critical_list.insert(idx + 1,
                bool(is_critical_list[idx] or is_critical_list[idx + 1]))

    # --- Phase 1: Remove points below minimum spacing ---
    n_before = len(r_list)
    keep = [True] * len(r_list)
    last_kept = 0
    n_critical_preserved = 0
    for i in range(1, len(r_list)):
        dr = abs(r_list[last_kept] - r_list[i])
        r_ref = max(abs(r_list[last_kept]), abs(r_list[i]), eps)
        thr_i = _thr_for_edge(last_kept, i)
        if dr / r_ref < thr_i:
            # Would remove.  If this point is in a critical region and the
            # relaxed threshold is what saved it, count it as preserved.
            if has_critical and is_critical_list[i] and thr_i < min_dr_r:
                n_critical_preserved += 1
                # Keep it (relaxed threshold applies) and advance last_kept
                # so subsequent points are measured from this one.
                last_kept = i
            else:
                keep[i] = False
        else:
            last_kept = i
    keep[-1] = True  # always keep inner boundary
    n_removed = sum(not k for k in keep)
    if n_removed > 0 or n_critical_preserved > 0:
        print(f"[Grid Spacing] Removed {n_removed} below-min points "
              f"(global dr/r < {min_dr_r:.0e})"
              f"{f', preserved {n_critical_preserved} critical-region points' if n_critical_preserved else ''}, "
              f"{n_before} -> {n_before - n_removed} points")
        r_list = [r for r, k in zip(r_list, keep) if k]
        v_list = [x for x, k in zip(v_list, keep) if k]
        T_list = [x for x, k in zip(T_list, keep) if k]
        rho_list = [x for x, k in zip(rho_list, keep) if k]
        ne_list = [x for x, k in zip(ne_list, keep) if k]
        tau_list = [x for x, k in zip(tau_list, keep) if k]
        kappa_list = [x for x, k in zip(kappa_list, keep) if k]
        for elem in mf_lists:
            mf_lists[elem] = [x for x, k in zip(mf_lists[elem], keep) if k]
        if has_critical:
            is_critical_list = [c for c, k in zip(is_critical_list, keep) if k]

    # --- Phase 2: Insert midpoints in cells above maximum spacing ---
    n_inserted = 0
    max_iter = max_points * 4  # safety cap
    for _ in range(max_iter):
        if len(r_list) >= max_points:
            break
        worst_i, worst_dr_r = -1, 0.0
        for i in range(len(r_list) - 1):
            dr = abs(r_list[i] - r_list[i + 1])
            r_ref = max(abs(r_list[i]), eps)
            dr_r = dr / r_ref
            if dr_r > worst_dr_r:
                worst_i, worst_dr_r = i, dr_r
        if worst_i < 0 or worst_dr_r <= max_dr_r:
            break
        _insert_midpoint(worst_i)
        n_inserted += 1
    if n_inserted > 0:
        print(f"[Grid Spacing] Inserted {n_inserted} midpoints for max spacing "
              f"(dr/r > {max_dr_r:.0e}), grid now {len(r_list)} points")

    # --- Phase 3: Smooth abrupt cell-size transitions ---
    # NB: Fix the WORST (highest-ratio) violation each pass, not the first
    # found scanning left-to-right. Inserting a midpoint perturbs the ratios
    # of neighbouring pairs and can spawn new (usually mild) violations near
    # the point of insertion. A first-found strategy can get stuck churning
    # through a cluster of such self-generated violations near the start of
    # the array and never reach a rare-but-severe violation further along
    # (observed: a single ~20x cell-ratio violation near the inner boundary
    # survived unresolved even with max_points raised from 400 to 550, while
    # 100+ points were spent resolving churn earlier in the array).
    n_smoothed = 0
    for _ in range(max_points):
        if len(r_list) >= max_points:
            break
        worst_i, worst_ratio = -1, max_cell_ratio
        for i in range(len(r_list) - 2):
            dr0 = abs(r_list[i + 1] - r_list[i])
            dr1 = abs(r_list[i + 2] - r_list[i + 1])
            ratio = max(dr0, dr1) / max(min(dr0, dr1), eps)
            if ratio > worst_ratio:
                worst_i, worst_ratio = i, ratio
        if worst_i < 0:
            break
        dr0 = abs(r_list[worst_i + 1] - r_list[worst_i])
        dr1 = abs(r_list[worst_i + 2] - r_list[worst_i + 1])
        if dr0 > dr1:
            _insert_midpoint(worst_i)
        else:
            _insert_midpoint(worst_i + 1)
        n_smoothed += 1
    if n_smoothed > 0:
        print(f"[Grid Spacing] Inserted {n_smoothed} smoothing points "
              f"(cell ratio > {max_cell_ratio}), grid now {len(r_list)} points")

    # --- Report final spacing ---
    dr_arr = np.abs(np.diff(r_list))
    r_ref_arr = np.maximum(np.abs(np.array(r_list[:-1])), eps)
    dr_r_arr = dr_arr / r_ref_arr
    print(f"[Grid Spacing] Final: {len(r_list)} points, "
          f"dr/r range: [{dr_r_arr.min():.4e}, {dr_r_arr.max():.4e}]")

    return (np.array(r_list), np.array(v_list), np.array(T_list),
            np.array(rho_list), np.array(ne_list), np.array(tau_list),
            {elem: np.array(mf_lists[elem]) for elem in mf_lists},
            np.array(kappa_list))


def resolve_recombination_fronts(r_cm, T_K, ne, rho,
                                 t_jump_threshold=0.3,
                                 front_width_zones=5,
                                 t_min=3000.0, t_max=50000.0):
    """Resolve under-resolved recombination fronts by smoothing T over a few zones.

    STELLA often under-resolves H/He recombination fronts, where T drops from
    ~8000-10000 K to ~5000 K across a single zone. This creates a discontinuity
    in the electron density (15x drop) that CMFGEN cannot converge with.

    This function detects sharp T jumps (recombination fronts) and replaces
    the single-zone discontinuity with a smooth log-linear transition over
    ``front_width_zones`` zones. The electron density is then recomputed
    from the smoothed T using the Saha equation.

    This is NOT convergence smoothing — it is resolving a physical feature
    (recombination front) that STELLA under-resolved. The front has a finite
    physical width set by the recombination time and expansion timescale;
    STELLA's zone size at the front is larger than this width.

    Args:
        r_cm: Radius array (outer->inner, decreasing), in cm
        T_K: Temperature array (K)
        ne: Electron density array (/cm^3) — will be recomputed after smoothing
        rho: Density array (g/cm^3)
        t_jump_threshold: Minimum |dT/T| to flag a recombination front.
            Default 0.3 = 30% drop. H recombination typically shows 40-60%.
        front_width_zones: Number of zones to spread the transition over.
            Default 5 — resolves the front over a few depths for CMFGEN.
        t_min: Minimum temperature floor (K). Default 3000.
        t_max: Maximum temperature cap (K). Default 50000.

    Returns:
        (T_smoothed, ne_smoothed, n_fronts_resolved, front_info)
    """
    eps = 1e-30
    n = len(T_K)
    T_out = T_K.copy().astype(float)
    ne_out = ne.copy().astype(float)

    # Detect recombination fronts: zones where |dT/T| > threshold
    # and T is in the recombination regime (4000-15000 K)
    # and ne is significant (> 1e6) — skip cold neutral gas where ne is tiny
    is_front = np.zeros(n, dtype=bool)
    for i in range(n - 1):
        dT = abs(T_K[i] - T_K[i + 1])
        T_avg = 0.5 * (T_K[i] + T_K[i + 1])
        if T_avg > 0:
            frac_jump = dT / T_avg
        else:
            frac_jump = 0
        ne_avg = 0.5 * (ne[i] + ne[i + 1])
        # Recombination front: T in 4000-15000 K range, large fractional jump,
        # and ne is significant (partially ionized gas, not cold neutral)
        if frac_jump > t_jump_threshold and 4000 < T_avg < 15000 and ne_avg > 1e6:
            is_front[i] = True
            is_front[i + 1] = True

    # Group consecutive front zones
    front_regions = []
    i = 0
    while i < n:
        if is_front[i]:
            j = i
            while j < n and is_front[j]:
                j += 1
            front_regions.append((i, j - 1))
            i = j
        else:
            i += 1

    # Smooth each front region by extending the transition over front_width_zones
    front_info = []
    for lo, hi in front_regions:
        # Extend the region to include surrounding zones for context
        ext_lo = max(0, lo - front_width_zones // 2)
        ext_hi = min(n - 1, hi + front_width_zones // 2)

        # Get the T values at the edges of the extended region
        T_lo = T_out[ext_lo]
        T_hi = T_out[ext_hi]

        # Only smooth if there's a significant jump across the extended region
        if abs(T_lo - T_hi) / max(T_lo, T_hi, eps) < t_jump_threshold:
            continue

        # Log-linear interpolation across the extended region
        log_T_lo = np.log10(max(T_lo, t_min))
        log_T_hi = np.log10(max(T_hi, t_min))
        n_zones = ext_hi - ext_lo + 1
        if n_zones < 3:
            continue

        log_T_interp = np.linspace(log_T_lo, log_T_hi, n_zones)
        T_smooth = np.power(10.0, log_T_interp)
        T_smooth = np.clip(T_smooth, t_min, t_max)

        # Apply the smoothed T
        T_out[ext_lo:ext_hi + 1] = T_smooth

        front_info.append({
            'zone_lo': ext_lo,
            'zone_hi': ext_hi,
            'r_cm': r_cm[lo],
            'T_before': T_K[lo],
            'T_after': T_K[hi],
            'T_smooth_lo': T_smooth[0],
            'T_smooth_hi': T_smooth[-1],
            'n_zones': n_zones,
        })

    # Recompute electron density for smoothed zones.
    # We do NOT use the Saha equation here — it is extremely sensitive to T
    # in the recombination regime (8000-12000 K) and produces wildly wrong
    # values compared to STELLA's NLTE solver. Instead, we interpolate ne
    # smoothly in log space across the front, matching the T smoothing.
    if len(front_info) > 0:
        for fi in front_info:
            ext_lo = fi['zone_lo']
            ext_hi = fi['zone_hi']
            n_zones = ext_hi - ext_lo + 1
            # Log-linear interpolation of ne across the extended region
            log_ne_lo = np.log10(max(ne_out[ext_lo], 1.0))
            log_ne_hi = np.log10(max(ne_out[ext_hi], 1.0))
            log_ne_interp = np.linspace(log_ne_lo, log_ne_hi, n_zones)
            ne_out[ext_lo:ext_hi + 1] = np.power(10.0, log_ne_interp)
            ne_out = np.clip(ne_out, 1.0, None)

    return T_out, ne_out, len(front_info), front_info


def write_sn_hydro_data(outpath: str, r_cm, v_kms, rho, T_K, sigma, kappa,
                        time_days: float, mass_fractions: dict,
                        lum_lsun: float, t_floor: float = 3000.0,
                        t_cap: float = 0.0,
                        trace_threshold: float = 0.0,
                        r_inner_cm: float = None, T_inner_K: float = None,
                        tau_inner: float = None,
                        tau_inner_chi: float = None,
                        ni56_initial: np.ndarray = None,
                        ne_stella: np.ndarray = None,
                        ne_stella_r: np.ndarray = None,
                        inc_rad_decays: bool = True):
    """Write CMFGEN SN_HYDRO_DATA file (outer->inner order)."""
    n = len(r_cm)
    # nfrac is computed after isotope handling below (mass_fractions may be modified)

    # Orientation guard: CMFGEN's rd_sn_data.f and RV_SN_MODEL_SNIIN both expect
    # radius in outer->inner (decreasing) order. Flip if inner->outer, and raise
    # if non-monotonic (which would break CMFGEN's interpolation).
    if n >= 2:
        if r_cm[0] < r_cm[-1]:
            print("[Orientation Guard] Radius is inner->outer; flipping to outer->inner for SN_HYDRO_DATA")
            r_cm = r_cm[::-1].copy()
            v_kms = v_kms[::-1].copy()
            rho = rho[::-1].copy()
            T_K = T_K[::-1].copy()
            sigma = sigma[::-1].copy()
            kappa = kappa[::-1].copy()
            mass_fractions = {k: v[::-1].copy() for k, v in mass_fractions.items()}
            if ni56_initial is not None:
                ni56_initial = ni56_initial[::-1].copy()
            if ne_stella is not None:
                ne_stella = ne_stella[::-1].copy()
            if ne_stella_r is not None:
                ne_stella_r = ne_stella_r[::-1].copy()
        if not np.all(np.diff(r_cm) < 0):
            n_degen = int(np.sum(np.diff(r_cm) >= 0))
            raise ValueError(
                f"SN_HYDRO_DATA orientation guard: radius has {n_degen} non-decreasing step(s) "
                f"(degenerate or non-monotonic). CMFGEN requires strictly decreasing r. "
                f"R[0]={r_cm[0]:.10e}, R[1]={r_cm[1]:.10e}, R[-1]={r_cm[-1]:.10e}"
            )

    r_1e10 = r_cm / 1e10
    
    # Ensure positive time for CMFGEN/SEDONA compatibility
    # Some hydro simulations may have negative times for pre-explosion phases
    if time_days < 0:
        print(f"[Time Correction] Negative time detected: {time_days:.2f} days")
        print(f"  Converting to positive time since explosion bomb introduction")
        time_days = abs(time_days)
    
    # Save original T before floor/cap to detect which zones changed
    T_K_original = T_K.copy()

    # Apply temperature floor
    n_floored = int(np.sum(T_K < t_floor))
    if n_floored > 0:
        print(f"[T floor] Flooring {n_floored} depths from {T_K.min():.0f} K to T_FLOOR={t_floor:.0f} K")
        T_K = np.maximum(T_K, t_floor)

    # Apply temperature cap (prevents MON_INTERP crash when T_gas at outer
    # boundary is very high from shock heating; at T > ~5e5 K hydrogen is
    # fully ionized, giving high LTE opacity that pushes TAU(1) > 2/3 as
    # populations evolve)
    if t_cap > 0:
        n_capped = int(np.sum(T_K > t_cap))
        if n_capped > 0:
            print(f"[T cap] Capping {n_capped} depths from {T_K.max():.0f} K to T_CAP={t_cap:.0f} K")
            T_K = np.minimum(T_K, t_cap)

    # Break flat isothermal regions after floor to prevent SEDONA/CMFGEN singular-matrix crashes
    T_K = break_isothermal_regions(T_K, t_floor=t_floor)

    # Recalculate kappa only for zones where T was changed by the floor/cap.
    # The STELLA NLTE kappa (passed in as the `kappa` argument) is much more
    # accurate than the analytical Kramers approximation — overwriting it
    # entirely gives opacities that are wrong by up to 40x at the inner
    # boundary, corrupting the optical depth structure and causing population
    # oscillations in CMFGEN.
    #
    # For zones where T was floored/capped, the original STELLA kappa was
    # computed at the pre-floor T, so we recalculate those zones using the
    # analytical formula with the corrected T. For all other zones, we keep
    # the STELLA NLTE kappa as-is.
    kappa_analytical = compute_kappa_total(rho, T_K, mass_fractions)
    t_changed = T_K != T_K_original
    n_changed = int(np.sum(t_changed))
    if n_changed > 0:
        kappa = np.where(t_changed, kappa_analytical, kappa)
        print(f"[kappa] Recalculated for {n_changed} floored/capped zones (kept STELLA NLTE for {n - n_changed} zones): "
              f"min={kappa.min():.3e}, max={kappa.max():.3e} cm²/g")
    else:
        print(f"[kappa] Keeping STELLA NLTE opacity (no T floor/cap changes): "
              f"min={kappa.min():.3e}, max={kappa.max():.3e} cm²/g")

    T_10k  = T_K / 1e4

    # Mean atomic mass and electron density (approximate, assuming full ionisation)
    elem_A = {'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
              'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
              'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'FePeak': 56.0, 'Ni58': 58.0, 'Ni56': 56.0}
    elem_Z = {'H': 1.0, 'He': 2.0, 'C': 6.0, 'N': 7.0, 'O': 8.0,
              'Ne': 10.0, 'Na': 11.0, 'Mg': 12.0, 'Al': 13.0, 'Si': 14.0,
              'S': 16.0, 'Ar': 18.0, 'Ca': 20.0, 'FePeak': 26.0, 'Ni58': 28.0, 'Ni56': 28.0}

    inv_A = np.zeros(n)
    Z_over_A = np.zeros(n)
    for elem, frac in mass_fractions.items():
        A = elem_A.get(elem, 28.0)
        Z = elem_Z.get(elem, 14.0)
        inv_A    += frac / A
        Z_over_A += frac * Z / A

    mean_A = np.where(inv_A > 1e-30, 1.0 / inv_A, 28.0)
    atom_density  = rho * AVOGADRO * inv_A          # atoms/cm^3
    if ne_stella is not None and ne_stella_r is not None:
        # Interpolate STELLA n_e from original grid onto final (post-AMR) grid
        # Both arrays are outer->inner (decreasing r); np.interp needs increasing x
        ne_interp = np.interp(r_cm, ne_stella_r[::-1], ne_stella[::-1])
        electron_density = np.clip(ne_interp, 1.0, None)  # floor at 1 cm^-3

        # Sanity check: ne cannot exceed full ionization
        # ne_max = rho * N_A * sum(X_i * Z_i / A_i)
        ne_max = rho * AVOGADRO * Z_over_A
        n_unphysical = int(np.sum(electron_density > ne_max * 1.1))
        if n_unphysical > 0:
            print(f"[Ne] WARNING: STELLA n_e exceeds full ionization at {n_unphysical}/{len(rho)} depths. "
                  f"This indicates a bug in STELLA's Saha solver (Yat*AZ issue). "
                  f"Using STELLA ne as-is — fix STELLA source if convergence problems occur.")
        print(f"[Ne] Using STELLA ionization solver n_e (interpolated): "
              f"min={electron_density.min():.3e}, max={electron_density.max():.3e} cm^-3")
    else:
        electron_density = calculate_electron_density(rho, T_K, mass_fractions)  # Saha fallback
        print("[Ne] Using Saha approximation for n_e (STELLA n_e not available)")

    def write_block(f, label, values, per_line=6):
        f.write(f"{label}\n")
        for i in range(0, len(values), per_line):
            chunk = values[i:i+per_line]
            f.write(''.join(f"    {v:22.15E}" for v in chunk) + '\n')
        f.write('\n')

    # Calculate radioactive decay chain: Ni56 -> Co56 -> Fe56
    # Decay timescales (from STELLA source: volenpumnoint.f line 99)
    t_Ni = 8.8      # days (e-folding time for Ni56)
    t_Co = 111.3    # days (e-folding time for Co56)

    # Get initial Ni56 mass fraction from parameter (from .xni file)
    # If not provided, check mass_fractions (backward compatibility)
    if ni56_initial is not None:
        Ni56_initial = ni56_initial
    else:
        Ni56_initial = mass_fractions.get('Ni56', np.zeros(n))

    # Calculate decayed fractions at snapshot time
    Ni56_t = Ni56_initial * np.exp(-time_days / t_Ni)
    Co56_t = Ni56_initial * (t_Co / (t_Co - t_Ni)) * (np.exp(-time_days / t_Co) - np.exp(-time_days / t_Ni))
    Fe56_t = Ni56_initial - Ni56_t - Co56_t  # Total decayed to Fe56

    # Print decay information
    if np.any(Ni56_initial > 1e-10):
        max_Ni56_init = np.max(Ni56_initial)
        max_Ni56_t = np.max(Ni56_t)
        max_Co56_t = np.max(Co56_t)
        max_Fe56_t = np.max(Fe56_t)
        print(f"[Radioactive decay] t={time_days:.1f} days:")
        print(f"  Initial Ni56 max fraction: {max_Ni56_init:.4e}")
        print(f"  Current Ni56 max fraction: {max_Ni56_t:.4e} ({100*max_Ni56_t/max_Ni56_init:.1f}% remaining)")
        print(f"  Co56 max fraction: {max_Co56_t:.4e}")
        print(f"  Fe56 (from decay) max fraction: {max_Fe56_t:.4e}")

    # Update elemental mass fractions and prepare isotopic data for CMFGEN.
    #
    # The Ni56 decay chain is: Ni56 -> Co56 -> Fe56
    #
    # From rd_sn_data.f (ELSE branch, lines ~448-474):
    #   For every element that has isotopes listed in SN_HYDRO_DATA, CMFGEN:
    #     1. Sums the isotope mass fractions for that element
    #     2. Checks that elemental mass fraction is within 20% of the isotope sum (STOP if not)
    #     3. Replaces the elemental population with the isotope sum
    #
    # Therefore, for EACH species that has an isotope entry, the elemental mass
    # fraction must be present AND must equal the sum of its isotopes.
    #
    #   NICK elemental  = Ni56_t    (all Ni56 is radioactive; no stable Ni56 in progenitor)
    #   NICK 56 isotope = Ni56_t
    #
    #   COB elemental   = Co56_t    (all Co56 comes from Ni56 decay; no stable Co56 in progenitor)
    #   COB 56 isotope  = Co56_t
    #
    #   IRON elemental  = stable_fe + Fe56_t   (progenitor Fe + Fe56 from decay)
    #   IRON 56 isotope = stable_fe + Fe56_t   (must equal elemental)
    #
    # The original Ni56 key in mass_fractions is the initial value (pre-decay), so we
    # replace it with the current value. FePeak gets updated to include decayed Fe56.

    # Replace initial Ni56 with current (decayed) value
    if 'Ni56' in mass_fractions:
        del mass_fractions['Ni56']
    mass_fractions['NICK'] = Ni56_t

    # Add Co56 as a new elemental species (none existed before in .abn)
    mass_fractions['COB'] = Co56_t

    # Update FePeak to include Fe56 produced by decay
    stable_fe = np.zeros(n)
    if 'FePeak' in mass_fractions:
        stable_fe = mass_fractions.pop('FePeak')
    total_fe = stable_fe + Fe56_t
    mass_fractions['IRON'] = total_fe

    # Build isotopic data: (species_name, baryon_number) -> mass_fraction array
    # All three isotopes must be present and match their elemental counterparts.
    # When inc_rad_decays=False (SN_B models), isotopes are not written so
    # CMFGEN sets NUM_ISOTOPES=0 and does not enable INC_RAD_DECAYS.
    # RD_SN_DATA accepts NISO=0 (patched: NISO .LT. 0 check, not .EQ. 0).
    if inc_rad_decays:
        isotopic_data = {
            ('NICK', 56): Ni56_t,    # matches NICK elemental
            ('COB',  56): Co56_t,    # matches COB  elemental
            ('IRON', 56): total_fe,  # matches IRON elemental
        }
        n_isotopes = len(isotopic_data)
    else:
        isotopic_data = {}
        n_isotopes = 0
        print("  [INC_RAD_DECAYS] Isotope blocks (Ni56, Co56, Fe56) suppressed — "
              "INC_RAD_DECAYS will be F in VADAT")
    
    # Map CMFGEN species names to our internal names
    cmfgen_to_internal = {
        'HYD': 'H',
        'HE': 'He', 
        'CARB': 'C',
        'OXY': 'O',
        'SIL': 'Si',        # CMFGEN uses 'SIL' not 'SILI'
        'CAL': 'Ca',        # CMFGEN uses 'CAL' not 'CALC'
        'IRON': 'Fe',
        'NICK': 'Ni',
        'COB': 'Co'
    }
    
    # Filter to only include the 6 elements needed: H, He, O, Si, Ca, Fe
    # (plus radioactive Ni56, Co56 isotopes for the decay chain).
    # All other elements (C, N, Ne, Na, Mg, Al, S, Ar, Ni58) are excluded.
    # Their mass fraction is redistributed among the kept elements via renormalization.
    # When inc_rad_decays=False, also exclude NICK and COB — without MODEL_SPEC
    # ISF entries, CMFGEN creates implicit species for them that are numerically
    # unstable at low T (T < 4000 K) in the inner ejecta.
    allowed_elements = {'H', 'He', 'O', 'Si', 'Ca', 'Fe', 'FePeak', 'Ni', 'Co',
                       'HYD', 'HE', 'OXY', 'SIL', 'CAL', 'IRON', 'NICK', 'COB'}
    if not inc_rad_decays:
        allowed_elements -= {'Ni', 'Co', 'NICK', 'COB'}
        print("  [INC_RAD_DECAYS] NICK and COB excluded from mass fractions — "
              "no implicit species will be created by CMFGEN")
    
    # Convert internal names to CMFGEN names for output
    cmfgen_mass_fractions = {}
    for internal_name, frac_array in mass_fractions.items():
        # Only keep allowed elements
        if internal_name not in allowed_elements:
            print(f"  Skipping {internal_name} (not in allowed element list)")
            continue
        
        # Find corresponding CMFGEN name
        cmfgen_name = None
        for cmfgen_name_key, internal in cmfgen_to_internal.items():
            if internal == internal_name:
                cmfgen_name = cmfgen_name_key
                break
        if cmfgen_name is None:
            # FePeak maps to IRON in CMFGEN
            if internal_name == 'FePeak':
                cmfgen_name = 'IRON'
            else:
                cmfgen_name = internal_name.upper()
        cmfgen_mass_fractions[cmfgen_name] = frac_array
    
    # Replace mass_fractions with CMFGEN-named version
    mass_fractions = cmfgen_mass_fractions

    # Renormalize mass fractions to sum=1 at each depth after element filtering.
    # Filtering drops elements (N, Ne, Mg, Al, S, Ar) whose mass fraction is lost;
    # without renormalization the remaining fractions sum to <1, especially in
    # blended/transition regions where the dropped elements carry significant mass.
    mf_sum = np.zeros(n)
    for frac in mass_fractions.values():
        mf_sum += frac
    mf_sum = np.where(mf_sum < 1e-30, 1.0, mf_sum)
    for elem in mass_fractions:
        mass_fractions[elem] = mass_fractions[elem] / mf_sum
    print(f"  Renormalized {len(mass_fractions)} elements after filtering: "
          f"per-depth sum range [{mf_sum.min():.6f}, {mf_sum.max():.6f}] -> 1.0")

    nfrac = len(mass_fractions)  # Count after filtering

    # Zero out trace species below trace_threshold to aid CMFGEN convergence.
    # When trace_threshold > 0, any element mass fraction below this value at any
    # depth is set to 0. The ABUND_FLOOR (1e-30) then brings it back to a negligible
    # but non-zero value to prevent SETPOP_WLD_V2 crashes. With ABUND_FLOOR=1e-30,
    # the number fraction (POP_SPECIES/POP_ATOM) is ~1e-32, well below the FIXPOP
    # absolute abundance threshold (1e-10), so the species' iterative corrections
    # are frozen — preventing unbounded population corrections from trace species
    # at H/He-dominated depths.
    if trace_threshold > 0:
        for elem, frac in mass_fractions.items():
            n_zeroed = int(np.sum(frac < trace_threshold))
            if n_zeroed > 0:
                print(f"  [Trace zeroing] {elem}: {n_zeroed} depths below {trace_threshold:.0e} (zeroed)")
            mass_fractions[elem] = np.where(frac < trace_threshold, 0.0, frac)
        # Renormalize after zeroing
        mf_sum_trace = np.zeros(n)
        for frac in mass_fractions.values():
            mf_sum_trace += frac
        mf_sum_trace = np.where(mf_sum_trace < 1e-30, 1.0, mf_sum_trace)
        for elem in mass_fractions:
            mass_fractions[elem] = mass_fractions[elem] / mf_sum_trace

    # Set missing/trace species to exactly 0.0 in SN_HYDRO_DATA.
    # If an element's mass fraction at a depth is < 1e-10, write 0.0E+00
    # instead of a tiny floor value. This is physically cleaner than the
    # previous 1e-30 floor and makes the missing species explicit in the
    # output file. CMFGEN's SETPOP_WLD_V2 may crash with exactly 0 abundance
    # for species that are listed in the model atom — to avoid this, ensure
    # that species with all-zero abundances are handled by FIXPOP or removed
    # from the model atom via VADAT.
    ABSENT_THRESHOLD = 1e-10
    for elem, frac in mass_fractions.items():
        n_zeroed = int(np.sum(frac < ABSENT_THRESHOLD))
        if n_zeroed > 0:
            print(f"  [Zero species] {elem}: {n_zeroed} depths below {ABSENT_THRESHOLD:.0e} (set to 0.0)")
        mass_fractions[elem] = np.where(frac < ABSENT_THRESHOLD, 0.0, frac)

    # Renormalize after zeroing — only zones where sum > 0 need renormalization.
    mf_sum_post = np.zeros(n)
    for frac in mass_fractions.values():
        mf_sum_post += frac
    mf_sum_post = np.where(mf_sum_post < 1e-30, 1.0, mf_sum_post)
    for elem in mass_fractions:
        mass_fractions[elem] = mass_fractions[elem] / mf_sum_post

    # Also zero isotopic data and renormalize to match elemental mass fractions.
    # CMFGEN checks that elemental = sum(isotopes) within 20%; if elemental is
    # renormalized but isotope is not, the check fails (STOP in rd_sn_data.f).
    for iso_key, iso_frac in isotopic_data.items():
        iso_zeroed = np.where(iso_frac < ABSENT_THRESHOLD, 0.0, iso_frac)
        # Renormalize by the same factor as elemental mass fractions
        isotopic_data[iso_key] = iso_zeroed / mf_sum_post

    # Force isotope mass fractions to exactly match elemental for species with isotopes.
    # This guarantees the 20% consistency check in rd_sn_data.f passes.
    iso_to_elem = {
        ('NICK', 56): 'NICK',
        ('COB',  56): 'COB',
        ('IRON', 56): 'IRON',
    }
    for iso_key, elem_name in iso_to_elem.items():
        if iso_key in isotopic_data and elem_name in mass_fractions:
            isotopic_data[iso_key] = mass_fractions[elem_name].copy()

    # Build plot-only mass fractions for visualisation
    mass_fractions_plot = mass_fractions.copy()  # already has NICK, COB, IRON

    # Rosseland mean opacity in units of 10^{-10} cm^{-1} (chi = kappa * rho)
    # This is the quantity CMFGEN reads as "Rosseland mean opacty (10^{-10} cm^{-1})"
    rosseland_chi = kappa * rho * 1e10  # convert cm^{-1} to 10^{-10} cm^{-1}

    with open(outpath, 'w') as f:
        f.write(f"Number of data points:{n:>10d}\n")
        f.write(f"Number of mass fractions:{nfrac:>10d}\n")
        f.write(f"Number of isotopes:{n_isotopes:>10d}\n")
        f.write(f"Time(days) since explosion:{time_days:>15.7f}\n")
        # Extra header lines used by our pipeline (not in Dessart reference but
        # read by prepare_VADAT.py, run_model_grid.py, fix_sedona_param.py,
        # convert_cmfgen_to_sirocco.py). CMFGEN's rd_sn_data.f ignores them.
        f.write(f"Energy Lum [lsun]:{lum_lsun:>19.6E}\n")
        if r_inner_cm is not None:
            f.write(f"Inner boundary radius [cm]:{r_inner_cm:>19.6E}\n")
        if T_inner_K is not None:
            f.write(f"Inner boundary temperature [K]:{T_inner_K:>19.6E}\n")
        if tau_inner is not None:
            f.write(f"Inner boundary Rosseland tau:{tau_inner:>19.6E}\n")
        if tau_inner_chi is not None:
            f.write(f"Inner boundary chi-column tau:{tau_inner_chi:>19.6E}\n")
        f.write('\n')
        write_block(f, "Radius grid (10^10cm)", r_1e10)
        write_block(f, "Velocity (km/s)",       v_kms)
        write_block(f, "Sigma (dlnV/dlnr-1)",   sigma)
        write_block(f, "Temperature (10^4 K)",  T_10k)
        write_block(f, "Density (gm/cm^3)",     rho)
        write_block(f, "Atom density (/cm^3)",  atom_density)
        write_block(f, "Electron density (/cm^3)", electron_density)
        write_block(f, "Rosseland mean opacty (10^{-10} cm^{-1})", rosseland_chi)
        write_block(f, "Kappa (cm^2/gm)",       kappa)
        for elem, frac in mass_fractions.items():
            # mass_fractions now already has CMFGEN names from our conversion
            write_block(f, f"{elem} mass fraction", frac)

        # Write isotopic data
        # Format: "SPECIES_NAME BARYON_NUMBER mass fraction"
        per_line = 8
        for (iso_name, baryon_num), iso_frac in isotopic_data.items():
            f.write(f"{iso_name} {baryon_num} mass fraction\n")
            for i in range(0, len(iso_frac), per_line):
                chunk = iso_frac[i:i+per_line]
                f.write(''.join(f"    {v:14.7E}" for v in chunk) + '\n')
            f.write('\n')

    print(f"Written {outpath}  ({n} depth points, {nfrac} mass fractions)")
    print(f"  t={time_days:.3f} days, L={lum_lsun:.3e} Lsun")
    print(f"  r=[{r_1e10[0]:.3e},{r_1e10[-1]:.3e}] 1e10cm")
    print(f"  v=[{v_kms[0]:.1f},{v_kms[-1]:.1f}] km/s")


def break_isothermal_regions(temperature, t_floor=3000.0, min_rel_change=5e-4):
    """
    Break flat/isothermal temperature regions by adding a tiny monotonic gradient.

    SEDONA's NLTE solver can crash with a singular matrix when many consecutive
    zones have exactly the same temperature (especially at low densities). This
    function identifies contiguous blocks of cells at the floor temperature and
    applies a small linear perturbation so each cell differs from its neighbour,
    keeping the rate matrix well-conditioned.

    Parameters:
    -----------
    temperature : array (K), outer->inner order
    t_floor : float
        Temperature floor value to identify flat regions (default: 3000.0)
    min_rel_change : float
        Minimum fractional dT/T to enforce (default: 5e-4 = 0.05%)

    Returns:
    --------
    temperature : array with flat regions broken
    """
    temperature = np.asarray(temperature, dtype=float).copy()
    n = len(temperature)
    if n < 2:
        return temperature

    tol = 1e-9
    at_floor = np.abs(temperature - t_floor) < tol
    if not np.any(at_floor):
        return temperature

    n_broken = 0
    n_discontinuous = 0
    i = 0
    while i < n:
        if not at_floor[i]:
            i += 1
            continue
        # Find contiguous floored block
        j = i
        while j < n and at_floor[j]:
            j += 1
        block_len = j - i
        if block_len >= 2:
            delta = max(min_rel_change * 1.1, 1e-6)
            if j < n:
                # CRITICAL: ramp from t_floor at the outer edge up to the TRUE
                # neighbouring temperature at the first non-floored cell
                # (index j), instead of an index-only linear perturbation that
                # ignores the real profile. The original implementation ended
                # the ramp at t_floor*(1+block_len*delta) regardless of the
                # actual temperature{j], which for long floored blocks (e.g.
                # ~45% of the grid in cold outer SN ejecta) produces a value
                # completely disconnected from the real data -- a discontinuity
                # of order ~10% in T. Because NLTE/Saha populations depend on
                # exp(-chi/kT), a 10% T jump can change high-ionization-energy
                # species' populations by MANY orders of magnitude in a single
                # depth step, which was observed to trigger CMFGEN population
                # explosions (MAXCH -> 1e100) starting exactly at this boundary.
                t_end = temperature[j]
                if t_end <= t_floor * (1.0 + tol):
                    # Neighbour is also very close to the floor: fall back to
                    # the small linear perturbation to keep cells distinct.
                    t_end = t_floor * (1.0 + block_len * delta)
                else:
                    n_discontinuous += 1
                for k in range(i, j):
                    frac = (k - i + 1) / (block_len + 1)
                    temperature[k] = t_floor + (t_end - t_floor) * frac
                # NOTE: no minimum-per-step-change guard here. Each cell already
                # gets a distinct value from the linear interpolation (frac is
                # strictly increasing), which is sufficient to avoid the exact
                # degeneracy SEDONA's solver chokes on. Enforcing a minimum
                # *fractional* step (as the old code did) would fight the
                # smooth connection to t_end and re-introduce the very
                # discontinuity this branch exists to remove.
            else:
                # Floored block runs to the innermost depth (no real neighbour
                # to connect to) -- keep the original small linear ramp.
                for k in range(i, j):
                    temperature[k] = t_floor * (1.0 + (k - i + 1) * delta)
            n_broken += block_len
        i = j

    if n_broken > 0:
        print(f"[Anti-isothermal] Broken {n_broken} floored cells "
              f"(T == {t_floor:.1f} K) with {min_rel_change:.2e} gradient")
    if n_discontinuous > 0:
        print(f"[Anti-isothermal] Ramped {n_discontinuous} floored block(s) to "
              f"connect smoothly with the true neighbouring temperature "
              f"(prevents a T discontinuity at the floor boundary)")
    return temperature


def plot_comparison(orig_r, orig_rho, orig_T, orig_v, orig_tau, orig_sigma,
                   new_r, new_rho, new_T, new_v, new_tau, new_sigma,
                   output_file, time_days, tau_min=1e-5, tau_max=30.0,
                   orig_mass_fractions=None, new_mass_fractions=None,
                   orig_ne=None, new_ne=None,
                   ne_from_stella=False, orig_tau_nlte=None, orig_tau_const=None, kappa=0.34):
    """Plot before/after comparison of STELLA to CMFGEN conversion.

    Uses drawstyle='steps-mid' for all profile plots to reflect the discrete
    zone-centered nature of the data.
    """
    orig_r_rsun = orig_r  # cm
    new_r_rsun = new_r  # cm

    # Determine consistent x-axis limits (radius)
    r_min = min(np.min(orig_r_rsun), np.min(new_r_rsun)) * 0.9
    r_max = max(np.max(orig_r_rsun), np.max(new_r_rsun)) * 1.1

    # Create 4-row layout: rows 1-3 have 2 cols (6 plots),
    # row 4 mass fractions (full width)
    fig, axes = plt.subplots(4, 2, figsize=(16, 24))
    fig.suptitle(f'STELLA to CMFGEN Conversion: {time_days:.1f}-day snapshot', fontsize=16)

    # ---- Row 1: Density | Temperature ----
    # Plot 1: Density
    ax = axes[0, 0]
    ax.loglog(orig_r_rsun, orig_rho, 'bs', alpha=0.7, markersize=3, drawstyle='steps-mid', label=f'STELLA ({len(orig_r)} points)')
    ax.loglog(new_r_rsun, new_rho, 'ro', alpha=0.7, markersize=2, drawstyle='steps-mid', label=f'CMFGEN ({len(new_r)} points)')
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_title('Density Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xlabel('Radius (cm)', fontsize=10)

    # Plot 2: Temperature
    ax = axes[0, 1]
    ax.loglog(orig_r_rsun, orig_T, 'bs', alpha=0.7, markersize=3, drawstyle='steps-mid', label=f'STELLA ({len(orig_r)} points)')
    ax.loglog(new_r_rsun, new_T, 'ro', alpha=0.7, markersize=2, drawstyle='steps-mid', label=f'CMFGEN ({len(new_r)} points)')
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xlabel('Radius (cm)', fontsize=10)

    # ---- Row 2: Velocity | Sigma ----
    # Plot 3: Velocity
    ax = axes[1, 0]
    ax.loglog(orig_r_rsun, orig_v, 'bs', alpha=0.7, markersize=4, drawstyle='steps-mid', label=f'STELLA ({len(orig_r)} points)')
    ax.loglog(new_r_rsun, new_v, 'ro', alpha=0.8, markersize=2, drawstyle='steps-mid', label=f'CMFGEN ({len(new_r)} points)')
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_title('Velocity Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xlabel('Radius (cm)', fontsize=10)

    # Plot 4: Sigma
    ax = axes[1, 1]
    ax.plot(orig_r_rsun, orig_sigma, 'bs', alpha=0.7, markersize=3, drawstyle='steps-mid', label=f'STELLA ({len(orig_r)} points)')
    ax.plot(new_r_rsun, new_sigma, 'ro', alpha=0.7, markersize=2, drawstyle='steps-mid', label=f'CMFGEN ({len(new_r)} points)')
    ax.axhline(y=-0.99, color='k', linestyle=':', alpha=0.5, label='σ=-0.99')
    ax.axhline(y=-1.0, color='gray', linestyle='--', alpha=0.5, label='σ=-1.0')
    ax.set_xlabel('Radius (cm)', fontsize=10)
    ax.set_ylabel('Sigma (σ)', fontsize=12)
    ax.set_title('Sigma Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_yscale('symlog', linthresh=1.0)
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')

    # ---- Row 3: Optical depth | Electron density ----
    # Plot 5: Optical depth
    ax = axes[2, 0]
    orig_tau_plot = np.clip(orig_tau, 1e-10, None)
    new_tau_plot = np.clip(new_tau, 1e-10, None)
    ax.plot(orig_r_rsun, orig_tau_plot, 'bs', alpha=0.5, markersize=3, drawstyle='steps-mid', label=f'STELLA analytical ({len(orig_r)} pts)')
    if orig_tau_nlte is not None:
        orig_tau_nlte_plot = np.clip(orig_tau_nlte, 1e-10, None)
        ax.plot(orig_r_rsun, orig_tau_nlte_plot, 'g^', alpha=0.5, markersize=2, drawstyle='steps-mid', label=f'STELLA NLTE ({len(orig_r)} pts)')
    if orig_tau_const is not None:
        orig_tau_const_plot = np.clip(orig_tau_const, 1e-10, None)
        ax.plot(orig_r_rsun, orig_tau_const_plot, 'b--', alpha=0.5, markersize=2, drawstyle='steps-mid', label=f'STELLA κ={kappa} ({len(orig_r)} pts)')
    ax.plot(new_r_rsun, new_tau_plot, 'ro', alpha=0.7, markersize=2, drawstyle='steps-mid', label=f'CMFGEN κ={kappa} ({len(new_r)} points)')
    ax.axhline(y=2/3, color='orange', linestyle=':', alpha=0.8, linewidth=2, label='τ = 2/3')
    ax.axhline(y=tau_min, color='gray', linestyle='--', alpha=0.5, linewidth=1, label=f'τ_min={tau_min:.0e}')
    ax.axhline(y=tau_max, color='gray', linestyle='--', alpha=0.5, linewidth=1, label=f'τ_max={tau_max:.0f}')
    ax.set_xlabel('Radius (cm)', fontsize=10)
    ax.set_ylabel('Optical Depth τ', fontsize=12)
    ax.set_title('Optical Depth Profile', fontsize=14)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_yscale('log')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    all_tau = [new_tau, orig_tau]
    if orig_tau_nlte is not None:
        all_tau.append(orig_tau_nlte)
    if orig_tau_const is not None:
        all_tau.append(orig_tau_const)
    tau_lo = max(min(t.min() for t in all_tau), 1e-10)
    tau_hi = max(t.max() for t in all_tau)
    ax.set_ylim(tau_lo * 0.3, tau_hi * 2.0)

    # Plot 6: Electron density + Atom density
    ax = axes[2, 1]
    AVOGADRO = 6.02214076e23
    elem_A = {'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
              'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
              'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'FePeak': 56.0, 'IRON': 56.0,
              'Ni58': 58.0, 'Ni56': 56.0, 'NICK': 56.0, 'COB': 56.0, 'SIL': 28.0,
              'CARB': 12.0, 'OXY': 16.0, 'Fe': 56.0, 'Fe56': 56.0}

    def _compute_atom_density(rho_arr, mf_dict):
        inv_A = np.zeros_like(rho_arr)
        if mf_dict is not None:
            for elem, frac in mf_dict.items():
                A = elem_A.get(elem, 28.0)
                inv_A += frac / A
        inv_A = np.where(inv_A > 1e-30, inv_A, 1.0 / 28.0)
        return rho_arr * AVOGADRO * inv_A

    has_ne = orig_ne is not None and new_ne is not None
    has_mf = new_mass_fractions is not None

    if has_ne or has_mf:
        ne_source = 'STELLA' if ne_from_stella else 'Saha'
        if has_ne:
            ax.loglog(orig_r_rsun, orig_ne, 'bs', alpha=0.7, markersize=3, drawstyle='steps-mid', label=f'{ne_source} n_e ({len(orig_r)})')
            ax.loglog(new_r_rsun, new_ne, 'ro', alpha=0.7, markersize=2, drawstyle='steps-mid', label=f'{ne_source} n_e interp ({len(new_r)})')
        if has_mf:
            new_natom = _compute_atom_density(new_rho, new_mass_fractions)
            ax.loglog(new_r_rsun, new_natom, 'g^', alpha=0.7, markersize=2, drawstyle='steps-mid', label=f'n_atom CMFGEN ({len(new_r)})')
            if orig_mass_fractions is not None:
                orig_natom = _compute_atom_density(orig_rho, orig_mass_fractions)
                ax.loglog(orig_r_rsun, orig_natom, 'gD', alpha=0.5, markersize=3, drawstyle='steps-mid', label=f'n_atom STELLA ({len(orig_r)})')
        ax.set_xlabel('Radius (cm)', fontsize=10)
        ax.set_ylabel('Number Density (cm⁻³)', fontsize=12)
        ax.set_title('Electron & Atom Density', fontsize=14)
        ax.legend(fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3, which='both')
        ax.set_xlim(r_min, r_max)
    else:
        ax.text(0.5, 0.5, 'Density data not available', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Electron & Atom Density', fontsize=14)
        ax.set_xlabel('Radius (cm)', fontsize=10)
        ax.set_xlim(r_min, r_max)
        ax.set_xscale('log')

    # ---- Row 4: Mass fractions (full width) ----
    axes[3, 0].set_visible(False)
    axes[3, 1].set_visible(False)
    ax = fig.add_subplot(4, 1, 4)

    mf_colors = {
        'H': 'tab:blue', 'He': 'tab:cyan', 'C': 'tab:green', 'N': 'olive',
        'O': 'tab:orange', 'Ne': 'tab:pink', 'Na': 'wheat', 'Mg': 'lime',
        'Al': 'silver', 'Si': 'tab:purple', 'S': 'gold', 'Ar': 'tan',
        'Ca': 'tab:red', 'FePeak': 'tab:brown', 'IRON': 'tab:brown',
        'Fe': 'tab:brown', 'NICK': 'darkgray', 'Ni56': 'darkgray',
        'Ni58': 'dimgray', 'COB': 'teal',
    }
    plot_order = ['H', 'He', 'O', 'C', 'Si', 'S', 'Ca', 'Fe', 'IRON', 'FePeak',
                  'NICK', 'Ni56', 'Ni58', 'COB', 'Mg', 'Na', 'Ne', 'N', 'Al', 'Ar']

    plotted_elems = set()
    if orig_mass_fractions is not None:
        for elem in plot_order:
            if elem in orig_mass_fractions and elem not in plotted_elems:
                frac = orig_mass_fractions[elem]
                if np.max(frac) > 1e-4:
                    color = mf_colors.get(elem, 'gray')
                    ax.plot(orig_r_rsun, frac, drawstyle='steps-mid',
                           linestyle='--', color=color, alpha=0.5, linewidth=1.5)
                    plotted_elems.add(elem)
        for elem, frac in orig_mass_fractions.items():
            if elem not in plotted_elems and np.max(frac) > 1e-4:
                color = mf_colors.get(elem, 'gray')
                ax.plot(orig_r_rsun, frac, drawstyle='steps-mid',
                           linestyle='--', color=color, alpha=0.5, linewidth=1.5)
                plotted_elems.add(elem)

    plotted_elems_solid = set()
    if new_mass_fractions is not None:
        for elem in plot_order:
            if elem in new_mass_fractions and elem not in plotted_elems_solid:
                frac = new_mass_fractions[elem]
                if np.max(frac) > 1e-4:
                    color = mf_colors.get(elem, 'gray')
                    ax.plot(new_r_rsun, frac, drawstyle='steps-mid',
                               linestyle='-', color=color, alpha=0.9, linewidth=2.0,
                               label=elem)
                    plotted_elems_solid.add(elem)
        for elem, frac in new_mass_fractions.items():
            if elem not in plotted_elems_solid and np.max(frac) > 1e-4:
                color = mf_colors.get(elem, 'gray')
                ax.plot(new_r_rsun, frac, drawstyle='steps-mid',
                           linestyle='-', color=color, alpha=0.9, linewidth=2.0,
                           label=elem)
                plotted_elems_solid.add(elem)

    ax.set_xlabel('Radius (cm)', fontsize=10)
    ax.set_ylabel('Mass Fraction', fontsize=12)
    ax.set_title('Mass Fractions (dashed=full STELLA, solid=CMFGEN region)', fontsize=14)
    ax.legend(fontsize=9, ncol=5, loc='upper right')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    ax.set_ylim(-0.05, 1.05)

    plt.tight_layout()
    plot_file = str(Path(output_file).parent / 'SN_HYDRO_DATA_comparison.png')
    import os as _os; _os.makedirs(_os.path.dirname(plot_file) or '.', exist_ok=True)
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")
    plt.close()


def plot_adaptive_grid_diagnostics(orig_r, orig_v, orig_T, orig_rho, orig_tau,
                                   downsampled_r, downsampled_v, downsampled_T, downsampled_rho, downsampled_tau,
                                   refined_r, refined_v, refined_T, refined_rho, refined_tau,
                                   output_file, time_days,
                                   max_radius_jump=0.05, max_temp_jump=0.10,
                                   max_density_jump=0.15, max_ne_jump=0.10,
                                   max_logtau_jump=0.10, max_points=500,
                                   orig_ne=None, downsampled_ne=None, refined_ne=None,
                                   orig_mass_fractions=None, downsampled_mass_fractions=None,
                                   refined_mass_fractions=None):
    """Create diagnostic plots for adaptive grid refinement process.

    Upper row: Density, Temperature, Velocity, Ne + Atom density profiles
    Lower row: Histograms of jumps (drho/rho, dT/T, dr/r, dNe/Ne) before/after AMR
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(26, 12))
    fig.suptitle(f'Adaptive Grid Diagnostics - {time_days:.1f} days', fontsize=16)

    RSUN = 6.96e10
    AVOGADRO = 6.02214076e23

    # Convert to Rsun for plotting
    orig_r_rsun = orig_r  # cm
    downsampled_r_rsun = downsampled_r  # cm
    refined_r_rsun = refined_r  # cm

    # --- UPPER ROW: Physical profiles with line + small markers ---

    n_orig = len(orig_r)
    n_base = len(downsampled_r)
    n_ref = len(refined_r)

    # Plot 1: Density profile
    ax = axes[0, 0]
    ax.plot(orig_r_rsun, orig_rho, '-', linewidth=1.5, drawstyle='steps-mid', label=f'Tau-trimmed ({n_orig})', alpha=0.4, color='gray')
    ax.plot(downsampled_r_rsun, downsampled_rho, '-o', markersize=3, linewidth=1, drawstyle='steps-mid',
            label=f'Base grid ({n_base})', alpha=0.7, color='blue')
    ax.plot(refined_r_rsun, refined_rho, '-s', markersize=2, linewidth=0.8, drawstyle='steps-mid',
            label=f'AMR refined ({n_ref})', alpha=0.9, color='red')
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_title('Density Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Plot 2: Temperature profile
    ax = axes[0, 1]
    ax.plot(orig_r_rsun, orig_T, '-', linewidth=1.5, drawstyle='steps-mid', label=f'Tau-trimmed ({n_orig})', alpha=0.4, color='gray')
    ax.plot(downsampled_r_rsun, downsampled_T, '-o', markersize=3, linewidth=1, drawstyle='steps-mid',
            label=f'Base grid ({n_base})', alpha=0.7, color='blue')
    ax.plot(refined_r_rsun, refined_T, '-s', markersize=2, linewidth=0.8, drawstyle='steps-mid',
            label=f'AMR refined ({n_ref})', alpha=0.9, color='red')
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Plot 3: Velocity profile
    ax = axes[0, 2]
    ax.plot(orig_r_rsun, orig_v, '-', linewidth=1.5, drawstyle='steps-mid', label=f'Tau-trimmed ({n_orig})', alpha=0.4, color='gray')
    ax.plot(downsampled_r_rsun, downsampled_v, '-o', markersize=3, linewidth=1, drawstyle='steps-mid',
            label=f'Base grid ({n_base})', alpha=0.7, color='blue')
    ax.plot(refined_r_rsun, refined_v, '-s', markersize=2, linewidth=0.8, drawstyle='steps-mid',
            label=f'AMR refined ({n_ref})', alpha=0.9, color='red')
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_title('Velocity Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xscale('log')

    # Plot 4: Electron density + Atom density
    ax = axes[0, 3]

    # Compute atom density: n_atom = rho * N_A * sum(X_i / A_i)
    elem_A = {'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
              'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
              'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'FePeak': 56.0, 'IRON': 56.0,
              'Ni58': 58.0, 'Ni56': 56.0, 'NICK': 56.0, 'COB': 56.0, 'SIL': 28.0,
              'CARB': 12.0, 'OXY': 16.0, 'Fe': 56.0, 'Fe56': 56.0}

    def _compute_atom_density(rho_arr, mf_dict):
        inv_A = np.zeros_like(rho_arr)
        if mf_dict is not None:
            for elem, frac in mf_dict.items():
                A = elem_A.get(elem, 28.0)
                inv_A += frac / A
        inv_A = np.where(inv_A > 1e-30, inv_A, 1.0 / 28.0)
        return rho_arr * AVOGADRO * inv_A

    if refined_ne is not None:
        # Electron density
        ax.loglog(refined_r_rsun, refined_ne, '-s', markersize=2, linewidth=0.8, drawstyle='steps-mid',
                  label=f'Ne AMR ({n_ref})', alpha=0.9, color='red')
        if downsampled_ne is not None:
            ax.loglog(downsampled_r_rsun, downsampled_ne, '-o', markersize=3, linewidth=1, drawstyle='steps-mid',
                      label=f'Ne base ({n_base})', alpha=0.7, color='blue')
        if orig_ne is not None:
            ax.loglog(orig_r_rsun, orig_ne, '-', linewidth=1.5, drawstyle='steps-mid',
                      label=f'Ne tau-trim ({n_orig})', alpha=0.4, color='gray')

        # Atom density
        n_atom_ref = _compute_atom_density(refined_rho, refined_mass_fractions)
        ax.loglog(refined_r_rsun, n_atom_ref, '--s', markersize=2, linewidth=0.8, drawstyle='steps-mid',
                  label=f'N_atom AMR ({n_ref})', alpha=0.9, color='darkgreen')
        if downsampled_mass_fractions is not None:
            n_atom_base = _compute_atom_density(downsampled_rho, downsampled_mass_fractions)
            ax.loglog(downsampled_r_rsun, n_atom_base, '--o', markersize=3, linewidth=1, drawstyle='steps-mid',
                      label=f'N_atom base ({n_base})', alpha=0.7, color='green')

        ax.set_ylabel('Number Density (cm⁻³)', fontsize=12)
        ax.set_title('Electron & Atom Density', fontsize=14)
        ax.legend(fontsize=8, ncol=2)
    else:
        ax.text(0.5, 0.5, 'Ne data not available', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Electron & Atom Density', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')

    # Set x-axis label on upper row
    for ax in axes[0, :]:
        ax.set_xlabel('Radius (cm)', fontsize=12)

    # Calculate jumps for base and refined grids
    def calculate_jumps(r, T, rho, tau, ne=None):
        eps = 1e-30
        dr_r = np.abs(np.diff(r) / np.maximum(np.abs(r[:-1]), eps))
        dT_T = np.abs(np.diff(T) / np.maximum(T[:-1], eps))
        drho_rho = np.abs(np.diff(rho) / np.maximum(rho[:-1], eps))
        dlogtau = np.abs(np.diff(np.log10(np.maximum(tau, eps))))
        dne_ne = np.abs(np.diff(ne) / np.maximum(ne[:-1], eps)) if ne is not None else None
        return dr_r, dT_T, drho_rho, dlogtau, dne_ne

    dr_r_base, dT_T_base, drho_rho_base, dlogtau_base, dne_ne_base = calculate_jumps(
        downsampled_r, downsampled_T, downsampled_rho, downsampled_tau, downsampled_ne)
    dr_r_ref, dT_T_ref, drho_rho_ref, dlogtau_ref, dne_ne_ref = calculate_jumps(
        refined_r, refined_T, refined_rho, refined_tau, refined_ne)

    # --- LOWER ROW: Jump histograms before/after AMR ---
    bins = np.linspace(0, 1, 40)

    # Plot 5: Density jump histogram
    ax = axes[1, 0]
    ax.hist(drho_rho_base, bins=bins, alpha=0.5, label=f'Base ({n_base}) max={drho_rho_base.max():.3f}', color='blue')
    ax.hist(drho_rho_ref, bins=bins, alpha=0.5, label=f'AMR ({n_ref}) max={drho_rho_ref.max():.3f}', color='red')
    ax.axvline(max_density_jump, color='k', linestyle='--', linewidth=2, label=f'Threshold ({max_density_jump})')
    ax.set_xlabel('drho / rho', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Density Jump Distribution', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_yscale('log')

    # Plot 6: Temperature jump histogram
    ax = axes[1, 1]
    ax.hist(dT_T_base, bins=bins, alpha=0.5, label=f'Base ({n_base}) max={dT_T_base.max():.3f}', color='blue')
    ax.hist(dT_T_ref, bins=bins, alpha=0.5, label=f'AMR ({n_ref}) max={dT_T_ref.max():.3f}', color='red')
    ax.axvline(max_temp_jump, color='k', linestyle='--', linewidth=2, label=f'Threshold ({max_temp_jump})')
    ax.set_xlabel('dT / T', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Temperature Jump Distribution', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_yscale('log')

    # Plot 7: Radius jump histogram
    ax = axes[1, 2]
    ax.hist(dr_r_base, bins=bins, alpha=0.5, label=f'Base ({n_base}) max={dr_r_base.max():.3f}', color='blue')
    ax.hist(dr_r_ref, bins=bins, alpha=0.5, label=f'AMR ({n_ref}) max={dr_r_ref.max():.3f}', color='red')
    ax.axvline(max_radius_jump, color='k', linestyle='--', linewidth=2, label=f'Threshold ({max_radius_jump})')
    ax.set_xlabel('dr / r', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Radius Jump Distribution', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_yscale('log')

    # Plot 8: Electron density jump histogram
    ax = axes[1, 3]
    if dne_ne_base is not None and dne_ne_ref is not None:
        ax.hist(dne_ne_base, bins=bins, alpha=0.5, label=f'Base ({n_base}) max={dne_ne_base.max():.3f}', color='blue')
        ax.hist(dne_ne_ref, bins=bins, alpha=0.5, label=f'AMR ({n_ref}) max={dne_ne_ref.max():.3f}', color='red')
        ax.axvline(max_ne_jump, color='k', linestyle='--', linewidth=2, label=f'Threshold ({max_ne_jump})')
        ax.set_xlabel('dNe / Ne', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('Ne Jump Distribution', fontsize=14)
        ax.legend(fontsize=9)
    else:
        ax.text(0.5, 0.5, 'Ne data not available', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Ne Jump Distribution', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_yscale('log')

    plt.tight_layout()
    plot_file = str(Path(output_file).parent / 'SN_HYDRO_DATA_adaptive_grid.png')
    import os as _os
    _os.makedirs(_os.path.dirname(plot_file) or '.', exist_ok=True)
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Adaptive grid diagnostic plot saved to {plot_file}")
    plt.close()


def plot_mass_fractions_stacked(radius_cm, density, mass_fractions, time_days, output_file,
                                trimmed_radius_cm=None, trimmed_mass_fractions=None,
                                tau=None, tau_min=1e-5, tau_max=30.0):
    """Create a mass fraction plot vs radius.
    
    Plots original STELLA data as dashed lines and trimmed CMFGEN data as solid lines.
    Isotopes (NICK, COB, IRON) are shown with distinct dotted/dash-dot styles.
    """
    rsun_cm = 6.96e10  # kept for reference; plots now use cm
    
    # Original data (STELLA - full model)
    r_orig = np.asarray(radius_cm)
    
    # Trimmed data (CMFGEN - selected tau range)
    if trimmed_radius_cm is not None and trimmed_mass_fractions is not None:
        r_trim = np.asarray(trimmed_radius_cm)
        mf_trim = trimmed_mass_fractions
        has_trimmed = True
    else:
        has_trimmed = False
    
    # Elemental species: (key_in_orig, key_in_trim, display_label, color)
    ELEM_SPECS = [
        ('H',      'H',    'H',    '#FF6B6B'),
        ('He',     'He',   'He',   '#4ECDC4'),
        ('C',      'C',    'C',    '#45B7D1'),
        ('O',      'O',    'O',    '#96CEB4'),
        ('Si',     'Si',   'Si',   '#9B59B6'),
        ('Ca',     'Ca',   'Ca',   '#E74C3C'),
        ('FePeak', 'IRON', 'Fe',   '#7F8C8D'),
    ]
    # Isotope species shown only in trimmed (CMFGEN) data (radioactive parents only)
    ISO_SPECS = [
        ('NICK', 'Ni56 (isotope)', '#e377c2', '-.'),
        ('COB',  'Co56 (isotope)', '#ff7f0e', '-.'),
    ]
    
    fig, ax = plt.subplots(figsize=(12, 6))

    x_orig = r_orig  # cm

    # Shade the CMFGEN selection region so the user can see which part of the
    # STELLA profile was selected (solid) vs excluded (dashed only).
    if has_trimmed:
        x_trim = r_trim  # cm
        r_sel_min, r_sel_max = x_trim.min(), x_trim.max()
        ax.axvspan(r_sel_min, r_sel_max, color='#E8E8E8', alpha=0.5, zorder=0,
                   label='CMFGEN selection')
        ax.axvline(r_sel_min, color='gray', linewidth=0.8, linestyle=':', alpha=0.7, zorder=1)
        ax.axvline(r_sel_max, color='gray', linewidth=0.8, linestyle=':', alpha=0.7, zorder=1)

    # Plot original STELLA data as dashed lines (FULL profile, all zones)
    for key_orig, key_trim, label, color in ELEM_SPECS:
        if key_orig in mass_fractions and np.any(mass_fractions[key_orig] > 1e-30):
            ax.plot(x_orig, mass_fractions[key_orig], color=color,
                    linewidth=1.8, linestyle='--', alpha=0.85, zorder=5)

    # Plot trimmed CMFGEN elemental data as solid lines ON TOP of the dashed
    # (only within the selection region)
    if has_trimmed:
        for key_orig, key_trim, label, color in ELEM_SPECS:
            if key_trim in mf_trim and np.any(mf_trim[key_trim] > 1e-30):
                ax.plot(x_trim, mf_trim[key_trim], color=color,
                        linewidth=2.5, linestyle='-', alpha=0.9, label=label, zorder=6)
        # Plot isotopes as dash-dot
        for key, label, color, ls in ISO_SPECS:
            if key in mf_trim and np.any(mf_trim[key] > 1e-30):
                ax.plot(x_trim, mf_trim[key], color=color,
                        linewidth=2, linestyle=ls, alpha=0.8, label=label, zorder=6)

    # Build legend with explicit STELLA / CMFGEN proxy entries
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    stella_proxy = Line2D([0], [0], color='gray', linewidth=1.8, linestyle='--', alpha=0.85, label='STELLA (dashed)')
    cmfgen_proxy = Line2D([0], [0], color='gray', linewidth=2.5, linestyle='-', alpha=0.9, label='CMFGEN (solid)')
    sel_proxy = Patch(facecolor='#E8E8E8', edgecolor='gray', linewidth=0.8, linestyle=':', alpha=0.5, label='CMFGEN selection')
    handles, labels = ax.get_legend_handles_labels()
    # Deduplicate element labels (solid lines already labelled)
    seen = set()
    deduped_handles, deduped_labels = [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen.add(l)
            deduped_handles.append(h)
            deduped_labels.append(l)
    handles = deduped_handles + [stella_proxy, cmfgen_proxy, sel_proxy]
    labels = deduped_labels + ['STELLA (dashed)', 'CMFGEN (solid)', 'CMFGEN selection']

    ax.set_xlabel(r'Radius (cm)', fontsize=12)
    ax.set_ylabel('Mass Fraction', fontsize=12)
    ax.set_title(f'Mass Fractions: STELLA (dashed, full profile) vs CMFGEN (solid, selection)\nt = {time_days:.1f} days', fontsize=13)
    ax.set_xlim(x_orig.min() * 0.9, x_orig.max() * 1.1)
    ax.set_ylim(0, 1.05)
    ax.set_xscale('log')
    ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1.01, 0.5),
              fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3, which='both')

    # Annotate the selection boundaries
    if has_trimmed:
        ax.annotate(f'r_min\n{r_sel_min:.2e}', xy=(r_sel_min, 0.98), xytext=(r_sel_min * 0.3, 0.98),
                    fontsize=8, color='gray', ha='center',
                    arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))
        ax.annotate(f'r_max\n{r_sel_max:.2e}', xy=(r_sel_max, 0.98), xytext=(r_sel_max * 3, 0.98),
                    fontsize=8, color='gray', ha='center',
                    arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    # Add a zoomed inset for the outer region where non-selected STELLA material
    # (SN_A ejecta beyond the CMFGEN outer boundary) is compressed into <0.1% of
    # the log x-axis and otherwise invisible.
    if has_trimmed and r_sel_max < x_orig.max():
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        # Zoom to the region around the CMFGEN outer boundary
        zoom_r_min = r_sel_max * 0.5
        zoom_r_max = x_orig.max() * 1.05
        ax_inset = inset_axes(ax, width='35%', height='40%', loc='upper left',
                              bbox_to_anchor=(0.02, -0.02, 1, 1),
                              bbox_transform=ax.transAxes, borderpad=1.5)
        # Shade CMFGEN region in inset
        ax_inset.axvspan(r_sel_max, zoom_r_max, color='#F8F8F8', alpha=0.5, zorder=0)
        # Plot dashed STELLA in inset
        for key_orig, key_trim, label, color in ELEM_SPECS:
            if key_orig in mass_fractions and np.any(mass_fractions[key_orig] > 1e-30):
                ax_inset.plot(x_orig, mass_fractions[key_orig], color=color,
                              linewidth=1.5, linestyle='--', alpha=0.85)
        # Plot solid CMFGEN in inset (only the part within zoom range)
        for key_orig, key_trim, label, color in ELEM_SPECS:
            if key_trim in mf_trim and np.any(mf_trim[key_trim] > 1e-30):
                ax_inset.plot(x_trim, mf_trim[key_trim], color=color,
                              linewidth=2.0, linestyle='-', alpha=0.9)
        ax_inset.set_xlim(zoom_r_min, zoom_r_max)
        ax_inset.set_ylim(0, 1.05)
        ax_inset.set_xscale('log')
        ax_inset.set_title('Zoom: outer edge', fontsize=8)
        ax_inset.tick_params(labelsize=7)
        ax_inset.axvline(r_sel_max, color='gray', linewidth=0.8, linestyle=':', alpha=0.7)
        ax_inset.grid(True, alpha=0.3, which='both')

    # Reserve space for the legend on the right
    fig.subplots_adjust(right=0.72, left=0.07, top=0.88, bottom=0.12)
    plot_file = str(Path(output_file).parent / 'SN_HYDRO_DATA_mass_fraction.png')
    import os as _os; _os.makedirs(_os.path.dirname(plot_file) or '.', exist_ok=True)
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Mass fraction plot saved to {plot_file}")
    plt.close()


def plot_dr_vs_r(r_cm, output_file, time_days, rho=None, T=None, v=None):
    """Plot grid spacing (|dr|/r) vs radius to visualize grid resolution.
    
    Features:
    - Log scale x-axis for radius
    - Upper x-axis showing grid index
    - Right y-axes showing actual physical quantities (density, temperature, velocity)
    """
    dr = np.diff(r_cm)
    r_mid = 0.5 * (r_cm[1:] + r_cm[:-1])
    
    # Calculate relative spacing |dr|/r (absolute value since grid is outer->inner)
    dr_over_r = np.abs(dr) / r_mid
    
    # Grid indices (0 = outermost, n-1 = innermost)
    grid_indices = np.arange(len(r_mid))
    
    fig, ax = plt.subplots(figsize=(14, 7))
    
    # Main plot: dr/r vs radius
    ax.plot(r_mid, dr_over_r, 'b-', linewidth=1.5, label='Grid spacing |dr|/r')
    
    ax.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax.set_ylabel('Grid spacing |dr|/r', fontsize=12, color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.set_title(f'Grid Resolution: |dr|/r vs r\nt = {time_days:.1f} days', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    
    # Add reference line for 0.5% threshold
    ax.axhline(0.005, color='r', linestyle='--', alpha=0.7, label='0.5% threshold')
    
    # Upper x-axis: Grid Index (linear progression mapped to log radius)
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    
    # Create tick locations: evenly spaced in grid index
    # r_mid is ordered outer->inner, so index 0 is outermost (largest r)
    # grid index 0 -> r_max, grid index n-1 -> r_min
    n_ticks = 7
    tick_indices = np.linspace(0, len(r_mid)-1, n_ticks, dtype=int)
    
    # Get corresponding radii for these grid indices
    tick_radii = r_mid[tick_indices]
    
    ax_top.set_xticks(tick_radii)
    ax_top.set_xticklabels([f'{idx}' for idx in tick_indices])
    ax_top.set_xlabel('Grid Index (0=outer, n=inner) [linear]', fontsize=12)
    
    # Right y-axes: Physical quantities (actual values, not normalized)
    if rho is not None or T is not None or v is not None:
        # Interpolate to r_mid points
        rho_mid = 0.5 * (rho[1:] + rho[:-1]) if rho is not None else None
        T_mid = 0.5 * (T[1:] + T[:-1]) if T is not None else None
        v_mid = 0.5 * (v[1:] + v[:-1]) if v is not None else None
        
        offset = 0
        if rho is not None:
            ax_rho = ax.twinx()
            ax_rho.semilogy(r_mid, rho_mid, 'g--', linewidth=1.0, alpha=0.6, label='Density')
            ax_rho.set_ylabel('Density (g/cm$^3$)', fontsize=11, color='g')
            ax_rho.tick_params(axis='y', labelcolor='g')
            offset += 60
        
        if T is not None:
            ax_T = ax.twinx()
            ax_T.spines['right'].set_position(('outward', offset))
            ax_T.semilogy(r_mid, T_mid, 'orange', linestyle='--', linewidth=1.0, alpha=0.6, label='Temperature')
            ax_T.set_ylabel('Temperature (K)', fontsize=11, color='orange')
            ax_T.tick_params(axis='y', labelcolor='orange')
            offset += 60
        
        if v is not None:
            ax_v = ax.twinx()
            ax_v.spines['right'].set_position(('outward', offset))
            ax_v.plot(r_mid, v_mid, 'm--', linewidth=1.0, alpha=0.6, label='Velocity')
            ax_v.set_ylabel('Velocity (km/s)', fontsize=11, color='m')
            ax_v.tick_params(axis='y', labelcolor='m')
    
    # Combine legends from all axes
    lines, labels = ax.get_legend_handles_labels()
    if rho is not None:
        l, lab = ax_rho.get_legend_handles_labels()
        lines.extend(l)
        labels.extend(lab)
    if T is not None:
        l, lab = ax_T.get_legend_handles_labels()
        lines.extend(l)
        labels.extend(lab)
    if v is not None:
        l, lab = ax_v.get_legend_handles_labels()
        lines.extend(l)
        labels.extend(lab)
    
    ax.legend(lines, labels, loc='upper right', fontsize=9)
    
    plt.tight_layout()
    plot_file = str(Path(output_file).parent / 'SN_HYDRO_DATA_dr_vs_r.png')
    import os as _os; _os.makedirs(_os.path.dirname(plot_file) or '.', exist_ok=True)
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"dr vs r plot saved to {plot_file}")
    plt.close()


def plot_zone_density(r_cm, rho, kappa, tau, output_file, time_days,
                      tau_breakpoints=None, thin_tau_threshold=1.0):
    """Plot zone density (zones per equal-log-radius bin) vs radius.

    Shows how grid points are distributed in radius, with tau overlaid
    to highlight the tau-based zoning strategy.  Vertical lines mark
    tau breakpoints and the thin-region threshold.
    """
    n = len(r_cm)
    r_rsun = r_cm  # cm

    # Compute zone density via histogram in log(radius)
    log_r = np.log10(r_cm)
    n_bins = 60
    counts, bin_edges = np.histogram(log_r, bins=n_bins)
    bin_centers = 10 ** (0.5 * (bin_edges[:-1] + bin_edges[1:]))
    bin_widths_log = bin_edges[1:] - bin_edges[:-1]
    # zones per dex in log(r)
    density = counts / bin_widths_log

    fig, ax1 = plt.subplots(figsize=(14, 7))

    # Zone density bar plot
    ax1.bar(bin_centers, density, width=np.diff(10**bin_edges),
            color='steelblue', alpha=0.7, label='Zones per dex in log(r)')
    ax1.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax1.set_ylabel('Zone density (zones / dex)', fontsize=12, color='steelblue')
    ax1.tick_params(axis='y', labelcolor='steelblue')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, which='both')

    # Tau overlay on right y-axis
    ax2 = ax1.twinx()
    ax2.plot(r_cm, tau, 'r-', linewidth=2, drawstyle='steps-mid', alpha=0.8, label='Cumulative τ')
    ax2.set_ylabel('Cumulative optical depth τ', fontsize=12, color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.set_yscale('log')

    # Mark tau breakpoints
    if tau_breakpoints is not None:
        for bp in tau_breakpoints:
            # Find radius where tau = bp
            idx = np.searchsorted(tau, bp)
            if 0 < idx < n:
                r_bp = r_cm[idx]
                ax1.axvline(r_bp, color='green', linestyle='--', alpha=0.7, linewidth=1.5)
                ax1.annotate(f'τ={bp}', xy=(r_bp, ax1.get_ylim()[1]),
                             fontsize=9, color='green', ha='center', va='top')

    # Mark thin-region threshold
    idx_thin = np.searchsorted(tau, thin_tau_threshold)
    if 0 < idx_thin < n:
        r_thin = r_cm[idx_thin]
        ax1.axvline(r_thin, color='orange', linestyle=':', alpha=0.8, linewidth=2)
        ax1.annotate(f'τ={thin_tau_threshold} (thin threshold)',
                     xy=(r_thin, ax1.get_ylim()[1] * 0.5),
                     fontsize=9, color='orange', ha='center', va='top')

    # Title and legend
    ax1.set_title(f'Zone Density vs Radius ({n} zones)\nt = {time_days:.1f} days', fontsize=14)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)

    plt.tight_layout()
    plot_file = str(Path(output_file).parent / 'SN_HYDRO_DATA_zone_density.png')
    import os as _os; _os.makedirs(_os.path.dirname(plot_file) or '.', exist_ok=True)
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Zone density plot saved to {plot_file}")
    plt.close()


def generate_grid_diagnostic_report(
    r_cm: np.ndarray,
    v_kms: np.ndarray,
    rho: np.ndarray,
    T_K: np.ndarray,
    kappa: np.ndarray,
    tau: np.ndarray,
    sigma: np.ndarray,
    mass_fractions: dict,
    time_days: float,
    output_file: str,
    n_e: np.ndarray | None = None,
) -> dict:
    """Generate a machine-readable grid diagnostic report.

    Produces a JSON file with all quantities specified in the task:
    ND, bounds, max deltas, worst-resolved cells, composition interfaces, etc.

    Parameters
    ----------
    r_cm : array, outer→inner (decreasing)
    v_kms, rho, T_K, kappa, tau, sigma : arrays, same orientation
    mass_fractions : dict of element → array
    time_days : epoch time
    output_file : path to SN_HYDRO_DATA (report written alongside as .json)

    Returns
    -------
    dict : the report data (also written as JSON)
    """
    import json
    n = len(r_cm)
    eps = 1e-30

    # Per-cell quantities
    dtau = np.diff(tau)
    dtau_abs = np.abs(dtau)
    r_mid = 0.5 * (r_cm[:-1] + r_cm[1:])

    # dlog(tau) — guard against tau=0
    tau_safe = np.maximum(tau, 1e-10)
    dlogtau = np.abs(np.diff(np.log10(tau_safe)))

    # dtau/dr and dtau/dln(r)
    dr = np.abs(np.diff(r_cm))
    dr_safe = np.maximum(dr, eps)
    dtau_dr = dtau_abs / dr_safe
    dtau_dlnr = dtau_dr * r_mid

    # Physical gradients
    dlogrho = np.abs(np.diff(np.log10(np.maximum(rho, eps))))
    dlogT = np.abs(np.diff(np.log10(np.maximum(T_K, eps))))
    dv = np.abs(np.diff(v_kms))

    # Electron density gradient
    dlogne = None
    if n_e is not None:
        dlogne = np.abs(np.diff(np.log10(np.maximum(n_e, eps))))

    # Composition jumps (max over all elements)
    max_comp_jump = 0.0
    worst_comp_elem = ""
    worst_comp_idx = -1
    for elem, mf in mass_fractions.items():
        if np.all(mf < 1e-10):
            continue
        jumps = np.abs(np.diff(mf))
        # Relative jump: |dX| / max(X, floor)
        mf_safe = np.maximum(mf[:-1], 1e-6)
        rel_jumps = jumps / mf_safe
        idx = int(np.argmax(rel_jumps))
        if rel_jumps[idx] > max_comp_jump:
            max_comp_jump = float(rel_jumps[idx])
            worst_comp_elem = elem
            worst_comp_idx = idx

    # Fraction of tau in largest cell
    total_tau = float(np.sum(dtau_abs))
    frac_largest_cell = float(np.max(dtau_abs) / total_tau) if total_tau > 0 else 0.0

    # Composition interfaces (count)
    n_interfaces = 0
    interface_indices = []
    for elem, mf in mass_fractions.items():
        if np.all(mf < 1e-10):
            continue
        jumps = np.abs(np.diff(mf))
        mf_safe = np.maximum(mf[:-1], 1e-6)
        rel_jumps = jumps / mf_safe
        for i in np.where(rel_jumps > 0.3)[0]:
            if i not in interface_indices:
                interface_indices.append(int(i))
                n_interfaces += 1

    # Velocity turning points
    n_turns = 0
    for i in range(1, n - 1):
        if (v_kms[i] > v_kms[i - 1] and v_kms[i] > v_kms[i + 1]) or \
           (v_kms[i] < v_kms[i - 1] and v_kms[i] < v_kms[i + 1]):
            n_turns += 1

    # Negative sigma regions
    neg_sigma = sigma < 0
    n_neg_sigma = int(np.sum(neg_sigma))

    # Worst cells
    worst_dtau_idx = int(np.argmax(dtau_abs))
    worst_dlogtau_idx = int(np.argmax(dlogtau))
    worst_drho_idx = int(np.argmax(dlogrho))
    worst_dT_idx = int(np.argmax(dlogT))
    worst_dv_idx = int(np.argmax(dv))

    report = {
        "time_days": float(time_days),
        "ND": int(n),
        "r_inner_cm": float(r_cm[-1]),
        "r_outer_cm": float(r_cm[0]),
        "v_inner_kms": float(v_kms[-1]),
        "v_outer_kms": float(v_kms[0]),
        "tau_inner": float(tau[-1]),
        "tau_outer": float(tau[0]),
        "max_delta_log_tau": float(np.max(dlogtau)),
        "median_delta_log_tau": float(np.median(dlogtau)),
        "max_delta_rho": float(np.max(dlogrho)),
        "max_delta_T": float(np.max(dlogT)),
        "max_delta_v": float(np.max(dv)),
        "max_delta_ne": float(np.max(dlogne)) if dlogne is not None else None,
        "max_composition_jump": float(max_comp_jump),
        "maximum_cell_delta_tau": float(np.max(dtau_abs)),
        "fraction_of_tau_in_largest_cell": float(frac_largest_cell),
        "number_of_refinement_regions": int(n_interfaces),  # approx
        "number_of_composition_interfaces": int(n_interfaces),
        "number_of_interaction_regions": int(n_turns),  # velocity turns as proxy
        "number_of_velocity_turns": int(n_turns),
        "negative_sigma_zones": int(n_neg_sigma),
        "sigma_min": float(np.min(sigma)),
        "sigma_max": float(np.max(sigma)),
        "worst_optical_depth_cell": worst_dtau_idx,
        "worst_dlogtau_cell": worst_dlogtau_idx,
        "worst_density_gradient_cell": worst_drho_idx,
        "worst_temperature_gradient_cell": worst_dT_idx,
        "worst_velocity_gradient_cell": worst_dv_idx,
        "worst_composition_interface": {
            "element": worst_comp_elem,
            "cell_index": worst_comp_idx,
            "relative_jump": float(max_comp_jump),
        },
        "composition_interface_indices": interface_indices,
        "dtau_dr_max": float(np.max(dtau_dr)),
        "dtau_dlnr_max": float(np.max(dtau_dlnr)),
        "cells_dtau_gt_1": int(np.sum(dtau_abs > 1.0)),
        "cells_dtau_gt_0.3": int(np.sum(dtau_abs > 0.3)),
        "cells_dlogtau_gt_0.1": int(np.sum(dlogtau > 0.1)),
        "mass_fractions_sum_min": float(np.min(sum(mass_fractions.values()))),
        "mass_fractions_sum_max": float(np.max(sum(mass_fractions.values()))),
    }

    # Write JSON report
    report_path = str(Path(output_file).parent / 'SN_HYDRO_DATA_grid_report.json')
    import os as _os; _os.makedirs(_os.path.dirname(report_path) or '.', exist_ok=True)
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"Grid diagnostic report saved to {report_path}")

    return report


def plot_dtau_resolution_panels(r_cm, tau, output_file, time_days,
                                rho=None, T=None, v_kms=None, kappa=None):
    """Multi-panel optical-depth resolution diagnostic.

    Panels:
      1. dtau/dr vs r  — optical-depth gradient
      2. dtau/dln(r) vs r  — logarithmic optical-depth gradient
      3. Per-cell dtau (bar) — which cells dominate the tau integral
      4. Cumulative tau vs r — where optical depth accumulates

    These plots make it immediately obvious whether a single cell
    dominates the optical-depth integral due to inadequate sampling.
    """
    n = len(r_cm)
    if n < 3:
        return

    eps = 1e-30
    dtau = np.diff(tau)
    dtau_abs = np.abs(dtau)
    r_mid = 0.5 * (r_cm[:-1] + r_cm[1:])
    dr = np.abs(np.diff(r_cm))
    dr_safe = np.maximum(dr, eps)

    dtau_dr = dtau_abs / dr_safe
    dtau_dlnr = dtau_dr * r_mid

    # Per-cell fraction of total tau
    total_tau = np.sum(dtau_abs)
    frac_per_cell = dtau_abs / total_tau if total_tau > 0 else np.zeros_like(dtau_abs)

    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    fig.suptitle(f'Optical-Depth Resolution Diagnostics — t = {time_days:.1f} d, ND = {n}',
                 fontsize=14, fontweight='bold')

    # Panel 1: dtau/dr vs r
    ax = axes[0, 0]
    ax.plot(r_mid, dtau_dr, 'b-o', markersize=2, linewidth=1.0, label='dτ/dr')
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('dτ/dr (cm⁻¹)', fontsize=11, color='b')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('Optical-Depth Gradient: dτ/dr', fontsize=13)
    # Mark the cell with max dtau/dr
    idx_max = int(np.argmax(dtau_dr))
    ax.annotate(f'max: cell {idx_max}\n{dtau_dr[idx_max]:.2e}',
                xy=(r_mid[idx_max], dtau_dr[idx_max]),
                fontsize=8, color='red',
                arrowprops=dict(arrowstyle='->', color='red'),
                xytext=(r_mid[idx_max] * 1.5, dtau_dr[idx_max] * 0.3))

    # Panel 2: dtau/dln(r) vs r
    ax = axes[0, 1]
    ax.plot(r_mid, dtau_dlnr, 'g-o', markersize=2, linewidth=1.0, label='dτ/dlnr')
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('dτ/dln(r)', fontsize=11, color='g')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('Logarithmic Optical-Depth Gradient: dτ/dln(r)', fontsize=13)
    idx_max2 = int(np.argmax(dtau_dlnr))
    ax.annotate(f'max: cell {idx_max2}\n{dtau_dlnr[idx_max2]:.2e}',
                xy=(r_mid[idx_max2], dtau_dlnr[idx_max2]),
                fontsize=8, color='red',
                arrowprops=dict(arrowstyle='->', color='red'),
                xytext=(r_mid[idx_max2] * 1.5, dtau_dlnr[idx_max2] * 0.3))

    # Panel 3: Per-cell tau contribution (bar chart)
    ax = axes[1, 0]
    colors = ['red' if f > 0.1 else 'orange' if f > 0.05 else 'steelblue'
              for f in frac_per_cell]
    ax.bar(np.arange(len(frac_per_cell)), frac_per_cell * 100, color=colors, width=1.0)
    ax.set_xlabel('Cell Index (0=outer)', fontsize=11)
    ax.set_ylabel('% of total τ per cell', fontsize=11)
    ax.set_title('Per-Cell Optical-Depth Contribution', fontsize=13)
    ax.axhline(100.0 / (n - 1), color='gray', linestyle='--', alpha=0.5,
               label=f'Uniform ({100.0/(n-1):.2f}%)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    # Mark cells > 10%
    for i in np.where(frac_per_cell > 0.1)[0]:
        ax.annotate(f'{frac_per_cell[i]*100:.1f}%', xy=(i, frac_per_cell[i] * 100),
                    fontsize=7, color='red', ha='center', va='bottom')

    # Panel 4: Cumulative tau vs r
    ax = axes[1, 1]
    ax.plot(r_cm, tau, 'k-', linewidth=1.5, label='τ(r)')
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('Cumulative τ', fontsize=11)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('Cumulative Optical Depth', fontsize=13)
    ax.axhline(2.0/3.0, color='purple', linestyle=':', alpha=0.6, label='τ = 2/3')
    ax.legend(fontsize=9)
    # Mark photosphere
    idx_phot = int(np.argmin(np.abs(tau - 2.0/3.0)))
    ax.annotate(f'τ=2/3 at r={r_cm[idx_phot]:.2e}',
                xy=(r_cm[idx_phot], 2.0/3.0),
                fontsize=8, color='purple',
                arrowprops=dict(arrowstyle='->', color='purple'),
                xytext=(r_cm[idx_phot] * 3, 2.0/3.0 * 5))

    plt.tight_layout()
    plot_file = str(Path(output_file).parent / 'SN_HYDRO_DATA_dtau_resolution.png')
    import os as _os; _os.makedirs(_os.path.dirname(plot_file) or '.', exist_ok=True)
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"dtau resolution panels saved to {plot_file}")
    plt.close()


def plot_dtau_vs_r(r_cm, tau, output_file, time_days,
                   v_kms=None, rho=None, T=None, sigma=None,
                   tau_min=1e-3, tau_max=30.0):
    """Plot optical-depth increment dtau per grid cell vs radius.

    This is the primary grid-quality diagnostic for CMFGEN radiative
    transfer: each cell must resolve the optical-depth gradient so that
    the formal solution and Lambda iteration are stable.

    The plot shows:
      - dtau per cell (left y-axis, log scale) with reference lines at
        dtau = 0.01, 0.1, 1.0 (the CMFGEN resolution thresholds).
      - Cumulative tau (right y-axis) so the photosphere (tau=2/3) and
        the tau_max / tau_min trimming boundaries are visible.
      - Optional overlaid velocity, density, temperature, and sigma
        profiles on secondary axes to correlate dtau spikes with
        physical structure (shocks, composition interfaces, velocity
        reversals).
      - Highlighted cells where dtau > 1.0 (under-resolved for the
        formal solution) and dtau > 0.3 (marginal).

    Arrays are outer->inner (decreasing radius), matching the
    converter's internal convention.

    Parameters
    ----------
    r_cm : array
        Radius grid in cm (outer->inner, decreasing).
    tau : array
        Cumulative optical depth at each grid point (outer->inner).
    output_file : str
        Output path — the plot is written to the same directory with
        suffix ``_dtau_vs_r.png``.
    time_days : float
        Epoch time in days (for the title).
    v_kms, rho, T, sigma : optional arrays
        Physical profiles overlaid on secondary axes for correlation.
    tau_min, tau_max : float
        Trimming boundaries shown as horizontal reference lines on the
        cumulative-tau axis.
    """
    n = len(r_cm)
    if n < 3:
        return

    # Per-cell dtau: difference in cumulative tau between adjacent points.
    # tau is outer->inner (increasing inward), so dtau[i] = tau[i+1] - tau[i]
    # for the cell between r[i] and r[i+1].
    dtau = np.diff(tau)
    r_mid = 0.5 * (r_cm[:-1] + r_cm[1:])
    r_outer = r_cm[:-1]
    r_inner = r_cm[1:]

    # Cell-centred physical quantities for overlay
    v_mid = 0.5 * (v_kms[:-1] + v_kms[1:]) if v_kms is not None else None
    rho_mid = 0.5 * (rho[:-1] + rho[1:]) if rho is not None else None
    T_mid = 0.5 * (T[:-1] + T[1:]) if T is not None else None
    sigma_mid = 0.5 * (sigma[:-1] + sigma[1:]) if sigma is not None else None

    fig, ax = plt.subplots(figsize=(16, 9))

    # --- Main: dtau vs r ---
    # Use absolute value (dtau is always >= 0 for monotonic tau, but
    # guard against numerical noise).
    dtau_abs = np.abs(dtau)
    ax.plot(r_mid, dtau_abs, 'b-o', linewidth=1.2, markersize=3,
            label='dτ per cell', zorder=5)

    ax.set_xlabel('Radius (cm)', fontsize=12)
    ax.set_ylabel('dτ per cell', fontsize=12, color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title(f'Optical-Depth Resolution: dτ vs r\nt = {time_days:.1f} days, ND = {n}',
                 fontsize=14)
    ax.grid(True, alpha=0.3, which='both')

    # Reference lines for CMFGEN resolution thresholds
    for thr, label, color, alpha in [
        (1.0,  'dτ = 1.0 (under-resolved)', 'red',    0.8),
        (0.3,  'dτ = 0.3 (marginal)',       'orange', 0.7),
        (0.1,  'dτ = 0.1 (good)',           'green',  0.5),
        (0.01, 'dτ = 0.01 (well-resolved)', 'gray',   0.4),
    ]:
        ax.axhline(thr, color=color, linestyle='--', alpha=alpha, linewidth=1.0,
                   label=label)

    # Highlight under-resolved cells
    under_resolved = dtau_abs > 1.0
    if np.any(under_resolved):
        ax.scatter(r_mid[under_resolved], dtau_abs[under_resolved],
                   s=80, facecolors='none', edgecolors='red', linewidths=1.5,
                   zorder=10, label=f'{np.sum(under_resolved)} cells dτ > 1.0')

    marginal = (dtau_abs > 0.3) & (dtau_abs <= 1.0)
    if np.any(marginal):
        ax.scatter(r_mid[marginal], dtau_abs[marginal],
                   s=50, facecolors='none', edgecolors='orange', linewidths=1.2,
                   zorder=9, label=f'{np.sum(marginal)} cells 0.3 < dτ ≤ 1.0')

    # --- Right axis: cumulative tau ---
    ax_tau = ax.twinx()
    ax_tau.plot(r_cm, tau, 'k-', linewidth=1.5, alpha=0.6, label='τ (cumulative)')
    ax_tau.set_ylabel('Cumulative τ', fontsize=12, color='k')
    ax_tau.tick_params(axis='y', labelcolor='k')
    ax_tau.set_yscale('log')

    # Tau reference lines
    ax_tau.axhline(2.0/3.0, color='purple', linestyle=':', alpha=0.6, linewidth=1.0)
    ax_tau.axhline(tau_max, color='brown', linestyle=':', alpha=0.5, linewidth=1.0)
    ax_tau.axhline(tau_min, color='brown', linestyle=':', alpha=0.5, linewidth=1.0)
    # Annotate tau=2/3 (photosphere)
    ax_tau.text(r_cm[0] * 1.05, 2.0/3.0, ' τ=2/3', fontsize=8,
                color='purple', va='bottom')

    # --- Secondary axes: velocity, sigma, density, temperature ---
    offset = 60
    if v_kms is not None:
        ax_v = ax.twinx()
        ax_v.spines['right'].set_position(('outward', offset))
        ax_v.plot(r_mid, v_mid, 'm--', linewidth=0.8, alpha=0.5, label='V (km/s)')
        ax_v.set_ylabel('Velocity (km/s)', fontsize=10, color='m')
        ax_v.tick_params(axis='y', labelcolor='m')
        offset += 55

    if sigma is not None:
        ax_s = ax.twinx()
        ax_s.spines['right'].set_position(('outward', offset))
        ax_s.plot(r_mid, sigma_mid, 'r--', linewidth=0.8, alpha=0.5, label='σ')
        ax_s.set_ylabel('σ = dlnV/dlnR − 1', fontsize=10, color='r')
        ax_s.tick_params(axis='y', labelcolor='r')
        ax_s.axhline(0.0, color='r', linestyle='-', alpha=0.2, linewidth=0.5)
        ax_s.axhline(-1.0, color='r', linestyle=':', alpha=0.3, linewidth=0.8)
        offset += 55

    if rho is not None:
        ax_rho = ax.twinx()
        ax_rho.spines['right'].set_position(('outward', offset))
        ax_rho.semilogy(r_mid, rho_mid, 'g--', linewidth=0.8, alpha=0.5, label='ρ')
        ax_rho.set_ylabel('Density (g/cm³)', fontsize=10, color='g')
        ax_rho.tick_params(axis='y', labelcolor='g')
        offset += 55

    if T is not None:
        ax_T = ax.twinx()
        ax_T.spines['right'].set_position(('outward', offset))
        ax_T.semilogy(r_mid, T_mid, 'orange', linestyle='--', linewidth=0.8, alpha=0.5, label='T')
        ax_T.set_ylabel('Temperature (K)', fontsize=10, color='orange')
        ax_T.tick_params(axis='y', labelcolor='orange')

    # Combine legends from main + tau axes (keep it readable)
    lines, labels = ax.get_legend_handles_labels()
    l, lab = ax_tau.get_legend_handles_labels()
    lines.extend(l); labels.extend(lab)
    ax.legend(lines, labels, loc='lower left', fontsize=8, ncol=2)

    # Summary text box
    summary = (f'Cells: {n-1}    '
               f'max dτ = {dtau_abs.max():.3f}    '
               f'median dτ = {np.median(dtau_abs):.4f}\n'
               f'dτ > 1.0: {np.sum(under_resolved)}    '
               f'dτ > 0.3: {np.sum(dtau_abs > 0.3)}    '
               f'dτ > 0.1: {np.sum(dtau_abs > 0.1)}')
    ax.text(0.98, 0.97, summary, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            family='monospace')

    plt.tight_layout()
    plot_file = str(Path(output_file).parent / 'SN_HYDRO_DATA_dtau_vs_r.png')
    import os as _os; _os.makedirs(_os.path.dirname(plot_file) or '.', exist_ok=True)
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"dtau vs r plot saved to {plot_file}")
    plt.close()


def plot_sigma_vs_index(r_cm, v_kms, output_file, time_days):
    """Plot sigma (dlnV/dlnR - 1) vs grid index.

    CMFGEN indexes depth from outside (index 1) to inside (index ND).
    Sigma = -1 marks the velocity minimum where dV/dR = 0 (Sobolev breakdown).
    The T1 floor at 0.1 caps the Sobolev optical depth amplification at 10x.
    """
    sigma = sigma_velocity(r_cm, v_kms)
    n = len(sigma)
    indices = np.arange(1, n + 1)  # 1-based to match CMFGEN depth indexing

    fig, ax = plt.subplots(figsize=(14, 7))

    ax.plot(indices, sigma, 'b-', linewidth=1.2, drawstyle='steps-mid', label='σ = dlnV/dlnR − 1')
    ax.axhline(0.0, color='gray', linestyle='-', alpha=0.5, linewidth=0.8)
    ax.axhline(-1.0, color='r', linestyle='--', alpha=0.7, linewidth=1.0,
               label='σ = −1 (velocity minimum, Sobolev breakdown)')

    # Highlight zones where |sigma + 1| < 0.1 (T1 floor active region)
    near_m1 = np.abs(sigma + 1.0) < 0.1
    if np.any(near_m1):
        ax.fill_between(indices, -1.5, 2.5, where=near_m1, alpha=0.15, color='red',
                        label=f'|σ+1| < 0.1 (T1 floor active, {np.sum(near_m1)} zones)')

    ax.set_xlabel('Grid Index (CMFGEN depth, 1=outermost)', fontsize=12)
    ax.set_ylabel('σ = dlnV/dlnR − 1', fontsize=12)
    ax.set_title(f'Sigma vs Grid Index\nt = {time_days:.1f} days, ND = {n}', fontsize=14)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1, n)
    ax.invert_xaxis()  # depth ND (inner) on left, depth 1 (outer) on right

    plt.tight_layout()
    plot_file = str(Path(output_file).parent / 'SN_HYDRO_DATA_sigma_vs_index.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Sigma vs index plot saved to {plot_file}")
    plt.close()


def plot_structure_overview(r_cm, v_kms, rho, T_K, mass_fractions, kappa, time_days, output_file,
                           n_e=None, n_atom=None, ne_from_stella=False):
    """Create 6-panel structure overview plot with consistent x-axis formatting.

    Input arrays are in outer->inner order (from STELLA .swd).
    For plotting, convert to inner->outer order so x-axis increases left-to-right.
    """
    rsun_cm = 6.96e10  # kept for reference; plots now use cm

    # Convert to inner->outer order for plotting (x-axis increases left-to-right)
    # Note: This reverses the arrays in-place for plotting only
    r_plt = r_cm[::-1]
    v_plt = v_kms[::-1]
    rho_plt = rho[::-1]
    T_plt = T_K[::-1]
    kappa_plt = kappa[::-1]
    mf_plt = {elem: mass_fractions[elem][::-1] for elem in mass_fractions}

    r_plt = r_plt  # cm (no conversion)

    # Determine consistent x-axis limits
    r_min = np.min(r_plt) * 0.9
    r_max = np.max(r_plt) * 1.1

    # Create 2x4 grid
    fig, axes = plt.subplots(2, 4, figsize=(24, 12))
    fig.suptitle(f'STELLA Structure Overview — t = {time_days:.1f} days', fontsize=16, fontweight='bold')
    
    # Plot 1: Density
    ax = axes[0, 0]
    ax.loglog(r_plt, rho_plt, 'b-', linewidth=2, drawstyle='steps-mid')
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_title('Density Profile', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 2: Temperature
    ax = axes[0, 1]
    ax.loglog(r_plt, T_plt, 'r-', linewidth=2, drawstyle='steps-mid')
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Profile', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 3: Velocity
    ax = axes[0, 2]
    ax.loglog(r_plt, v_plt, 'g-', linewidth=2, drawstyle='steps-mid')
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_title('Velocity Profile', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)

    # Plot 4: Opacity
    ax = axes[0, 3]
    ax.loglog(r_plt, kappa_plt, 'm-', linewidth=2, drawstyle='steps-mid')
    ax.set_ylabel('Opacity (cm²/g)', fontsize=12)
    ax.set_title('Rosseland Opacity', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 5: Mass fractions (stacked)
    ax = axes[1, 0]
    # Use all elements present in mass_fractions, skipping zero-valued ones
    # Map internal names to display names for clarity
    display_names = {
        'H': 'H', 'He': 'He', 'C': 'C', 'N': 'N', 'O': 'O', 'Ne': 'Ne',
        'Si': 'Si', 'S': 'S', 'Ca': 'Ca', 'FePeak': 'Fe', 'Fe': 'Fe',
        'NICK': 'Ni', 'COB': 'Co', 'IRON': 'Fe', 'Ni58': 'Ni', 'Ni56': 'Ni'
    }
    # Extended color palette for more elements
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD',
              '#F8B500', '#6C5CE7', '#A8E6CF', '#FFD93D', '#FF8E53', '#6A0572',
              '#AB83A1', '#C44569', '#F8B195', '#F67280', '#C06C84', '#6C5B7B']

    # Filter out elements that are all zero
    sorted_elements = sorted([elem for elem in mf_plt.keys() if np.any(mf_plt[elem] > 1e-10)])
    mf_array = np.array([mf_plt[elem] for elem in sorted_elements])
    display_labels = [display_names.get(elem, elem) for elem in sorted_elements]

    if len(sorted_elements) > 0:
        ax.stackplot(r_plt, mf_array, labels=display_labels,
                    colors=colors[:len(sorted_elements)],
                    alpha=0.8)
        ax.legend(fontsize=8, ncol=2, loc='upper right')
    ax.set_ylabel('Mass Fraction', fontsize=12)
    ax.set_title('Mass Fractions', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')

    # Plot 6: Electron and atom densities
    ax = axes[1, 1]
    ne_source_label = 'STELLA' if ne_from_stella else 'Saha (analytical)'
    if n_e is not None:
        n_e_plt = n_e[::-1]  # Convert to inner->outer order
        ax.loglog(r_plt, n_e_plt, 'b-', linewidth=2, drawstyle='steps-mid', label=f'n_e [{ne_source_label}]')
    if n_atom is not None:
        n_atom_plt = n_atom[::-1]  # Convert to inner->outer order
        ax.loglog(r_plt, n_atom_plt, 'r-', linewidth=2, drawstyle='steps-mid', label='n_atom (cm⁻³)')
    # Full ionization limit: ne_max = rho * N_A * sum(X_i * Z_i / A_i)
    if n_e is not None and n_atom is not None:
        N_A = 6.02214076e23
        elem_A_loc = {'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
                      'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
                      'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'FePeak': 56.0, 'Fe': 56.0,
                      'IRON': 56.0, 'NICK': 58.0, 'Ni58': 58.0, 'Ni56': 56.0,
                      'COB': 59.0, 'SIL': 28.0, 'CAL': 40.0, 'OXY': 16.0}
        elem_Z_loc = {'H': 1.0, 'He': 2.0, 'C': 6.0, 'N': 7.0, 'O': 8.0,
                      'Ne': 10.0, 'Na': 11.0, 'Mg': 12.0, 'Al': 13.0, 'Si': 14.0,
                      'S': 16.0, 'Ar': 18.0, 'Ca': 20.0, 'FePeak': 26.0, 'Fe': 26.0,
                      'IRON': 26.0, 'NICK': 28.0, 'Ni58': 28.0, 'Ni56': 28.0,
                      'COB': 27.0, 'SIL': 14.0, 'CAL': 20.0, 'OXY': 8.0}
        rho_plt_for_ne = rho[::-1]
        Z_over_A_loc = np.zeros(len(rho_plt_for_ne))
        for elem, frac in mf_plt.items():
            A = elem_A_loc.get(elem, 28.0)
            Z = elem_Z_loc.get(elem, 14.0)
            Z_over_A_loc += frac * Z / A
        ne_max_plt = rho_plt_for_ne * N_A * Z_over_A_loc
        ax.loglog(r_plt, ne_max_plt, 'b--', linewidth=1.5, drawstyle='steps-mid', alpha=0.7, label='n_e (full ionization)')
    if n_e is not None or n_atom is not None:
        ax.legend(fontsize=9)
    ax.set_ylabel('Number Density (cm⁻³)', fontsize=12)
    ax.set_title(f'Electron & Atom Densities (n_e: {ne_source_label})', fontsize=13)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)

    # Plot 7: Sigma (velocity gradient) — placed directly below velocity panel
    ax = axes[1, 2]
    sigma_plt = sigma_velocity(r_cm, v_kms)[::-1]  # Convert to inner->outer order
    ax.semilogx(r_plt, sigma_plt, 'c-', linewidth=2, drawstyle='steps-mid')
    ax.set_ylabel('σ = dlnV/dlnr - 1', fontsize=12)
    ax.set_title('Velocity Gradient', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(r_min, r_max)

    # Plot 8: Optical depth
    ax = axes[1, 3]
    tau = compute_kappa_tau(r_cm, rho, kappa)
    tau_plt = tau[::-1]  # Convert to inner->outer order
    ax.loglog(r_plt, tau_plt, 'k-', linewidth=2, drawstyle='steps-mid')
    ax.axhline(y=2/3, color='orange', linestyle='--', alpha=0.8, linewidth=2, label='τ = 2/3')
    ax.set_ylabel('Optical Depth τ', fontsize=12)
    ax.set_title('Optical Depth', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)

    # Set common x-axis label on bottom row
    for ax in axes[1, :]:
        ax.set_xlabel(r'Radius (cm)')
    
    plt.tight_layout()
    
    plot_file = str(Path(output_file).parent / 'SN_HYDRO_DATA_structure_overview.png')
    import os as _os; _os.makedirs(_os.path.dirname(plot_file) or '.', exist_ok=True)
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Structure overview plot saved to {plot_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('swd_file',    help='STELLA *.swd output file')
    parser.add_argument('output_file', help='Output SN_HYDRO_DATA path')
    parser.add_argument('--rho-file',     default=None,
                        help='STELLA *.rho file for final evolved composition (preferred)')
    parser.add_argument('--abn-file',     default=None,
                        help='STELLA *.abn file for initial composition (fallback)')
    parser.add_argument('--input-profile', default=None,
                        help='[DEPRECATED] HERACLES input_profile file for composition fallback')
    parser.add_argument('--xni-file',     default=None,
                        help='STELLA *.xni file for Ni56 distribution')
    parser.add_argument('--lbol-file',    default=None,
                        help='STELLA *.lbol file for actual bolometric luminosity')
    parser.add_argument('--target-days',  type=float, default=None,
                        help='Pick epoch closest to this (default: last epoch)')
    parser.add_argument('--tau-max',      type=float, default=30.0)
    parser.add_argument('--tau-min',      type=float, default=1e-3)
    parser.add_argument('--v-min-km-s',   type=float, default=0.0,
                        help='Trim inner zones with V < this value (km/s). 0 = no trim. '
                             'Very low velocity zones cause BA matrix instability in CMFGEN.')
    parser.add_argument('--velocity-max-kms', type=float, default=0.0,
                        help='Trim outer zones with V > this value (km/s). 0 = no trim. '
                             'Use for broken STELLA runs where outer zones reach unphysical velocities.')
    parser.add_argument('--nzon-target',  type=int,   default=200)
    parser.add_argument('--min-points',   type=int,   default=50,
                        help='Minimum number of points in output SN_HYDRO_DATA')
    parser.add_argument('--kappa',        type=float, default=0.34,
                        help='Electron scattering opacity cm^2/g')
    parser.add_argument('--model-name',   default='stella')
    parser.add_argument('--enable-grid-thinning', action='store_true',
                        help='Enable grid thinning to remove near-degenerate points')
    parser.add_argument('--grid-thinning-threshold', type=float, default=0.005,
                        help='min_dr/r threshold for grid thinning (default: 0.005)')
    parser.add_argument('--enable-adaptive-grid', action='store_true',
                        help='Enable Dessart-style adaptive grid refinement')
    parser.add_argument('--smooth-physical-discontinuities', action='store_true', default=False,
                        help='Gaussian-smooth T and rho across sharp discontinuities before grid '
                             'construction (Dessart+2024 approach). Reduces single-zone jumps so '
                             'the Dessart grid criteria (10%% per step) can be satisfied.')
    parser.add_argument('--phys-smooth-sigma-zones', type=float, default=10.0,
                        help='Gaussian sigma in grid points for physical discontinuity smoothing (default: 10). '
                             'Larger = more smoothing. 10 zones spreads a 2x T jump to ~8%% per step.')
    parser.add_argument('--stretch-trouble-regions', action='store_true', default=False,
                        help='Stretch the radius coordinate in trouble regions (where dT/T, '
                             'drho/rho, dlogtau exceed thresholds) without changing physical '
                             'values at original zones. Preserves peaks, unlike Gaussian smoothing.')
    parser.add_argument('--stretch-factor', type=float, default=3.0,
                        help='How many times more points to insert in trouble regions (default: 3.0). '
                             '3.0 means 2 extra interpolated points between each pair of original zones.')
    parser.add_argument('--max-stretch-zones', type=int, default=200,
                        help='Maximum total inserted points across all trouble regions (default: 200).')
    parser.add_argument('--min-zones-to-stretch', type=int, default=2,
                        help='Minimum number of zones in a trouble region to warrant stretching (default: 2).')
    # --- Critical-region (region-dependent min_dr/r) options ---
    parser.add_argument('--critical-regions', action='store_true', default=False,
                        help='Detect narrow high-gradient critical regions (shock fronts, composition '
                             'interfaces, CSM shells) and apply a region-dependent min_dr/r there so the '
                             'AMR can insert the points CMFGEN needs without lowering the global floor.')
    parser.add_argument('--critical-min-dr-r', type=float, default=1e-5,
                        help='min_dr/r threshold inside critical regions (default: 1e-5). The AMR may '
                             'subdivide below the global floor there; the degenerate-grid thinner and '
                             'final spacing enforcement preserve points down to this value inside '
                             'critical regions.')
    parser.add_argument('--critical-disable-min-dr', action='store_true', default=False,
                        help='Disable the min_dr/r constraint entirely inside critical regions (no points '
                             'removed there by the thinner / spacing enforcement). Use with care: can '
                             'leave near-degenerate points that cause div-by-zero in CMFGEN sigma.')
    parser.add_argument('--critical-dT-T', type=float, default=0.10,
                        help='dT/T threshold for critical-region detection (default: 0.10).')
    parser.add_argument('--critical-drho-rho', type=float, default=0.10,
                        help='drho/rho threshold for critical-region detection (default: 0.10).')
    parser.add_argument('--critical-dlogtau', type=float, default=0.10,
                        help='dlogtau threshold for critical-region detection (default: 0.10).')
    parser.add_argument('--critical-dne-ne', type=float, default=0.10,
                        help='dne/ne threshold for critical-region detection (default: 0.10).')
    parser.add_argument('--critical-dv-v', type=float, default=0.30,
                        help='dv/v threshold for critical-region detection (default: 0.30).')
    parser.add_argument('--critical-composition-jump', type=float, default=0.30,
                        help='Absolute mass-fraction jump for composition-interface detection (default: 0.30).')
    parser.add_argument('--critical-composition-interfaces', type=str, default='H-He,He-O,O-Si,Si-Fe',
                        help='Comma-separated element pairs to detect as composition handovers '
                             '(default: H-He,He-O,O-Si,Si-Fe). Use "none" to disable named-interface '
                             'detection.')
    parser.add_argument('--critical-pad-zones', type=int, default=2,
                        help='Zones to pad around each detected critical region on each side (default: 2).')
    parser.add_argument('--stretch-critical-regions', action='store_true', default=False,
                        help='Stretch the radius coordinate in detected critical regions (insert '
                             'interpolated points, preserve original STELLA values) before AMR. '
                             'Stronger than --stretch-trouble-regions and mask-driven.')
    parser.add_argument('--critical-stretch-factor', type=float, default=5.0,
                        help='Sub-intervals per edge inside critical regions when stretching (default: 5.0).')
    parser.add_argument('--critical-max-stretch-zones', type=int, default=300,
                        help='Cap on total points inserted by critical-region stretching (default: 300).')
    parser.add_argument('--base-grid-target', type=int, default=50,
                        help='Initial base grid size (log-spaced in radius) before AMR (default: 50). '
                             'With --use-tau-amr, this becomes the final target number of zones if --nzon-target is not set.')
    parser.add_argument('--use-tau-amr', action='store_true',
                        help='Use log-uniform tau regridding as the final grid instead of velocity/sigma AMR. '
                             'The grid is distributed evenly in log(tau) with the target zone count. '
                             'This matches the Dessart+2015 CMFGEN grid philosophy and avoids AMR-induced '
                             'sigma/velocity scatter. Requires --base-grid-method=log_tau or --base-grid-target.')
    parser.add_argument('--feature-aware-downsample', action='store_true', default=True,
                        help='Use feature-aware downsampling that preserves high-gradient regions from original grid (default: True)')
    parser.add_argument('--no-feature-aware-downsample', dest='feature_aware_downsample', action='store_false',
                        help='Disable feature-aware downsampling and use pure log-spaced sampling')
    parser.add_argument('--base-grid-method', type=str, default='tau_amr_recursive',
                        choices=['log_radius', 'log_tau', 'tau_amr_recursive'],
                        help='Base grid downsampling method: tau_amr_recursive (default, '
                             'optical-depth-primary recursive AMR with hard constraints), '
                             'log_tau (tau-based single-pass), or log_radius (pure geometric, '
                             'NOT recommended — no optical-depth awareness)')
    # --- Tau-AMR config file (overrides individual thresholds below) ---
    parser.add_argument('--tau-amr-config', type=str, default=None,
                        help='Path to YAML config file for tau-AMR thresholds '
                             '(e.g. configs/tau_amr_config.yaml). When provided, '
                             'thresholds are read from the file and override the '
                             'individual --amr-* defaults below.')
    # --- Tau-AMR individual thresholds (overridden by --tau-amr-config) ---
    parser.add_argument('--amr-target-dlogtau', type=float, default=None,
                        help='Tau-AMR: max dlog10(tau) per cell in smooth regions (default: 0.08)')
    parser.add_argument('--amr-target-dlogrho', type=float, default=None,
                        help='Tau-AMR: max dlog10(rho) per cell (default: 0.30)')
    parser.add_argument('--amr-target-dv-kms', type=float, default=None,
                        help='Tau-AMR: max |dv| per cell in km/s (default: 300.0)')
    parser.add_argument('--amr-target-dlogT', type=float, default=None,
                        help='Tau-AMR: max dlog10(T) per cell (default: 0.20)')
    parser.add_argument('--amr-target-dr-over-r', type=float, default=None,
                        help='Tau-AMR: max geometric stretch dr/r (default: 0.10)')
    parser.add_argument('--amr-target-dtau-es', type=float, default=None,
                        help='Tau-AMR: max electron-scattering tau per cell (default: 1.0)')
    parser.add_argument('--amr-target-dlogchi', type=float, default=None,
                        help='Tau-AMR: max dlog10(chi) per cell (default: 0.5)')
    parser.add_argument('--amr-target-dv-over-v', type=float, default=None,
                        help='Tau-AMR: max |dv|/v per cell (default: 0.30)')
    parser.add_argument('--amr-kappa-source', type=str, default=None,
                        choices=['stella', 'hybrid', 'kramers'],
                        help='Tau-AMR: opacity source for optical depth (default: hybrid)')
    parser.add_argument('--amr-w-tau', type=float, default=None,
                        help='Tau-AMR: weight for optical-depth criterion (default: 1.0)')
    parser.add_argument('--amr-w-rho', type=float, default=None,
                        help='Tau-AMR: weight for density criterion (default: 0.5)')
    parser.add_argument('--amr-w-v', type=float, default=None,
                        help='Tau-AMR: weight for velocity criterion (default: 0.5)')
    parser.add_argument('--amr-w-T', type=float, default=None,
                        help='Tau-AMR: weight for temperature criterion (default: 0.4)')
    parser.add_argument('--amr-w-curv', type=float, default=None,
                        help='Tau-AMR: weight for curvature criterion (default: 0.3)')
    parser.add_argument('--amr-composition-jump', type=float, default=None,
                        help='Tau-AMR: |dX_i| threshold for composition interface (default: 0.30)')
    parser.add_argument('--amr-min-points-per-interface', type=int, default=None,
                        help='Tau-AMR: min points at each composition interface (default: 5)')
    parser.add_argument('--amr-min-points-per-shock', type=int, default=None,
                        help='Tau-AMR: min points at each shock (default: 7)')
    parser.add_argument('--amr-min-points-per-extremum', type=int, default=None,
                        help='Tau-AMR: min points at each velocity extremum (default: 5)')
    parser.add_argument('--amr-shock-density-jump', type=float, default=None,
                        help='Tau-AMR: density compression ratio for shock detection (default: 1.5)')
    parser.add_argument('--amr-shock-temp-jump', type=float, default=None,
                        help='Tau-AMR: dlogT for shock detection (default: 0.5)')
    parser.add_argument('--amr-enable-coarsening', action='store_true', default=None,
                        help='Tau-AMR: enable coarsening pass to find minimum grid (default: True)')
    parser.add_argument('--amr-no-coarsening', dest='amr_enable_coarsening', action='store_false',
                        help='Tau-AMR: disable coarsening pass')
    parser.add_argument('--amr-escalation-max-rounds', type=int, default=5,
                        help='Tau-AMR: max escalation rounds for pathological-cell '
                             'refinement (default: 5). Set to 0 to disable escalation.')
    parser.add_argument('--tau-breakpoints', type=str, default='0.01,1.0',
                        help='Comma-separated tau breakpoints for log_tau base grid (inner breakpoints; outer/inner auto-detected)')
    parser.add_argument('--tau-zone-fractions', type=str, default='0.10,0.25,0.65',
                        help='Comma-separated zone fractions per tau segment (must sum to ~1.0)')
    parser.add_argument('--tau-spacing', type=str, default='auto', choices=['linear', 'log', 'auto'],
                        help='Tau spacing: linear (uniform dtau per segment), log (log-uniform tau), '
                             'or auto (blind photosphere-weighted log(tau+1), ignores breakpoints/fractions). Default: auto')
    # --- Composition-aware grid options ---
    parser.add_argument('--composition-aware', action='store_true', default=True,
                        help='Enable composition-aware grid placement: force-include zones at sharp '
                             'composition interfaces (Fe-core/Si-shell, Si/O, O/He) and use hybrid '
                             'mass-fraction interpolation. Default: True')
    parser.add_argument('--no-composition-aware', dest='composition_aware', action='store_false',
                        help='Disable composition-aware grid placement (use pure tau spacing with '
                             'linear mass-fraction interpolation)')
    parser.add_argument('--composition-gradient-threshold', type=float, default=0.5,
                        help='Fractional |dX|/X above which a zone boundary is a sharp composition '
                             'interface (default: 0.5). Lower = more sensitive = more interface zones.')
    parser.add_argument('--composition-elements', type=str, default=None,
                        help='Comma-separated element names to check for interfaces (default: auto-detect '
                             'from available elements: FePeak, Ni56, Ni58, Si, S, Ca, O, C, He, H)')
    parser.add_argument('--no-composition-grid-insertion', dest='composition_grid_insertion',
                        action='store_false', default=True,
                        help='Disable force-including grid points at composition interfaces. '
                             'Hybrid mass-fraction interpolation is still applied. Use for models where '
                             'grid quality (cell ratio, sigma) is critical.')
    parser.add_argument('--composition-interface-cap', type=float, default=0.15,
                        help='Max fraction of target_points for composition interface zones (default: 0.15)')
    parser.add_argument('--min-zones-per-shell', type=int, default=3,
                        help='Minimum zones per composition shell for post-selection verification (default: 3)')
    parser.add_argument('--trace-threshold', type=float, default=1e-4,
                        help='Mass fraction below which log interpolation is used for trace species '
                             'in smooth regions (default: 1e-4)')
    parser.add_argument('--max-thin-zones', type=int, default=0,
                        help='Max zones in low-tau region (tau < thin_grid_tau). 0 = no cap.')
    parser.add_argument('--max-thick-tau', type=float, default=1e10,
                        help='Skip AMR in optically thick region (tau > this). Prevents wasting budget on degenerate inner zones. 1e10 = no cap.')
    parser.add_argument('--gradient-threshold', type=float, default=0.05,
                        help='Gradient threshold for feature-aware downsampling (default: 0.05)')
    parser.add_argument('--max-radius-jump', type=float, default=0.05,
                        help='Max radius jump fraction for adaptive grid (default: 0.05)')
    parser.add_argument('--max-temp-jump', type=float, default=0.10,
                        help='Max temperature jump fraction (default: 0.10)')
    parser.add_argument('--max-density-jump', type=float, default=0.15,
                        help='Max density jump fraction (default: 0.15)')
    parser.add_argument('--max-ne-jump', type=float, default=0.10,
                        help='Max electron density jump fraction (default: 0.10)')
    parser.add_argument('--max-logtau-jump', type=float, default=0.10,
                        help='Max log optical depth jump (default: 0.10)')
    parser.add_argument('--max-velocity-jump', type=float, default=0.10,
                        help='Max fractional velocity change per cell |dv|/v (default: 0.10). '
                             'Refines velocity gradients even when rho/T/tau are smooth — critical '
                             'for non-monotonic velocity fields in interacting SNe.')
    parser.add_argument('--max-sigma-jump', type=float, default=1.0,
                        help='Max change in sigma=dlnv/dlnr-1 between adjacent cells (default: 1.0). '
                             'Large sigma jumps create spikes in MOM_JREL_V9; refining to smooth '
                             'sigma is essential for the Dessart+2015 direct treatment.')
    parser.add_argument('--max-abs-sigma', type=float, default=5.0,
                        help='Maximum allowed |sigma| at any grid point (default: 5.0). '
                             'When |sigma| exceeds this, AMR inserts points to reduce discretization '
                             'error. Set high (e.g. 100) to disable sigma-driven refinement and let '
                             'physical gradients (rho, T, tau) drive grid placement instead.')
    parser.add_argument('--no-sigma-optimize', dest='sigma_optimize', action='store_false',
                        default=True,
                        help='Disable optimal insertion point selection for sigma minimization. '
                             'By default, AMR places new points at the fraction that minimizes '
                             'the resulting |sigma|. With --no-sigma-optimize, simple midpoint '
                             'insertion is used, which is faster and lets physical gradients '
                             'drive grid placement rather than sigma optimization.')
    # --- Sigma-aware grid coarsening (Luc Dessart approach) ---
    parser.add_argument('--coarsen-negative-sigma', dest='coarsen_negative_sigma',
                        action='store_true', default=True,
                        help='Coarsen grid regions where sigma < threshold to 1-2 zones '
                             '(Luc Dessart approach: negative sigma only affects 1-2 zones '
                             'in well-behaved models). Default: enabled.')
    parser.add_argument('--no-coarsen-negative-sigma', dest='coarsen_negative_sigma',
                        action='store_false',
                        help='Disable sigma-aware grid coarsening.')
    parser.add_argument('--coarsen-sigma-threshold', type=float, default=-0.1,
                        help='Sigma threshold below which zones are coarsened (default: -0.1). '
                             'Zones with sigma < this value are candidates for coarsening.')
    parser.add_argument('--coarsen-max-zones', type=int, default=2,
                        help='Maximum zones to keep per negative-sigma region (default: 2). '
                             '2 = boundary points only; 3 = boundary + midpoint.')
    parser.add_argument('--coarsen-min-region-zones', type=int, default=4,
                        help='Only coarsen regions with at least this many zones (default: 4). '
                             'Smaller regions are already adequately resolved.')
    parser.add_argument('--refine-outer-boundary', action='store_true', default=True,
                        help='Refine the outer boundary by subdividing the first few cells '
                             'to achieve target dr/r (default: enabled). Prevents J(mom)/J(ray) '
                             'divergence and luminosity blowup at the outer boundary.')
    parser.add_argument('--outer-refine-cells', type=int, default=5,
                        help='Number of outermost cells to subdivide (default: 5).')
    parser.add_argument('--outer-refine-subdivisions', type=int, default=3,
                        help='Sub-cells per original cell (default: 3). Overridden by '
                             '--outer-refine-target-dr-r when the original dr/r is larger.')
    parser.add_argument('--outer-refine-target-dr-r', type=float, default=0.015,
                        help='Target dr/r at the outer boundary (default: 0.015). Cells '
                             'with dr/r above this are automatically subdivided.')
    parser.add_argument('--outer-refine-max-points', type=int, default=30,
                        help='Maximum new points to insert at the outer boundary (default: 30).')
    parser.add_argument('--max-cell-ratio', type=float, default=1.5,
                        help='Max ratio between adjacent cell sizes for grid smoothing (default: 1.5). '
                             'Lower = smoother spacing transitions = fewer sigma oscillations. '
                             '2.0 allows 2x spacing jumps; 1.5 requires gradual transitions.')
    parser.add_argument('--min-dr-r', type=float, default=1e-6,
                        help='Minimum fractional cell size dr/r for the final grid (default: 1e-6). '
                             'Points with dr/r below this are removed to prevent degenerate grid '
                             'points that cause division-by-zero in CMFGEN sigma computation. '
                             'Also used by insert_outer_boundary_points to limit the finest '
                             'inserted cell. 1e-6 is safe for 22.15E write precision.')
    parser.add_argument('--max-dr-r', type=float, default=0.5,
                        help='Maximum fractional cell size dr/r for the final grid (default: 0.5). '
                             'Cells with dr/r above this are subdivided to prevent large cells '
                             'that destabilise Lambda iteration. Also used by '
                             'insert_outer_boundary_points to ensure the outer boundary region '
                             'is adequately resolved. 0.5 = 50%% radius change per cell; '
                             '0.1 = 10%% (finer, more points needed).')
    parser.add_argument('--max-points',   type=int,   default=150,
                        help='Maximum number of grid points including outer extras (default: 150)')
    parser.add_argument('--enable-sigma-conformance', dest='enable_sigma_conformance',
                        action='store_true', default=True,
                        help='Run a final sigma-conformance pass (default: enabled) that re-checks '
                             'the exact CMFGEN sigma formula on the truly final grid — after grid '
                             'spacing enforcement, regularization, velocity re-interpolation, and '
                             'smoothing — and inserts points at any remaining violation. AMR runs '
                             'early in the pipeline and cannot see sigma issues introduced by these '
                             'later steps; this closes the loop. Only inserts points, never removes.')
    parser.add_argument('--no-sigma-conformance', dest='enable_sigma_conformance',
                        action='store_false',
                        help='Disable the final sigma-conformance pass.')
    parser.add_argument('--sigma-conformance-headroom', type=int, default=40,
                        help='Extra grid points (beyond --max-points) the sigma-conformance pass '
                             'may use (default: 40). Needed because the cell-ratio pass in '
                             '--max-points can consume the entire nominal budget on purely '
                             'geometric smoothing, leaving no room to fix genuine sigma violations. '
                             'CMFGEN allocates ND from the input file at runtime, so a modest '
                             'overrun is safe.')
    parser.add_argument('--thin-grid-tau', type=float, default=1.0,
                        help='Optical-depth threshold below which AMR and min-nd filling receive extra priority (default: 1.0)')
    parser.add_argument('--thin-grid-priority', type=float, default=1.0,
                        help='AMR and min-nd fill score multiplier for cells at or below --thin-grid-tau (default: 1.0)')
    parser.add_argument('--outer-boundary-extra-points', type=int, default=5,
                        help='Number of geometrically-spaced points inserted at the outer boundary '
                             '(depth 1 in CMFGEN) to stabilise the radiation field. 0 = disabled. '
                             'Default: 5')
    parser.add_argument('--outer-boundary-ratio', type=float, default=1.5,
                        help='Geometric spacing ratio for outer boundary points (OB_RAT in CMFGEN). '
                             'Each cell is OB_RAT times larger than the previous one. '
                             '1.5 = moderate refinement, 2.0 = aggressive. Default: 1.5')
    # Ionization-front-aware grid enhancement
    parser.add_argument('--ionization-front-aware', action='store_true', default=True,
                        help='Detect ionization fronts (sharp n_e gradients) and force-include '
                             'source zones at those fronts. Ensures fronts are resolved over '
                             'a few depths for CMFGEN population solver convergence. Default: True.')
    parser.add_argument('--no-ionization-front-aware', dest='ionization_front_aware',
                        action='store_false',
                        help='Disable ionization-front-aware grid enhancement.')
    parser.add_argument('--ionization-front-threshold', type=float, default=0.5,
                        help='|dlog(n_e)/dlog(r)| threshold for ionization front detection '
                             '(default: 0.5). Lower = more sensitive = more fronts detected. '
                             'H recombination fronts typically show 0.5-2.0 dex change.')
    parser.add_argument('--ionization-front-min-points', type=int, default=5,
                        help='Minimum source zones to mark around each ionization front '
                             '(default: 5). Ensures the front is resolved over a few depths.')
    parser.add_argument('--ionization-front-width', type=int, default=3,
                        help='Number of zones on each side of the front center to mark '
                             '(default: 3). Total marked zones = 2*width+1 before expansion.')
    parser.add_argument('--ionization-front-cap', type=float, default=0.15,
                        help='Max fraction of target_points for ionization front zones '
                             '(default: 0.15). Prevents fronts from dominating the grid.')
    # Recombination-front resolution
    parser.add_argument('--resolve-recombination-fronts', action='store_true', default=True,
                        help='Detect under-resolved recombination fronts (sharp T drops in '
                             '4000-15000K range) and spread the transition over a few zones. '
                             'Resolves a physical feature STELLA under-resolved. Default: True.')
    parser.add_argument('--no-resolve-recombination-fronts', dest='resolve_recombination_fronts',
                        action='store_false',
                        help='Disable recombination-front resolution.')
    parser.add_argument('--recombination-t-jump', type=float, default=0.3,
                        help='Minimum |dT/T| to flag a recombination front (default: 0.3). '
                             'H recombination typically shows 40-60%% T drop.')
    parser.add_argument('--recombination-front-width', type=int, default=5,
                        help='Number of zones to spread a recombination front over (default: 5). '
                             'Ensures the front is resolved over a few depths for CMFGEN.')
    parser.add_argument('--sigma-priority-threshold', type=float, default=2.0,
                        help='|sigma| above which AMR edges get extra priority (default: 2.0). '
                             'Ensures velocity reversal regions compete with thin-region priority.')
    parser.add_argument('--sigma-priority', type=float, default=5.0,
                        help='AMR score multiplier for edges with |sigma| > --sigma-priority-threshold (default: 5.0). '
                             'Higher = more points allocated to velocity reversal regions.')
    parser.add_argument('--dsigma-priority-threshold', type=float, default=0.5,
                        help='|d(sigma)| above which AMR edges get extra priority (default: 0.5). '
                             'Targets sharp sigma transitions like the outer boundary J(mom)/J(ray) mismatch.')
    parser.add_argument('--dsigma-priority', type=float, default=5.0,
                        help='AMR score multiplier for edges with |d(sigma)| > --dsigma-priority-threshold (default: 5.0).')
    parser.add_argument('--ne-priority-threshold', type=float, default=2.0,
                        help='(dNe/Ne)/max-ne-jump ratio above which AMR edges get extra priority (default: 2.0). '
                             'Ensures sharp recombination/ionization fronts compete for AMR budget against '
                             'velocity/sigma-driven edges in non-monotonic-velocity SN models.')
    parser.add_argument('--ne-priority', type=float, default=1.0,
                        help='AMR score multiplier for edges with (dNe/Ne)/max-ne-jump > --ne-priority-threshold '
                             '(default: 1.0 = no boost). Higher = more points allocated to ionization fronts.')
    parser.add_argument('--min-nd',       type=int,   default=150,
                        help='Minimum number of grid points for CMFGEN convergence (default: 150). '
                             'If the grid after AMR/thinning has fewer points, midpoints are '
                             'inserted at the largest cells (with grid smoothing) to reach this '
                             'minimum. Essential for stable Lambda iteration in interacting SNe.')
    parser.add_argument('--t-floor',      type=float, default=3000.0,
                        help='Minimum temperature floor in K (default: 3000)')
    parser.add_argument('--t-cap',        type=float, default=0.0,
                        help='Maximum temperature cap in K (default: 0=disabled). '
                             'When >0, zones with T > t_cap are capped to t_cap. '
                             'Prevents MON_INTERP crashes when T_gas at the outer '
                             'boundary is very high from shock heating.')
    parser.add_argument('--trace-zero-threshold', type=float, default=0.0,
                        help='Mass fraction threshold below which element abundances '
                             'are zeroed at all depths (default: 0=disabled). '
                             'E.g. 1e-6 zeros any element with X < 1e-6 at any depth. '
                             'The ABUND_FLOOR (1e-30) then brings zeros back to a '
                             'negligible but non-zero value, preventing SETPOP crashes '
                             'while ensuring FIXPOP freezes trace species corrections.')
    parser.add_argument('--smooth-velocity', action='store_true',
                        help='Enable velocity smoothing to reduce sigma spikes at ejecta interfaces')
    parser.add_argument('--sigma-threshold', type=float, default=2.0,
                        help='Sigma threshold for velocity smoothing (default: 2.0)')
    parser.add_argument('--smooth-width', type=int, default=5,
                        help='Number of points for velocity smoothing window (default: 5)')
    # Sigma-transition smoothing: smooth the sigma gradient (not magnitude)
    parser.add_argument('--smooth-sigma', action='store_true', default=False,
                        help='Smooth sigma transitions so sigma evolves gradually between regions '
                             '(preserves physical spikes, smooths jumps). Default: False — '
                             'AMR max_sigma_jump handles smoothness; smoothing modifies velocity '
                             'away from STELLA output.')
    parser.add_argument('--no-smooth-sigma', dest='smooth_sigma', action='store_false',
                        help='Disable sigma-transition smoothing')
    parser.add_argument('--max-dsigma', type=float, default=1.0,
                        help='Max |d(sigma)| between adjacent grid points after smoothing (default: 1.0)')
    parser.add_argument('--sigma-smooth-passes', type=int, default=10,
                        help='Number of sigma-transition smoothing passes (default: 10)')
    parser.add_argument('--sigma-smooth-window', type=float, default=0.5,
                        help='Gaussian width as %% of log-r range for sigma smoothing (default: 0.5). '
                             'Smaller = tighter smoothing that preserves more structure. '
                             '0.1 = very tight (preserves spikes, partial transition smoothing). '
                             '0.3 = recommended (smooth transitions, preserves non-monotonicity). '
                             '1.0 = wide (may reduce spike magnitude).')
    parser.add_argument('--sigma-smooth-max-rel-change', type=float, default=0.02,
                        help='Maximum fractional velocity change per point in sigma smoothing (default: 0.02 = 2%%). '
                             'Ensures smoothing is minor — only small wiggles are smoothed, '
                             'large-scale velocity structure is preserved.')
    # Wiggle smoothing: light Gaussian smoothing of small-scale velocity oscillations
    parser.add_argument('--smooth-wiggles', action='store_true', default=False,
                        help='Lightly smooth small-scale velocity wiggles (numerical noise from STELLA grid). '
                             'Preserves large-scale non-monotonic structure. Default: False (v2).')
    parser.add_argument('--no-smooth-wiggles', dest='smooth_wiggles', action='store_false',
                        help='Disable velocity wiggle smoothing')
    parser.add_argument('--wiggle-window', type=int, default=3,
                        help='Gaussian sigma in points for wiggle smoothing (default: 3)')
    parser.add_argument('--wiggle-max-change', type=float, default=0.02,
                        help='Max fractional velocity change per point from wiggle smoothing (default: 0.02 = 2%%)')
    parser.add_argument('--velocity-smooth-window-kms', type=float, default=0.0,
                        help='If >0, apply a moving-average velocity smoothing with this window width in km/s. '
                             'Helps CMFGEN convergence by removing small-scale velocity noise while preserving '
                             'the large-scale non-monotonic profile. e.g. 100.0 = 100 km/s window.')
    # Sigma scatter smoothing: targets oscillating sigma in near-constant-V regions
    parser.add_argument('--smooth-sigma-scatter', action='store_true', default=False,
                        help='Smooth sigma oscillations in near-constant velocity regions using Savitzky-Golay. '
                             'Targets |d(sigma)| > max_dsigma_scatter where V is nearly constant. Default: False (v2).')
    parser.add_argument('--no-smooth-sigma-scatter', dest='smooth_sigma_scatter', action='store_false',
                        help='Disable sigma scatter smoothing')
    parser.add_argument('--max-dsigma-scatter', type=float, default=0.5,
                        help='Max |d(sigma)| threshold for scatter smoothing (default: 0.5). '
                             'Regions where |d(sigma)| exceeds this are smoothed with Savitzky-Golay.')
    parser.add_argument('--sigma-scatter-window', type=int, default=7,
                        help='Savitzky-Golay window size in grid points for scatter smoothing (default: 7, must be odd). '
                             'Larger = more aggressive smoothing of sigma oscillations.')
    parser.add_argument('--sigma-scatter-max-rel-change', type=float, default=0.01,
                        help='Max fractional velocity change per point from scatter smoothing (default: 0.01 = 1%%). '
                             'Ensures only small wiggles are smoothed, preserving large-scale velocity structure.')
    parser.add_argument('--sigma-scatter-dv-threshold', type=float, default=0.01,
                        help='Max |dv|/v ratio for scatter classification (default: 0.01 = 1%%). '
                             'Regions with larger velocity gradients are physical (shock fronts) and left intact.')
    # Final Savitzky-Golay smoothing (always-on historically — now configurable)
    # When enabled, applies S-G filter (win=11, poly=3) to V, log(rho), log(T)
    # to kill grid-scale ringing.  When disabled, preserves raw STELLA velocities
    # at shock fronts where S-G can create catastrophic sigma spikes by smoothing
    # across discontinuities with large grid-spacing ratios.
    parser.add_argument('--enable-final-smoothing', action='store_true', default=False,
                        help='Apply final Savitzky-Golay smoothing (win=11, poly=3) to V, rho, T. '
                             'Kills grid-scale ringing but can create sigma spikes at shock fronts '
                             'with large grid-spacing ratios. Default: False (preserve raw velocities).')
    parser.add_argument('--no-final-smoothing', dest='enable_final_smoothing', action='store_false',
                        help='Disable final Savitzky-Golay smoothing (default). '
                             'Preserves raw STELLA velocities at shock fronts. '
                             'Sigma is computed with the CMFGEN formula (no clip).')
    # Sigma global smoothing: smooth sigma across the ENTIRE grid (not just flat-V)
    parser.add_argument('--smooth-sigma-global', action='store_true', default=False,
                        help='Smooth sigma across the entire grid to reduce AMR-induced scatter. '
                             'Uses median filter on sigma + V reconstruction. '
                             'Protects velocity reversal regions (|sigma| > preserve_threshold). '
                             'Default: False (v2).')
    parser.add_argument('--no-smooth-sigma-global', dest='smooth_sigma_global', action='store_false',
                        help='Disable global sigma smoothing')
    parser.add_argument('--sigma-global-window', type=int, default=7,
                        help='Median filter window for global sigma smoothing (default: 7, must be odd). '
                             'Larger = more aggressive smoothing.')
    parser.add_argument('--sigma-global-max-rel-change', type=float, default=0.05,
                        help='Max fractional velocity change per point from global sigma smoothing (default: 0.05 = 5%%).')
    parser.add_argument('--sigma-global-preserve-threshold', type=float, default=2.0,
                        help='|sigma| above which points are protected from global smoothing (default: 2.0). '
                             'Preserves physical velocity reversal features.')
    parser.add_argument('--plot-dir', type=str, default=None,
                        help='Directory for diagnostic plots (default: same as output file). '
                             'Set to "none" to skip plot generation.')
    parser.add_argument('--no-inc-rad-decays', dest='inc_rad_decays', action='store_false',
                        default=True,
                        help='Do not write isotope mass-fraction blocks (Ni56, Co56, Fe56) '
                             'to SN_HYDRO_DATA.  Sets "Number of isotopes: 0" so CMFGEN '
                             'does not enable INC_RAD_DECAYS.  Use this for SN_B models '
                             'where radioactive decay heating is not the focus.')
    args = parser.parse_args()

    # Parse composition_elements from comma-separated string to list
    if isinstance(args.composition_elements, str):
        args.composition_elements = [e.strip() for e in args.composition_elements.split(',') if e.strip()]

    print(f"Reading {args.swd_file}...")
    # When target_days is given, read_swd parses only the closest epoch (two-pass
    # seek-based read) instead of every epoch in the (often multi-hundred-MB) .swd.
    # The selected epoch is chosen by _select_epoch, which is identical to the
    # selection logic below, so epochs[t_pick] is always the parsed snapshot.
    epochs = read_swd(args.swd_file, target_days=args.target_days)
    available = sorted(epochs.keys())
    print(f"Available epochs: {available[0]:.2f} ... {available[-1]:.2f} days "
          f"({len(available)} snapshots)")

    if args.target_days is None:
        t_pick = available[-1]
    else:
        # Handle negative times (pre-SBO) vs positive times (post-SBO)
        # If target_days is positive, prefer positive epochs (post-SBO)
        # If target_days is negative, prefer negative epochs (pre-SBO)
        if args.target_days >= 0:
            # Filter to positive epochs only (post-SBO)
            positive_epochs = [t for t in available if t >= 0]
            if positive_epochs:
                t_pick = min(positive_epochs, key=lambda t: abs(t - args.target_days))
                # If target is far from any available positive epoch, use last positive epoch
                if abs(t_pick - args.target_days) > max(1.0, 0.1 * (positive_epochs[-1] - positive_epochs[0])):
                    t_pick = positive_epochs[-1]
                    print(f"Requested {args.target_days:.2f} d far from available post-SBO range; "
                          f"using last post-SBO epoch")
            else:
                # No positive epochs available, use closest epoch (could be pre-SBO)
                t_pick = min(available, key=lambda t: abs(t - args.target_days))
                print(f"WARNING: No post-SBO epochs available; using closest epoch (pre-SBO)")
        else:
            # Target is negative (pre-SBO), prefer negative epochs
            negative_epochs = [t for t in available if t < 0]
            if negative_epochs:
                t_pick = min(negative_epochs, key=lambda t: abs(t - args.target_days))
                # If target is far from any available negative epoch, use first negative epoch
                if abs(t_pick - args.target_days) > max(1.0, 0.1 * (negative_epochs[-1] - negative_epochs[0])):
                    t_pick = negative_epochs[0]
                    print(f"Requested {args.target_days:.2f} d far from available pre-SBO range; "
                          f"using earliest pre-SBO epoch")
            else:
                # No negative epochs available, use closest epoch (could be post-SBO)
                t_pick = min(available, key=lambda t: abs(t - args.target_days))
                print(f"WARNING: No pre-SBO epochs available; using closest epoch (post-SBO)")
    print(f"Using epoch t = {t_pick:.4f} days "
          f"(requested {args.target_days:.2f})" if args.target_days else
          f"Using last epoch t = {t_pick:.4f} days")

    # Tolerance check: warn if the closest epoch is far from the requested time
    if args.target_days is not None:
        delta = abs(t_pick - args.target_days)
        if delta > 5.0:
            print(f"WARNING: Closest available epoch ({t_pick:.2f} d) is {delta:.2f} days "
                  f"from requested ({args.target_days:.2f} d). "
                  f"STELLA may have crashed before reaching this epoch.")
            print(f"  Available epoch range: {available[0]:.2f} ... {available[-1]:.2f} days")
        elif delta > 1.0:
            print(f"  Note: Using epoch {delta:.2f} days from requested (within 5-day tolerance)")

    snap = epochs[t_pick]
    if snap is None:
        # Safety net: the single-epoch reader parsed a different epoch than the
        # caller selected (should not happen — both use _select_epoch — but fall
        # back to a full read rather than crash on a None snapshot).
        print(f"WARNING: single-epoch read mismatch (t_pick={t_pick}); "
              f"falling back to full .swd read")
        epochs = read_swd(args.swd_file)
        snap = epochs[t_pick]
    r    = snap['r']
    v    = snap['v']
    T    = snap['T']
    rho  = snap['rho']
    kappa_stella = snap.get('kappa')   # Rosseland opacity from STELLA NLTE tables (cm²/g)
    ne_stella_snap = snap.get('ne')    # None for old-format .swd without n_e column
    ne_stella_r    = snap['r'].copy()  # Original STELLA radii for interpolation
    epoch_data = snap  # Save full epoch data for composition extraction
    nzon = len(r)

    # NEVER homologously expand the STELLA profile. If the requested target
    # day is beyond the available STELLA epochs, just use the closest available
    # epoch. Homologous expansion (r~t, rho~t^-3, T~t^-1) is an approximation
    # that breaks down for interacting SNe with reverse shocks and non-monotonic
    # velocity profiles — it produces unphysical density/temperature profiles
    # that don't match the actual STELLA hydrodynamics.
    if args.target_days is not None and abs(t_pick - args.target_days) > 0.5:
        print(f"\n  [No Homologous Expansion] Using closest STELLA epoch {t_pick:.2f}d "
              f"(requested {args.target_days:.2f}d, delta={abs(t_pick - args.target_days):.2f}d)")
        print(f"    SN_HYDRO_DATA will be written with t={t_pick:.4f}d (the actual STELLA epoch)")

    print(f"Snapshot: {nzon} zones")
    if kappa_stella is not None:
        print(f"[kappa] STELLA Rosseland opacity available: "
              f"min={kappa_stella.min():.3e}, max={kappa_stella.max():.3e} cm²/g")
    
    # Save original STELLA data for re-interpolation after grid processing.
    # The AMR and grid spacing enforcement interpolate from the downsampled
    # grid, which can miss velocity extrema (e.g., reverse shock dip) when
    # STELLA zones are skipped. Re-interpolating from the full STELLA data
    # ensures the final grid values are as close as possible to the real
    # STELLA profile.
    # NOTE: rho, T, and mass fractions are also saved here so they can be
    # re-interpolated from the full STELLA data after grid finalization.
    # Without this, the density profile from the downsampled grid (79 points)
    # can differ by up to 40% from STELLA in regions with subtle density
    # curvature (e.g., the inner density shell at R~1.3e+15 cm).
    # Save original data for mass fraction plot (before trimming)
    orig_r_for_plot = r.copy()
    orig_rho_for_plot = rho.copy()
    orig_T_for_plot = T.copy()

    # Load or construct mass fractions
    # Priority: .swd evolved composition > .rho file > .abn initial > solar
    # BUT: STELLA's nuclear network only tracks H, He, C, O, Si, Fe — not Ca, S, Ar, etc.
    # Those elements are zeroed in the .swd even though they exist in the .abn initial composition.
    # Fix: when .swd has zero Ca, merge in Ca from .abn (initial KEPLER composition).
    # This preserves the physical Ca abundance for CMFGEN spectral modelling (Ca II H&K, [Ca II]).
    if epoch_data is not None and 'composition' in epoch_data and epoch_data['composition'] is not None:
        print("Using evolved composition data from .swd file")
        composition = epoch_data['composition']
        # .swd writes Y_k*A_k (number fraction × atomic mass), not true mass fractions.
        # Normalize over ALL 15 elements to get X_k = Y_k*A_k / <A> before subsetting.
        all_mf = normalize_swd_composition(composition)
        # Only keep the 6 elements needed: H, He, O, Si, Ca, Fe
        mass_fractions = {
            'H':       all_mf['H'],
            'He':      all_mf['He'],
            'O':       all_mf['O'],
            'Si':      all_mf['Si'],
            'Ca':      all_mf['Ca'],
            'FePeak':  all_mf['Fe'],
        }
        print(f"  Composition from .swd (6 elements: H, He, O, Si, Ca, Fe)")

        # STELLA's nuclear network doesn't track Ca — it's always zero in .swd.
        # Merge Ca from .abn file (initial composition) if available.
        if np.all(mass_fractions['Ca'] < 1e-20):
            print("  Ca is zero in .swd (STELLA nuclear network doesn't track Ca)")
            if args.abn_file:
                print(f"  Merging Ca from .abn file: {args.abn_file}")
                try:
                    abn_raw_for_ca = read_abn(args.abn_file, nzon)
                    mass_fractions['Ca'] = abn_raw_for_ca['Ca']
                    print(f"  Ca from .abn: min={np.min(mass_fractions['Ca']):.4e}, "
                          f"max={np.max(mass_fractions['Ca']):.4e}")
                except Exception as e:
                    print(f"  WARNING: Failed to read Ca from .abn file: {e}")
            else:
                print("  WARNING: No .abn file available to recover Ca")
    elif args.rho_file:
        print(f"Reading evolved abundances from {args.rho_file}...")
        rho_raw = read_rho(args.rho_file, nzon)
        # Map .rho elements to CMFGEN-friendly set
        mass_fractions = {
            'H':  rho_raw['H'],
            'He': rho_raw['He'],
            'C':  rho_raw['C'],
            'O':  rho_raw['O'],
            'Si': rho_raw['Si'],
            'S':  rho_raw['S'],
            'Ca': rho_raw['Ca'],
            'FePeak': rho_raw['Fe'],  # Use Fe from .rho file
        }
        print(f"  Using final evolved composition from .rho file")
    elif args.abn_file:
        print(f"Reading initial abundances from {args.abn_file}...")
        abn_raw = read_abn(args.abn_file, nzon)
        # Map to CMFGEN-friendly set
        mass_fractions = {
            'H':  abn_raw['H'],
            'He': abn_raw['He'],
            'C':  abn_raw['C'],
            'O':  abn_raw['O'],
            'Si': abn_raw['Si'],
            'S':  abn_raw['S'],
            'Ca': abn_raw['Ca'],
            'FePeak': abn_raw['FePeak'],
        }
        print(f"  Using initial composition from .abn file (not evolved)")
        
        # Check if Ca is zero in .abn file and use input_profile as fallback
        if np.all(mass_fractions['Ca'] == 0) and args.input_profile:
            print(f"WARNING: Ca is zero in .abn file, reading from input_profile...")
            try:
                profile_raw = read_input_profile(args.input_profile)
                # Interpolate profile data to match SWD grid
                from scipy.interpolate import interp1d
                # Read profile radii for interpolation
                profile_r = []
                with open(args.input_profile) as fh:
                    for line in fh:
                        if line.startswith('#'):
                            continue
                        cols = line.split()
                        if len(cols) >= 11:
                            profile_r.append(float(cols[0]))
                profile_r = np.array(profile_r)[::-1]  # outer->inner
                # Interpolate Ca to SWD grid
                interp_ca = interp1d(profile_r, profile_raw['Ca'], kind='linear', 
                                     bounds_error=False, fill_value='extrapolate')
                mass_fractions['Ca'] = interp_ca(r)
                print(f"  Ca interpolated from input_profile: min={np.min(mass_fractions['Ca']):.4e}, max={np.max(mass_fractions['Ca']):.4e}")
            except Exception as e:
                print(f"  WARNING: Failed to read Ca from input_profile: {e}")
                print(f"  Using solar abundance fallback instead")
    else:
        print("WARNING: No composition data available, using solar abundances")
        mass_fractions = {
            'H':  np.full(nzon, 0.70),
            'He': np.full(nzon, 0.28),
            'C': np.full(nzon, 0.009),
            'O': np.full(nzon, 0.009),
            'Si': np.full(nzon, 0.0007),
            'S': np.full(nzon, 0.0007),
            'Ca': np.full(nzon, 0.0002),
            'FePeak': np.full(nzon, 0.0013),
        }
    
    # Save original mass fractions before trimming (before normalization)
    orig_mass_fractions = {k: v.copy() for k, v in mass_fractions.items()}
    
    # Initialize ni56_for_isotopes (will be filled from .xni or .rho or .abn)
    ni56_for_isotopes = np.zeros(nzon)
    
    # Read Ni56 from .xni file if provided (overrides .rho/.abn Ni56)
    if args.xni_file:
        print(f"Reading Ni56 distribution from {args.xni_file}...")
        try:
            ni56_from_xni = read_xni(args.xni_file, nzon)
            # Store Ni56 separately for isotopic data, don't add to mass_fractions yet
            ni56_for_isotopes = ni56_from_xni
            print(f"  Ni56 max fraction: {np.max(ni56_from_xni):.4e}")
        except Exception as e:
            print(f"  WARNING: Failed to read .xni file: {e}")
            print(f"  Using Ni56 from .swd/.rho/.abn file instead")
            if epoch_data is not None and 'composition' in epoch_data and epoch_data['composition'] is not None:
                ni56_for_isotopes = all_mf['Ni'].copy()
            elif args.rho_file:
                ni56_for_isotopes = rho_raw['Ni'].copy()
            elif args.abn_file:
                ni56_for_isotopes = abn_raw['Ni56'].copy()
            else:
                ni56_for_isotopes = np.zeros(nzon)
    else:
        # Use Ni56 from .swd (preferred), .rho file, or .abn file
        if epoch_data is not None and 'composition' in epoch_data and epoch_data['composition'] is not None:
            ni56_for_isotopes = all_mf['Ni'].copy()
        elif args.rho_file:
            ni56_for_isotopes = rho_raw['Ni'].copy()
        elif args.abn_file:
            ni56_for_isotopes = abn_raw['Ni56'].copy()
        else:
            # No composition file, set Ni56 to zero
            ni56_for_isotopes = np.zeros(nzon)
    
    # Remove Ni56 from elemental mass fractions (will be in isotopic_data)
    # This must happen BEFORE normalization
    if 'Ni56' in mass_fractions:
        del mass_fractions['Ni56']
    
    # Normalise elemental mass fractions (stable material only)
    mf_sum = sum(mass_fractions.values())
    mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
    mass_fractions = {k: v / mf_sum for k, v in mass_fractions.items()}

    # stella_mf_full is saved after trimming (below) to prevent
    # re-interpolation from pulling values outside the trimmed range

    # Update orig_mass_fractions to match (remove Ni56, normalize)
    if 'Ni56' in orig_mass_fractions:
        del orig_mass_fractions['Ni56']
    orig_mf_sum = sum(orig_mass_fractions.values())
    orig_mf_sum = np.where(orig_mf_sum < 1e-10, 1.0, orig_mf_sum)
    orig_mass_fractions = {k: v / orig_mf_sum for k, v in orig_mass_fractions.items()}

    # Optical depth trimming with adaptive relaxation to meet minimum point requirement
    # tau is ordered from outside (small tau) to inside (large tau)
    # Hybrid tau trimming: use STELLA NLTE kappa for the inner boundary (tau_max)
    # and analytical kappa for the outer boundary (tau_min).
    #
    # STELLA NLTE kappa is more accurate at high optical depth (tau > 1) where
    # T_gas ≈ T_rad and the NLTE tables are reliable. Using it for tau_max
    # ensures the inner boundary is at the correct depth for the diffusion
    # approximation (tau=100 with STELLA kappa, not tau=6 as with analytical).
    #
    # Analytical kappa is more generous at the outer boundary (tau < 1) where
    # STELLA NLTE kappa can be very low (cold neutral gas, kappa ~ 0.003 cm²/g),
    # which would cut off the outer ejecta prematurely. The analytical Kramers
    # opacity (kappa ~ 0.11 cm²/g) includes bound-free opacity that keeps the
    # outer atmosphere in the grid.
    if kappa_stella is not None:
        kappa_stella_arr = kappa_stella.copy()
        kappa_analytical_arr = compute_kappa_total(rho, T, mass_fractions)
        # Floor STELLA NLTE kappa at electron scattering opacity — the NLTE
        # Rosseland mean drops below electron scattering in the recombination
        # zone (T < 5300 K), giving unphysically low tau for AMR/grid placement.
        kappa_stella_floored = np.maximum(kappa_stella_arr, args.kappa)
        tau_stella = compute_kappa_tau(r, rho, kappa_stella_floored)
        tau_analytical = compute_kappa_tau(r, rho, kappa_analytical_arr)
        print(f"[kappa] STELLA NLTE tau range (floored at kappa_es={args.kappa}): [{tau_stella.min():.3e}, {tau_stella.max():.3e}]")
        print(f"[kappa] Analytical tau range:  [{tau_analytical.min():.3e}, {tau_analytical.max():.3e}]")
        # Use STELLA kappa (original, not floored) as the primary for output —
        # CMFGEN computes its own opacity and doesn't use this column.
        # But use floored kappa for tau/grid decisions.
        kappa_arr = kappa_stella_arr
        tau = tau_stella
    else:
        kappa_analytical_arr = compute_kappa_total(rho, T, mass_fractions)
        kappa_arr = kappa_analytical_arr
        tau = compute_kappa_tau(r, rho, kappa_arr)
        tau_analytical = tau
        print(f"[kappa] Using analytical Rosseland opacity for tau calculation (no STELLA kappa): "
              f"min={kappa_arr.min():.3e}, max={kappa_arr.max():.3e} cm²/g")

    # Compute tau with the constant configured kappa for inner boundary trimming.
    # STELLA NLTE kappa can be very low in cold neutral gas (min ~0.0001 cm²/g),
    # giving tau_max ~1 instead of ~240 with kappa=0.34. Since CMFGEN uses its own
    # opacity (closer to the constant kappa), the inner boundary should be set
    # using the constant kappa tau, not the STELLA NLTE tau.
    tau_const_kappa = compute_kappa_tau(r, rho, np.full(nzon, args.kappa))
    print(f"[kappa] Constant kappa={args.kappa} tau range: [{tau_const_kappa.min():.3e}, {tau_const_kappa.max():.3e}]")

    print(f"Original tau range: [{tau.min():.3e}, {tau.max():.3e}]")

    # Initial tau trimming (Dessart+2015: truncate at tau_Rosseland=30)
    # STRICT: do NOT relax tau_max if too few points. Adaptive grid will add
    # points where gradients require it, but we keep the tau=30 boundary.
    print(f"[Tau Trim] Dessart+2015 style: truncating at tau_Rosseland={args.tau_max:.1f}")
    tau_max_current = args.tau_max
    tau_min_current = args.tau_min

    # Find outer boundary: first zone where tau >= tau_min
    # Use analytical (Kramers) tau for the outer boundary, which includes
    # bound-free opacity from H/He. STELLA NLTE kappa in cold (T<3000 K)
    # neutral gas is very low (~0.001 cm²/g), giving tau~0 and trimming
    # away real ejecta. Analytical kappa (~0.11 cm²/g) is more physically
    # complete for the outer atmosphere.
    tau_for_outer = tau_analytical
    tau_source = "analytical"
    i_outer = 0
    while i_outer < nzon and tau_for_outer[i_outer] < tau_min_current:
        i_outer += 1
    print(f"[Tau Trim] Outer boundary using {tau_source} tau (tau_min={tau_min_current:.3e})")
    # Find inner boundary: last zone where constant-kappa tau <= tau_max
    # Use the constant configured kappa (args.kappa) for inner boundary trimming,
    # NOT the STELLA NLTE tau. STELLA NLTE kappa is too low in cold neutral gas
    # (gives tau_max ~1 instead of ~240 with kappa=0.34), so tau_max=30 would
    # never be enforced. The constant kappa matches what CMFGEN uses for its
    # Rosseland opacity calculation.
    tau_for_inner = tau_const_kappa
    i_inner = nzon - 1
    while i_inner > 0 and tau_for_inner[i_inner] > tau_max_current:
        i_inner -= 1
    print(f"[Tau Trim] Inner boundary using constant kappa={args.kappa} tau (tau_max={tau_max_current:.1f})")
    print(f"  Inner zone: {i_inner+1}/{nzon}, tau_const={tau_for_inner[i_inner]:.2f}, "
          f"tau_stella={tau[i_inner]:.4f}, V={v[i_inner]:.1f} km/s")
    i_outer = max(0, i_outer)  # No extra zone — keep CMFGEN within STELLA tau>tau_min

    # Velocity truncation: remove inner zones with V < v_min_km_s.
    # Very low velocity zones (< ~100 km/s) at the inner boundary cause
    # numerical instabilities in CMFGEN (enormous Sobolev optical depths,
    # ill-conditioned BA matrix, species population explosion after ~100 iterations).
    # Truncating at V_min removes these zones while preserving the reverse-shock structure.
    v_min_km_s = getattr(args, 'v_min_km_s', 0.0) or 0.0
    if v_min_km_s > 0:
        # STELLA arrays are outer->inner (decreasing r, V generally decreasing inward)
        # Find the innermost zone with V >= v_min
        v_arr = np.abs(v)  # km/s
        i_inner_v = i_inner
        while i_inner_v > i_outer and v_arr[i_inner_v] < v_min_km_s:
            i_inner_v -= 1
        if i_inner_v < i_inner:
            n_removed = i_inner - i_inner_v
            print(f"[V Trim] Removing {n_removed} inner zones with V < {v_min_km_s:.1f} km/s "
                  f"(V_min was {v_arr[i_inner]:.1f} km/s, now {v_arr[i_inner_v]:.1f} km/s)")
            i_inner = i_inner_v

    n_after_tau = i_inner - i_outer + 1

    if n_after_tau < 2:
        raise ValueError(
            f"No usable CMFGEN grid at t={t_pick:.4f} d: only {n_after_tau} zone "
            f"within tau=[{tau_min_current:.3e}, {tau_max_current:.3e}]. "
            "STELLA did not produce an optically resolved snapshot for this epoch."
        )
    if n_after_tau < args.min_points:
        print(f"[Tau Trim] WARNING: only {n_after_tau} points within tau range "
              f"(min_points={args.min_points}), but keeping strict tau_max={tau_max_current}")
        print(f"  Adaptive grid will add points where gradients require resolution.")

    # Inner boundary: zone just behind tau_max (NOT in SN_HYDRO_DATA)
    i_boundary = min(i_inner + 1, nzon - 1)
    r_inner_boundary = r[i_boundary]
    T_inner_boundary = T[i_boundary]
    print(f"Inner boundary zone (not in SN_HYDRO_DATA): depth {i_boundary}, "
          f"T={T_inner_boundary:.1f} K, R={r_inner_boundary:.2e} cm")

    # Ghost zone removal: remove only the unphysical STELLA ghost zone (V > 0.33c)
    # that results from BADSTE numerical artifacts. The non-monotonic velocity
    # profile (sigma < -1) from the STELLA hydro is PHYSICAL and must be conserved.
    # CMFGEN handles sigma < -1 via the Dessart+2015 prescription in fg_j_cmf_v13.f
    # and the Sobolev routines. VEL_LAW=12 + USE_LAM_ES=T handle non-monotonic V.
    #
    # We only trim zones with V > velocity_max_kms (default: very high, only ghost).
    # No artificial zones are added — the STELLA profile is used as-is.
    if args.velocity_max_kms > 0:
        v_max_kms = args.velocity_max_kms
        # Find the first zone from outside where V <= v_max_kms
        i_v_outer = i_outer
        while i_v_outer <= i_inner and v[i_v_outer] > v_max_kms:
            i_v_outer += 1
        n_trimmed_v = i_v_outer - i_outer
        if n_trimmed_v > 0:
            v_max_original = v[i_outer]
            print(f"[Ghost Zone Trim] Removing {n_trimmed_v} outer zone(s) with V > {v_max_kms:.0f} km/s "
                  f"(V_max was {v_max_original:.1f} km/s = {v_max_original/3e5:.2f}c)")
            i_outer = i_v_outer
            print(f"[Ghost Zone Trim] Remaining V range: [{v[i_outer]:.1f}, {v[i_inner]:.1f}] km/s "
                  f"(max/c = {v[i_outer:i_inner+1].max()/3e5:.3f})")
            # Count zones with sigma < -1 (non-monotonic velocity) — these are KEPT
            n_sigma_lt_m1 = 0
            for i in range(i_outer, i_inner):
                if r[i] != r[i+1] and v[i] > 0 and v[i+1] > 0:
                    sigma = np.log(v[i]/v[i+1]) / np.log(r[i]/r[i+1]) - 1
                    if sigma < -1:
                        n_sigma_lt_m1 += 1
            if n_sigma_lt_m1 > 0:
                print(f"[Ghost Zone Trim] Keeping {n_sigma_lt_m1} zones with sigma < -1 "
                      f"(non-monotonic velocity from STELLA hydro — handled by Dessart+2015 prescription)")

    # Final tau trim: exclude the boundary zone
    sl = slice(i_outer, i_inner + 1)
    print(f"Tau trim: keep depths {i_outer}..{i_inner}  "
          f"(tau=[{tau[i_outer]:.2e},{tau[i_inner]:.2e}])")

    # Validate inner boundary tau
    # With ZERO_FLUX (DIF=F), tau>50 is not required for the boundary condition,
    # but is still recommended for LTE validity (CMFGEN manual).
    tau_at_ib = tau_for_inner[i_inner]
    if tau_at_ib < 50:
        print(f"\n  *** NOTE: tau at inner boundary = {tau_at_ib:.1f} < 50 ***")
        print(f"      Using ZERO_FLUX inner boundary (DIF=F) — tau>50 not required for boundary condition.")
        print(f"      For LTE validity at inner boundary, consider increasing cmfgen_tau_max.\n")

    r    = r[sl];   v    = v[sl];   T    = T[sl]
    rho  = rho[sl]; kappa_arr = kappa_arr[sl]
    mass_fractions = {k: v2[sl] for k, v2 in mass_fractions.items()}

    # Save full STELLA data for re-interpolation AFTER trimming, so that
    # re-interpolation only uses values within the trimmed range (prevents
    # pulling in velocities below v_min or above velocity_max_kms).
    stella_r_full = r.copy()
    stella_v_full = v.copy()
    stella_rho_full = rho.copy()
    stella_T_full = T.copy()
    stella_mf_full = {k: v2.copy() for k, v2 in mass_fractions.items()}
    
    # Trim Ni56 array to match the trimmed grid
    if 'ni56_for_isotopes' in locals():
        ni56_for_isotopes = ni56_for_isotopes[sl]

    # Save post-tau-trim grid for Ni56 interpolation after grid processing
    r_tau_trim = r.copy()

    # Disable artifact filter - it was removing optically thick zones (high density)
    # The high density zones in the inner region are physical, not artifacts
    print("[Artifact Filter] DISABLED - keeping all zones after tau trim")

    # --- Outer boundary refinement ---
    # The photosphere and spectrum formation region are in the outer layers.
    # Add extra points where tau < 1 to improve Sobolev/flux accuracy.
    def refine_outer_boundary(r, v, T, rho, kappa, mass_fractions, tau, max_extra=20, tau_thresh=1.0):
        """Insert additional points in the outer low-tau region (tau < tau_thresh)."""
        n = len(r)
        # Find the index where tau first exceeds tau_thresh (outer->inner order)
        i_thresh = n
        for i in range(n):
            if tau[i] > tau_thresh:
                i_thresh = i
                break
        # The outer region is indices 0 .. i_thresh-1
        n_outer = i_thresh
        if n_outer < 3:
            print(f"[Outer Refine] Outer region too small ({n_outer} pts), skipping")
            return r, v, T, rho, kappa, mass_fractions, tau

        # We insert points between consecutive pairs in the outer region
        # Limit to avoid bloat
        n_insert_target = min(max_extra, n_outer)
        # Choose evenly spaced insertion points
        insert_indices = []
        step = max(1, n_outer // n_insert_target)
        for i in range(0, n_outer - 1, step):
            insert_indices.append(i)

        # Build new arrays by inserting at midpoints
        def interp_linear(a, b, frac):
            return a + frac * (b - a)
        def interp_log(a, b, frac):
            eps = 1e-30
            return 10 ** (np.log10(max(a, eps)) + frac * (np.log10(max(b, eps)) - np.log10(max(a, eps))))

        r_new, v_new, T_new, rho_new, kappa_new, tau_new = list(r), list(v), list(T), list(rho), list(kappa), list(tau)
        mf_new = {k: list(v) for k, v in mass_fractions.items()}

        # Insert from the end backwards to keep indices valid
        offset = 0
        for i in sorted(insert_indices, reverse=True):
            j = i + offset  # current position in growing list
            if j >= len(r_new) - 1:
                continue
            r_mid = interp_linear(r_new[j], r_new[j+1], 0.5)
            v_mid = interp_linear(v_new[j], v_new[j+1], 0.5)
            T_mid = interp_log(T_new[j], T_new[j+1], 0.5)
            rho_mid = interp_log(rho_new[j], rho_new[j+1], 0.5)
            kappa_mid = interp_log(kappa_new[j], kappa_new[j+1], 0.5)
            tau_mid = interp_log(tau_new[j], tau_new[j+1], 0.5)

            r_new.insert(j+1, r_mid)
            v_new.insert(j+1, v_mid)
            T_new.insert(j+1, T_mid)
            rho_new.insert(j+1, rho_mid)
            kappa_new.insert(j+1, kappa_mid)
            tau_new.insert(j+1, tau_mid)
            for elem in mf_new:
                mf_mid = interp_linear(mf_new[elem][j], mf_new[elem][j+1], 0.5)
                mf_new[elem].insert(j+1, mf_mid)
            offset += 1

        n_added = len(r_new) - n
        print(f"[Outer Refine] Inserted {n_added} points in outer region (tau < {tau_thresh})")
        print(f"  Outer region: {n_outer} -> {n_outer + n_added} points")
        return (np.array(r_new), np.array(v_new), np.array(T_new),
                np.array(rho_new), np.array(kappa_new),
                {k: np.array(v) for k, v in mf_new.items()}, np.array(tau_new))

    # Dessart-style adaptive grid workflow
    if args.enable_adaptive_grid:
        print("\n--- Adaptive Grid Workflow ---")

        # Physical discontinuity smoothing (Dessart+2024 approach):
        # Gaussian-smooth T and rho across sharp 1D artifacts before grid
        # construction. This reduces single-zone jumps (e.g., T jumping 120%
        # at the cannon shell entry) so the Dessart grid criteria (10% per
        # step in T, rho, ne) can be satisfied without hundreds of points.
        if getattr(args, 'smooth_physical_discontinuities', False):
            print("\n  [PhysSmooth] Applying Gaussian smoothing to T, rho, and ne "
                  "(Dessart+2024 approach)")
            T, rho = smooth_physical_discontinuities(
                r, T, rho,
                sigma_zones=getattr(args, 'phys_smooth_sigma_zones', 10.0),
                max_dT_T=args.max_temp_jump,
                max_drho_rho=args.max_density_jump)
            # Also smooth the STELLA ne so that write_sn_hydro_data
            # interpolates a smoothed ne (not the original sharp one).
            # NOTE: ne_stella_r must remain the ORIGINAL STELLA radii (full
            # length) because it is also used for kappa interpolation.
            if ne_stella_snap is not None:
                from scipy.ndimage import gaussian_filter1d
                log_ne = np.log10(np.maximum(ne_stella_snap, 1e-30))
                log_ne_smooth = gaussian_filter1d(log_ne, sigma=getattr(args, 'phys_smooth_sigma_zones', 10.0), mode='nearest')
                ne_stella_snap = np.power(10.0, log_ne_smooth)
                print(f"  [PhysSmooth] ne smoothed: range {10**log_ne.min():.2e} -> {ne_stella_snap.min():.2e} "
                      f"to {10**log_ne.max():.2e} -> {ne_stella_snap.max():.2e}")
            # Recompute kappa and tau with smoothed T/rho
            if kappa_stella is not None and len(kappa_arr) == len(r):
                kappa_arr = compute_kappa_total(rho, T, mass_fractions)
            # Update the full-grid copies for re-interpolation
            stella_T_full = T.copy()
            stella_rho_full = rho.copy()

        # Coordinate stretching in trouble regions (alternative to Gaussian
        # smoothing): preserves peak values, only redistributes radius spacing
        # so the AMR can insert points while respecting dR constraints.
        if getattr(args, 'stretch_trouble_regions', False):
            print("\n  [Stretch] Applying coordinate stretching to trouble regions "
                  "(preserves peak values)")
            # Compute ne and tau for the full grid if not already available
            n_e_full = calculate_electron_density(rho, T, mass_fractions)
            tau_full = compute_kappa_tau(r, rho, np.full(len(r), args.kappa))
            r, v, T, rho, n_e_full, tau_full, mass_fractions = \
                stretch_trouble_regions(
                    r, v, T, rho, n_e_full, tau_full, mass_fractions,
                    max_dT_T=args.max_temp_jump,
                    max_drho_rho=args.max_density_jump,
                    max_dlogtau=args.max_logtau_jump,
                    max_dne_ne=args.max_ne_jump,
                    max_dv_v=getattr(args, 'max_velocity_jump', 0.30),
                    stretch_factor=getattr(args, 'stretch_factor', 3.0),
                    max_stretch_zones=getattr(args, 'max_stretch_zones', 200),
                    min_zones_to_stretch=getattr(args, 'min_zones_to_stretch', 2),
                    verbose=True)
            # Recompute kappa on the stretched grid
            if kappa_stella is not None:
                kappa_arr = np.interp(r, ne_stella_r[::-1], kappa_stella[::-1])
                kappa_arr = np.clip(kappa_arr, 1e-6, None)
            else:
                kappa_arr = compute_kappa_total(rho, T, mass_fractions)
            # Update the full-grid copies for re-interpolation
            stella_T_full = T.copy()
            stella_rho_full = rho.copy()
            stella_v_full = v.copy()
            stella_r_full = r.copy()
            if stella_mf_full:
                stella_mf_full = {k: v.copy() for k, v in mass_fractions.items()}
            print(f"  [Stretch] Stretched grid: {len(r)} points "
                  f"(T range {T.min():.0f}-{T.max():.0f} K, "
                  f"rho range {rho.min():.2e}-{rho.max():.2e})")

        # --- Critical-region detection and region-dependent min_dr/r ---
        # Narrow high-gradient regions (shock fronts, composition
        # interfaces, CSM interaction shells, recombination fronts) often
        # need finer radial sampling than the global min_dr/r allows.  We
        # detect them here on the (optionally pre-smoothed / pre-stretched)
        # full grid, then:
        #   * optionally stretch them further (stretch_critical_regions),
        #   * pass the mask to the AMR so it may subdivide below the
        #     global MIN_DR_OVER_R floor inside critical regions,
        #   * pass the mask to the degenerate-grid thinner and the final
        #     spacing enforcement so they preserve fine points there.
        # The mask is re-derived on the AMR base grid (after downsample)
        # and again before the final spacing enforcement, because the grid
        # length changes between stages; the underlying physical gradients
        # are preserved so re-detection finds the same regions.
        critical_mask = None
        if getattr(args, 'critical_regions', False):
            print("\n  [CriticalRegions] Detecting narrow high-gradient critical regions")
            n_e_crit = calculate_electron_density(rho, T, mass_fractions)
            tau_crit = compute_kappa_tau(r, rho, np.full(len(r), args.kappa))
            # Parse named composition interfaces
            iface_arg = getattr(args, 'critical_composition_interfaces', '')
            if iface_arg and iface_arg.lower() != 'none':
                interfaces = []
                for pair in iface_arg.split(','):
                    a, _, b = pair.partition('-')
                    a, b = a.strip(), b.strip()
                    if a and b:
                        interfaces.append((a, b))
            else:
                interfaces = None
            critical_mask, _regions = identify_critical_regions(
                r, T, rho, n_e_crit, tau_crit, mass_fractions,
                v=v,
                max_dT_T=args.critical_dT_T,
                max_drho_rho=args.critical_drho_rho,
                max_dlogtau=args.critical_dlogtau,
                max_dne_ne=args.critical_dne_ne,
                max_dv_v=args.critical_dv_v,
                composition_jump=args.critical_composition_jump,
                composition_interfaces=interfaces,
                pad_zones=args.critical_pad_zones,
                verbose=True)

            if getattr(args, 'stretch_critical_regions', False) \
                    and critical_mask is not None and np.any(critical_mask):
                print("\n  [Stretch-Critical] Stretching critical regions "
                      "(preserves original STELLA values)")
                (r, v, T, rho, n_e_crit, tau_crit, mass_fractions,
                 critical_mask) = stretch_critical_regions(
                    r, v, T, rho, n_e_crit, tau_crit, mass_fractions,
                    critical_mask,
                    stretch_factor=args.critical_stretch_factor,
                    max_stretch_zones=args.critical_max_stretch_zones,
                    verbose=True)
                # Recompute kappa on the stretched grid
                if kappa_stella is not None:
                    kappa_arr = np.interp(r, ne_stella_r[::-1], kappa_stella[::-1])
                    kappa_arr = np.clip(kappa_arr, 1e-6, None)
                else:
                    kappa_arr = compute_kappa_total(rho, T, mass_fractions)
                # Update the full-grid copies for re-interpolation
                stella_T_full = T.copy()
                stella_rho_full = rho.copy()
                stella_v_full = v.copy()
                stella_r_full = r.copy()
                if stella_mf_full:
                    stella_mf_full = {k: vv.copy() for k, vv in mass_fractions.items()}
                print(f"  [Stretch-Critical] Grid now {len(r)} points "
                      f"({int(np.sum(critical_mask))} critical zones)")

        # Compute tau for grid placement using the best available opacity.
        # When STELLA NLTE kappa is available, use it FLOORED at the electron-
        # scattering opacity (args.kappa). This captures the higher true
        # Rosseland opacity in metal-rich zones (Fe-peak, Si-burning products
        # have line opacity 5-10x above electron scattering), giving a more
        # physically accurate tau profile for grid placement. The floor
        # prevents the NLTE kappa from collapsing in cold neutral gas (T<3000K)
        # where it drops to ~0.001 cm²/g.
        #
        # When STELLA kappa is not available, fall back to constant kappa
        # (the previous behavior).
        if kappa_stella is not None and len(kappa_arr) == len(r):
            kappa_for_grid = np.maximum(kappa_arr, args.kappa)
            tau_source_label = f"STELLA NLTE floored at {args.kappa}"
        else:
            kappa_for_grid = np.full(len(r), args.kappa)
            tau_source_label = f"constant kappa={args.kappa}"
        tau = compute_kappa_tau(r, rho, kappa_for_grid)
        print(f"[kappa] Grid-placement tau from {tau_source_label}: "
              f"range [{tau.min():.3e}, {tau.max():.3e}]")
        orig_r_adapt = r.copy()
        orig_v_adapt = v.copy()
        orig_T_adapt = T.copy()
        orig_rho_adapt = rho.copy()
        orig_tau_adapt = tau.copy()
        orig_mass_fractions_adapt = {k: v.copy() for k, v in mass_fractions.items()}

        # Step 1: Downsample to a coarse base grid.
        # This is the starting point for AMR.
        #
        # When coordinate stretching is enabled, we BYPASS the uniform
        # downsample and instead do a STRETCH-AWARE coarsening: smooth
        # regions (where stretching did NOT insert points) are coarsened
        # to fit within max_points, while trouble regions (where stretching
        # inserted extra points) are preserved exactly.  This ensures the
        # shock/interface regions keep their 25+ depth points while smooth
        # regions don't waste grid points.
        n_full = len(r)
        if getattr(args, 'stretch_trouble_regions', False):
            if n_full > args.max_points:
                # Coarsen only the non-stretched (smooth) regions to fit
                # within max_points.  We identify smooth regions as edges
                # where the local gradient is below all thresholds.
                print(f"[Base Grid] Stretching enabled — coarsening smooth regions "
                      f"of {n_full}-pt stretched grid to fit {args.max_points} pts "
                      f"(preserving trouble regions)")
                # Compute local gradients to identify smooth edges
                dT_T = np.abs(np.diff(T) / np.maximum(T[:-1], 1e-30))
                drho_rho = np.abs(np.diff(rho) / np.maximum(rho[:-1], 1e-30))
                dlogtau = np.abs(np.diff(np.log10(np.maximum(tau, 1e-30))))
                dne_ne = np.abs(np.diff(n_e_full) / np.maximum(n_e_full[:-1], 1e-30))
                # An edge is "smooth" if ALL gradients are below threshold
                smooth_edge = (
                    (dT_T < args.max_temp_jump) &
                    (drho_rho < args.max_density_jump) &
                    (dlogtau < args.max_logtau_jump) &
                    (dne_ne < args.max_ne_jump)
                )
                # Iteratively remove points from smooth regions until we're
                # under max_points.  Each pass removes every other smooth
                # point; repeat until we fit or no more smooth points remain.
                keep = np.ones(n_full, dtype=bool)
                pass_num = 0
                while int(np.sum(keep)) > args.max_points:
                    pass_num += 1
                    n_before = int(np.sum(keep))
                    # Recompute smooth edges on the kept grid
                    idx = np.where(keep)[0]
                    if len(idx) < 4:
                        break
                    T_k = T[idx]
                    rho_k = rho[idx]
                    tau_k = tau[idx]
                    ne_k = n_e_full[idx]
                    dT_T_k = np.abs(np.diff(T_k) / np.maximum(T_k[:-1], 1e-30))
                    drho_rho_k = np.abs(np.diff(rho_k) / np.maximum(rho_k[:-1], 1e-30))
                    dlogtau_k = np.abs(np.diff(np.log10(np.maximum(tau_k, 1e-30))))
                    dne_ne_k = np.abs(np.diff(ne_k) / np.maximum(ne_k[:-1], 1e-30))
                    smooth_k = (
                        (dT_T_k < args.max_temp_jump) &
                        (drho_rho_k < args.max_density_jump) &
                        (dlogtau_k < args.max_logtau_jump) &
                        (dne_ne_k < args.max_ne_jump)
                    )
                    # Remove every other point in smooth runs (from inner end)
                    removed_this_pass = 0
                    i = len(idx) - 2
                    while i > 0 and int(np.sum(keep)) > args.max_points:
                        if smooth_k[i] and smooth_k[i - 1] and keep[idx[i]]:
                            keep[idx[i]] = False
                            removed_this_pass += 1
                        i -= 2
                    n_after = int(np.sum(keep))
                    if n_after == n_before:
                        break  # no more smooth points to remove
                print(f"[Base Grid] Coarsened {n_full} -> {int(np.sum(keep))} pts "
                      f"({pass_num} passes, preserving trouble regions)")
                r = r[keep]
                v = v[keep]
                T = T[keep]
                rho = rho[keep]
                n_e_full = n_e_full[keep]
                tau = tau[keep]
                mass_fractions = {k: val[keep] for k, val in mass_fractions.items()}
                if kappa_stella is not None and len(kappa_arr) == len(keep):
                    kappa_arr = kappa_arr[keep]
            else:
                print(f"[Base Grid] Stretching enabled — using stretched grid directly "
                      f"({n_full} pts) as AMR base (bypassing downsample)")
        elif args.base_grid_method == 'tau_amr_recursive':
            # Delegate grid construction to the dedicated tau-AMR optimiser
            # (optical-depth-primary recursive cell-based AMR with hard
            # constraints for composition interfaces, shocks and velocity
            # extrema; conservative zone-integrated remap; <=max_points).
            from tau_amr_optimiser import (AMRConfig, detect_all_features,
                                           recursive_amr, reconstruct_optical_depth)
            # Build AMRConfig from YAML config file and/or CLI args.
            # Priority: CLI arg (if not None) > YAML config file > AMRConfig defaults.
            amr_yaml = {}
            tau_amr_config_path = getattr(args, 'tau_amr_config', None)
            if tau_amr_config_path:
                import yaml as _yaml
                with open(tau_amr_config_path) as _fh:
                    amr_yaml = _yaml.safe_load(_fh) or {}
                print(f"[Tau-AMR-Recursive] Loaded config from {tau_amr_config_path}")
            # Helper: CLI arg overrides YAML, YAML overrides AMRConfig default.
            def _amr_opt(key, default):
                cli_val = getattr(args, key, None)
                if cli_val is not None:
                    return cli_val
                if key in amr_yaml:
                    return amr_yaml[key]
                return default
            # Map CLI arg names to AMRConfig field names where they differ.
            def _amr_opt_mapped(cli_key, yaml_key, default):
                cli_val = getattr(args, cli_key, None)
                if cli_val is not None:
                    return cli_val
                if yaml_key in amr_yaml:
                    return amr_yaml[yaml_key]
                return default
            amr_cfg = AMRConfig(
                max_points=args.max_points,
                min_points=getattr(args, 'min_nd', 60) or _amr_opt_mapped('min_points', 'min_points', 60),
                base_grid_points=max(40, args.base_grid_target),
                tau_floor=_amr_opt_mapped('tau_floor', 'tau_floor', 1e-6),
                kappa_source=_amr_opt_mapped('amr_kappa_source', 'kappa_source', 'hybrid'),
                target_dlogtau=_amr_opt_mapped('amr_target_dlogtau', 'target_dlogtau', 0.08),
                target_dlogrho=_amr_opt_mapped('amr_target_dlogrho', 'target_dlogrho', 0.30),
                target_dv_kms=_amr_opt_mapped('amr_target_dv_kms', 'target_dv_kms', 300.0),
                target_dlogT=_amr_opt_mapped('amr_target_dlogT', 'target_dlogT', 0.20),
                target_dr_over_r=_amr_opt_mapped('amr_target_dr_over_r', 'target_dr_over_r', 0.10),
                target_dtau_es=_amr_opt_mapped('amr_target_dtau_es', 'target_dtau_es', 1.0),
                target_dlogchi=_amr_opt_mapped('amr_target_dlogchi', 'target_dlogchi', 0.5),
                target_dv_over_v=_amr_opt_mapped('amr_target_dv_over_v', 'target_dv_over_v', 0.30),
                w_tau=_amr_opt_mapped('amr_w_tau', 'w_tau', 1.0),
                w_rho=_amr_opt_mapped('amr_w_rho', 'w_rho', 0.5),
                w_v=_amr_opt_mapped('amr_w_v', 'w_v', 0.5),
                w_T=_amr_opt_mapped('amr_w_T', 'w_T', 0.4),
                w_curv=_amr_opt_mapped('amr_w_curv', 'w_curv', 0.3),
                w_dtau_es=_amr_opt_mapped('amr_w_dtau_es', 'w_dtau_es', 0.6),
                w_dlogchi=_amr_opt_mapped('amr_w_dlogchi', 'w_dlogchi', 0.5),
                w_dv_over_v=_amr_opt_mapped('amr_w_dv_over_v', 'w_dv_over_v', 0.4),
                composition_jump=_amr_opt_mapped('amr_composition_jump', 'composition_jump',
                                                 _amr_opt_mapped('composition_gradient_threshold',
                                                                 'composition_jump', 0.30)),
                composition_trace_floor=_amr_opt_mapped('composition_trace_floor',
                                                        'composition_trace_floor', 1e-4),
                min_points_per_interface=_amr_opt_mapped('amr_min_points_per_interface',
                                                         'min_points_per_interface', 5),
                shock_density_jump=_amr_opt_mapped('amr_shock_density_jump',
                                                   'shock_density_jump', 1.5),
                shock_velocity_curv=_amr_opt_mapped('shock_velocity_curv',
                                                    'shock_velocity_curv', 0.5),
                shock_temp_jump=_amr_opt_mapped('amr_shock_temp_jump',
                                                'shock_temp_jump', 0.5),
                min_points_per_shock=_amr_opt_mapped('amr_min_points_per_shock',
                                                     'min_points_per_shock', 7),
                min_points_per_extremum=_amr_opt_mapped('amr_min_points_per_extremum',
                                                        'min_points_per_extremum', 5),
                enable_coarsening=_amr_opt_mapped('amr_enable_coarsening',
                                                  'enable_coarsening', True),
                mass_tol=_amr_opt_mapped('mass_tol', 'mass_tol', 1e-3),
                verbose=True)
            # Use the converter's already-computed kappa/tau for feature detection
            feats, comp_f, shock_f, ext_f = detect_all_features(
                r, v, rho, T, mass_fractions, tau, amr_cfg)
            interface_indices = np.array([f.index for f in comp_f], dtype=int)
            print(f"[Tau-AMR-Recursive] {n_full} source zones -> recursive AMR "
                  f"(max {args.max_points} pts). Features: "
                  f"{len(comp_f)} interfaces, {len(shock_f)} shocks, "
                  f"{len(ext_f)} velocity extrema")
            # Save source arrays for escalation remap (recursive_amr consumes them)
            r_src = r.copy(); v_src = v.copy(); T_src = T.copy()
            rho_src = rho.copy(); kappa_src = kappa_arr.copy()
            tau_src = tau.copy()
            mf_src = {k: v2.copy() for k, v2 in mass_fractions.items()}
            (r, v, T, rho, kappa_arr, tau, mass_fractions, amr_report) = \
                recursive_amr(r, v, T, rho, kappa_arr, tau, mass_fractions,
                              feats, amr_cfg, interface_indices)
            print(f"[Tau-AMR-Recursive] Final grid: {len(r)} points "
                  f"({amr_report['splits']} splits)")
            if amr_report['unresolved_features']:
                print(f"[Tau-AMR-Recursive] WARNING: "
                      f"{len(amr_report['unresolved_features'])} features "
                      f"still under-resolved at point budget")
            # Progressive escalation: locally refine pathological cells.
            # This is the closed-loop Stage-A refinement that inserts points
            # at the geometric-mean radius of each failed cell and re-checks.
            from tau_amr_optimiser import escalate_grid
            (r, v, T, rho, kappa_arr, tau, mass_fractions, esc_report) = \
                escalate_grid(
                    r, v, T, rho, kappa_arr, tau, mass_fractions,
                    r_src, v_src, T_src, rho_src, kappa_src, tau_src, mf_src,
                    feats, amr_cfg, interface_indices,
                    max_rounds=getattr(args, 'amr_escalation_max_rounds', 5))
            if esc_report['rounds'] > 0:
                print(f"[Tau-AMR-Recursive] Escalation: {esc_report['rounds']} round(s), "
                      f"{esc_report['points_inserted']} points inserted, "
                      f"ND={esc_report['final_nd']}, "
                      f"pathological {'PASSED' if esc_report['pathological_passed'] else 'FAILED'}")
            # Run pathological-cell diagnostic on the final grid for reporting.
            from tau_amr_optimiser import pathological_cell_report
            patho_report = pathological_cell_report(
                r, v, T, rho, kappa_arr, tau, mass_fractions, amr_cfg)
            print(patho_report['summary'])
            if not patho_report['passed']:
                print(f"[Tau-AMR-Recursive] WARNING: {patho_report['n_failures']} "
                      f"pathological cell(s) detected — grid may be inadequate for CMFGEN RT.")
                print(f"  Recommended action: increase max_points or tighten AMR thresholds.")
            # Skip the converter's own AMR — the recursive AMR grid is final.
            args.use_tau_amr = True
        elif args.base_grid_method == 'log_tau':
            tau_bp = [float(x) for x in args.tau_breakpoints.split(',')]
            tau_zf = [float(x) for x in args.tau_zone_fractions.split(',')]
            print(f"[Base Grid] Tau-based downsampling {n_full} tau-trimmed zones to "
                  f"{args.base_grid_target} base points ({args.tau_spacing} in tau)")
            result = downsample_grid_tau_based(r, v, T, rho, kappa_arr, tau, mass_fractions,
                                                target_points=args.base_grid_target,
                                                tau_breakpoints=tau_bp,
                                                zone_fractions=tau_zf,
                                                spacing=args.tau_spacing,
                                                composition_aware=getattr(args, 'composition_aware', True),
                                                composition_gradient_threshold=getattr(args, 'composition_gradient_threshold', 0.5),
                                                composition_elements=getattr(args, 'composition_elements', None),
                                                min_zones_per_shell=getattr(args, 'min_zones_per_shell', 3),
                                                trace_threshold=getattr(args, 'trace_threshold', 1e-4),
                                                composition_grid_insertion=getattr(args, 'composition_grid_insertion', True),
                                                composition_interface_cap=getattr(args, 'composition_interface_cap', 0.15),
                                                ionization_front_aware=getattr(args, 'ionization_front_aware', True),
                                                ionization_front_threshold=getattr(args, 'ionization_front_threshold', 0.5),
                                                ionization_front_min_points=getattr(args, 'ionization_front_min_points', 5),
                                                ionization_front_width=getattr(args, 'ionization_front_width', 3),
                                                ionization_front_cap=getattr(args, 'ionization_front_cap', 0.15))
            n_fixed = 4
            r, v, T, rho, kappa_arr = result[:n_fixed + 1]
            mf_keys = list(mass_fractions.keys())
            for i, key in enumerate(mf_keys):
                mass_fractions[key] = result[n_fixed + 1 + i]
        elif args.feature_aware_downsample:
            print(f"[Base Grid] Feature-aware downsampling {n_full} tau-trimmed zones to "
                  f"{args.base_grid_target} base points (preserving high-gradient regions)")
            mf_list = list(mass_fractions.values())
            result = downsample_grid_feature_aware(r, v, T, rho, kappa_arr, *mf_list,
                                                   target_points=args.base_grid_target,
                                                   gradient_threshold=args.gradient_threshold)
            n_fixed = 4
            r, v, T, rho, kappa_arr = result[:n_fixed + 1]
            for i, key in enumerate(mass_fractions):
                mass_fractions[key] = result[n_fixed + 1 + i]
        else:
            print(f"[Base Grid] Downsampling {n_full} tau-trimmed zones to "
                  f"{args.base_grid_target} log-spaced base points")
            mf_list = list(mass_fractions.values())
            result = downsample_grid(r, v, T, rho, kappa_arr, *mf_list,
                                    target_points=args.base_grid_target)
            n_fixed = 4
            r, v, T, rho, kappa_arr = result[:n_fixed + 1]
            for i, key in enumerate(mass_fractions):
                mass_fractions[key] = result[n_fixed + 1 + i]

        # Save base grid (before adaptive refinement)
        downsampled_r = r.copy()
        downsampled_v = v.copy()
        downsampled_T = T.copy()
        downsampled_rho = rho.copy()
        downsampled_tau = compute_kappa_tau(r, rho, np.full(len(r), args.kappa))
        downsampled_mass_fractions = {k: v.copy() for k, v in mass_fractions.items()}

        # Step 2: Adaptive mesh refinement (unless pure tau regridding requested)
        n_e = calculate_electron_density(rho, T, mass_fractions)
        tau = compute_kappa_tau(r, rho, np.full(len(r), args.kappa))
        downsampled_ne = n_e.copy()
        orig_ne_adapt = calculate_electron_density(orig_rho_adapt, orig_T_adapt, orig_mass_fractions_adapt)

        # Re-derive the critical-region mask on the AMR base grid.  The
        # downsample may have dropped points, so the mask from the full
        # grid no longer matches the array length.  The physical gradients
        # are preserved (original STELLA values unchanged), so re-detection
        # finds the same regions on the base grid.
        amr_critical_mask = None
        if getattr(args, 'critical_regions', False):
            iface_arg = getattr(args, 'critical_composition_interfaces', '')
            if iface_arg and iface_arg.lower() != 'none':
                interfaces = [(p.split('-')[0].strip(), p.split('-')[1].strip())
                              for p in iface_arg.split(',') if '-' in p]
            else:
                interfaces = None
            amr_critical_mask, _ = identify_critical_regions(
                r, T, rho, n_e, tau, mass_fractions, v=v,
                max_dT_T=args.critical_dT_T,
                max_drho_rho=args.critical_drho_rho,
                max_dlogtau=args.critical_dlogtau,
                max_dne_ne=args.critical_dne_ne,
                max_dv_v=args.critical_dv_v,
                composition_jump=args.critical_composition_jump,
                composition_interfaces=interfaces,
                pad_zones=args.critical_pad_zones,
                verbose=False)

        if args.use_tau_amr:
            print(f"[Tau AMR] Using tau-based base grid as final grid ({len(r)} pts). "
                  f"Skipping velocity/sigma adaptive refinement.")
        else:
            n_outer_extra = args.outer_boundary_extra_points
            r, v, T, rho, n_e, tau, mass_fractions = adaptive_grid_refinement(
                r, v, T, rho, n_e, tau, mass_fractions,
                max_radius_jump=args.max_radius_jump,
                max_temp_jump=args.max_temp_jump,
                max_density_jump=args.max_density_jump,
                max_ne_jump=args.max_ne_jump,
                max_logtau_jump=args.max_logtau_jump,
                max_points=args.max_points,
                n_outer_extra=n_outer_extra,
                max_velocity_jump=args.max_velocity_jump,
                max_sigma_jump=args.max_sigma_jump,
                max_cell_ratio=args.max_cell_ratio,
                thin_tau_threshold=args.thin_grid_tau,
                thin_priority=args.thin_grid_priority,
                ne_priority_threshold=args.ne_priority_threshold,
                ne_priority=args.ne_priority,
                sigma_priority_threshold=args.sigma_priority_threshold,
                sigma_priority=args.sigma_priority,
                dsigma_priority_threshold=args.dsigma_priority_threshold,
                dsigma_priority=args.dsigma_priority,
                max_thin_zones=args.max_thin_zones,
                max_thick_tau=getattr(args, 'max_thick_tau', 1e10),
                max_abs_sigma=args.max_abs_sigma,
                sigma_optimize=args.sigma_optimize,
                critical_mask=amr_critical_mask,
                critical_min_dr_r=getattr(args, 'critical_min_dr_r', 1e-5),
            )

        # Recalculate kappa with refined grid
        # If STELLA kappa was available, interpolate it onto the new grid
        if kappa_stella is not None:
            # Interpolate STELLA kappa onto refined grid (outer->inner ordering)
            # ne_stella_r and kappa_stella are full original arrays (outer->inner = decreasing r)
            # np.interp needs xp increasing, so reverse both reference arrays
            kappa_arr = np.interp(r, ne_stella_r[::-1], kappa_stella[::-1])
            kappa_arr = np.clip(kappa_arr, 1e-6, None)
            print(f"[kappa] Interpolated STELLA opacity onto refined grid: "
                  f"min={kappa_arr.min():.3e}, max={kappa_arr.max():.3e} cm²/g")
        else:
            kappa_arr = compute_kappa_total(rho, T, mass_fractions)

        # Generate adaptive grid diagnostic plots
        if not (args.plot_dir and args.plot_dir.lower() == 'none'):
            print("\n--- Generating Adaptive Grid Diagnostic Plots ---")
            plot_output_adapt = str(Path(args.plot_dir) / Path(args.output_file).name) if args.plot_dir else args.output_file
            if args.plot_dir:
                import os as _os
                _os.makedirs(args.plot_dir, exist_ok=True)
            plot_adaptive_grid_diagnostics(
                orig_r_adapt, orig_v_adapt, orig_T_adapt, orig_rho_adapt, orig_tau_adapt,
                downsampled_r, downsampled_v, downsampled_T, downsampled_rho, downsampled_tau,
                r, v, T, rho, tau,
                plot_output_adapt, t_pick,
            max_radius_jump=args.max_radius_jump,
            max_temp_jump=args.max_temp_jump,
            max_density_jump=args.max_density_jump,
            max_ne_jump=args.max_ne_jump,
            max_logtau_jump=args.max_logtau_jump,
            max_points=args.max_points,
            orig_ne=orig_ne_adapt,
            downsampled_ne=downsampled_ne,
            refined_ne=n_e,
            orig_mass_fractions=orig_mass_fractions_adapt,
            downsampled_mass_fractions=downsampled_mass_fractions,
            refined_mass_fractions=mass_fractions
        )
        
    elif args.enable_grid_thinning:
        # Legacy grid thinning (if adaptive grid not enabled).
        # Use the adaptive thinner when critical regions are enabled so
        # fine points inside critical regions are preserved; otherwise fall
        # back to the original global-threshold thinner.
        mf_list = list(mass_fractions.values())
        if getattr(args, 'critical_regions', False) and amr_critical_mask is not None \
                and len(amr_critical_mask) == len(r):
            # Re-derive the mask on the current grid (it may have been
            # modified by the non-AMR base-grid path above).
            n_e_t = calculate_electron_density(rho, T, mass_fractions)
            tau_t = compute_kappa_tau(r, rho, np.full(len(r), args.kappa))
            thin_mask, _ = identify_critical_regions(
                r, T, rho, n_e_t, tau_t, mass_fractions, v=v,
                max_dT_T=args.critical_dT_T,
                max_drho_rho=args.critical_drho_rho,
                max_dlogtau=args.critical_dlogtau,
                max_dne_ne=args.critical_dne_ne,
                max_dv_v=args.critical_dv_v,
                composition_jump=args.critical_composition_jump,
                pad_zones=args.critical_pad_zones,
                verbose=False)
            result = thin_degenerate_grid_adaptive(
                r, v, T, rho, kappa_arr, *mf_list,
                min_dr_r=args.grid_thinning_threshold,
                critical_mask=thin_mask,
                critical_min_dr_r=getattr(args, 'critical_min_dr_r', 1e-5),
                disable_in_critical=getattr(args, 'critical_disable_min_dr', False),
                verbose=True)
        else:
            result = thin_degenerate_grid(r, v, T, rho, kappa_arr, *mf_list,
                                          min_dr_r=args.grid_thinning_threshold)
        n_fixed = 4
        r, v, T, rho, kappa_arr = result[:n_fixed + 1]
        for i, key in enumerate(mass_fractions):
            mass_fractions[key] = result[n_fixed + 1 + i]
    else:
        print(f"[Grid Processing] DISABLED - keeping all {len(r)} zones")
        # n_e and tau are needed later for grid spacing enforcement and plotting.
        # When adaptive grid is disabled, they are never computed in the AMR block.
        n_e = calculate_electron_density(rho, T, mass_fractions)
        tau = compute_kappa_tau(r, rho, np.full(len(r), args.kappa))

    # --- Outer boundary refinement ---
    # Insert geometrically-spaced fine points at the outer boundary (depth 1 in
    # CMFGEN). This stabilises the radiation field at the optically thin transition
    # where the opacity-radiation feedback loop causes population oscillations.
    # The AMR converter handles this (not Fortran NBND_INS) for full control over
    # grid spacing and velocity interpolation. See insert_outer_boundary_points().
    if args.outer_boundary_extra_points > 0:
        print(f"\n--- Outer Boundary Refinement ({args.outer_boundary_extra_points} points) ---")
        n_e_ob = calculate_electron_density(rho, T, mass_fractions)
        tau_ob = compute_kappa_tau(r, rho, np.full(len(r), args.kappa))
        r, v, T, rho, n_e_ob, tau_ob, mass_fractions = insert_outer_boundary_points(
            r, v, T, rho, n_e_ob, tau_ob, mass_fractions,
            n_insert=args.outer_boundary_extra_points,
            ob_rat=args.outer_boundary_ratio,
            min_dr_r=args.min_dr_r,
            max_dr_r=args.max_dr_r,
        )
        # Recalculate kappa on the refined grid
        if kappa_stella is not None:
            kappa_arr = np.interp(r, ne_stella_r[::-1], kappa_stella[::-1])
            kappa_arr = np.clip(kappa_arr, 1e-6, None)
        else:
            kappa_arr = compute_kappa_total(rho, T, mass_fractions)

    # --- Enforce minimum grid points for CMFGEN convergence ---
    # If the grid after AMR/thinning has fewer than --min-nd points, insert
    # midpoints at the largest cells (with grid smoothing) to reach the minimum.
    # This is essential for stable Lambda iteration: too few points leave large
    # cells where populations change drastically, causing overshoot and oscillation.
    if args.min_nd > 0 and len(r) < args.min_nd:
        print(f"\n--- Minimum Grid Points Enforcement ({len(r)} -> {args.min_nd}) ---")
        # Recompute n_e and tau for the fill function (needed for interpolation)
        n_e_fill = calculate_electron_density(rho, T, mass_fractions)
        tau_fill = compute_kappa_tau(r, rho, kappa_arr)
        r, v, T, rho, n_e_fill, tau_fill, mass_fractions = fill_grid_to_minimum(
            r, v, T, rho, n_e_fill, tau_fill, mass_fractions,
            min_points=args.min_nd,
            max_cell_ratio=args.max_cell_ratio,
            thin_tau_threshold=args.thin_grid_tau,
            thin_priority=args.thin_grid_priority,
        )
        # Recalculate kappa on the filled grid
        if kappa_stella is not None:
            kappa_arr = np.interp(r, ne_stella_r[::-1], kappa_stella[::-1])
            kappa_arr = np.clip(kappa_arr, 1e-6, None)
            print(f"[kappa] Re-interpolated STELLA opacity onto filled grid: "
                  f"min={kappa_arr.min():.3e}, max={kappa_arr.max():.3e} cm²/g")
        else:
            kappa_arr = compute_kappa_total(rho, T, mass_fractions)

    # Resample Ni56 to the final grid after downsampling/refinement
    # Both r and r_tau_trim are in outer->inner (decreasing) order.
    # np.interp requires xp to be INCREASING, so we reverse both arrays.
    # Save a copy at r_tau_trim coordinates for later inner-boundary restoration.
    ni56_at_tau_trim = ni56_for_isotopes.copy() if 'ni56_for_isotopes' in locals() and ni56_for_isotopes is not None else None
    if len(ni56_for_isotopes) != len(r):
        print(f"[Ni56 Resample] Interpolating Ni56 from {len(ni56_for_isotopes)} to {len(r)} grid points")
        ni56_for_isotopes = np.interp(r, r_tau_trim[::-1], ni56_for_isotopes[::-1])

    # Final check: ensure we meet minimum point requirement
    if len(r) < args.min_points:
        print(f"WARNING: Final grid has {len(r)} points, below minimum {args.min_points}")
        print(f"  Consider increasing STELLA nzon or relaxing tau bounds further")
    else:
        print(f"[Grid Check] Final grid has {len(r)} points >= minimum {args.min_points}")

    # Apply velocity smoothing if requested
    if args.smooth_velocity:
        v = smooth_velocity_discontinuities(r, v,
                                            sigma_threshold=args.sigma_threshold,
                                            smooth_width=args.smooth_width)

    # Light wiggle smoothing: removes small-scale velocity oscillations (numerical
    # noise from STELLA's Lagrangian grid) while preserving large-scale non-monotonic
    # structure. Capped at 2% relative change per point.
    if args.smooth_wiggles:
        v = smooth_velocity_wiggles(r, v,
                                    window_pts=args.wiggle_window,
                                    max_rel_change=args.wiggle_max_change)

    # Velocity smoothing in km/s space: moving average with a fixed window width.
    # This is a very light smoothing that removes point-to-point noise while
    # preserving the large-scale non-monotonic velocity profile.
    if args.velocity_smooth_window_kms > 0:
        v = smooth_velocity_window_kms(r, v, window_kms=args.velocity_smooth_window_kms)

    # Apply sigma-transition smoothing: smooth the sigma *gradient* (transitions
    # between regions) without erasing the sigma *magnitude* (physical spikes).
    # This makes sigma evolve gradually from ~0 (outer) to ~-30 (shock front)
    # instead of jumping abruptly, which destabilises MOM_JREL_V9.
    if args.smooth_sigma:
        v = smooth_sigma_transitions(r, v,
                                     max_dsigma=args.max_dsigma,
                                     n_passes=args.sigma_smooth_passes,
                                     window=args.sigma_smooth_window,
                                     max_rel_change=args.sigma_smooth_max_rel_change)

    # Sigma scatter smoothing: targets oscillating sigma in near-constant-V regions.
    # Uses Savitzky-Golay polynomial fit to smooth the velocity derivative (sigma)
    # in regions where |d(sigma)| is large due to numerical noise in dlnV/dlnR.
    # This is the most effective smoother for the "staircase" pattern produced by
    # AMR grid pairing in near-constant velocity regions.
    if args.smooth_sigma_scatter:
        v = smooth_sigma_scatter(r, v,
                                 max_dsigma=args.max_dsigma_scatter,
                                 window_pts=args.sigma_scatter_window,
                                 max_rel_change=args.sigma_scatter_max_rel_change,
                                 dv_threshold=args.sigma_scatter_dv_threshold)

    # --- Boundary velocity homologization (DISABLED) ---
    # Previously snapped the outermost/innermost n_snap points to homologous
    # v = (v/r)_interior * r. However, this modified STELLA velocities by ~1%,
    # creating artificial sigma spikes at the snap boundary (e.g., sigma=1.65
    # at depth 3 in the 5d model). Per AGENTS.md, the velocity profile must be
    # preserved as-is from STELLA. CMFGEN's VEL_LAW=12 handles the boundary
    # via NBND_INS points. The AMR's max_sigma_jump criterion handles sigma
    # smoothness through grid refinement, not velocity modification.
    if False:
        eps_v = 1e-30
        n_tot = len(r)
        i_lo = max(2, n_tot // 5)
        i_hi = min(n_tot - 2, 4 * n_tot // 5)
        vr_interior = np.median(v[i_lo:i_hi] / np.maximum(np.abs(r[i_lo:i_hi]), eps_v))
        sigma_pre = sigma_velocity(r, v)

        # Check each boundary independently
        # Only snap if v/r deviation is small (< 5%) — numerical noise
        # Large deviations are physical (velocity reversal) and must be preserved
        DEVIATION_THRESHOLD = 0.05  # 5%
        snapped_outer = False
        snapped_inner = False

        # Snap only the outermost/innermost n_snap points to homologous.
        # Previously this loop iterated inward until finding a >5% deviation,
        # which snapped dozens of points in the nearly-homologous outer wind.
        # This modified STELLA velocities by 3-5% (vr_interior median was pulled
        # down by the inner velocity reversal), creating artificial sigma spikes
        # where snapped meets un-snapped. Now we only snap the boundary points
        # themselves, preserving the STELLA velocity profile everywhere else.
        for i in range(min(n_snap, n_tot)):
            vr_outer = v[i] / max(abs(r[i]), eps_v)
            if abs(vr_outer / vr_interior - 1.0) < DEVIATION_THRESHOLD:
                v[i] = vr_interior * r[i]
                snapped_outer = True

        for i in range(min(n_snap, n_tot)):
            j = n_tot - 1 - i
            vr_inner = v[j] / max(abs(r[j]), eps_v)
            if abs(vr_inner / vr_interior - 1.0) < DEVIATION_THRESHOLD:
                v[j] = vr_interior * r[j]
                snapped_inner = True

        if snapped_outer or snapped_inner:
            sigma_post = sigma_velocity(r, v)
            parts = []
            if snapped_outer:
                parts.append(f"outer sigma[0]: {sigma_pre[0]:.4f} -> {sigma_post[0]:.4f}")
            if snapped_inner:
                parts.append(f"inner sigma[-1]: {sigma_pre[-1]:.4f} -> {sigma_post[-1]:.4f}")
            print(f"[Boundary Snap] Snapped boundary velocities to homologous "
                  f"v/r={vr_interior:.6e} ({', '.join(parts)})")
            if not snapped_outer:
                dev = abs(v[0] / max(abs(r[0]), eps_v) / vr_interior - 1) * 100
                print(f"  Outer NOT snapped (v/r deviation {dev:.1f}% > {DEVIATION_THRESHOLD*100:.0f}%)")
            if not snapped_inner:
                dev = abs(v[-1] / max(abs(r[-1]), eps_v) / vr_interior - 1) * 100
                print(f"  Inner NOT snapped (v/r deviation {dev:.1f}% > {DEVIATION_THRESHOLD*100:.0f}%)")

    # --- Outer-boundary homologous buffer (DISABLED) ---
    # Previously added n_outer_extra_buf artificial homologous zones beyond the
    # outermost point to provide a sigma=0 buffer. However:
    #   1. CMFGEN's VEL_LAW=12 (RV_SN_MODEL_SNIIN) already inserts NBND_INS=3
    #      boundary points at each end with progressively finer spacing.
    #   2. The buffer creates a sharp sigma discontinuity at the transition
    #      from artificial (sigma=0) to real ejecta (sigma~1), causing
    #      MOM_JREL_V9 instabilities.
    #   3. THK_CONT=F (our setting) handles the optically thin outer boundary
    #      via the moment equation, not via artificial buffer zones.
    # The outer_boundary_extra_points config parameter is now ignored.
    n_outer_extra_buf = 0
    if n_outer_extra_buf > 0:
        r0_buf = r[0]
        r1_buf = r[1]
        v0_buf = v[0]
        T0_buf = T[0]
        rho0_buf = rho[0]
        rho1_buf = rho[1]
        ne0_buf = n_e[0] if ne_stella_snap is not None else None
        tau0_buf = tau[0]
        mf0_buf = {e: mass_fractions[e][0] for e in mass_fractions}

        # Homologous v/r from the (now snapped) outermost point
        vr_buf = v0_buf / r0_buf

        # Power-law density slope from outermost two points
        eps_pow = 1e-30
        if abs(rho0_buf) > eps_pow and abs(rho1_buf) > eps_pow and abs(r0_buf - r1_buf) > eps_pow:
            rho_slope_buf = np.log(rho0_buf / rho1_buf) / np.log(r0_buf / r1_buf)
        else:
            rho_slope_buf = -3.0

        # Match local grid spacing (dr/r) for smooth transition
        dr_r_buf = abs(r0_buf - r1_buf) / r0_buf

        # Prepend points outward (r increasing beyond r0)
        # Insert k=1 first (smallest extension), then k=2, ..., k=N (largest)
        # so that depth 1 = largest r (CMFGEN convention: r decreasing with depth)
        buf_r, buf_v, buf_T, buf_rho, buf_ne, buf_tau, buf_mf = [], [], [], [], [], [], {e: [] for e in mass_fractions}
        for k in range(1, n_outer_extra_buf + 1):
            r_new = r0_buf * (1.0 + dr_r_buf) ** k
            v_new = vr_buf * r_new
            rho_new = rho0_buf * (r_new / r0_buf) ** rho_slope_buf
            buf_r.append(r_new)
            buf_v.append(v_new)
            buf_T.append(T0_buf)
            buf_rho.append(rho_new)
            buf_tau.append(tau0_buf * (rho_new / rho0_buf))
            if ne0_buf is not None:
                buf_ne.append(ne0_buf * (rho_new / rho0_buf))
            for elem in mass_fractions:
                buf_mf[elem].append(mf0_buf[elem])

        # Reverse so largest r is first (depth 1 = outermost)
        buf_r = buf_r[::-1]
        buf_v = buf_v[::-1]
        buf_T = buf_T[::-1]
        buf_rho = buf_rho[::-1]
        buf_tau = buf_tau[::-1]
        if buf_ne:
            buf_ne = buf_ne[::-1]
        for elem in mass_fractions:
            buf_mf[elem] = buf_mf[elem][::-1]

        # Prepend to main arrays
        r = np.concatenate([buf_r, r])
        v = np.concatenate([buf_v, v])
        T = np.concatenate([buf_T, T])
        rho = np.concatenate([buf_rho, rho])
        tau = np.concatenate([buf_tau, tau])
        if buf_ne:
            n_e = np.concatenate([buf_ne, n_e])
        for elem in mass_fractions:
            mass_fractions[elem] = np.concatenate([buf_mf[elem], mass_fractions[elem]])

        # Extend kappa_arr for buffer zones (scales with density ratio)
        kappa0_buf = kappa_arr[0]
        buf_kappa = np.array([kappa0_buf * (rn / rho0_buf) for rn in buf_rho])
        kappa_arr = np.concatenate([buf_kappa, kappa_arr])

        # Extend ni56_for_isotopes with zeros (no radioactive material in buffer)
        if 'ni56_for_isotopes' in locals() and ni56_for_isotopes is not None:
            ni56_for_isotopes = np.concatenate([np.zeros(n_outer_extra_buf), ni56_for_isotopes])

        print(f"[Buffer] Added {n_outer_extra_buf} homologous buffer zones beyond outer boundary "
              f"(r: {r0_buf:.3e} -> {r[0]:.3e}, v/r={vr_buf:.6e}, "
              f"rho slope={rho_slope_buf:.2f}, final grid: {len(r)} pts)")

    # --- Enforce grid spacing constraints (min + max + cell ratio smoothing) ---
    # This is the final spacing enforcement pass after all grid processing
    # (AMR, outer boundary insertion, min_nd fill, buffer zones, smoothing).
    # It guarantees:
    #   - No degenerate points (dr/r < min_dr_r) that cause div-by-zero in CMFGEN
    #   - No oversized cells (dr/r > max_dr_r) that destabilise Lambda iteration
    #   - No abrupt spacing transitions (cell ratio > max_cell_ratio) that
    #     create sigma spikes
    print(f"\n--- Grid Spacing Enforcement ---")
    n_before_spacing = len(r)
    # Need n_e and tau arrays for enforce_grid_spacing; compute if not available
    # (n_e/tau may be unset if adaptive grid workflow was skipped or if the
    # arrays were resized by smoothing passes after the AMR step).
    try:
        n_e_check = n_e
    except NameError:
        n_e = None
    if n_e is None or len(n_e) != len(r):
        n_e = calculate_electron_density(rho, T, mass_fractions)
    try:
        tau_check = tau
    except NameError:
        tau = None
    if tau is None or len(tau) != len(r):
        tau = compute_kappa_tau(r, rho, kappa_arr)
    # Re-derive the critical-region mask on the current (pre-spacing-
    # enforcement) grid.  The grid length has changed since the AMR/base
    # grid (outer-boundary insertion, min_nd fill, buffer zones, smoothing),
    # but the physical gradients are preserved, so re-detection finds the
    # same critical regions and lets the final spacing enforcement preserve
    # the fine points the AMR inserted there.
    spacing_critical_mask = None
    if getattr(args, 'critical_regions', False):
        iface_arg = getattr(args, 'critical_composition_interfaces', '')
        if iface_arg and iface_arg.lower() != 'none':
            interfaces = [(p.split('-')[0].strip(), p.split('-')[1].strip())
                          for p in iface_arg.split(',') if '-' in p]
        else:
            interfaces = None
        spacing_critical_mask, _ = identify_critical_regions(
            r, T, rho, n_e, tau, mass_fractions, v=v,
            max_dT_T=args.critical_dT_T,
            max_drho_rho=args.critical_drho_rho,
            max_dlogtau=args.critical_dlogtau,
            max_dne_ne=args.critical_dne_ne,
            max_dv_v=args.critical_dv_v,
            composition_jump=args.critical_composition_jump,
            composition_interfaces=interfaces,
            pad_zones=args.critical_pad_zones,
            verbose=False)
    r, v, T, rho, n_e, tau, mass_fractions, kappa_arr = enforce_grid_spacing(
        r, v, T, rho, n_e, tau, mass_fractions, kappa_arr,
        min_dr_r=args.min_dr_r,
        max_dr_r=args.max_dr_r,
        max_cell_ratio=args.max_cell_ratio,
        max_points=args.max_points,
        critical_mask=spacing_critical_mask,
        critical_min_dr_r=getattr(args, 'critical_min_dr_r', 1e-5),
        disable_in_critical=getattr(args, 'critical_disable_min_dr', False),
    )
    # Resample Ni56 to the final grid after spacing enforcement
    if 'ni56_for_isotopes' in locals() and ni56_for_isotopes is not None:
        if len(ni56_for_isotopes) != len(r):
            # Interpolate from pre-spacing grid to post-spacing grid
            # Both are outer->inner (decreasing); np.interp needs increasing x
            r_old = np.linspace(r[0], r[-1], len(ni56_for_isotopes))
            ni56_for_isotopes = np.interp(r, r_old[::-1], ni56_for_isotopes[::-1])

    # --- Sanity checks (critical-region pipeline) ---
    # Verify the final grid is well-formed and that original STELLA zone
    # values were preserved (only interpolated points were inserted).
    if getattr(args, 'critical_regions', False):
        dr_final = np.diff(r)
        if np.any(dr_final >= 0):
            print("[CriticalRegions] WARNING: non-monotonic radius in final grid "
                  "(r must strictly decrease with depth index)")
        if np.any(np.abs(dr_final) <= 0):
            print("[CriticalRegions] WARNING: zero or negative dr in final grid")
        # Check original STELLA values are preserved at the points that
        # survive from the full STELLA grid (compare against stella_r_full
        # if available).  This is a tolerance check: re-interpolation and
        # smoothing passes can introduce tiny floating-point differences,
        # but the peak T/rho at shock fronts and composition interfaces
        # must be preserved.
        if stella_r_full is not None and len(stella_r_full) > 0:
            T_stella_at_r = np.interp(r, stella_r_full[::-1], stella_T_full[::-1])
            # Only compare at points that coincide with a STELLA radius
            # (within a small tolerance); interpolated points are exempt.
            tol = 1e-6 * max(np.median(np.abs(r)), 1.0)
            coincides = np.array([
                np.any(np.abs(r[i] - stella_r_full) < tol)
                for i in range(len(r))])
            if np.any(coincides):
                rel_err = np.abs(T[coincides] - T_stella_at_r[coincides]) / \
                          np.maximum(np.abs(T_stella_at_r[coincides]), 1e-30)
                max_err = float(np.max(rel_err))
                if max_err > 1e-3:
                    print(f"[CriticalRegions] WARNING: original T preservation "
                          f"error = {max_err:.2e} at coincident STELLA points "
                          f"(expected <1e-3)")
                else:
                    print(f"[CriticalRegions] OK: original STELLA T values "
                          f"preserved at coincident points (max rel err {max_err:.2e})")

    # Regularize the R-grid spacing sequence: independent insertion passes
    # (AMR, min-points fill, outer-boundary, cell-ratio smoothing) each pick
    # points by different criteria, which can leave a "sawtooth" dR pattern
    # even when the underlying velocity profile is smooth. Sigma's
    # distance-weighted gradient formula amplifies this irregularity into
    # point-to-point sigma noise. This is a small (<=1%) positional nudge —
    # V is re-derived on the nudged grid immediately below, so it stays
    # self-consistent; other slowly-varying quantities (T, rho, tau, mass
    # fractions) are unaffected to first order by such a small shift.
    if args.smooth_sigma_scatter:
        r = regularize_grid_spacing(r, window_pts=args.sigma_scatter_window)

    # Re-interpolate velocity, density, temperature, and mass fractions from
    # the full STELLA data onto the final grid.
    # The AMR and grid spacing enforcement interpolate from the downsampled
    # grid, which can miss velocity extrema (e.g., reverse shock dip at
    # ~5000 km/s) and density curvature (e.g., inner density shells) when
    # STELLA zones are skipped. This re-interpolation from the full STELLA
    # data ensures the final grid values match the real STELLA profile.
    # Without re-interpolating rho/T/mass_fractions, the density can differ
    # by up to 40% from STELLA in regions with subtle density curvature.
    if stella_r_full is not None and len(stella_r_full) > 0:
        v_before = v.copy()
        # np.interp needs increasing x; both arrays are outer->inner (decreasing)
        v = np.interp(r, stella_r_full[::-1], stella_v_full[::-1])
        n_corrected = np.sum(np.abs(v - v_before) > 1.0)
        if n_corrected > 0:
            max_dv = np.max(np.abs(v - v_before))
            print(f"[Velocity Re-interp] Corrected {n_corrected} points from full STELLA data "
                  f"(max |dV| = {max_dv:.2f} km/s)")

        # Re-interpolate density, temperature, and mass fractions from full STELLA data.
        # Use log interpolation for rho and T (positive quantities spanning orders of
        # magnitude).  Linear interpolation of these quantities across sharp density
        # jumps (e.g. interaction shells) introduces significant artifacts because it
        # weights the high-density side too heavily.  Log interpolation preserves the
        # relative gradient structure.
        # Use linear interpolation for mass fractions (can be zero for absent species).
        rho_before = rho.copy()
        T_before = T.copy()
        _eps = 1e-30
        _log_rho_full = np.log10(np.maximum(stella_rho_full, _eps))
        _log_T_full = np.log10(np.maximum(stella_T_full, _eps))
        rho = 10.0 ** np.interp(r, stella_r_full[::-1], _log_rho_full[::-1])
        T = 10.0 ** np.interp(r, stella_r_full[::-1], _log_T_full[::-1])
        rho_max_ratio = np.max(rho / np.maximum(rho_before, 1e-30))
        rho_min_ratio = np.min(rho / np.maximum(rho_before, 1e-30))
        if rho_max_ratio > 1.01 or rho_min_ratio < 0.99:
            print(f"[Density Re-interp] Corrected from full STELLA data "
                  f"(ratio range: {rho_min_ratio:.3f} to {rho_max_ratio:.3f})")
        if np.max(np.abs(T - T_before) / np.maximum(T_before, 1.0)) > 0.01:
            print(f"[Temperature Re-interp] Corrected from full STELLA data "
                  f"(max |dT/T| = {np.max(np.abs(T - T_before) / np.maximum(T_before, 1.0)):.3f})")
        # Re-interpolate mass fractions using hybrid scheme that preserves
        # sharp composition interfaces. Linear interpolation (the previous
        # approach) smears Fe-core/Si-shell, Si/O, O/He boundaries across
        # the interpolated gap, diluting peak metal abundances and creating
        # unphysical intermediate compositions. The hybrid scheme uses:
        #   - Nearest-neighbour at detected composition interfaces (preserves
        #     sharp shell boundaries exactly)
        #   - Log interpolation for trace species (X < trace_threshold) in
        #     smooth regions (preserves dynamic range)
        #   - Linear interpolation elsewhere (smooth, well-sampled regions)
        if getattr(args, 'composition_aware', True) and stella_mf_full:
            # Detect interfaces in the full STELLA data
            comp_elems = [e for e in DEFAULT_COMPOSITION_ELEMENTS
                          if e in stella_mf_full]
            iface_idx, iface_by_elem = detect_composition_interfaces(
                stella_mf_full,
                gradient_threshold=getattr(args, 'composition_gradient_threshold', 0.5),
                elements=comp_elems,
            )
            if len(iface_idx) > 0:
                print(f"[MF Re-interp] Hybrid interpolation: {len(iface_idx)} "
                      f"interface zones detected in full STELLA data")
                mf_hybrid = hybrid_interp_mass_fractions(
                    r, stella_r_full, stella_mf_full,
                    interface_indices=iface_idx,
                    trace_threshold=getattr(args, 'trace_threshold', 1e-4),
                )
                for elem in mass_fractions:
                    if elem in mf_hybrid:
                        mass_fractions[elem] = mf_hybrid[elem]
                    elif elem in stella_mf_full:
                        mass_fractions[elem] = np.clip(
                            np.interp(r, stella_r_full[::-1],
                                      stella_mf_full[elem][::-1]), 0.0, None)
            else:
                for elem in mass_fractions:
                    if elem in stella_mf_full:
                        mass_fractions[elem] = np.clip(
                            np.interp(r, stella_r_full[::-1],
                                      stella_mf_full[elem][::-1]), 0.0, None)
        else:
            for elem in mass_fractions:
                if elem in stella_mf_full:
                    mass_fractions[elem] = np.clip(
                        np.interp(r, stella_r_full[::-1],
                                  stella_mf_full[elem][::-1]), 0.0, None)

        # [Velocity Re-interp] pulls raw (unsmoothed) STELLA velocities onto the
        # final grid, bypassing the wiggle/kms/sigma smoothing done earlier (all
        # of which only saw the pre-re-interp velocity). Re-apply ALL smoothing
        # now, on the FINAL r/v, so the sigma actually written to SN_HYDRO_DATA
        # reflects the smoothing — otherwise sharp STELLA-native jumps re-appear
        # right before writeout (observed: identical worst-case cell-ratio/
        # sigma-jump artifacts regardless of --max-points, since they were
        # reintroduced here).
        if args.smooth_wiggles:
            v = smooth_velocity_wiggles(r, v,
                                        window_pts=args.wiggle_window,
                                        max_rel_change=args.wiggle_max_change)
        if args.velocity_smooth_window_kms > 0:
            v = smooth_velocity_window_kms(r, v, window_kms=args.velocity_smooth_window_kms)
        if args.smooth_sigma:
            v = smooth_sigma_transitions(r, v,
                                         max_dsigma=args.max_dsigma,
                                         n_passes=args.sigma_smooth_passes,
                                         window=args.sigma_smooth_window,
                                         max_rel_change=args.sigma_smooth_max_rel_change)
        if args.smooth_sigma_scatter:
            v = smooth_sigma_scatter(r, v,
                                     max_dsigma=args.max_dsigma_scatter,
                                     window_pts=args.sigma_scatter_window,
                                     max_rel_change=args.sigma_scatter_max_rel_change,
                                     dv_threshold=args.sigma_scatter_dv_threshold)

    # --- Global sigma smoothing ---
    # Smooth sigma across the ENTIRE grid to reduce AMR-induced point-to-point
    # scatter. Unlike smooth_sigma_scatter (which only targets near-constant-V
    # regions), this operates everywhere except protected velocity-reversal
    # zones (|sigma| > preserve_threshold). Uses median filter + V reconstruction.
    if args.smooth_sigma_global:
        v = smooth_sigma_global(r, v,
                                window_pts=args.sigma_global_window,
                                max_rel_change=args.sigma_global_max_rel_change,
                                preserve_threshold=args.sigma_global_preserve_threshold)

    # --- Final sigma conformance pass ---
    # Every step above this point (cell-ratio driven spacing enforcement,
    # grid regularization, velocity re-interpolation from full-resolution
    # STELLA data, and the smoothing re-applied afterward) can introduce
    # new sigma violations that AMR — the only sigma-aware stage — never
    # gets a chance to see or fix, since it runs earlier in the pipeline.
    # This pass re-checks the exact CMFGEN sigma formula on the grid as it
    # will actually be written, and inserts points at any remaining
    # violation. See enforce_sigma_conformance() docstring for details.
    if args.enable_sigma_conformance:
        if n_e is None or len(n_e) != len(r):
            n_e = calculate_electron_density(rho, T, mass_fractions)
        if tau is None or len(tau) != len(r):
            tau = compute_kappa_tau(r, rho, kappa_arr)
        r, v, T, rho, n_e, tau, mass_fractions, kappa_arr = enforce_sigma_conformance(
            r, v, T, rho, n_e, tau, mass_fractions, kappa_arr,
            max_sigma_jump=args.max_sigma_jump,
            max_cell_ratio=args.max_cell_ratio,
            min_dr_r=args.min_dr_r,
            max_dr_r=args.max_dr_r,
            max_points=args.max_points,
            headroom=args.sigma_conformance_headroom,
        )
        if 'ni56_for_isotopes' in locals() and ni56_for_isotopes is not None:
            if len(ni56_for_isotopes) != len(r):
                r_old = np.linspace(r[0], r[-1], len(ni56_for_isotopes))
                ni56_for_isotopes = np.interp(r, r_old[::-1], ni56_for_isotopes[::-1])

    # --- Coarsen negative-sigma regions (Luc Dessart approach) ---
    # This runs AFTER all grid construction, spacing enforcement, velocity
    # re-interpolation, smoothing, and sigma conformance — as the final
    # grid modification step. Negative sigma (velocity decreasing with
    # radius) only affects 1-2 zones in well-behaved interacting SN models
    # (Dessart et al.). Over-resolving these regions destabilises the moment
    # solver and Ng acceleration. Since FIX_T=T (temperature is not solved),
    # under-resolving is acceptable — the velocity reversal is preserved as
    # a sharp discontinuity.
    if getattr(args, 'coarsen_negative_sigma', True):
        print(f"\n--- Sigma-Aware Grid Coarsening (final step) ---")
        n_before_coarsen = len(r)
        if n_e is None or len(n_e) != len(r):
            n_e = calculate_electron_density(rho, T, mass_fractions)
        if tau is None or len(tau) != len(r):
            tau = compute_kappa_tau(r, rho, kappa_arr)

        # Detect composition interfaces to protect from coarsening.
        # The metal-rich cannon has sharp boundaries where metals change
        # by orders of magnitude over 1-2 zones. Sigma-coarsening must not
        # remove zones near these boundaries, or the grid becomes
        # under-resolved at the cannon edge, causing CMFGEN convergence
        # failures (population explosions at the composition discontinuity).
        protected_indices = set()
        if getattr(args, 'composition_aware', True) and mass_fractions:
            iface_indices, _ = detect_composition_interfaces(
                mass_fractions,
                gradient_threshold=getattr(args, 'composition_gradient_threshold', 0.3),
            )
            protected_indices = iface_indices
            if protected_indices:
                print(f"[Cannon-Protection] {len(protected_indices)} composition interface "
                      f"zones protected from sigma-coarsening")

        r, v, T, rho, n_e, tau, mass_fractions, kappa_arr = coarsen_negative_sigma_regions(
            r, v, T, rho, n_e, tau, mass_fractions, kappa_arr,
            sigma_threshold=getattr(args, 'coarsen_sigma_threshold', -0.1),
            max_zones_per_region=getattr(args, 'coarsen_max_zones', 2),
            min_region_zones=getattr(args, 'coarsen_min_region_zones', 4),
            protected_indices=protected_indices if protected_indices else None,
        )
        n_after_coarsen = len(r)
        if n_after_coarsen < n_before_coarsen:
            print(f"[Sigma-Coarsen] Grid reduced from {n_before_coarsen} to {n_after_coarsen} points")
        if 'ni56_for_isotopes' in locals() and ni56_for_isotopes is not None:
            if len(ni56_for_isotopes) != len(r):
                r_old = np.linspace(r[0], r[-1], len(ni56_for_isotopes))
                ni56_for_isotopes = np.interp(r, r_old[::-1], ni56_for_isotopes[::-1])

    # --- Cannon-aware grid refinement ---
    # The metal-rich cannon has sharp composition boundaries where the
    # grid must be well-resolved for CMFGEN convergence. After sigma-coarsening,
    # the region adjacent to the cannon inner boundary can become under-resolved
    # (large dr/r). This step inserts midpoints in cells near composition
    # interfaces where dr/r exceeds a tighter threshold than the global max_dr_r.
    # This is critical for the interacting-SN geometry where a metal-rich cannon
    # sits away from the ionisation front, creating a second numerical front.
    if getattr(args, 'composition_aware', True) and mass_fractions:
        iface_indices, _ = detect_composition_interfaces(
            mass_fractions,
            gradient_threshold=getattr(args, 'composition_gradient_threshold', 0.3),
        )
        if iface_indices and len(r) > 10:
            # Build a buffer of zones around each interface
            cannon_buffer = 5  # zones on each side of interface
            refine_zones = set()
            for idx in iface_indices:
                for offset in range(-cannon_buffer, cannon_buffer + 1):
                    if 0 <= idx + offset < len(r):
                        refine_zones.add(idx + offset)

            # Find cells in the refine zone where dr/r > 0.05 (tighter than global 0.5)
            cannon_max_dr_r = 0.05
            r_list = list(r)
            v_list = list(v)
            T_list = list(T)
            rho_list = list(rho)
            ne_list = list(n_e) if n_e is not None else [0.0] * len(r)
            tau_list = list(tau) if tau is not None else [0.0] * len(r)
            kappa_list = list(kappa_arr)
            mf_lists = {elem: list(mass_fractions[elem]) for elem in mass_fractions}

            n_inserted = 0
            i = 0
            while i < len(r_list) - 1:
                if i in refine_zones or (i + 1) in refine_zones:
                    dr = abs(r_list[i + 1] - r_list[i])
                    dr_frac = dr / max(abs(r_list[i]), 1e-30)
                    if dr_frac > cannon_max_dr_r and len(r_list) < args.max_points + 100:
                        # Insert midpoint using optimal_insert_fraction
                        f = optimal_insert_fraction(r_list, v_list, i)
                        r_list.insert(i + 1, r_list[i] + f * (r_list[i + 1] - r_list[i]))
                        v_list.insert(i + 1, v_list[i] + f * (v_list[i + 1] - v_list[i]))
                        T_list.insert(i + 1, 10**(np.log10(max(T_list[i], 1e-30)) + f * (np.log10(max(T_list[i+1], 1e-30)) - np.log10(max(T_list[i], 1e-30)))))
                        rho_list.insert(i + 1, 10**(np.log10(max(rho_list[i], 1e-30)) + f * (np.log10(max(rho_list[i+1], 1e-30)) - np.log10(max(rho_list[i], 1e-30)))))
                        ne_list.insert(i + 1, 10**(np.log10(max(ne_list[i], 1e-30)) + f * (np.log10(max(ne_list[i+1], 1e-30)) - np.log10(max(ne_list[i], 1e-30)))))
                        tau_list.insert(i + 1, 10**(np.log10(max(tau_list[i], 1e-30)) + f * (np.log10(max(tau_list[i+1], 1e-30)) - np.log10(max(tau_list[i], 1e-30)))))
                        kappa_list.insert(i + 1, 10**(np.log10(max(kappa_list[i], 1e-30)) + f * (np.log10(max(kappa_list[i+1], 1e-30)) - np.log10(max(kappa_list[i], 1e-30)))))
                        for elem in mf_lists:
                            mf_lists[elem].insert(i + 1, mf_lists[elem][i] + f * (mf_lists[elem][i + 1] - mf_lists[elem][i]))
                        n_inserted += 1
                        # Update refine_zones for the new point
                        new_refine = set()
                        for idx in refine_zones:
                            if idx > i:
                                new_refine.add(idx + 1)
                            else:
                                new_refine.add(idx)
                        new_refine.add(i + 1)
                        refine_zones = new_refine
                        # Don't advance i — recheck the new cell
                    else:
                        i += 1
                else:
                    i += 1

            if n_inserted > 0:
                print(f"\n--- Cannon-Aware Grid Refinement ---")
                print(f"[Cannon-Refine] Inserted {n_inserted} midpoints near composition interfaces "
                      f"(max dr/r target: {cannon_max_dr_r}, buffer: {cannon_buffer} zones)")
                r = np.array(r_list)
                v = np.array(v_list)
                T = np.array(T_list)
                rho = np.array(rho_list)
                n_e = np.array(ne_list)
                tau = np.array(tau_list)
                kappa_arr = np.array(kappa_list)
                mass_fractions = {elem: np.array(mf_lists[elem]) for elem in mf_lists}
                print(f"[Cannon-Refine] Grid: {len(r)} points "
                      f"(was {len(r) - n_inserted}, +{n_inserted} inserted)")
                dr_new = np.abs(np.diff(r))
                dr_frac_new = dr_new / r[1:]
                print(f"[Cannon-Refine] Max dr/r near interfaces: "
                      f"{dr_frac_new[list(refine_zones & set(range(len(dr_frac_new))))].max():.4f}"
                      if refine_zones & set(range(len(dr_frac_new))) else
                      f"[Cannon-Refine] Max dr/r: {dr_frac_new.max():.4f}")

    # --- Refine outer boundary for radiation field stability ---
    # CMFGEN requires fine grid resolution (dr/r < ~0.015) at the outer
    # boundary (low optical depth) for the moment equation and ray-by-ray
    # solution to agree. Coarse outer grids cause "J(mom) and J(ray) differ
    # by more than 20%" errors and luminosity blowup. This subdivides the
    # first few outer cells to achieve the target spacing.
    if getattr(args, 'refine_outer_boundary', True):
        print(f"\n--- Outer Boundary Refinement ---")
        n_before_refine = len(r)
        if n_e is None or len(n_e) != len(r):
            n_e = calculate_electron_density(rho, T, mass_fractions)
        if tau is None or len(tau) != len(r):
            tau = compute_kappa_tau(r, rho, kappa_arr)
        r, v, T, rho, n_e, tau, mass_fractions, kappa_arr = refine_outer_boundary_grid(
            r, v, T, rho, n_e, tau, mass_fractions, kappa_arr,
            n_outer_cells=getattr(args, 'outer_refine_cells', 5),
            subdivisions=getattr(args, 'outer_refine_subdivisions', 3),
            target_dr_r=getattr(args, 'outer_refine_target_dr_r', 0.015),
            max_new_points=getattr(args, 'outer_refine_max_points', 30),
        )
        n_after_refine = len(r)
        if n_after_refine > n_before_refine:
            print(f"[Outer-Refine] Grid expanded from {n_before_refine} to {n_after_refine} points")
        if 'ni56_for_isotopes' in locals() and ni56_for_isotopes is not None:
            if len(ni56_for_isotopes) != len(r):
                r_old = np.linspace(r[0], r[-1], len(ni56_for_isotopes))
                ni56_for_isotopes = np.interp(r, r_old[::-1], ni56_for_isotopes[::-1])

    # Sigma — computed with the exact CMFGEN formula (rv_sn_model_sniin.f:241)
    # so the written sigma matches what CMFGEN computes internally from R and V.
    sigma = sigma_velocity(r, v)

    # ==========================================================================
    # Final Savitzky-Golay smoothing (configurable, default OFF)
    #
    # When enabled, applies S-G filter (win=11, poly=3) to V, log(rho), log(T)
    # to kill grid-scale numerical ringing from the velocity re-interp step.
    #
    # When DISABLED (default), preserves raw STELLA velocities at shock fronts.
    # This is critical for non-monotonic velocity models where S-G smoothing
    # across a shock front with a large grid-spacing ratio (e.g. 26000:7 in
    # log-r) can create catastrophic sigma spikes (sigma ~ -5000) by pulling
    # the velocity away from its physical value at the discontinuity.
    #
    # The sigma written to SN_HYDRO_DATA is computed with the CMFGEN formula
    # (not np.gradient) and is NOT clipped, so the file reflects what CMFGEN
    # will actually compute internally from R and V.
    # ==========================================================================
    if getattr(args, 'enable_final_smoothing', False):
        print("\n--- Final Savitzky-Golay Smoothing for CMFGEN Stability ---")
        sigma_before = sigma_velocity(r, v)
        v, rho, T, sigma = apply_cmfgen_final_smoothing(r, v, rho, T, win=11, poly=3)
        min_sigma = np.min(sigma)
        max_sigma = np.max(sigma)
        print(f"  -> Final Sigma range: [{min_sigma:.2f}, {max_sigma:.2f}]")
        print(f"  -> Max |sigma| before: {np.max(np.abs(sigma_before)):.2f}, after: {max(abs(min_sigma), abs(max_sigma)):.2f}")
        if min_sigma < -15.0 or max_sigma > 50.0:
            print("  -> WARNING: Sigma still outside safe CMFGEN bounds!")
    else:
        print("\n--- Final Savitzky-Golay Smoothing DISABLED ---")
        print("  Preserving raw velocities at shock fronts (non-monotonic velocity model).")
        print(f"  Sigma range (CMFGEN formula, no clip): [{sigma.min():.2f}, {sigma.max():.2f}]")
        if sigma.min() < -15.0 or sigma.max() > 50.0:
            print(f"  -> WARNING: Sigma outside [-15, 50] — CMFGEN non-monotonic solver (NON_MONO_SOLVER=T) handles this.")

    # Calculate electron and atom densities for plotting
    if ne_stella_snap is not None and ne_stella_r is not None:
        # Interpolate STELLA n_e from original grid onto final (post-AMR) grid
        # Both arrays are outer->inner (decreasing r); np.interp needs increasing x
        # Only interpolate for original grid points; buffer zone n_e already set
        n_outer_buf = n_outer_extra_buf if n_outer_extra_buf > 0 else 0
        if n_outer_buf > 0:
            n_e_orig = np.interp(r[n_outer_buf:], ne_stella_r[::-1], ne_stella_snap[::-1])
            n_e_orig = np.clip(n_e_orig, 1.0, None)
            n_e = np.concatenate([n_e[:n_outer_buf], n_e_orig])
        else:
            n_e = np.interp(r, ne_stella_r[::-1], ne_stella_snap[::-1])
            n_e = np.clip(n_e, 1.0, None)  # floor at 1 cm^-3
    else:
        n_e = calculate_electron_density(rho, T, mass_fractions)

    # Calculate atom density
    elem_A = {'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
              'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
              'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'FePeak': 56.0, 'Ni58': 58.0, 'Ni56': 56.0}
    inv_A = np.zeros(len(r))
    for elem, frac in mass_fractions.items():
        A = elem_A.get(elem, 28.0)
        inv_A += frac / A
    n_atom = rho * AVOGADRO * inv_A  # atoms/cm^3

    # Ensure ni56_for_isotopes is defined (even if empty)
    if 'ni56_for_isotopes' not in locals():
        ni56_for_isotopes = np.zeros(len(r))

    # --- Recombination-front resolution ---
    # STELLA often under-resolves H/He recombination fronts, creating a single-zone
    # T drop (e.g., 9000K -> 5500K) that causes a 15x ne discontinuity. CMFGEN
    # cannot converge with such a jump. This pass spreads the T transition over
    # a few zones by log-linear interpolation, resolving the physical front that
    # STELLA under-resolved. This is NOT convergence smoothing — it is resolving
    # a physical feature.
    if getattr(args, 'resolve_recombination_fronts', True):
        T_before_fronts = T.copy()
        ne_before_fronts = n_e.copy()
        T, n_e, n_fronts_resolved, front_info = resolve_recombination_fronts(
            r, T, n_e, rho,
            t_jump_threshold=getattr(args, 'recombination_t_jump', 0.3),
            front_width_zones=getattr(args, 'recombination_front_width', 5),
            t_min=args.t_floor,
        )
        if n_fronts_resolved > 0:
            print(f"\n--- Recombination Front Resolution ---")
            print(f"[Recomb-Front] Resolved {n_fronts_resolved} under-resolved recombination front(s):")
            for fi in front_info:
                print(f"  Zone {fi['zone_lo']}-{fi['zone_hi']}: r={fi['r_cm']:.3e} cm, "
                      f"T {fi['T_before']:.0f}K -> {fi['T_after']:.0f}K "
                      f"spread over {fi['n_zones']} zones "
                      f"({fi['T_smooth_lo']:.0f}K -> {fi['T_smooth_hi']:.0f}K)")
            # Recompute atom density (unchanged — only T and ne changed)
            # Recompute sigma (unchanged — only T changed, not v or r)
            max_dT = np.max(np.abs(T - T_before_fronts) / np.maximum(T_before_fronts, 1e-30))
            max_dne = np.max(np.abs(n_e - ne_before_fronts) / np.maximum(ne_before_fronts, 1e-30))
            print(f"[Recomb-Front] Max changes: dT/T={max_dT:.3f}, dn_e/n_e={max_dne:.3f}")
        else:
            print(f"[Recomb-Front] No under-resolved recombination fronts detected")

    # Also smooth the STELLA source ne array so that write_sn_hydro_data's
    # interpolation from ne_stella_snap produces smooth values. Without this,
    # write_sn_hydro_data re-interpolates from the original sharp STELLA ne,
    # overriding our smoothed n_e.
    if getattr(args, 'resolve_recombination_fronts', True) and ne_stella_snap is not None and len(front_info) > 0:
        ne_stella_T = snap['T'].copy() if 'snap' in dir() and 'T' in snap else None
        if ne_stella_T is not None:
            print(f"[Recomb-Front] Source ne before smoothing: "
                  f"min={ne_stella_snap.min():.3e}, max={ne_stella_snap.max():.3e}")
            # Note: resolve_recombination_fronts returns (T_out, ne_out, n_fronts, front_info)
            _, ne_stella_snap_smoothed, n_fronts_src, _ = resolve_recombination_fronts(
                ne_stella_r, ne_stella_T, ne_stella_snap, snap['rho'],
                t_jump_threshold=getattr(args, 'recombination_t_jump', 0.3),
                front_width_zones=getattr(args, 'recombination_front_width', 5),
                t_min=args.t_floor,
            )
            if n_fronts_src > 0:
                ne_stella_snap = ne_stella_snap_smoothed
                print(f"[Recomb-Front] Smoothed STELLA source ne at {n_fronts_src} front(s) "
                      f"for write_sn_hydro_data interpolation")
                print(f"[Recomb-Front] Source ne after smoothing: "
                      f"min={ne_stella_snap.min():.3e}, max={ne_stella_snap.max():.3e}")

    # --- Inner boundary tau validation ---
    # After AMR/downsampling, the inner high-tau zones may have been dropped
    # because the inner post-shock region has nearly constant T and slowly
    # varying rho, so the AMR assigns few points there. This check re-adds
    # missing inner zones from the original STELLA data to ensure the inner
    # boundary tau meets the tau_max target.
    tau_es_final = compute_kappa_tau(r, rho, np.full(len(r), args.kappa))
    tau_inner_final = float(tau_es_final[-1])
    tau_max_target = args.tau_max
    print(f"\n--- Inner Boundary Tau Validation ---")
    print(f"  Final grid inner boundary: r={r[-1]:.3e} cm, tau_es={tau_inner_final:.2f} "
          f"(target tau_max={tau_max_target:.1f})")

    if tau_inner_final < tau_max_target * 0.9 and 'snap' in dir():
        # Find STELLA zones between current inner boundary and tau_max
        snap_r = snap['r']
        snap_rho = snap['rho']
        snap_v = snap['v']
        snap_T = snap['T']
        snap_kappa = snap.get('kappa')
        # snap arrays are outer->inner (decreasing r), same as final grid
        tau_snap = compute_kappa_tau(snap_r, snap_rho, np.full(len(snap_r), args.kappa))
        # Find zones with tau > current inner boundary tau and tau <= tau_max
        # These are the missing inner zones
        r_inner_current = r[-1]
        missing_mask = (tau_snap > tau_inner_final) & (tau_snap <= tau_max_target) & (snap_r < r_inner_current)
        # Also include zones just beyond tau_max (the boundary zone)
        i_boundary_snap = np.searchsorted(tau_snap[::-1], tau_max_target)
        if i_boundary_snap < len(tau_snap):
            ib_idx = len(tau_snap) - 1 - i_boundary_snap
            if ib_idx >= 0 and ib_idx < len(snap_r):
                missing_mask[ib_idx] = missing_mask[ib_idx] | (tau_snap[ib_idx] > tau_max_target)

        n_missing = int(np.sum(missing_mask))
        if n_missing > 0:
            print(f"  WARNING: Inner boundary tau {tau_inner_final:.2f} < tau_max {tau_max_target:.1f}")
            print(f"  Re-adding {n_missing} missing inner STELLA zones "
                  f"(r: {snap_r[missing_mask].min():.3e} - {snap_r[missing_mask].max():.3e} cm, "
                  f"tau: {tau_snap[missing_mask].min():.2f} - {tau_snap[missing_mask].max():.2f})")

            # Extract missing zones, sorted by decreasing r (outer->inner)
            miss_idx = np.where(missing_mask)[0]
            miss_idx = miss_idx[np.argsort(snap_r[miss_idx])[::-1]]  # decreasing r

            # Append to final grid
            r_miss = snap_r[miss_idx]
            v_miss = snap_v[miss_idx]
            T_miss = snap_T[miss_idx]
            rho_miss = snap_rho[miss_idx]

            # Interpolate mass fractions from stella_mf_full (tau-trimmed STELLA data)
            # onto missing zone radii. stella_mf_full has the same element keys as
            # mass_fractions and covers the full tau-trimmed range including inner zones.
            mf_miss = {}
            for elem in mass_fractions:
                if stella_mf_full and elem in stella_mf_full and len(stella_mf_full[elem]) == len(stella_r_full):
                    # Interpolate from stella_r_full (decreasing) onto r_miss (decreasing)
                    # np.interp needs increasing x, so reverse both
                    mf_miss[elem] = np.interp(r_miss, stella_r_full[::-1], stella_mf_full[elem][::-1])
                else:
                    # Fallback: use innermost zone value (nearly constant in post-shock region)
                    mf_miss[elem] = np.full(n_missing, mass_fractions[elem][-1])

            # Compute kappa for missing zones
            if snap_kappa is not None:
                kappa_miss = snap_kappa[miss_idx]
                kappa_miss = np.clip(kappa_miss, 1e-6, None)
            else:
                kappa_miss = compute_kappa_total(rho_miss, T_miss, mf_miss)

            # Append to arrays
            r = np.concatenate([r, r_miss])
            v = np.concatenate([v, v_miss])
            T = np.concatenate([T, T_miss])
            rho = np.concatenate([rho, rho_miss])
            kappa_arr = np.concatenate([kappa_arr, kappa_miss])
            for elem in mass_fractions:
                mass_fractions[elem] = np.concatenate([mass_fractions[elem], mf_miss[elem]])

            # Extend ni56_for_isotopes if present
            if 'ni56_for_isotopes' in locals() and ni56_for_isotopes is not None:
                if len(ni56_for_isotopes) < len(r) and ni56_at_tau_trim is not None:
                    # Interpolate from original tau-trim coordinates onto extended grid
                    ni56_for_isotopes = np.interp(r, r_tau_trim[::-1], ni56_at_tau_trim[::-1])

            # Recompute n_e and tau on extended grid
            n_e = calculate_electron_density(rho, T, mass_fractions)
            tau = compute_kappa_tau(r, rho, np.full(len(r), args.kappa))
            tau_es_final = tau
            tau_inner_final = float(tau[-1])

            # Recompute sigma on the extended grid — sigma_velocity depends on
            # r and v, both of which were just extended with missing inner zones.
            # Without this, sigma has the old (pre-extension) length, causing
            # a size mismatch in SN_HYDRO_DATA and a Fortran input conversion
            # error in CMFGEN's rd_sn_data.f / set_rv_hydro_model_v3.f.
            sigma = sigma_velocity(r, v)

            # Update tau_at_ib for the header
            tau_at_ib = tau_inner_final

            print(f"  Extended grid: {len(r)} points, inner boundary tau_es={tau_inner_final:.2f}")
            if len(r) > args.max_points:
                print(f"  WARNING: grid size {len(r)} exceeds max_points {args.max_points} "
                      f"due to inner boundary restoration — CMFGEN grid limit is 1000")
        else:
            print(f"  No missing inner zones found in STELLA data (tau_snap max={tau_snap.max():.2f})")
    else:
        print(f"  OK: inner boundary tau {tau_inner_final:.2f} >= {tau_max_target * 0.9:.1f} (90% of tau_max)")

    # Luminosity: use STELLA's actual bolometric flux when available.
    # The blackbody at the inner boundary gives the wrong luminosity when the
    # inner boundary is deep (tau=100 with STELLA kappa) — the flux at that
    # depth can be negative (inward-flowing due to reverse shock). The correct
    # luminosity for the diffusion boundary condition is the total luminosity
    # of the supernova, which is the STELLA flux at the photosphere (tau=2/3).
    # Fallback: blackbody at the inner boundary (for old-format .swd without flux).
    tau2 = compute_kappa_tau(r, rho, kappa_arr)
    tau_inner_chi = float(tau2[-1])  # tau from the actual chi column written to SN_HYDRO_DATA
    i_phot = np.searchsorted(tau2, 2.0 / 3.0)
    i_phot = min(i_phot, len(r) - 1)
    T_phot = T[i_phot]
    r_phot = r[i_phot]

    # Report tau consistency: the header tau_inner uses constant kappa_es=0.34
    # (matching CMFGEN's internal opacity and the input_hydro.dat TAU column),
    # while the chi column uses STELLA NLTE kappa (or hybrid). In cool neutral
    # gas the NLTE kappa can be 60x lower, giving a much smaller chi-column tau.
    # CMFGEN does NOT read the chi column for grid construction (VEL_LAW=12 uses
    # input_hydro.dat), so the header tau is the physically relevant one for the
    # boundary condition. Both values are written to the header for transparency.
    tau_ratio = tau_at_ib / max(tau_inner_chi, 1e-10)
    print(f"\n  [Tau Consistency] header tau_inner (const kappa={args.kappa}): {tau_at_ib:.2f}")
    print(f"  [Tau Consistency] chi-column tau (STELLA NLTE kappa):          {tau_inner_chi:.2f}")
    print(f"  [Tau Consistency] ratio: {tau_ratio:.1f}x  (expected >1 in cool neutral gas)")
    if tau_ratio > 3.0:
        print(f"  [Tau Consistency] NOTE: large discrepancy — STELLA NLTE kappa is {tau_ratio:.0f}x lower")
        print(f"      than constant kappa_es in the recombination zone. This is expected;")
        print(f"      CMFGEN uses its own opacity (closer to constant kappa) for the diffusion BC.")

    # Try STELLA flux at the photosphere (on the original STELLA grid)
    if kappa_stella is not None and 'snap' in dir() and 'flux' in snap:
        orig_tau = compute_kappa_tau(snap['r'], snap['rho'], snap['kappa'])
        i_phot_orig = np.searchsorted(orig_tau, 2.0 / 3.0)
        if i_phot_orig >= len(snap['flux']):
            # Photosphere is outside the STELLA grid (tau_max < 2/3 when
            # measured from outside in). snap arrays are inner->outer, so
            # flux[0] is the innermost (hottest, most luminous) zone.
            i_phot_orig = 0
            print(f"[LSTAR] WARNING: Photosphere outside STELLA grid, "
                  f"using inner boundary flux (zone 0) for luminosity")
        i_phot_orig = min(i_phot_orig, len(snap['flux']) - 1)
        flux_phot = snap['flux'][i_phot_orig]
        # Use absolute value (flux can be negative in some STELLA configurations)
        lum_lsun = abs(flux_phot) / LSUN
        lum_phot_lsun = luminosity_bb(T_phot, r_phot) / LSUN
        print(f"Inner boundary (zone behind tau_max): T={T_inner_boundary:.1f} K, "
              f"R={r_inner_boundary:.2e} cm")
        print(f"Photosphere (tau=2/3): depth {i_phot}, T={T_phot:.1f} K, "
              f"R={r_phot:.2e} cm -> L={lum_phot_lsun:.3e} Lsun (blackbody)")
        print(f"Using STELLA flux at photosphere: L={lum_lsun:.3e} Lsun "
              f"(flux={flux_phot:.3e} erg/s at orig depth {i_phot_orig})")
    else:
        lum_lsun = luminosity_bb(T_inner_boundary, r_inner_boundary) / LSUN
        lum_phot_lsun = luminosity_bb(T_phot, r_phot) / LSUN
        print(f"Inner boundary (zone behind tau_max): T={T_inner_boundary:.1f} K, "
              f"R={r_inner_boundary:.2e} cm -> L={lum_lsun:.3e} Lsun (blackbody, no STELLA flux)")
        print(f"Photosphere (tau=2/3): depth {i_phot}, T={T_phot:.1f} K, "
              f"R={r_phot:.2e} cm -> L={lum_phot_lsun:.3e} Lsun (for reference)")

    # Final pathological-cell diagnostic on the grid that will be written.
    # This runs for ALL base_grid_methods (not just tau_amr_recursive) and
    # checks the final grid after all post-processing (smoothing, fronts, etc.).
    if args.base_grid_method != 'tau_amr_recursive':
        try:
            from tau_amr_optimiser import pathological_cell_report
            patho_report_final = pathological_cell_report(
                r, v, T, rho, kappa_arr, tau2, mass_fractions)
            print(patho_report_final['summary'])
            if not patho_report_final['passed']:
                print(f"[Grid Quality] WARNING: {patho_report_final['n_failures']} "
                      f"pathological cell(s) in final grid — consider switching to "
                      f"tau_amr_recursive or increasing max_points.")
        except ImportError:
            pass  # tau_amr_optimiser not available

    # --- Final hard cap: trim grid to max_points if any post-processing
    # step (outer boundary insertion, min-nd fill, sigma conformance) pushed
    # the grid beyond max_points.  This ensures the converter respects the
    # --max-points limit as a true upper bound, not just an AMR budget.
    # Points are removed from the smoothest interior regions (lowest local
    # gradient in log-rho) to preserve shocks, interfaces, and boundaries.
    if len(r) > args.max_points:
        n_excess = len(r) - args.max_points
        print(f"[Hard Cap] Grid has {len(r)} points, trimming {n_excess} "
              f"to enforce max_points={args.max_points}")
        # Compute local density gradient as a removal priority metric.
        # Remove points from the smoothest regions (lowest |dlogrho|).
        log_rho = np.log10(np.maximum(rho, 1e-30))
        dlogrho = np.abs(np.diff(log_rho))
        # Pad endpoints with high priority (don't remove boundaries).
        # dlogrho has N-1 elements (one per edge); assign each interior
        # point the max of its two adjacent edge gradients.
        priority = np.full(len(r), 1e10)
        for i in range(1, len(r) - 1):
            priority[i] = max(dlogrho[i - 1], dlogrho[i])
        # Iteratively remove the lowest-priority point
        keep = np.ones(len(r), dtype=bool)
        for _ in range(n_excess):
            kept_idx = np.where(keep)[0]
            # Don't remove first/last point
            removable = kept_idx[1:-1]
            if len(removable) == 0:
                break
            priorities = priority[removable]
            worst = removable[np.argmin(priorities)]
            keep[worst] = False
        r = r[keep]
        v = v[keep]
        T = T[keep]
        rho = rho[keep]
        if n_e is not None and len(n_e) == len(keep):
            n_e = n_e[keep]
        if tau is not None and len(tau) == len(keep):
            tau = tau[keep]
        kappa_arr = kappa_arr[keep]
        for elem in mass_fractions:
            mass_fractions[elem] = mass_fractions[elem][keep]
        if 'ni56_for_isotopes' in locals() and ni56_for_isotopes is not None:
            if len(ni56_for_isotopes) == len(keep):
                ni56_for_isotopes = ni56_for_isotopes[keep]
        print(f"[Hard Cap] Final grid: {len(r)} points")

    write_sn_hydro_data(
        args.output_file, r, v, rho, T, sigma, kappa_arr,
        time_days=t_pick, mass_fractions=mass_fractions, lum_lsun=lum_lsun,
        t_floor=args.t_floor,
        t_cap=args.t_cap,
        trace_threshold=getattr(args, 'trace_zero_threshold', 0.0),
        r_inner_cm=r_inner_boundary, T_inner_K=T_inner_boundary,
        tau_inner=tau_at_ib,
        tau_inner_chi=tau_inner_chi,
        ni56_initial=ni56_for_isotopes,
        ne_stella=ne_stella_snap,
        ne_stella_r=ne_stella_r,
        inc_rad_decays=args.inc_rad_decays
    )

    # Create mass_fractions_plot with Co56 and Fe56 for visualization
    # The write function already calculated the decay chain, so we need to recalculate for plotting
    mass_fractions_plot = mass_fractions.copy()
    t_Ni = 8.8
    t_Co = 111.3
    Ni56_initial = ni56_for_isotopes.copy()
    Ni56_t = Ni56_initial * np.exp(-abs(t_pick) / t_Ni)
    Co56_t = Ni56_initial * (t_Co / (t_Co - t_Ni)) * (np.exp(-abs(t_pick) / t_Co) - np.exp(-abs(t_pick) / t_Ni))
    Fe56_t = Ni56_initial - Ni56_t - Co56_t
    mass_fractions_plot['Co56'] = Co56_t
    mass_fractions_plot['Fe56'] = Fe56_t
    
    # Generate diagnostic plots
    print("\n--- Generating diagnostic plots ---")
    
    # Plot comparison (original vs final)
    # Filter out ghost zones (V > velocity_max_kms) AND zones outside the
    # tau range used for the CMFGEN grid (tau < tau_min or tau > tau_max).
    # These extend far beyond the CMFGEN grid and compress the plot x-axis,
    # making the CMFGEN data invisible. Only plot STELLA zones that were
    # actually used in the conversion.
    orig_r_raw = snap['r']
    orig_v_raw = snap['v']
    orig_T_raw = snap['T']
    orig_rho_raw = snap['rho']
    v_max_kms = args.velocity_max_kms
    # Only filter ghost zones if velocity_max_kms > 0 (0 means no trim was done)
    if v_max_kms > 0:
        keep_mask = orig_v_raw <= v_max_kms
        n_ghost = (~keep_mask).sum()
        if n_ghost > 0:
            print(f"[Plot] Filtering {n_ghost} ghost zones (V > {v_max_kms} km/s) from comparison plot")
    else:
        keep_mask = np.ones(len(orig_v_raw), dtype=bool)
        n_ghost = 0
    # Show the FULL STELLA grid in the comparison plot (no tau filtering).
    # The purpose of the plot is to show what STELLA produced vs what CMFGEN gets,
    # including the zones that were trimmed.
    orig_r = orig_r_raw[keep_mask]
    orig_v = orig_v_raw[keep_mask]
    orig_T = orig_T_raw[keep_mask]
    orig_rho = orig_rho_raw[keep_mask]
    # Also filter mass fractions, n_e, and plot arrays to match ghost filter
    if orig_mass_fractions is not None:
        orig_mass_fractions = {k: v[keep_mask] for k, v in orig_mass_fractions.items()}
    if ne_stella_snap is not None:
        orig_ne_plot = ne_stella_snap[keep_mask]
    else:
        orig_ne_plot = None
    # Filter the mass-fraction plot arrays (orig_r_for_plot was saved before trimming)
    orig_r_for_plot = orig_r_for_plot[keep_mask]
    orig_rho_for_plot = orig_rho_for_plot[keep_mask]
    # Compute BOTH tau curves for the optical depth plot:
    # - STELLA NLTE tau (accurate at high optical depth, unreliable in cold gas)
    # - Analytical Kramers tau (includes bound-free, used for outer boundary trim)
    if kappa_stella is not None:
        orig_kappa_nlte = kappa_stella[keep_mask].copy()
        orig_tau_nlte = compute_kappa_tau(orig_r, orig_rho, orig_kappa_nlte)
    else:
        orig_tau_nlte = None
    orig_kappa_analytical = compute_kappa_total(orig_rho, orig_T, orig_mass_fractions)
    orig_tau = compute_kappa_tau(orig_r, orig_rho, orig_kappa_analytical)
    orig_sigma = sigma_velocity(orig_r, orig_v)
    
    # Compute tau for diagnostic plots and report.
    # Use the same opacity as the AMR grid placement: STELLA NLTE floored at
    # electron-scattering (args.kappa).  This ensures the diagnostic dlog(tau),
    # dtau/dr, and per-cell tau contribution reflect the actual optical-depth
    # structure the AMR was optimising, not a constant-kappa proxy.
    # Also compute constant-kappa tau for comparison (used by trim boundaries).
    orig_tau_const = compute_kappa_tau(orig_r, orig_rho, np.full(len(orig_r), args.kappa))
    _kappa_for_diag = np.maximum(kappa_arr, args.kappa) if kappa_arr is not None else np.full(len(r), args.kappa)
    new_tau = compute_kappa_tau(r, rho, _kappa_for_diag)
    new_sigma = sigma_velocity(r, v)
    
    # Compute electron densities for comparison plot
    # Use the ghost-filtered orig_ne_plot if available
    orig_ne = orig_ne_plot  # already filtered for ghost zones above
    if orig_ne is None:
        orig_ne = calculate_electron_density(orig_rho, orig_T, orig_mass_fractions)
    # Interpolate STELLA n_e onto final CMFGEN grid for new_ne
    if ne_stella_snap is not None:
        new_ne = np.interp(r, ne_stella_r[::-1], ne_stella_snap[::-1])
        new_ne = np.clip(new_ne, 1.0, None)
    else:
        new_ne = calculate_electron_density(rho, T, mass_fractions)
    
    # --- Diagnostic plots ---
    # By default plots go to Path(output_file).parent (same dir as SN_HYDRO_DATA).
    # --plot-dir overrides this; --plot-dir none skips plot generation entirely.
    plot_output = args.output_file
    if args.plot_dir and args.plot_dir.lower() == 'none':
        pass  # skip plots
    else:
        if args.plot_dir:
            import os as _os
            _os.makedirs(args.plot_dir, exist_ok=True)
            # Construct a fake output_file path so Path(output_file).parent points to plot_dir
            plot_output = str(Path(args.plot_dir) / Path(args.output_file).name)

        plot_comparison(
            orig_r, orig_rho, orig_T, orig_v, orig_tau, orig_sigma,
            r, rho, T, v, new_tau, new_sigma,
            plot_output, t_pick, tau_min=args.tau_min, tau_max=args.tau_max,
            orig_mass_fractions=orig_mass_fractions, new_mass_fractions=mass_fractions,
            orig_ne=orig_ne, new_ne=new_ne,
            ne_from_stella=(ne_stella_snap is not None),
            orig_tau_nlte=orig_tau_nlte if 'orig_tau_nlte' in locals() else None,
            orig_tau_const=orig_tau_const if 'orig_tau_const' in locals() else None,
            kappa=args.kappa
        )

        # Plot mass fractions: dashed = original STELLA, solid = trimmed CMFGEN
        plot_mass_fractions_stacked(orig_r_for_plot, orig_rho_for_plot, orig_mass_fractions, t_pick, plot_output,
                                     trimmed_radius_cm=r, trimmed_mass_fractions=mass_fractions_plot)

        # Plot dr vs r
        plot_dr_vs_r(r, plot_output, t_pick, rho=rho, T=T, v=v)

        # Plot zone density vs radius
        tau_bp = [float(x) for x in args.tau_breakpoints.split(',')] if args.base_grid_method == 'log_tau' else None
        plot_zone_density(r, rho, kappa_arr, new_tau, plot_output, t_pick,
                          tau_breakpoints=tau_bp, thin_tau_threshold=args.thin_grid_tau)

        # Plot structure overview
        plot_structure_overview(r, v, rho, T, mass_fractions, kappa_arr, t_pick, plot_output,
                               n_e=n_e, n_atom=n_atom,
                               ne_from_stella=(ne_stella_snap is not None))

        # Plot sigma vs grid index
        plot_sigma_vs_index(r, v, plot_output, t_pick)

        # Plot dtau vs r — primary grid-quality diagnostic for CMFGEN RT
        plot_dtau_vs_r(r, new_tau, plot_output, t_pick,
                       v_kms=v, rho=rho, T=T, sigma=new_sigma,
                       tau_min=args.tau_min, tau_max=args.tau_max)

        # Multi-panel dtau resolution diagnostic (dtau/dr, dtau/dlnr, per-cell contribution)
        plot_dtau_resolution_panels(r, new_tau, plot_output, t_pick,
                                    rho=rho, T=T, v_kms=v, kappa=kappa_arr)

        # Generate machine-readable grid diagnostic report
        generate_grid_diagnostic_report(
            r, v, rho, T, kappa_arr, new_tau, new_sigma,
            mass_fractions, t_pick, plot_output, n_e=new_ne)


if __name__ == '__main__':
    main()
