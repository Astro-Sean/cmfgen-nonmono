#!/usr/bin/env python3
"""
Generate a physics-motivated EDDFACTOR file for CMFGEN FIXED_J initialization.

The J-feedback loop (population changes -> opacity changes -> J changes ->
population changes at other depths) causes oscillatory divergence in
interacting supernova models. The CMFGEN manual recommends USE_FIXED_J=T
to break this loop: populations converge to consistency with a fixed J,
then J is released for final convergence.

USE_FIXED_J=T requires a valid EDDFACTOR file with the correct frequency
grid (generated internally by CMFGEN based on atomic data and MODEL_SPEC).
This script generates one in two phases:

Phase 1 (J_INIT): Run CMFGEN for 1 iteration with USE_FIXED_J=F, FIX_T=T.
  This generates the correct frequency grid and a rough J from the moment
  equations with the initial proxy populations. The job stops after 1
  iteration (NUM_ITS=1 in IN_ITS), leaving FINISH_REC=1.0 in EDDFACTOR.

Phase 2 (J_REPLACE): Read the EDDFACTOR from phase 1, replace the J values
  with a physics-motivated diluted Planck function:

    J(nu, r) = W(r) * B(T(r), nu)

  where W(r) is the geometric dilution factor and B(T, nu) is the Planck
  function. This removes dependence on the initial proxy populations and
  provides a cleaner physics-motivated J for the FIXED_J run.

Usage:
  python scripts/generate_eddfactor.py <cmfgen_dir> [--phase both|init|replace]
  python scripts/generate_eddfactor.py <cmfgen_dir> --phase both
  python scripts/generate_eddfactor.py <cmfgen_dir> --phase replace --source-eddfactor <path>

The generated EDDFACTOR enables USE_FIXED_J=T in the subsequent CMFGEN run.
After populations converge with fixed J, rerun with USE_FIXED_J=F for final
self-consistent J.
"""

import argparse
import os
import struct
import sys
import shutil
import subprocess
from pathlib import Path

import numpy as np

# Physical constants (CGS)
H_PLANCK = 6.62607015e-27    # erg s
K_BOLTZ = 1.380649e-16       # erg/K
C_LIGHT = 2.99792458e10      # cm/s

# EDDFACTOR record structure (from eddfac_rec_defs_mod.f)
DIE_CONT_REC = 1
RV_REC = 2
EDD_CONT_REC = 3
LINE_CONT_REC = 4
FINISH_REC = 5
INITIAL_RV_REC = 10
INITIAL_ACCESS_REC = 15


def read_eddfactor_info(info_path):
    """Read EDDFACTOR_INFO and return ND, RECL, WORD_SIZE, UNIT_SIZE."""
    with open(info_path) as f:
        lines = f.readlines()
    parts = lines[2].split()
    nd = int(parts[0])
    recl = int(parts[1])
    word_size = int(parts[2])
    unit_size = int(parts[3])
    return nd, recl, word_size, unit_size


def read_eddfactor(edd_path, info_path):
    """Read an EDDFACTOR file and return frequencies, J values, R, V, and metadata.

    Reads both continuum J records (from EDD_CONT_REC) and line J records
    (from LINE_CONT_REC). The line J records are preserved as-is for use
    when USE_FIXED_J=T.
    """
    nd, recl, word_size, unit_size = read_eddfactor_info(info_path)
    rec_size = word_size * (nd + 1)  # bytes per record
    filesize = os.path.getsize(edd_path)
    n_records = filesize // rec_size

    with open(edd_path, 'rb') as f:
        # Read FINISH_REC (record 5, 1-based)
        f.seek((FINISH_REC - 1) * rec_size)
        finish_val = struct.unpack('d', f.read(8))[0]

        if finish_val < 0.5:
            raise ValueError(
                f"EDDFACTOR is incomplete (FINISH_REC={finish_val}). "
                f"Run phase 1 (J_INIT) first to generate a complete EDDFACTOR."
            )

        # Read EDD_CONT_REC (record 3): ACCESS_F, NCF, NDEXT as INTEGERS
        f.seek((EDD_CONT_REC - 1) * rec_size)
        access_f, ncf, ndext = struct.unpack('iii', f.read(12))

        # If ACCESS_F is 0 (old format or not written), scan for first
        # data record with plausible frequency
        if access_f == 0:
            access_f = INITIAL_ACCESS_REC
            for rec in range(INITIAL_ACCESS_REC - 1, n_records):
                f.seek(rec * rec_size)
                data = np.frombuffer(f.read(rec_size), dtype=np.float64)
                freq = data[-1]
                if 1e-3 < freq < 1e3:
                    access_f = rec + 1
                    break

        # Read LINE_CONT_REC (record 4): ACCESS_F for line data
        f.seek((LINE_CONT_REC - 1) * rec_size)
        line_access_f = struct.unpack('i', f.read(4))[0]

        # Read R, V, LANG_COORD from records 10-12
        f.seek((INITIAL_RV_REC - 1) * rec_size)
        r = np.frombuffer(f.read(nd * 8), dtype=np.float64).copy()
        f.seek(INITIAL_RV_REC * rec_size)
        v = np.frombuffer(f.read(nd * 8), dtype=np.float64).copy()
        f.seek((INITIAL_RV_REC + 1) * rec_size)
        lang = np.frombuffer(f.read(nd * 8), dtype=np.float64).copy()

        # Read continuum J data records
        n_cont = ncf if ncf > 0 else 0
        cont_frequencies = np.zeros(n_cont)
        cont_j = np.zeros((n_cont, nd))
        for i in range(n_cont):
            rec = access_f - 1 + i
            if rec >= n_records:
                break
            f.seek(rec * rec_size)
            data = np.frombuffer(f.read(rec_size), dtype=np.float64)
            cont_frequencies[i] = data[-1]
            cont_j[i, :] = data[:nd]

        # Read line J data records (if LINE_CONT_REC is valid)
        line_frequencies = np.zeros(0)
        line_j = np.zeros((0, nd))
        if line_access_f > 0 and line_access_f < n_records:
            n_line = n_records - line_access_f + 1
            line_frequencies = np.zeros(n_line)
            line_j = np.zeros((n_line, nd))
            for i in range(n_line):
                rec = line_access_f - 1 + i
                if rec >= n_records:
                    break
                f.seek(rec * rec_size)
                data = np.frombuffer(f.read(rec_size), dtype=np.float64)
                line_frequencies[i] = data[-1]
                line_j[i, :] = data[:nd]

    return {
        'nd': nd,
        'rec_size': rec_size,
        'access_f': access_f,
        'ncf': n_cont,
        'ndext': ndext,
        'line_access_f': line_access_f,
        'ncf_line': len(line_frequencies),
        'cont_frequencies': cont_frequencies,
        'cont_j': cont_j,
        'line_frequencies': line_frequencies,
        'line_j': line_j,
        'frequencies': cont_frequencies,  # backward compat
        'j_values': cont_j,              # backward compat
        'r': r,
        'v': v,
        'lang': lang,
        'finish_val': finish_val,
    }


def read_sn_hydro_data(hydro_path):
    """Read SN_HYDRO_DATA and return R, V, T, density arrays.

    SN_HYDRO_DATA format (outer-to-inner, matching CMFGEN depth convention):
      Header lines with metadata (ND, time, luminosity, etc.)
      Then labeled sections:
        "Radius grid (10^10cm)"   -> R in units of 10^10 cm
        "Velocity (km/s)"         -> V in km/s
        "Temperature (10^4 K)"    -> T in units of 10^4 K
        "Density (gm/cm^3)"       -> rho in g/cm^3
        "Electron density (/cm^3)" -> ne in cm^-3
      Each section has ND values (possibly spanning multiple lines).
    """
    with open(hydro_path) as f:
        content = f.read()

    # Parse ND from header
    nd = None
    for line in content.splitlines():
        if "Number of data points" in line:
            nd = int(line.split(":")[-1].strip())
            break
    if nd is None:
        raise ValueError("Could not parse ND from SN_HYDRO_DATA header")

    # Helper to extract a named section
    def extract_section(label):
        idx = content.find(label)
        if idx < 0:
            return None
        # Find end of label line
        eol = content.find('\n', idx)
        # Read until next non-numeric line or end of file
        vals = []
        for line in content[eol+1:].splitlines():
            parts = line.split()
            if not parts:
                continue
            try:
                for p in parts:
                    vals.append(float(p))
            except ValueError:
                break
            if len(vals) >= nd:
                break
        return np.array(vals[:nd])

    r = extract_section("Radius grid")
    v = extract_section("Velocity")
    t = extract_section("Temperature")
    rho = extract_section("Density")
    # Rosseland mean opacity (10^{-10} cm^{-1}) — used for tau-aware W
    kappa_ross = extract_section("Rosseland mean opacty")

    # Unit conversions
    if r is not None:
        r = r * 1e10       # 10^10 cm -> cm
    if t is not None:
        t = t * 1e4        # 10^4 K -> K
    # v already in km/s, rho in g/cm^3
    # kappa_ross in 10^{-10} cm^{-1} -> cm^{-1}
    if kappa_ross is not None:
        kappa_ross = kappa_ross * 1e-10

    # Compute Rosseland optical depth at each depth (outside-in, cumulative)
    tau = None
    if r is not None and kappa_ross is not None:
        tau = np.zeros(nd)
        for i in range(1, nd):
            dr = abs(r[i] - r[i - 1])  # cm (R decreases with depth)
            tau[i] = tau[i - 1] + kappa_ross[i] * dr

    return {'nd': nd, 'r': r, 'v': v, 't': t, 'rho': rho, 'tau': tau}


def compute_dilution_factor(r, tau=None, tau_trans=1.0):
    """Compute dilution factor W at each depth.

    With tau (Rosseland optical depth, outside-in cumulative):
      W_eff(i) = W_geom(i) + (1 - W_geom(i)) * (1 - exp(-tau(i)/tau_trans))
      -> W_geom in optically thin (tau << 1)
      -> 1.0 in optically thick (tau >> 1)
    This matches grey_w_eff() in proxy_populations_v2.py, ensuring the
    EDDFACTOR J is consistent with the proxy populations' ionization balance.

    Without tau (fallback): pure geometric dilution.
    """
    r_outer = r[0]
    w_geom = np.zeros_like(r)
    for i in range(len(r)):
        ratio = r[i] / r_outer
        if ratio >= 1.0:
            w_geom[i] = 0.5
        else:
            w_geom[i] = 0.5 * (1.0 - np.sqrt(max(0.0, 1.0 - ratio**2)))

    if tau is None:
        return w_geom

    w_eff = w_geom + (1.0 - w_geom) * (1.0 - np.exp(-tau / tau_trans))
    return np.minimum(w_eff, 1.0)


def planck_function(nu_hz, t_k):
    """Compute Planck function B(T, nu) in erg/cm^2/s/Hz/sr.

    B(T, nu) = (2*h*nu^3/c^2) / (exp(h*nu/(k*T)) - 1)

    Uses Wien approximation for h*nu/(k*T) > 500 to avoid overflow.
    Uses Rayleigh-Jeans for h*nu/(k*T) < 1e-10.
    """
    x = H_PLANCK * nu_hz / (K_BOLTZ * t_k)
    if np.isscalar(x):
        if x > 500:
            return (2 * H_PLANCK * nu_hz**3 / C_LIGHT**2) * np.exp(-x)
        elif x < 1e-10:
            return 2 * K_BOLTZ * t_k * nu_hz**2 / C_LIGHT**2
        else:
            return (2 * H_PLANCK * nu_hz**3 / C_LIGHT**2) / (np.exp(x) - 1)
    else:
        b = np.zeros_like(x, dtype=np.float64)
        wien = x > 500
        rj = x < 1e-10
        mid = ~(wien | rj)
        b[wien] = (2 * H_PLANCK * nu_hz[wien]**3 / C_LIGHT**2) * np.exp(-x[wien])
        b[rj] = 2 * K_BOLTZ * t_k * nu_hz[rj]**2 / C_LIGHT**2
        b[mid] = (2 * H_PLANCK * nu_hz[mid]**3 / C_LIGHT**2) / (np.exp(x[mid]) - 1)
        return b


def compute_grey_j(frequencies_1e15, r, t, tau=None):
    """Compute J = W * B(T, nu) at each frequency and depth.

    Args:
        frequencies_1e15: array of frequencies in units of 10^15 Hz (NCF,)
        r: array of radii in cm (ND,)
        t: array of temperatures in K (ND,)
        tau: optional Rosseland optical depth (ND,) for tau-aware W

    Returns:
        j_grey: array of shape (NCF, ND) with J in erg/cm^2/s/Hz/sr
    """
    w = compute_dilution_factor(r, tau=tau)
    ncf = len(frequencies_1e15)
    nd = len(r)
    j_grey = np.zeros((ncf, nd), dtype=np.float64)

    for i in range(ncf):
        nu_hz = frequencies_1e15[i] * 1e15
        for k in range(nd):
            if t[k] > 0 and w[k] > 0:
                j_grey[i, k] = w[k] * planck_function(nu_hz, t[k])

    return j_grey


def write_eddfactor(edd_path, info_path, nd, frequencies_1e15, j_values, r, v, lang,
                    line_frequencies=None, line_j=None):
    """Write a complete EDDFACTOR file with the given J values.

    Writes both continuum J records (starting at record 15) and line J records
    (starting at record 15+NCF) if line data is provided.
    """
    ncf = len(frequencies_1e15)
    ncf_line = len(line_frequencies) if line_frequencies is not None else 0
    line_access = INITIAL_ACCESS_REC + ncf if ncf_line > 0 else 0
    rec_size = 8 * (nd + 1)  # bytes
    recl_fortran = rec_size // 4  # Fortran RECL in 4-byte units

    # Write EDDFACTOR_INFO
    with open(info_path, 'w') as f:
        f.write(" 12-Apr-2017                                       !INFO format date\n")
        f.write(" 23-Jan-2017                                       !File format date\n")
        f.write(f"         {nd:5d}  {recl_fortran:10d}           8           4           4            T\n")
        f.write("          ND        RECL   WORD_SIZE   UNIT_SIZE    INT_SIZE     LIT_END\n")

    # Write EDDFACTOR binary
    # Control records use Fortran INTEGER (4 bytes) for ACCESS_F/NCF/NDEXT,
    # and REAL(KIND=LDP) (8 bytes, double) for FINISH_REC.
    # J data records use REAL(KIND=LDP) (8 bytes) for J values and frequency.
    with open(edd_path, 'wb') as f:
        # Records 1-5: control records
        # Rec 1 (DIE_CONT_REC): ACCESS_F as INTEGER (0 = no dielectronic data)
        f.write(struct.pack('i', 0))
        f.write(b'\x00' * (rec_size - 4))
        # Rec 2 (RV_REC): IREC as INTEGER (INITIAL_RV_REC = 10)
        f.write(struct.pack('i', INITIAL_RV_REC))
        f.write(b'\x00' * (rec_size - 4))
        # Rec 3 (EDD_CONT_REC): ACCESS_F, NCF, NDEXT as INTEGERS
        f.write(struct.pack('iii', INITIAL_ACCESS_REC, ncf, nd))
        f.write(b'\x00' * (rec_size - 12))
        # Rec 4 (LINE_CONT_REC): ACCESS_F for line data
        f.write(struct.pack('i', line_access))
        f.write(b'\x00' * (rec_size - 4))
        # Rec 5 (FINISH_REC): T1 as REAL(KIND=LDP) (1.0 = complete)
        f.write(struct.pack('d', 1.0))
        f.write(b'\x00' * (rec_size - 8))

        # Records 6-9: reserved
        for _ in range(4):
            f.write(b'\x00' * rec_size)

        # Records 10-12: R, V, LANG_COORD
        f.write(r.astype(np.float64).tobytes())
        f.write(b'\x00' * (rec_size - nd * 8))
        f.write(v.astype(np.float64).tobytes())
        f.write(b'\x00' * (rec_size - nd * 8))
        f.write(lang.astype(np.float64).tobytes())
        f.write(b'\x00' * (rec_size - nd * 8))

        # Records 13-14: gap
        for _ in range(2):
            f.write(b'\x00' * rec_size)

        # Records 15+: continuum J data (ND doubles + 1 double for frequency)
        for i in range(ncf):
            record = np.zeros(nd + 1, dtype=np.float64)
            record[:nd] = j_values[i]
            record[nd] = frequencies_1e15[i]
            f.write(record.tobytes())

        # Line J data (if present): preserve from phase 1
        if ncf_line > 0 and line_j is not None:
            for i in range(ncf_line):
                record = np.zeros(nd + 1, dtype=np.float64)
                record[:nd] = line_j[i]
                record[nd] = line_frequencies[i]
                f.write(record.tobytes())

    total_recs = INITIAL_ACCESS_REC + ncf + ncf_line
    print(f"  Wrote EDDFACTOR: {ncf} cont + {ncf_line} line = {ncf+ncf_line} freqs, "
          f"ND={nd}, {total_recs} records, {os.path.getsize(edd_path)} bytes")
    print(f"  Wrote EDDFACTOR_INFO: ND={nd}, RECL={recl_fortran}")


def phase_init(cmfgen_dir, timeout_minutes=30):
    """Phase 1: Run CMFGEN for 1 iteration to generate the frequency grid."""
    cmfgen_dir = Path(cmfgen_dir)
    in_its_path = cmfgen_dir / "IN_ITS"
    vadat_path = cmfgen_dir / "VADAT"
    edd_path = cmfgen_dir / "EDDFACTOR"
    edd_info_path = cmfgen_dir / "EDDFACTOR_INFO"

    print("\n=== Phase 1: J_INIT (1-iteration CMFGEN run) ===")
    print(f"  Directory: {cmfgen_dir}")

    if not (cmfgen_dir / "SN_HYDRO_DATA").exists():
        print("  ERROR: SN_HYDRO_DATA not found", file=sys.stderr)
        return False
    if not in_its_path.exists():
        print("  ERROR: IN_ITS not found", file=sys.stderr)
        return False
    if not vadat_path.exists():
        print("  ERROR: VADAT not found", file=sys.stderr)
        return False

    # Back up IN_ITS and VADAT
    in_its_backup = in_its_path.read_text()
    vadat_backup = vadat_path.read_text()

    # Modify IN_ITS: NUM_ITS=1
    lines = in_its_backup.splitlines()
    modified = []
    for line in lines:
        if "[NUM_ITS]" in line:
            modified.append(f"1 [NUM_ITS] ")
        else:
            modified.append(line)
    in_its_path.write_text("\n".join(modified) + "\n")
    print("  Set NUM_ITS=1 in IN_ITS")

    # Modify VADAT: ensure USE_FIXED_J=F and COMP_F=T
    # (COMP_F=T forces CMFGEN to compute a new EDDFACTOR with the correct
    # frequency grid; USE_FIXED_J=F allows J to be computed self-consistently)
    vadat_lines = vadat_backup.splitlines()
    modified_v = []
    for line in vadat_lines:
        if "[USE_FIXED_J]" in line:
            modified_v.append("F [USE_FIXED_J]")
        elif "[COMP_F]" in line:
            modified_v.append("T [COMP_F]")
        else:
            modified_v.append(line)
    vadat_path.write_text("\n".join(modified_v) + "\n")
    print("  Set USE_FIXED_J=F, COMP_F=T in VADAT")

    # Remove existing EDDFACTOR (fresh start for COMPUTE_EDDFAC=T)
    if edd_path.exists():
        edd_path.unlink()
    if edd_info_path.exists():
        edd_info_path.unlink()
    print("  Removed existing EDDFACTOR (fresh start)")

    # Run CMFGEN
    batch_sh = cmfgen_dir.resolve() / "batch.sh"
    if not batch_sh.exists():
        print("  ERROR: batch.sh not found", file=sys.stderr)
        in_its_path.write_text(in_its_backup)
        vadat_path.write_text(vadat_backup)
        return False

    print(f"  Running CMFGEN (timeout: {timeout_minutes} min)...")
    try:
        result = subprocess.run(
            ["tcsh", "-f", str(batch_sh)],
            cwd=str(cmfgen_dir.resolve()),
            timeout=timeout_minutes * 60,
            capture_output=True,
            text=True
        )
        print(f"  CMFGEN exit code: {result.returncode}")
        if result.returncode != 0 and result.stderr:
            print(f"  STDERR (last 500 chars): {result.stderr[-500:]}", file=sys.stderr)
    except subprocess.TimeoutExpired:
        print("  ERROR: CMFGEN timed out", file=sys.stderr)
        in_its_path.write_text(in_its_backup)
        vadat_path.write_text(vadat_backup)
        return False
    finally:
        # Restore original IN_ITS and VADAT
        in_its_path.write_text(in_its_backup)
        vadat_path.write_text(vadat_backup)
        print("  Restored original IN_ITS and VADAT")

    # Check FINISH_REC
    if not edd_path.exists() or not edd_info_path.exists():
        print("  ERROR: EDDFACTOR not generated", file=sys.stderr)
        return False

    try:
        nd, recl, word_size, unit_size = read_eddfactor_info(edd_info_path)
        rec_size = word_size * (nd + 1)
        with open(edd_path, 'rb') as f:
            f.seek((FINISH_REC - 1) * rec_size)
            finish_val = struct.unpack('d', f.read(8))[0]
        if finish_val < 0.5:
            print(f"  ERROR: EDDFACTOR incomplete (FINISH_REC={finish_val})", file=sys.stderr)
            return False
        print(f"  EDDFACTOR valid: FINISH_REC={finish_val}, ND={nd}")
    except Exception as e:
        print(f"  ERROR reading EDDFACTOR: {e}", file=sys.stderr)
        return False

    return True


def phase_replace(cmfgen_dir, source_eddfactor=None):
    """Phase 2: Replace J values with W * B(T, nu) from SN_HYDRO_DATA."""
    cmfgen_dir = Path(cmfgen_dir)
    edd_path = cmfgen_dir / "EDDFACTOR"
    edd_info_path = cmfgen_dir / "EDDFACTOR_INFO"
    hydro_path = cmfgen_dir / "SN_HYDRO_DATA"

    src_edd = Path(source_eddfactor) if source_eddfactor else edd_path
    src_info = src_edd.parent / (src_edd.name + "_INFO")

    print("\n=== Phase 2: J_REPLACE (physics-motivated J) ===")
    print(f"  Source EDDFACTOR: {src_edd}")
    print(f"  Target directory: {cmfgen_dir}")

    # Read source EDDFACTOR
    print("  Reading source EDDFACTOR...")
    edd = read_eddfactor(src_edd, src_info)
    nd = edd['nd']
    frequencies = edd['cont_frequencies']
    ncf = edd['ncf']
    ncf_line = edd['ncf_line']
    line_freqs = edd['line_frequencies']
    line_j = edd['line_j']
    print(f"  ND={nd}, NCF={ncf} (continuum), NCF_LINE={ncf_line} (line)")
    print(f"  Cont frequency range: {frequencies[0]:.6e} to {frequencies[-1]:.6e} (10^15 Hz)")
    if ncf_line > 0:
        print(f"  Line frequency range: {line_freqs[0]:.6e} to {line_freqs[-1]:.6e} (10^15 Hz)")

    # Read SN_HYDRO_DATA for T (and V if needed)
    print("  Reading SN_HYDRO_DATA...")
    hydro = read_sn_hydro_data(hydro_path)
    if hydro['nd'] != nd:
        print(f"  WARNING: ND mismatch (EDDFACTOR={nd}, SN_HYDRO_DATA={hydro['nd']})")

    # The EDDFACTOR R/V arrays are on CMFGEN's regridded grid (RV_SN_MODEL_SNIIN
    # redistributes to a log-uniform tau grid).  We must compute J on THIS grid,
    # not the original SN_HYDRO_DATA grid.  EDDFACTOR R is in 1e10 cm units.
    r_regrid = edd['r'].copy()        # 1e10 cm, outside-in
    v_regrid = edd['v'].copy()        # km/s, outside-in
    r_regrid_cm = r_regrid * 1e10     # cm, for W computation
    r_hydro = hydro['r']              # cm, outside-in
    t_hydro = hydro['t']              # K, outside-in
    lang = np.zeros(nd)

    # Read tau from input_hydro.dat (total tau including scattering).
    # The Rosseland opacity in SN_HYDRO_DATA gives absorption-only tau
    # (~3.6 vs ~122 total), so we use input_hydro.dat which has the
    # correct total tau that CMFGEN sees.
    tau_hydro = None
    input_hydro_path = cmfgen_dir / "input_hydro.dat"
    if input_hydro_path.exists():
        try:
            tau_vals = []
            with open(input_hydro_path) as fih:
                for line in fih:
                    parts = line.split()
                    if len(parts) >= 7:
                        try:
                            float(parts[0])
                            tau_vals.append(float(parts[6]))
                        except ValueError:
                            continue
                    if len(tau_vals) >= hydro['nd']:
                        break
            if len(tau_vals) >= hydro['nd']:
                tau_hydro = np.array(tau_vals[:hydro['nd']])
                print(f"  Read tau from input_hydro.dat ({len(tau_hydro)} values)")
        except Exception as e:
            print(f"  WARNING: could not read tau from input_hydro.dat: {e}")

    if tau_hydro is None:
        tau_hydro = hydro.get('tau')
        if tau_hydro is not None:
            print(f"  WARNING: using Rosseland absorption tau (may underestimate total tau)")

    # Interpolate T and tau from SN_HYDRO_DATA grid onto the regridded grid.
    # Both grids are outside-in (R decreasing with depth).  np.interp requires
    # increasing x, so we reverse, interpolate, then reverse back.
    r_hydro_inc = r_hydro[::-1]
    t_hydro_inc = t_hydro[::-1]
    r_regrid_inc = r_regrid_cm[::-1]
    t_regrid = np.interp(r_regrid_inc, r_hydro_inc, t_hydro_inc)[::-1]

    tau_regrid = None
    if tau_hydro is not None:
        tau_hydro_inc = tau_hydro[::-1]
        tau_regrid = np.interp(r_regrid_inc, r_hydro_inc, tau_hydro_inc)[::-1]

    r = r_regrid_cm   # cm, on regridded grid
    t = t_regrid      # K, on regridded grid
    v = v_regrid      # km/s, on regridded grid
    tau = tau_regrid

    print(f"  Regridded R range: {r[0]:.4e} to {r[-1]:.4e} cm")
    print(f"  Regridded T range: {t[0]:.1f} to {t[-1]:.1f} K")
    print(f"  Regridded V range: {v[0]:.1f} to {v[-1]:.1f} km/s")
    if tau is not None:
        print(f"  Regridded Tau range: {tau[0]:.4e} to {tau[-1]:.4e}")

    # Compute dilution factor (tau-aware if tau available)
    w = compute_dilution_factor(r, tau=tau)
    print(f"  W range: {w[0]:.6f} to {w[-1]:.6f}")

    # Compute J = W * B(T, nu) for continuum frequencies on the regridded grid
    print("  Computing J = W * B(T, nu) for continuum...")
    j_grey = compute_grey_j(frequencies, r, t, tau=tau)

    # Show sample J values
    mid_freq = ncf // 2
    print(f"  J at mid-frequency ({frequencies[mid_freq]:.4e} 1e15 Hz):")
    print(f"    J[0]    = {j_grey[mid_freq, 0]:.4e} (surface, W={w[0]:.4f}, T={t[0]:.1f})")
    print(f"    J[nd/2] = {j_grey[mid_freq, nd//2]:.4e} (mid, W={w[nd//2]:.4f}, T={t[nd//2]:.1f})")
    print(f"    J[nd-1] = {j_grey[mid_freq, nd-1]:.4e} (inner, W={w[nd-1]:.4f}, T={t[nd-1]:.1f})")

    # Write modified EDDFACTOR (continuum J replaced, line J preserved)
    # Use the regridded R/V (in 1e10 cm units for R, km/s for V)
    print("  Writing modified EDDFACTOR...")
    write_eddfactor(edd_path, edd_info_path, nd, frequencies, j_grey, r_regrid, v_regrid, lang,
                    line_frequencies=line_freqs if ncf_line > 0 else None,
                    line_j=line_j if ncf_line > 0 else None)

    # Verify
    print("  Verifying...")
    edd_check = read_eddfactor(edd_path, edd_info_path)
    assert edd_check['finish_val'] > 0.5, "FINISH_REC not set"
    assert edd_check['ncf'] == ncf, f"NCF mismatch: {edd_check['ncf']} vs {ncf}"
    assert np.allclose(edd_check['cont_frequencies'], frequencies), "Cont freq mismatch"
    if ncf_line > 0:
        assert edd_check['ncf_line'] == ncf_line, f"Line NCF mismatch: {edd_check['ncf_line']} vs {ncf_line}"
        assert np.allclose(edd_check['line_frequencies'], line_freqs), "Line freq mismatch"
    print(f"  Verification passed: FINISH_REC={edd_check['finish_val']}, "
          f"NCF={edd_check['ncf']}, NCF_LINE={edd_check['ncf_line']}")

    return True


def phase_stella_replace(cmfgen_dir, swd_rad_path, target_days=None, plot_dir=None):
    """Phase 2 (STELLA): Replace J values with STELLA .swd_rad radiation moments.

    Reads the Phase 1 EDDFACTOR (from a 1-iteration CMFGEN run) to get the
    correct CMFGEN frequency grid and depth grid, then interpolates STELLA's
    comoving-frame J_nu onto that grid and writes the modified EDDFACTOR.

    Args:
        cmfgen_dir: CMFGEN run directory with Phase 1 EDDFACTOR
        swd_rad_path: Path to STELLA .swd_rad file
        target_days: Epoch to select from .swd_rad (default: last)
        plot_dir: Optional directory for diagnostic plots
    """
    cmfgen_dir = Path(cmfgen_dir)
    edd_path = cmfgen_dir / "EDDFACTOR"
    edd_info_path = cmfgen_dir / "EDDFACTOR_INFO"

    print("\n=== Phase 2 (STELLA): J_REPLACE from .swd_rad ===")
    print(f"  CMFGEN directory: {cmfgen_dir}")
    print(f"  STELLA .swd_rad: {swd_rad_path}")

    if not edd_path.exists() or not edd_info_path.exists():
        print("  ERROR: Phase 1 EDDFACTOR not found. Run --phase init first.", file=sys.stderr)
        return False

    if not Path(swd_rad_path).exists():
        print(f"  ERROR: .swd_rad file not found: {swd_rad_path}", file=sys.stderr)
        return False

    # Read Phase 1 EDDFACTOR for frequency grid and R/V
    print("  Reading Phase 1 EDDFACTOR...")
    edd = read_eddfactor(edd_path, edd_info_path)
    nd = edd['nd']
    frequencies = edd['cont_frequencies']  # in 1e15 Hz
    ncf = edd['ncf']
    ncf_line = edd['ncf_line']
    line_freqs = edd['line_frequencies']
    line_j = edd['line_j']
    r_edd = edd['r']  # in 1e10 cm, outside-in
    v_edd = edd['v']  # in km/s, outside-in
    print(f"  ND={nd}, NCF={ncf} (continuum), NCF_LINE={ncf_line} (line)")
    print(f"  Cont freq range: {frequencies[0]:.6e} to {frequencies[-1]:.6e} (1e15 Hz)")
    print(f"  R range: {r_edd[0]:.4e} to {r_edd[-1]:.4e} (1e10 cm)")
    print(f"  V range: {v_edd[0]:.1f} to {v_edd[-1]:.1f} km/s")

    # Read STELLA .swd_rad
    print("  Reading STELLA .swd_rad...")
    try:
        from convert_stella_swd_to_edd_factor import read_swd_rad, interpolate_to_cmfgen_grid, interpolate_frequency_grid
    except ImportError:
        sys.path.insert(0, str(Path(__file__).parent))
        from convert_stella_swd_to_edd_factor import read_swd_rad, interpolate_to_cmfgen_grid, interpolate_frequency_grid

    epochs = read_swd_rad(swd_rad_path)
    if not epochs:
        print("  ERROR: No epochs found in .swd_rad", file=sys.stderr)
        return False

    if target_days is not None:
        available = sorted(epochs.keys())
        t_best = min(available, key=lambda t: abs(t - target_days))
    else:
        t_best = max(epochs.keys())
    epoch = epochs[t_best]
    print(f"  Selected epoch: t={epoch.t_day:.4f} days ({epoch.nzon} zones, {epoch.nfrus} freq groups)")

    # Convert STELLA frequencies from Hz to 1e15 Hz for comparison
    stella_freq_1e15 = epoch.freq_hz / 1.0e15
    cmfgen_freq_hz = frequencies * 1.0e15  # back to Hz for interpolation

    # Get STELLA J in physical units (erg/cm^2/s/sr/Hz)
    j_nu_stella = epoch.j_nu_comoving  # (NZON, NFRUS)
    if np.all(j_nu_stella == 0):
        print("  WARNING: FradJ is all zeros — falling back to lab-frame Y_J")
        j_nu_stella = epoch.j_nu_lab

    # Convert EDDFACTOR R from 1e10 cm to cm for interpolation
    r_edd_cm = r_edd * 1.0e10

    # Interpolate STELLA J to CMFGEN depth grid (outer-to-inner)
    j_on_cmfgen_depth = interpolate_to_cmfgen_grid(
        j_nu_stella, epoch.r_cm, r_edd_cm
    )  # (ND, NFRUS)

    # Interpolate to CMFGEN frequency grid
    j_stella_cmfgen = interpolate_frequency_grid(
        j_on_cmfgen_depth, epoch.freq_hz, cmfgen_freq_hz
    )  # (ND, NCF)

    # Transpose to (NCF, ND) for write_eddfactor
    j_for_write = j_stella_cmfgen.T  # (NCF, ND)

    # Clamp negative J to zero
    n_neg = np.sum(j_for_write < 0)
    if n_neg > 0:
        print(f"  WARNING: clamping {n_neg} negative J values to zero")
        j_for_write = np.maximum(j_for_write, 0.0)

    # Report sample J values
    mid_freq = ncf // 2
    print(f"  J at mid-frequency ({frequencies[mid_freq]:.4e} 1e15 Hz):")
    print(f"    J[0]    = {j_for_write[mid_freq, 0]:.4e} (surface)")
    print(f"    J[nd/2] = {j_for_write[mid_freq, nd//2]:.4e} (mid)")
    print(f"    J[nd-1] = {j_for_write[mid_freq, nd-1]:.4e} (inner)")

    # Write modified EDDFACTOR (continuum J replaced, line J preserved)
    lang = np.zeros(nd)
    print("  Writing modified EDDFACTOR with STELLA J...")
    write_eddfactor(edd_path, edd_info_path, nd, frequencies, j_for_write,
                    r_edd, v_edd, lang,
                    line_frequencies=line_freqs if ncf_line > 0 else None,
                    line_j=line_j if ncf_line > 0 else None)

    # Verify
    print("  Verifying...")
    edd_check = read_eddfactor(edd_path, edd_info_path)
    assert edd_check['finish_val'] > 0.5, "FINISH_REC not set"
    assert edd_check['ncf'] == ncf, f"NCF mismatch: {edd_check['ncf']} vs {ncf}"
    assert np.allclose(edd_check['cont_frequencies'], frequencies), "Cont freq mismatch"
    print(f"  Verification passed: FINISH_REC={edd_check['finish_val']}, "
          f"NCF={edd_check['ncf']}, NCF_LINE={edd_check['ncf_line']}")

    # Generate diagnostic plots if requested
    if plot_dir:
        plot_path = Path(plot_dir)
        plot_path.mkdir(parents=True, exist_ok=True)
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            # J_nu profiles at selected frequencies
            fig, ax = plt.subplots(figsize=(10, 6))
            for j in [0, ncf//4, ncf//2, 3*ncf//4, ncf-1]:
                if 0 <= j < j_for_write.shape[0]:
                    ax.plot(r_edd_cm / 1e15, np.log10(np.maximum(j_for_write[j], 1e-300)),
                            label=f"nu={frequencies[j]:.2e} 1e15 Hz")
            ax.set_xlabel("Radius (1e15 cm)")
            ax.set_ylabel("log10 J_nu (erg/cm^2/s/sr/Hz)")
            ax.set_title("STELLA J_nu on CMFGEN grid")
            ax.legend()
            fig.savefig(plot_path / "stella_jnu_on_cmfgen_grid.png", dpi=150)
            plt.close(fig)

            # Compare STELLA T_rad vs CMFGEN grid T
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(epoch.r_cm / 1e15, epoch.t_gas, "b-", label="T_gas (STELLA)")
            ax.plot(epoch.r_cm / 1e15, epoch.t_rad, "r-", label="T_rad (STELLA)")
            ax.set_xlabel("Radius (1e15 cm)")
            ax.set_ylabel("Temperature (K)")
            ax.set_title("STELLA Temperature Profiles")
            ax.legend()
            ax.set_yscale("log")
            fig.savefig(plot_path / "stella_temperature_profiles.png", dpi=150)
            plt.close(fig)

            print(f"  Diagnostic plots written to {plot_dir}")
        except ImportError:
            print("  WARNING: matplotlib not available, skipping plots")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Generate a physics-motivated EDDFACTOR for CMFGEN FIXED_J initialization."
    )
    parser.add_argument("cmfgen_dir", help="CMFGEN run directory (with SN_HYDRO_DATA, VADAT, etc.)")
    parser.add_argument("--phase", choices=["both", "init", "replace", "stella_replace", "stella_both"],
                        default="both",
                        help="Phase: 'init' (1-iter CMFGEN), 'replace' (W*B(T)), "
                             "'stella_replace' (STELLA J), 'stella_both' (init+stella), 'both' (init+replace)")
    parser.add_argument("--source-eddfactor", default=None,
                        help="Source EDDFACTOR for replace phase (default: <cmfgen_dir>/EDDFACTOR)")
    parser.add_argument("--swd-rad", default=None,
                        help="Path to STELLA .swd_rad file (for stella_replace/stella_both phases)")
    parser.add_argument("--target-days", type=float, default=None,
                        help="Epoch to select from .swd_rad (default: last)")
    parser.add_argument("--plot-dir", default=None,
                        help="Directory for diagnostic plots")
    parser.add_argument("--timeout", type=int, default=30,
                        help="Timeout for CMFGEN init run in minutes (default: 30)")
    args = parser.parse_args()

    cmfgen_dir = Path(args.cmfgen_dir)
    if not cmfgen_dir.is_dir():
        print(f"ERROR: {cmfgen_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    success = True

    if args.phase in ("both", "init", "stella_both"):
        success = phase_init(cmfgen_dir, timeout_minutes=args.timeout)
        if not success:
            print("\nPhase 1 failed.", file=sys.stderr)
            sys.exit(1)

    if args.phase in ("both", "replace"):
        success = phase_replace(cmfgen_dir, source_eddfactor=args.source_eddfactor)
        if not success:
            print("\nPhase 2 failed.", file=sys.stderr)
            sys.exit(1)

    if args.phase in ("stella_replace", "stella_both"):
        if not args.swd_rad:
            print("ERROR: --swd-rad required for stella_replace/stella_both phases", file=sys.stderr)
            sys.exit(1)
        success = phase_stella_replace(cmfgen_dir, args.swd_rad,
                                       target_days=args.target_days,
                                       plot_dir=args.plot_dir)
        if not success:
            print("\nPhase 2 (STELLA) failed.", file=sys.stderr)
            sys.exit(1)

    if success:
        if args.phase in ("stella_replace", "stella_both"):
            print("\n=== Success ===")
            print("EDDFACTOR generated with STELLA J_nu from .swd_rad.")
            print("Now set USE_FIXED_J=T in VADAT and run CMFGEN.")
            print("After populations converge, rerun with USE_FIXED_J=F for final J.")
        else:
            print("\n=== Success ===")
            print("EDDFACTOR generated with physics-motivated J = W * B(T, nu).")
            print("Now set USE_FIXED_J=T in VADAT and run CMFGEN.")
            print("After populations converge, rerun with USE_FIXED_J=F for final J.")


if __name__ == "__main__":
    main()
