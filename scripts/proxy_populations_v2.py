
#!/usr/bin/env python3
"""
proxy_populations_v2.py
Complete replacement script with advanced non-monotonic optimizations & diagnostics.
Generates proxy CMFGEN _IN population files.

NEW OPTIMIZATIONS FOR NON-MONOTONIC FLOWS:
1. Saha Inversion for n_e: Decouples ionization from T_gas by inverting the Saha 
   equation to find T_eff_ion that exactly reproduces the STELLA electron density.
2. Non-Local Radiation in Reversals: Replaces unphysical local T_rad and W in 
   velocity reversal zones with spatially averaged values from monotonic shells.
3. Two-Temperature Approach: Uses T_eff_ion for ionization balance, and a 
   tau-blended T_exc for level excitation.
4. Advanced Diagnostics: Generates targeted plots to visualize reversals, 
   temperature decoupling, n_e consistency, and DC caps.
"""
from __future__ import annotations
import argparse
import dataclasses
import json
import math
import os
import re
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence, Tuple
import numpy as np

try:
    from scipy.special import logsumexp
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

def logsumexp(x, axis=None):
    x = np.asarray(x, dtype=float)
    m = np.max(x, axis=axis, keepdims=True)
    m_safe = np.where(np.isfinite(m), m, 0.0)
    out = np.squeeze(m_safe, axis=axis) + np.log(
        np.squeeze(np.sum(np.exp(x - m_safe), axis=axis, keepdims=True), axis=axis)
    )
    if axis is None:
        return float(out)
    return np.asarray(out)

# ===========================================================================
# Constants
# ===========================================================================
HDKT = 4.7994145
SAHA_CONST = 4.828e15
LTEPOP_C = 2.07078e-22
SPEED_LIGHT = 2.9979248e10
M_H = 1.6733e-24
K_B = 1.3807e-16
H_PLANCK = 6.6261e-27
SIGMA_T = 6.6524e-25
SIGMA_SB = 5.6704e-5
A_RAD = 7.5657e-15
L_SUN = 3.828e33
R_SUN_CM = 6.957e10
K_B_OVER_HC_CM = 0.695

FLOOR_T_1E4K = 1e-3
FLOOR_NE = 1.0
FLOOR_DI = 1e-300
FLOOR_DC = 1e-300
FLOOR_POP = 1e-300
FLOOR_GION = 1.0
FLOOR_FRAC = 1e-100
LOG_CLIP = 700.0

# ===========================================================================
# Conventions
# ===========================================================================
@dataclass
class CmfgenConventions:
    """Single source of truth for CMFGEN file conventions and proxy options."""
    depth_order_outer_to_inner: bool = True
    di_is_next_higher_stage: bool = True
    dc_linear_in_file: bool = True
    radius_unit_1e10_cm: bool = True
    temperature_unit_1e4_K: bool = True
    ed_unit_cm_inv3: bool = True
    di_unit_cm_inv3: bool = True

    write_level_count: str = "NF"
    t_in_file_is_physical: bool = False
    strategy: str = "LTE"
    write_t_ion_to_in: bool = False
    smooth: bool = True
    min_ion_frac: float = 1.0e-100
    t_exc_min: float = 0.3
    tau_trans: float = 1.0
    population_tolerance: float = 1e-6
    require_vadat_for_nlte: bool = True

    ne_max_iter: int = 50
    ne_tol: float = 1e-3

    target_lte_est: Optional[bool] = None
    target_dc_meth: Optional[str] = None

    use_w_in_ltepop: bool = False
    use_t_ion: bool = False
    w_apply_mode: str = "first"

    kappa_fallback: str = "thomson"
    kappa_bf_scale: float = 1.0e-2
    tau_mode: str = "plane"
    velocity_vth_cm_s: float = 1.0e6
    velocity_dvdr_floor: float = 1.0e-12
    velocity_trapping_factor: float = 5.0
    velocity_min_factor: float = 0.05
    velocity_aware_smoothing: bool = False

    temperature_aware_smoothing: bool = False
    temperature_log_slope_threshold: float = 0.5

    t_exc_source: str = "gas"
    t_hybrid_tau_scale: float = 1.0

    smooth_mode: str = "uniform"
    front_log_slope_threshold: float = 2.0

    ne_damping_mode: str = "fixed"
    ne_max_log_step: float = 0.0
    floor_ne: float = 1.0

    t_ion_source: str = "gas"
    stella_file: Optional[str] = None
    stella_radius_column: Optional[str] = None
    stella_velocity_column: Optional[str] = None
    stella_tau_column: Optional[str] = None
    stella_t_rad_column: Optional[str] = None
    stella_t_ion_column: Optional[str] = None
    stella_luminosity_column: Optional[str] = None
    stella_energy_density_column: Optional[str] = None
    stella_radius_unit: str = "cm"
    stella_temperature_unit: str = "K"
    stella_luminosity_unit: str = "erg/s"
    stella_interpolation_coordinate: str = "radius"
    stella_fallback: str = "grey"
    stella_t_ion_blend_mode: str = "linear"
    stella_t_ion_floor: float = 1.0e-3
    stella_t_ion_ceiling: float = 0.0

    # NEW OPTIMIZATIONS FOR NON-MONOTONIC FLOWS
    saha_inversion_for_ne: bool = False
    non_local_reversals: bool = False

    def validate(self) -> None:
        if self.write_level_count not in ("NF", "NS"):
            raise ValueError(f"write_level_count must be NF or NS, got {self.write_level_count!r}")
        if self.strategy not in ("LTE", "GREY", "CONTINUATION"):
            raise ValueError(f"unknown strategy {self.strategy!r}")
        if self.kappa_fallback not in ("thomson", "kramers"):
            raise ValueError("kappa_fallback must be thomson or kramers")
        if self.tau_mode not in ("plane", "sobolev", "nonmonotonic"):
            raise ValueError("tau_mode must be plane, sobolev, or nonmonotonic")
        if self.smooth_mode not in ("uniform", "front_aware"):
            raise ValueError("smooth_mode must be uniform or front_aware")
        if self.ne_damping_mode not in ("fixed", "adaptive"):
            raise ValueError("ne_damping_mode must be fixed or adaptive")
        if self.t_exc_source not in ("gas", "t_ion", "rad", "hybrid"):
            raise ValueError("t_exc_source must be gas, t_ion, rad, or hybrid")
        if self.t_ion_source not in ("gas", "grey", "stella", "stella-direct"):
            raise ValueError("t_ion_source must be gas, grey, stella, or stella-direct")
        if self.stella_radius_unit not in ("cm", "1e10cm", "rsun", "km"):
            raise ValueError("stella_radius_unit must be cm, 1e10cm, rsun, or km")
        if self.stella_temperature_unit not in ("K", "1e4K", "eV"):
            raise ValueError("stella_temperature_unit must be K, 1e4K, or eV")
        if self.stella_luminosity_unit not in ("erg/s", "lsun"):
            raise ValueError("stella_luminosity_unit must be erg/s or lsun")
        if self.stella_interpolation_coordinate not in ("radius", "velocity"):
            raise ValueError("stella_interpolation_coordinate must be radius or velocity")
        if self.stella_fallback not in ("grey", "gas", "error"):
            raise ValueError("stella_fallback must be grey, gas, or error")
        if self.stella_t_ion_blend_mode not in ("linear", "gas_floor"):
            raise ValueError("stella_t_ion_blend_mode must be linear or gas_floor")
        if self.w_apply_mode not in ("first", "all", "sqrt", "none"):
            raise ValueError("w_apply_mode must be first, all, sqrt, or none")
        if self.floor_ne <= 0:
            raise ValueError("floor_ne must be > 0")
        if self.front_log_slope_threshold <= 0:
            raise ValueError("front_log_slope_threshold must be > 0")
        if self.temperature_log_slope_threshold <= 0:
            raise ValueError("temperature_log_slope_threshold must be > 0")
        if self.velocity_vth_cm_s <= 0:
            raise ValueError("velocity_vth_cm_s must be > 0")
        if self.velocity_dvdr_floor <= 0:
            raise ValueError("velocity_dvdr_floor must be > 0")
        if self.velocity_trapping_factor < 1:
            raise ValueError("velocity_trapping_factor must be >= 1")
        if not (0 < self.velocity_min_factor <= 1):
            raise ValueError("velocity_min_factor must be in (0, 1]")
        if self.t_hybrid_tau_scale <= 0:
            raise ValueError("t_hybrid_tau_scale must be > 0")

# ===========================================================================
# Element/species tables
# ===========================================================================
ELEMENT_INFO: Dict[str, Tuple[int, float]] = {
    'HYD': (1, 1.008), 'HE': (2, 4.003), 'CARB': (6, 12.011), 'NIT': (7, 14.007),
    'OXY': (8, 15.999), 'NEON': (10, 20.180), 'FLUO': (9, 18.998),
    'SOD': (11, 22.990), 'MASS': (12, 24.305), 'ALUM': (13, 26.982),
    'SIL': (14, 28.085), 'PHOS': (15, 30.974), 'SULF': (16, 32.065),
    'CHLO': (17, 35.453), 'ARGO': (18, 39.948), 'POTA': (19, 39.098),
    'CAL': (20, 40.078), 'SCAN': (21, 44.956), 'TITA': (22, 47.867),
    'VANA': (23, 50.942), 'CHRO': (24, 51.996), 'MANG': (25, 54.938),
    'IRON': (26, 55.845), 'COBA': (27, 58.933), 'NICK': (28, 58.693),
    'COPP': (29, 63.546), 'ZINC': (30, 65.380),
}

_PREFIX_TO_LABEL: Dict[str, str] = {
    'H': 'HYD', 'He': 'HE', 'C': 'CARB', 'N': 'NIT', 'O': 'OXY', 'Ne': 'NEON',
    'F': 'FLUO', 'Na': 'SOD', 'Mg': 'MASS', 'Al': 'ALUM', 'Sk': 'SIL', 'Si': 'SIL',
    'P': 'PHOS', 'S': 'SULF', 'Cl': 'CHLO', 'Ar': 'ARGO', 'K': 'POTA', 'Ca': 'CAL',
    'Sc': 'SCAN', 'Ti': 'TITA', 'V': 'VANA', 'Cr': 'CHRO', 'Mn': 'MANG', 'Fe': 'IRON',
    'Co': 'COBA', 'Ni': 'NICK', 'Cu': 'COPP', 'Zn': 'ZINC',
}

_ION_SUFFIX: Dict[str, int] = {
    'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8,
    'IX': 9, 'X': 10, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7,
    '8': 8, '9': 9, '10': 10,
}

ROMAN = {1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V', 6: 'VI', 7: 'VII', 8: 'VIII', 9: 'IX', 10: 'X'}

def parse_species_name(name: str) -> Tuple[Optional[str], Optional[int]]:
    m = re.match(r'^([A-Z][a-z]?)([IVX]+|\d{1,2})$', name)
    if not m:
        m = re.match(r'^([A-Z])([IVX]+|\d{1,2})$', name)
    if not m:
        return None, None
    suffix = m.group(2)
    if suffix not in _ION_SUFFIX:
        return None, None
    return m.group(1), _ION_SUFFIX[suffix]

def species_info(name: str) -> Tuple[str, float, int]:
    elem, ion_stage = parse_species_name(name)
    if elem is None or elem not in _PREFIX_TO_LABEL:
        raise KeyError(f"cannot parse element for species {name!r}")
    label = _PREFIX_TO_LABEL[elem]
    if label not in ELEMENT_INFO:
        raise KeyError(f"element label {label} for {name!r} not in ELEMENT_INFO")
    _, mass = ELEMENT_INFO[label]
    return label, mass, ion_stage

# ===========================================================================
# Dataclasses for inputs and diagnostics
# ===========================================================================
@dataclass
class HydroData:
    nd: int
    time_days: float
    lum_lsun: float
    r_inner_cm: float
    t_inner_k: float
    radius: np.ndarray
    velocity: np.ndarray
    sigma: np.ndarray
    temperature: np.ndarray
    density: np.ndarray
    ed: np.ndarray
    natom: np.ndarray
    kappa: np.ndarray
    mass_fracs: Dict[str, np.ndarray]
    raw_path: str = ""

@dataclass
class SpeciesRow:
    name: str
    n_f: int
    n_s: int
    z: float
    gion: float
    mass: float
    at_no: float
    n_v: int = 0

@dataclass
class AtomicData:
    name: str
    n_levels: int
    ionization_energy: float
    stat_wt: np.ndarray
    energy: np.ndarray
    fedge: np.ndarray
    transitions: List[Dict[str, float]] = field(default_factory=list)
    f_to_s: Optional[np.ndarray] = None
    int_seq: Optional[np.ndarray] = None

@dataclass
class ModelInputs:
    hydro: HydroData
    species: List[SpeciesRow]
    vadat_abundances: Dict[str, float]
    vadat_flags: Dict[str, str]
    atom_data: Dict[str, AtomicData]
    gammas_in_ed: Optional[np.ndarray] = None
    model_dir: str = ""

@dataclass
class StellaRadiationField:
    source: str
    time_days: Optional[float] = None
    radius: Optional[np.ndarray] = None
    velocity: Optional[np.ndarray] = None
    tau: Optional[np.ndarray] = None
    t_rad: Optional[np.ndarray] = None
    t_ion: Optional[np.ndarray] = None
    luminosity: Optional[np.ndarray] = None
    energy_density: Optional[np.ndarray] = None
    notes: Dict[str, Any] = field(default_factory=dict)

@dataclass
class VelocityFieldDiagnostics:
    nonmonotonic: bool
    n_extrema: int
    extrema_mask: np.ndarray
    small_gradient_mask: np.ndarray
    dvds_signed: np.ndarray
    dvds_abs: np.ndarray
    sobolev_length: np.ndarray

@dataclass
class TemperatureProfileDiagnostics:
    n_inversions: int
    n_extrema: int
    inversion_mask: np.ndarray
    extrema_mask: np.ndarray
    large_gradient_mask: np.ndarray
    dTds: np.ndarray
    log_gradient: np.ndarray
    temperature_break_mask: np.ndarray

# ===========================================================================
# IO: CMFGEN / hydro / atomic data
# ===========================================================================
def _read_array(lines: List[str], start: int, nd: int) -> np.ndarray:
    vals: List[float] = []
    i = start + 1
    while len(vals) < nd:
        if i >= len(lines):
            raise ValueError(f"unexpected EOF reading {nd} values from line {start}")
        for p in lines[i].split():
            vals.append(float(p))
        if len(vals) >= nd:
            break
        i += 1
    return np.asarray(vals, dtype=float)

def _find_section(lines: List[str], name: str) -> int:
    for i, line in enumerate(lines):
        if line.strip().lower().startswith(name.lower()):
            return i
    return -1

def parse_sn_hydro_data(filepath: str) -> HydroData:
    with open(filepath) as f:
        lines = f.readlines()
    if not lines:
        raise ValueError(f"empty SN_HYDRO_DATA: {filepath}")
    
    m0 = re.search(r'\d+', lines[0])
    if not m0:
        raise ValueError(f"cannot read ND from header: {lines[0]!r}")
    nd = int(m0.group())
    
    m1 = re.search(r'\d+', lines[1])
    if not m1:
        raise ValueError(f"cannot read NMFRAC from line 1: {lines[1]!r}")
    
    time_days = float(lines[3].split()[-1])
    lum_lsun = float(lines[4].split()[-1])
    r_inner_cm = float(lines[5].split()[-1])
    t_inner_k = float(lines[6].split()[-1])

    def arr(name: str) -> np.ndarray:
        idx = _find_section(lines, name)
        if idx < 0:
            raise ValueError(f"section {name!r} not found in {filepath}")
        return _read_array(lines, idx, nd)

    radius = arr('Radius grid')
    velocity = arr('Velocity')
    sigma = arr('Sigma')
    temperature = arr('Temperature')
    density = arr('Density')
    ed = arr('Electron density')
    natom = arr('Atom density')
    
    kappa_idx = _find_section(lines, 'Kappa')
    kappa = _read_array(lines, kappa_idx, nd) if kappa_idx >= 0 else np.zeros(nd)

    elem_labels = [
        'HYD', 'HE', 'CARB', 'NIT', 'OXY', 'NEON', 'FLUO', 'SOD',
        'MASS', 'ALUM', 'SIL', 'PHOS', 'SULF', 'CHLO', 'ARGO', 'POTA',
        'CAL', 'SCAN', 'TITA', 'VANA', 'CHRO', 'MANG', 'IRON', 'COBA',
        'NICK', 'COPP', 'ZINC', 'NICK 56', 'COB 56'
    ]
    mass_fracs: Dict[str, np.ndarray] = {}
    for label in elem_labels:
        idx = _find_section(lines, f'{label} mass fraction')
        if idx >= 0:
            mass_fracs[label] = _read_array(lines, idx, nd)

    return HydroData(
        nd=nd, time_days=time_days, lum_lsun=lum_lsun, r_inner_cm=r_inner_cm,
        t_inner_k=t_inner_k, radius=radius, velocity=velocity, sigma=sigma,
        temperature=temperature, density=density, ed=ed, natom=natom,
        kappa=kappa, mass_fracs=mass_fracs, raw_path=str(filepath),
    )

def parse_vadat_abundances(vadat_path: str) -> Dict[str, float]:
    abundances: Dict[str, float] = {}
    if not os.path.exists(vadat_path):
        return abundances
    with open(vadat_path) as f:
        for line in f:
            m = re.match(r'\s*(-?[\d.EDed+-]+)\s+\[(\w+)/X\]', line)
            if m:
                val = float(m.group(1).replace('D', 'E').replace('d', 'e'))
                abundances[m.group(2)] = -val
    return abundances

def parse_vadat_flags(vadat_path: str) -> Dict[str, str]:
    flags: Dict[str, str] = {}
    if not os.path.exists(vadat_path):
        return flags
    with open(vadat_path) as f:
        for line in f:
            s = line.strip()
            if '[' in s and ']' in s:
                key = s[s.index('[') + 1:s.index(']')].strip()
                val = s[:s.index('[')].strip().split()
                flags[key] = val[0] if val else ""
    return flags

def parse_model_species(filepath: str) -> List[SpeciesRow]:
    with open(filepath) as f:
        lines = f.readlines()
    out: List[SpeciesRow] = []
    for line in lines:
        m = re.match(
            r'\s+([A-Z][a-zA-Z0-9]*)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+'
            r'([\d.]+)\s+([\d.]+)', line
        )
        if not m:
            continue
        name = m.group(1)
        if not re.match(r'^[A-Z][a-z]*[IVX]*$', name) and not re.match(r'^[A-Z][a-z]*\d+$', name):
            continue
        out.append(SpeciesRow(
            name=name, n_f=int(m.group(2)), n_s=int(m.group(3)),
            z=float(m.group(4)), gion=float(m.group(5)),
            mass=float(m.group(6)), at_no=float(m.group(7)),
        ))
    return out

def parse_model_spec_species(filepath: str) -> List[SpeciesRow]:
    with open(filepath) as f:
        lines = f.readlines()
    out: List[SpeciesRow] = []
    for line in lines:
        m = re.match(r'(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s+\[(\w+)_ISF\]', line)
        if not m:
            continue
        name = m.group(4)
        try:
            label, mass, ion_stage = species_info(name)
            z = float(ELEMENT_INFO[label][0])
        except KeyError:
            continue
        out.append(SpeciesRow(
            name=name, n_f=int(m.group(3)), n_s=int(m.group(2)),
            z=float(ion_stage), gion=1.0, mass=mass, at_no=z, n_v=int(m.group(1)),
        ))
    return out

def parse_f_oscdat(filepath: str, name: str = "") -> AtomicData:
    with open(filepath) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if '!Date' in lines[i]:
            i += 1
            while i < len(lines) and '!Date' in lines[i]:
                i += 1
            break
        i += 1
    if i >= len(lines):
        raise ValueError(f"no !Date line found in {filepath}")
    
    n_levels = int(lines[i].split()[0])
    i += 1
    ionization_en = float(lines[i].split()[0].replace('D', 'E').replace('d', 'e'))
    i += 1
    i += 1  # screened nuclear charge
    i += 1  # number of transitions
    
    while i < len(lines) and ('!' in lines[i] or lines[i].strip() == ''):
        i += 1
        
    stat_wt = np.zeros(n_levels)
    energy = np.zeros(n_levels)
    fedge = np.zeros(n_levels)
    
    for j in range(n_levels):
        parts = lines[i].split()
        if len(parts) < 3:
            raise ValueError(f"malformed level {j+1} in {filepath}: {lines[i]!r}")
        stat_wt[j] = float(parts[1])
        energy[j] = float(parts[2].replace('D', 'E').replace('d', 'e'))
        fedge[j] = (ionization_en - energy[j]) * SPEED_LIGHT * 1.0e-15
        i += 1
        
    return AtomicData(
        name=name, n_levels=n_levels, ionization_energy=ionization_en,
        stat_wt=stat_wt, energy=energy, fedge=fedge, transitions=[],
    )

def parse_f_to_s(filepath: str, n_f_used: int) -> Tuple[np.ndarray, np.ndarray]:
    with open(filepath) as f:
        lines = f.readlines()
    i = 0
    dates = 0
    while i < len(lines) and dates < 2:
        if '!Date' in lines[i]:
            dates += 1
        i += 1
        
    f_to_s: List[int] = []
    int_seq: List[int] = []
    for line in lines[i:]:
        s = line.strip()
        if not s or s.startswith('!') or s.startswith('*'):
            continue
        parts = s.split()
        if len(parts) < 8:
            continue
        try:
            sl = int(parts[-3])
            iseq = int(parts[-2])
        except (ValueError, IndexError):
            continue
        f_to_s.append(sl)
        int_seq.append(iseq)
        if len(f_to_s) >= n_f_used:
            break
            
    if len(f_to_s) < n_f_used:
        raise ValueError(f"F_TO_S {filepath} has {len(f_to_s)} levels, need {n_f_used}")
        
    return np.asarray(f_to_s, dtype=int), np.asarray(int_seq, dtype=int)

def parse_gammas_in_ed(gammas_path: str, nd: int) -> Optional[np.ndarray]:
    if not os.path.exists(gammas_path):
        return None
    with open(gammas_path) as f:
        lines = f.readlines()
    start = -1
    for i, line in enumerate(lines):
        if 'Electron density' in line:
            start = i + 1
            break
    if start < 0:
        return None
        
    vals: List[float] = []
    i = start
    while len(vals) < nd and i < len(lines):
        s = lines[i].strip()
        if s.startswith('!'):
            break
        for p in s.split():
            try:
                vals.append(float(p))
            except ValueError:
                break
        if len(vals) >= nd:
            break
        i += 1
    if len(vals) < nd:
        return None
    return np.asarray(vals, dtype=float)


def update_gammas_in_ne(gammas_path: str, ne_new: np.ndarray, nd: int) -> bool:
    """Patch GAMMAS_IN electron density section with self-consistent Ne.

    Reads the existing GAMMAS_IN, replaces the electron density values with
    ne_new (computed from self-consistent Saha/proxy ionization), and rewrites
    the file. All other sections (Gamma rates) are preserved.
    """
    if not os.path.exists(gammas_path):
        return False
    with open(gammas_path) as f:
        lines = f.readlines()

    # Find electron density section boundaries
    ed_start = -1
    ed_end = -1
    for i, line in enumerate(lines):
        if 'Electron density' in line:
            ed_start = i + 1
            break
    if ed_start < 0:
        return False

    # Find end of ED section (next comment line starting with !)
    vals_count = 0
    ed_end = ed_start
    for i in range(ed_start, len(lines)):
        s = lines[i].strip()
        if s.startswith('!') and vals_count >= nd:
            ed_end = i
            break
        try:
            parts = s.split()
            for p in parts:
                float(p)
                vals_count += 1
        except ValueError:
            ed_end = i
            break
        ed_end = i + 1
        if vals_count >= nd:
            ed_end = i + 1
            break

    # Build new ED section (10 per line, same format as original)
    new_ed_lines = []
    for i in range(0, nd, 10):
        chunk = ne_new[i:i+10]
        new_ed_lines.append('   ' + '   '.join(f'{v:10.4E}' for v in chunk) + '\n')

    # Reassemble file
    out = lines[:ed_start] + new_ed_lines + lines[ed_end:]
    with open(gammas_path, 'w') as f:
        f.writelines(out)

    print(f"  Updated GAMMAS_IN Ne: range {ne_new.min():.4e} to {ne_new.max():.4e}")
    return True

def load_model(model_dir: str, hydro_path: Optional[str] = None,
               conv: Optional[CmfgenConventions] = None) -> ModelInputs:
    model_dir = os.path.abspath(model_dir)
    if hydro_path is None:
        hydro_path = os.path.join(model_dir, 'SN_HYDRO_DATA')
    if not os.path.exists(hydro_path):
        alt = os.path.join(os.path.dirname(os.path.dirname(model_dir)), 'SN_HYDRO_DATA')
        if os.path.exists(alt):
            hydro_path = alt
        else:
            raise FileNotFoundError(f"SN_HYDRO_DATA not found at {hydro_path}")
            
    hydro = parse_sn_hydro_data(hydro_path)
    model_path = os.path.join(model_dir, 'MODEL')
    if os.path.exists(model_path):
        species = parse_model_species(model_path)
    else:
        species = parse_model_spec_species(os.path.join(model_dir, 'MODEL_SPEC'))
        
    vadat_path = os.path.join(model_dir, 'VADAT')
    vadat_abundances = parse_vadat_abundances(vadat_path)
    vadat_flags = parse_vadat_flags(vadat_path)
    
    atom_data: Dict[str, AtomicData] = {}
    for sp in species:
        oscdat = os.path.join(model_dir, f"{sp.name}_F_OSCDAT")
        if os.path.exists(oscdat):
            ad = parse_f_oscdat(oscdat, name=sp.name)
            f_to_s_path = os.path.join(model_dir, f"{sp.name}_F_TO_S")
            if os.path.exists(f_to_s_path) and sp.n_f > 0:
                try:
                    f_to_s, int_seq = parse_f_to_s(f_to_s_path, sp.n_f)
                    ad.f_to_s = f_to_s
                    ad.int_seq = int_seq
                except Exception:
                    pass
            atom_data[sp.name] = ad
            
    gammas_ed = parse_gammas_in_ed(os.path.join(model_dir, 'GAMMAS_IN'), hydro.nd)
    
    if conv is not None:
        conv.target_lte_est = (vadat_flags.get('LTE_EST', 'F').upper() == 'T')
        conv.target_dc_meth = vadat_flags.get('DC_METH', '')
        
    return ModelInputs(
        hydro=hydro, species=species, vadat_abundances=vadat_abundances,
        vadat_flags=vadat_flags, atom_data=atom_data, gammas_in_ed=gammas_ed,
        model_dir=model_dir,
    )

# ===========================================================================
# STELLA radiation-field adapter
# ===========================================================================
def _stella_float(tok: str) -> float:
    return float(tok.replace('D', 'E').replace('d', 'e'))

def parse_stella_radiation_file(
    path: str, radius_column: Optional[str] = None, velocity_column: Optional[str] = None,
    tau_column: Optional[str] = None, t_rad_column: Optional[str] = None,
    t_ion_column: Optional[str] = None, luminosity_column: Optional[str] = None,
    energy_density_column: Optional[str] = None, radius_unit: str = "cm",
    temperature_unit: str = "K", luminosity_unit: str = "erg/s",
    time_days: Optional[float] = None,
) -> StellaRadiationField:
    if not path: raise ValueError("STELLA file path is empty")
    if not os.path.exists(path): raise FileNotFoundError(f"STELLA file not found: {path}")
    
    with open(path) as f:
        raw_lines = f.readlines()
    clean: List[str] = []
    for line in raw_lines:
        s = line.strip()
        if not s or s.startswith('#'): continue
        clean.append(s)
    if not clean: raise ValueError(f"STELLA file is empty: {path}")
    
    header = clean[0].split()
    colidx = {name: i for i, name in enumerate(header)}
    rows: List[List[float]] = []
    for s in clean[1:]:
        toks = s.split()
        try: rows.append([_stella_float(t) for t in toks])
        except ValueError: continue
    if not rows: raise ValueError(f"STELLA file has no numeric rows: {path}")
    
    def col(name: Optional[str]) -> Optional[np.ndarray]:
        if not name or name not in colidx: return None
        idx = colidx[name]
        vals = [r[idx] for r in rows if len(r) > idx]
        return np.asarray(vals, dtype=float) if vals else None

    radius = col(radius_column); velocity = col(velocity_column)
    tau = col(tau_column); t_rad = col(t_rad_column)
    t_ion = col(t_ion_column); luminosity = col(luminosity_column)
    energy_density = col(energy_density_column)
    
    arrays = [a for a in [radius, velocity, tau, t_rad, t_ion, luminosity, energy_density] if a is not None]
    if arrays:
        nmin = min(len(a) for a in arrays)
        radius = radius[:nmin] if radius is not None else None
        velocity = velocity[:nmin] if velocity is not None else None
        tau = tau[:nmin] if tau is not None else None
        t_rad = t_rad[:nmin] if t_rad is not None else None
        t_ion = t_ion[:nmin] if t_ion is not None else None
        luminosity = luminosity[:nmin] if luminosity is not None else None
        energy_density = energy_density[:nmin] if energy_density is not None else None

    if radius is not None:
        if radius_unit == "cm": radius = radius / 1.0e10
        elif radius_unit == "1e10cm": pass
        elif radius_unit == "rsun": radius = radius * R_SUN_CM / 1.0e10
        elif radius_unit == "km": radius = radius * 1.0e5 / 1.0e10

    for arr, unit in [(t_rad, temperature_unit), (t_ion, temperature_unit)]:
        if arr is not None:
            if unit == "K": arr /= 1.0e4
            elif unit == "eV": arr = (arr * 11604.5) / 1.0e4

    if luminosity is not None and luminosity_unit == "lsun":
        luminosity = luminosity * L_SUN

    return StellaRadiationField(
        source=path, time_days=time_days, radius=radius, velocity=velocity,
        tau=tau, t_rad=t_rad, t_ion=t_ion, luminosity=luminosity,
        energy_density=energy_density, notes={"header": header, "nrows": len(rows)}
    )

def _interp_to_grid(x_target: np.ndarray, x_src: np.ndarray, y_src: np.ndarray) -> np.ndarray:
    mask = np.isfinite(x_src) & np.isfinite(y_src)
    x = np.asarray(x_src[mask], dtype=float)
    y = np.asarray(y_src[mask], dtype=float)
    if x.size == 0: raise ValueError("no finite STELLA interpolation points")
    order = np.argsort(x)
    x, y = x[order], y[order]
    ux, idx = np.unique(x, return_index=True)
    return np.interp(x_target, ux, y[idx], left=y[idx[0]], right=y[idx[-1]])

def interpolate_stella_to_hydro(stella: StellaRadiationField, hydro: HydroData, coordinate: str = "radius") -> StellaRadiationField:
    if coordinate == "radius":
        if stella.radius is None: raise ValueError("STELLA radius column not available")
        x_src, x_target = stella.radius, hydro.radius
    elif coordinate == "velocity":
        if stella.velocity is None: raise ValueError("STELLA velocity column not available")
        x_src, x_target = stella.velocity, hydro.velocity
    else: raise ValueError("coordinate must be radius or velocity")

    def interp(arr): return _interp_to_grid(x_target, x_src, arr) if arr is not None else None
    
    return StellaRadiationField(
        source=stella.source, time_days=stella.time_days,
        radius=interp(stella.radius), velocity=interp(stella.velocity),
        tau=interp(stella.tau), t_rad=interp(stella.t_rad),
        t_ion=interp(stella.t_ion), luminosity=interp(stella.luminosity),
        energy_density=interp(stella.energy_density), notes=stella.notes
    )

def t_rad_from_energy_density(energy_density: np.ndarray) -> np.ndarray:
    return (np.maximum(energy_density, 0.0) / A_RAD) ** 0.25 / 1.0e4

def t_rad_from_luminosity_radius(luminosity: np.ndarray, radius_cm: np.ndarray) -> np.ndarray:
    r = np.maximum(radius_cm, 1.0)
    l = np.maximum(luminosity, 0.0)
    return (l / (4.0 * np.pi * r * r * SIGMA_SB)) ** 0.25 / 1.0e4

# ===========================================================================
# Velocity and temperature diagnostics
# ===========================================================================
def compute_velocity_diagnostics(hydro: HydroData, velocity_dvdr_floor: float = 1.0e-12, velocity_vth_cm_s: float = 1.0e6) -> VelocityFieldDiagnostics:
    nd = hydro.nd
    r_cm = hydro.radius * 1.0e10
    v_cm_s = hydro.velocity * 1.0e5
    s, ds, dv = np.zeros(nd), np.zeros(nd), np.zeros(nd)
    for i in range(1, nd):
        dr = max(abs(r_cm[i] - r_cm[i - 1]), 1.0)
        s[i] = s[i - 1] + dr
        ds[i] = dr
        dv[i] = v_cm_s[i] - v_cm_s[i - 1]
        
    dvds = np.zeros(nd)
    if nd > 1:
        interval = np.zeros(nd)
        interval[1:] = dv[1:] / np.maximum(ds[1:], 1.0)
        dvds[0], dvds[-1] = interval[1], interval[-1]
        if nd > 2: dvds[1:-1] = 0.5 * (interval[1:-1] + interval[2:])

    extrema = np.zeros(nd, dtype=bool)
    if nd > 2:
        for i in range(1, nd - 1):
            if dv[i] == 0.0 or dv[i + 1] == 0.0: extrema[i] = True
            elif np.sign(dv[i]) != np.sign(dv[i + 1]): extrema[i] = True

    small_gradient = np.abs(dvds) < velocity_dvdr_floor
    sobolev_length = np.full(nd, np.inf)
    large = np.abs(dvds) > velocity_dvdr_floor
    sobolev_length[large] = velocity_vth_cm_s / np.abs(dvds[large])

    return VelocityFieldDiagnostics(
        nonmonotonic=bool(np.any(extrema)), n_extrema=int(np.sum(extrema)),
        extrema_mask=extrema, small_gradient_mask=small_gradient,
        dvds_signed=dvds, dvds_abs=np.abs(dvds), sobolev_length=sobolev_length
    )

def compute_temperature_diagnostics(hydro: HydroData, temperature_log_slope_threshold: float = 0.5) -> TemperatureProfileDiagnostics:
    nd = hydro.nd
    T = hydro.temperature
    r_cm = hydro.radius * 1.0e10
    inversion_mask = np.zeros(nd, dtype=bool)
    extrema_mask = np.zeros(nd, dtype=bool)
    large_gradient_mask = np.zeros(nd, dtype=bool)
    dTds = np.zeros(nd)
    log_gradient = np.zeros(nd)
    
    if nd < 2:
        return TemperatureProfileDiagnostics(0, 0, inversion_mask, extrema_mask, large_gradient_mask, dTds, log_gradient, np.zeros(nd, dtype=bool))

    s = np.zeros(nd)
    for i in range(1, nd): s[i] = s[i - 1] + max(abs(r_cm[i] - r_cm[i - 1]), 1.0)
    
    dT = np.diff(T)
    ds = np.maximum(np.diff(s), 1.0)
    dT_interval = dT / ds
    dTds[:-1] = dT_interval
    dTds[-1] = dT_interval[-1]
    
    logT = np.log10(np.maximum(T, 1e-10))
    log_grad_interval = np.abs(np.diff(logT))
    log_gradient[:-1] = log_grad_interval
    log_gradient[-1] = log_grad_interval[-1]
    
    for i in range(len(dT_interval)):
        if log_grad_interval[i] > temperature_log_slope_threshold:
            large_gradient_mask[i] = True
            large_gradient_mask[i + 1] = True
        if dT_interval[i] < -1.0e-8:
            inversion_mask[i] = True
            inversion_mask[i + 1] = True
            
    for i in range(1, len(dT_interval)):
        if dT_interval[i - 1] == 0.0 or dT_interval[i] == 0.0: extrema_mask[i] = True
        elif np.sign(dT_interval[i - 1]) != np.sign(dT_interval[i]): extrema_mask[i] = True

    return TemperatureProfileDiagnostics(
        n_inversions=int(np.sum(inversion_mask)), n_extrema=int(np.sum(extrema_mask)),
        inversion_mask=inversion_mask, extrema_mask=extrema_mask,
        large_gradient_mask=large_gradient_mask, dTds=dTds,
        log_gradient=log_gradient, temperature_break_mask=inversion_mask | extrema_mask | large_gradient_mask
    )

# ===========================================================================
# NEW OPTIMIZATIONS: Saha Inversion & Non-Local Reversals
# ===========================================================================
def compute_t_eff_for_ne_target(
    target_ne: np.ndarray,
    abundances: Dict[str, np.ndarray],
    hydro: HydroData,
    chi: Dict[str, List[float]],
    gratio: Dict[str, List[float]],
    bounds: Tuple[float, float] = (2e3, 2e5),
    max_iter: int = 60,
    tol: float = 1e-4
) -> np.ndarray:
    """
    Finds the effective temperature that yields the target electron density
    assuming Saha-Boltzmann ionization balance. Decouples ionization from T_gas.
    """
    nd = len(target_ne)
    t_eff = np.full(nd, 1e4)
    
    n_elem_dict = {}
    for elem, fracs in abundances.items():
        if elem in chi:
            A = ELEMENT_INFO.get(elem, (0, 1.0))[1]
            n_elem_dict[elem] = hydro.density * fracs / (A * M_H)
            
    for i in range(nd):
        ne_target = target_ne[i]
        if ne_target <= 0: continue
            
        t_low, t_high = bounds
        for _ in range(max_iter):
            t_mid = 0.5 * (t_low + t_high)
            ne_pred = 0.0
            
            for elem, n_elem_arr in n_elem_dict.items():
                n_elem = n_elem_arr[i]
                if n_elem <= 0: continue
                
                chi_list = chi[elem]
                gratio_list = gratio.get(elem, [1.0]*len(chi_list))
                
                cumlog = np.zeros(len(chi_list) + 1)
                for j in range(len(chi_list)):
                    kt_cm = K_B_OVER_HC_CM * t_mid
                    arg = -chi_list[j] / max(kt_cm, 1e-30)
                    lr = (math.log(SAHA_CONST) + 1.5 * math.log(t_mid) 
                          - math.log(ne_target) + math.log(max(gratio_list[j], 1e-300)) + arg)
                    cumlog[j+1] = cumlog[j] + lr
                
                max_c = np.max(cumlog)
                exp_c = np.exp(cumlog - max_c)
                sum_c = np.sum(exp_c)
                fracs_i = exp_c / sum_c
                
                charges = np.arange(len(fracs_i))
                ne_pred += n_elem * np.sum(fracs_i * charges)
                
            if ne_pred < ne_target: t_low = t_mid
            else: t_high = t_mid
                
            if (t_high - t_low) / t_mid < tol: break
                
        t_eff[i] = 0.5 * (t_low + t_high)
        
    return t_eff

def apply_non_local_radiation_in_reversals(
    t_rad: np.ndarray, 
    w: np.ndarray, 
    vel_diag: VelocityFieldDiagnostics,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    In velocity reversals, replace local T_rad and W with an average 
    from the surrounding monotonic regions.
    """
    t_rad_nl = t_rad.copy()
    w_nl = w.copy()
    
    if vel_diag is None: return t_rad_nl, w_nl
        
    reversal_mask = vel_diag.extrema_mask | (vel_diag.dvds_signed < 0)
    
    if np.any(reversal_mask):
        monotonic_mask = ~reversal_mask
        if np.any(monotonic_mask):
            monotonic_indices = np.where(monotonic_mask)[0]
            reversal_indices = np.where(reversal_mask)[0]
            
            for idx in reversal_indices:
                inner_m = monotonic_indices[monotonic_indices < idx]
                outer_m = monotonic_indices[monotonic_indices > idx]
                
                t_val, w_val = t_rad[idx], w[idx]
                
                if len(inner_m) > 0 and len(outer_m) > 0:
                    t_val = 0.5 * (t_rad[inner_m[-1]] + t_rad[outer_m[0]])
                    w_val = 0.5 * (w[inner_m[-1]] + w[outer_m[0]])
                elif len(inner_m) > 0:
                    t_val, w_val = t_rad[inner_m[-1]], w[inner_m[-1]]
                elif len(outer_m) > 0:
                    t_val, w_val = t_rad[outer_m[0]], w[outer_m[0]]
                    
                t_rad_nl[idx] = t_val
                w_nl[idx] = w_val
                
    return t_rad_nl, w_nl

# ===========================================================================
# Saha-chain construction
# ===========================================================================
def build_saha_chains(species: List[SpeciesRow], atom_data: Dict[str, AtomicData]) -> Tuple[Dict[str, List[float]], Dict[str, List[float]], Dict[str, Dict[int, float]], Dict[str, Dict[int, Tuple[float, float]]]]:
    by_elem: Dict[str, Dict[int, Tuple[str, SpeciesRow]]] = defaultdict(dict)
    for sp in species:
        try:
            label, _, ion_stage = species_info(sp.name)
        except KeyError: continue
        by_elem[label][ion_stage] = (sp.name, sp)

    chi: Dict[str, List[float]] = {}
    gratio: Dict[str, List[float]] = {}
    gion: Dict[str, Dict[int, float]] = {}
    chi_to_next: Dict[str, Dict[int, Tuple[float, float]]] = {}

    for label, stages in by_elem.items():
        chi_list, gr_list = [], []
        gion[label] = {}
        chi_to_next[label] = {}
        stage = 1
        while stage in stages:
            name, sp = stages[stage]
            ad = atom_data.get(name)
            if ad is None: break
            g_lower = float(ad.stat_wt[0]) if ad.stat_wt.size else 1.0
            if stage + 1 not in stages:
                g_upper = float(sp.gion) if sp.gion > 0 else 1.0
                gion[label][stage] = g_upper
                chi_to_next[label][stage] = (float(ad.ionization_energy), g_upper / g_lower if g_lower > 0 else 1.0)
                stage += 1
                continue
            next_name, next_sp = stages[stage + 1]
            next_ad = atom_data.get(next_name)
            g_upper = float(next_ad.stat_wt[0]) if next_ad is not None and next_ad.stat_wt.size else (float(next_sp.gion) if next_sp.gion > 0 else 1.0)
            chi_list.append(float(ad.ionization_energy))
            gr_list.append(g_upper / g_lower if g_lower > 0 else 1.0)
            gion[label][stage] = g_upper
            chi_to_next[label][stage] = (float(ad.ionization_energy), g_upper / g_lower if g_lower > 0 else 1.0)
            stage += 1
            
        highest_consecutive = stage - 1
        if highest_consecutive >= 1 and highest_consecutive in chi_to_next[label]:
            chi_next, gr_next = chi_to_next[label][highest_consecutive]
            chi_list.append(chi_next)
            gr_list.append(gr_next)
        if chi_list:
            chi[label] = chi_list
            gratio[label] = gr_list

        for stage, (name, sp) in stages.items():
            if stage in chi_to_next[label]: continue
            ad = atom_data.get(name)
            if ad is None: continue
            g_lower = float(ad.stat_wt[0]) if ad.stat_wt.size else 1.0
            g_upper = float(sp.gion) if sp.gion > 0 else 1.0
            if stage not in gion[label]: gion[label][stage] = g_upper
            chi_to_next[label][stage] = (float(ad.ionization_energy), g_upper / g_lower if g_lower > 0 else 1.0)

    return chi, gratio, gion, chi_to_next

# ===========================================================================
# Ionisation: log-space Saha
# ===========================================================================
@dataclass
class SahaResult:
    fracs: np.ndarray
    log_fracs: np.ndarray
    n_elem: np.ndarray
    n_stages: int
    omitted_stage_frac: np.ndarray

def _saha_log_ratios(chi_list, gratio_list, t_K, ne, w_dilute=None, w_apply_mode='first'):
    n_trans = len(chi_list)
    nd = t_K.size
    log_r = np.zeros((n_trans, nd))
    log_T = np.log(np.maximum(t_K, 1.0))
    log_ne = np.log(np.maximum(ne, FLOOR_NE))
    log_w = np.log(np.maximum(w_dilute, FLOOR_DC)) if w_dilute is not None else None

    for j in range(n_trans):
        kt_cm = K_B_OVER_HC_CM * t_K
        arg = -np.asarray(chi_list[j]) / np.maximum(kt_cm, 1e-30)
        log_r[j] = math.log(SAHA_CONST) + 1.5 * log_T - log_ne + math.log(max(gratio_list[j], FLOOR_DC)) + arg
        if log_w is not None:
            if w_apply_mode == 'first' and j == 0: log_r[j] = log_r[j] + log_w
            elif w_apply_mode == 'all': log_r[j] = log_r[j] + log_w
            elif w_apply_mode == 'sqrt': log_r[j] = log_r[j] + 0.5 * log_w
    return log_r

def solve_saha_chain(chi_list, gratio_list, t_K, ne, w_dilute=None, w_apply_mode='first'):
    n_trans = len(chi_list)
    n_stages = n_trans + 1
    nd = t_K.size
    if n_stages == 1:
        ones = np.ones((1, nd))
        return SahaResult(ones, np.zeros((1, nd)), np.zeros(nd), 1, np.zeros(nd))

    log_r = _saha_log_ratios(chi_list, gratio_list, t_K, ne, w_dilute, w_apply_mode)
    cumlog = np.zeros((n_stages, nd))
    for j in range(1, n_stages): cumlog[j] = cumlog[j - 1] + log_r[j - 1]
    log_denom = logsumexp(cumlog, axis=0)
    log_fracs = cumlog - log_denom[None, :]
    fracs = np.exp(log_fracs)
    return SahaResult(fracs, log_fracs, np.zeros(nd), n_stages, fracs[-1].copy())

@dataclass
class ElementIonization:
    elem_label: str
    stages: List[int]
    n_elem: np.ndarray
    fracs: np.ndarray
    n_stages: int
    omitted_next_frac: np.ndarray
    omitted_next_n: np.ndarray
    saha_result: SahaResult
    stage_to_idx: Dict[int, int] = field(default_factory=dict)

def group_species_by_element(species):
    by_elem = defaultdict(dict)
    for sp in species:
        try: label, _, ion_stage = species_info(sp.name)
        except KeyError: continue
        by_elem[label][ion_stage] = sp
    return by_elem

def compute_element_ionization(elem_label, stages_map, atom_data, hydro, chi, gratio, gion, abundances, t_override=None, w_dilute=None, t_floor_min=0.3, chi_to_next=None, w_apply_mode='first', floor_ne=FLOOR_NE):
    nd = hydro.nd
    t_local = t_override if t_override is not None else hydro.temperature
    chi_list = chi.get(elem_label, [])
    gratio_list = gratio.get(elem_label, [1.0] * len(chi_list))
    
    ref_chi = max(c[0] for c in chi_to_next[elem_label].values()) if chi_to_next and elem_label in chi_to_next else (chi_list[-1] if chi_list else 0.0)
    edge_equiv = ref_chi * SPEED_LIGHT / 1e15
    t_floor = t_floor_min * (1.0 + max(0.0, edge_equiv - 3.5) / 20.0)
    t_K = np.maximum(t_local, t_floor) * 1.0e4
    ne = np.maximum(hydro.ed, floor_ne)
    
    res = solve_saha_chain(chi_list, gratio_list, t_K, ne, w_dilute, w_apply_mode)
    A = ELEMENT_INFO.get(elem_label, (0, 1.0))[1]
    X = abundances.get(elem_label, np.zeros(nd))
    n_elem = hydro.density * X / (A * M_H)
    res.n_elem = n_elem
    
    stages_sorted = sorted(stages_map.keys())
    n_modelled = len(stages_sorted)
    n_chain = res.fracs.shape[0]
    fracs = np.zeros((n_modelled, nd))
    stage_to_idx = {}
    for i, stage in enumerate(stages_sorted):
        idx = min(stage - 1, n_chain - 1)
        fracs[i] = res.fracs[idx]
        stage_to_idx[stage] = i
        
    highest_modelled = max(stages_sorted) if stages_sorted else 1
    if n_chain > n_modelled:
        omitted_frac = res.fracs[-1].copy()
        omitted_n = n_elem * omitted_frac
    else:
        ctn = (chi_to_next or {}).get(elem_label, {}).get(highest_modelled)
        if ctn is not None:
            chi_next, gratio_next = ctn
            log_r = _saha_log_ratios([chi_next], [gratio_next], t_K, ne, w_dilute, w_apply_mode='first')
            log_1p = np.log1p(np.exp(-np.clip(log_r[0], -700, 700)))
            omitted_frac = np.exp(np.clip(-log_1p, -700.0, 0.0))
            omitted_n = n_elem * omitted_frac
        else:
            omitted_frac = fracs[-1].copy() if n_modelled else np.zeros(nd)
            omitted_n = np.zeros(nd)

    return ElementIonization(elem_label, stages_sorted, n_elem, fracs, n_modelled, omitted_frac, omitted_n, res, stage_to_idx)

def compute_ne_from_ionization(element_solutions):
    if not element_solutions: return np.zeros(0)
    nd = next(iter(element_solutions.values())).n_elem.size
    ne = np.zeros(nd)
    for sol in element_solutions.values():
        chain_fracs = sol.saha_result.fracs
        for s in range(chain_fracs.shape[0]):
            ne += sol.n_elem * chain_fracs[s] * s
    return ne

def compute_all_element_ionizations(species, atom_data, hydro, chi, gratio, gion, abundances, t_override=None, w_dilute=None, t_floor_min=0.3, chi_to_next=None, w_apply_mode='first', floor_ne=FLOOR_NE):
    by_elem = group_species_by_element(species)
    solutions = {}
    for label, stages_map in by_elem.items():
        if not any(atom_data.get(sp.name) is not None for sp in stages_map.values()): continue
        solutions[label] = compute_element_ionization(label, stages_map, atom_data, hydro, chi, gratio, gion, abundances, t_override, w_dilute, t_floor_min, chi_to_next, w_apply_mode, floor_ne)
    return solutions

# ===========================================================================
# Excitation / LTEPOP / DC
# ===========================================================================
def ltepop_wld_v2(ad, nlev, nd, t_1e4, ed, di_next, gion_next, w_occ=None):
    n = min(nlev, ad.n_levels)
    g = ad.stat_wt[:n]
    edge = ad.fedge[:n]
    T = np.maximum(t_1e4, FLOOR_T_1E4K)
    if w_occ is None: w_occ = np.ones((n, nd))
    
    log_rgu = math.log(LTEPOP_C) + np.log(np.maximum(ed, FLOOR_NE)) + np.log(np.maximum(di_next, FLOOR_DI)) - 1.5 * np.log(T) - math.log(max(gion_next, FLOOR_GION))
    log_pop = np.log(w_occ) + np.log(g)[:, None] + edge[:, None] * HDKT / T[None, :] + log_rgu[None, :]
    log_pop = np.clip(log_pop, -LOG_CLIP, LOG_CLIP)
    pop = np.exp(log_pop)
    
    out = np.zeros((nlev, nd))
    out[:n] = pop
    return out

def dc_from_target(target_pop, lte_pop, dc_min=FLOOR_DC, dc_max=1e300):
    safe_lte = np.maximum(lte_pop, 1e-300)
    uncapped = target_pop / safe_lte
    capped = (uncapped > dc_max) | (uncapped < dc_min) | ~np.isfinite(uncapped)
    dc = np.clip(uncapped, dc_min, dc_max)
    dc = np.where(np.isfinite(uncapped), dc, dc_min)
    return dc, capped, uncapped

# ===========================================================================
# Optical depth / grey radiation field
# ===========================================================================
def compute_tau_enhanced(hydro, kappa=None, conv=None):
    nd = hydro.nd
    k = kappa if kappa is not None else hydro.kappa
    k = np.where(np.isfinite(k), k, 0.0)
    
    kappa_fallback = conv.kappa_fallback if conv else "thomson"
    kappa_bf_scale = conv.kappa_bf_scale if conv else 1.0e-2
    tau_mode = conv.tau_mode if conv else "plane"
    v_th = conv.velocity_vth_cm_s if conv else 1.0e6
    dvdr_floor = conv.velocity_dvdr_floor if conv else 1.0e-12
    trapping_factor = conv.velocity_trapping_factor if conv else 5.0
    min_factor = conv.velocity_min_factor if conv else 0.05

    if np.all(k == 0.0):
        xh = hydro.mass_fracs.get('HYD', np.full(nd, 0.7))
        k_es = 0.2 * (1.0 + xh)
        if kappa_fallback == "kramers":
            t_safe = np.maximum(hydro.temperature, 0.05)
            k = k_es + kappa_bf_scale * t_safe ** (-3.5)
        else:
            k = k_es

    r = hydro.radius * 1.0e10
    tau_cont = np.zeros(nd)
    for i in range(1, nd):
        dr = abs(r[i] - r[i - 1])
        tau_cont[i] = tau_cont[i - 1] + 0.5 * (k[i] + k[i - 1]) * 0.5 * (hydro.density[i] + hydro.density[i - 1]) * dr

    tau_eff = tau_cont.copy()
    vel_diag = None
    diag = {"tau_mode": tau_mode, "kappa_fallback": kappa_fallback}

    if tau_mode == "plane": return tau_eff, tau_cont, vel_diag, diag

    vel_diag = compute_velocity_diagnostics(hydro, dvdr_floor, v_th)
    tau_eff = np.zeros(nd)
    tau_eff[0] = tau_cont[0]
    correction_factor = np.ones(nd)
    trapping_applied = np.zeros(nd, dtype=bool)
    v_cm_s = hydro.velocity * 1.0e5

    for i in range(1, nd):
        dtau_cont = max(tau_cont[i] - tau_cont[i - 1], 0.0)
        dr = max(abs(r[i] - r[i - 1]), 1.0)
        dv_abs = abs(v_cm_s[i] - v_cm_s[i - 1])
        dvdr = dv_abs / dr
        factor = 1.0
        
        if dvdr > dvdr_floor:
            L_sob = v_th / dvdr
            factor = min(1.0, L_sob / dr)
            factor = max(factor, min_factor)
        else:
            if tau_mode == "nonmonotonic":
                factor = trapping_factor
                trapping_applied[i] = True

        if tau_mode == "nonmonotonic":
            if vel_diag.extrema_mask[i] or vel_diag.extrema_mask[i - 1] or vel_diag.small_gradient_mask[i] or vel_diag.small_gradient_mask[i - 1]:
                factor = max(factor, trapping_factor)
                trapping_applied[i] = True

        correction_factor[i] = factor
        tau_eff[i] = tau_eff[i - 1] + factor * dtau_cont

    tau_eff = np.maximum.accumulate(tau_eff)
    diag.update({"velocity_tau_mode": tau_mode, "n_velocity_extrema": vel_diag.n_extrema, "n_trapping_zones": int(np.sum(trapping_applied))})
    return tau_eff, tau_cont, vel_diag, diag

def find_photosphere(hydro, tau, tau_target=2.0/3.0):
    nd = hydro.nd
    for i in range(nd - 1):
        if tau[i] >= tau_target:
            if i == 0: return float(hydro.radius[0]), 0
            t0, t1 = tau[i - 1], tau[i]
            if t1 == t0: return float(hydro.radius[i]), i
            frac = (tau_target - t0) / (t1 - t0)
            return float(hydro.radius[i - 1] + frac * (hydro.radius[i] - hydro.radius[i - 1])), i
    return float(hydro.radius[nd - 1]), nd - 1

def grey_w_eff(radius, r_ph, tau, tau_trans=1.0):
    w = np.empty_like(radius)
    for i in range(radius.size):
        if radius[i] <= r_ph:
            w_geom = 1.0
        else:
            ratio = r_ph / radius[i]
            w_geom = 0.5 * (1.0 - math.sqrt(max(0.0, 1.0 - ratio * ratio)))
        w[i] = min(w_geom + (1.0 - w_geom) * (1.0 - math.exp(-tau[i] / tau_trans)), 1.0)
    return w

def grey_t_rad(hydro, tau, r_ph, tau_trans=1.0):
    nd = hydro.nd
    t_local = hydro.temperature
    t_ph = t_local[0]
    for i in range(nd - 1):
        if tau[i] >= 2.0 / 3.0:
            if i == 0: t_ph = t_local[0]
            else:
                t0, t1 = tau[i - 1], tau[i]
                if t1 == t0: t_ph = t_local[i]
                else:
                    frac = (2.0 / 3.0 - t0) / (t1 - t0)
                    t_ph = t_local[i - 1] + frac * (t_local[i] - t_local[i - 1])
            break
    weight = np.exp(-tau / tau_trans)
    return np.maximum(t_local + (t_ph - t_local) * weight, 0.01)

def make_t_ion_from_radiation(t_gas, t_rad, tau, tau_trans=1.0, blend_mode="linear", t_floor_1e4k=1.0e-3, t_ceiling_1e4k=None):
    weight = np.exp(-np.maximum(tau, 0.0) / max(tau_trans, 1e-30))
    t_ion = t_gas + (t_rad - t_gas) * weight
    if blend_mode == "gas_floor": t_ion = np.maximum(t_ion, t_gas)
    t_ion = np.maximum(t_ion, t_floor_1e4k)
    if t_ceiling_1e4k is not None and t_ceiling_1e4k > 0: t_ion = np.minimum(t_ion, t_ceiling_1e4k)
    return t_ion

# ===========================================================================
# Smoothing
# ===========================================================================
def _median_filter(vals):
    n = vals.size
    if n < 3: return vals.copy()
    out = vals.copy()
    for i in range(n):
        a, b, c = vals[max(0, i - 1)], vals[i], vals[min(n - 1, i + 1)]
        out[i] = sorted((a, b, c))[1]
    return out

def _median_filter_front_aware(vals, log_slope_threshold=2.0, break_mask=None):
    n = vals.size
    if n < 3: return vals.copy()
    if break_mask is None: break_mask = np.zeros(n, dtype=bool)
    log_v = np.log10(np.maximum(vals, 1e-300))
    out = vals.copy()
    for i in range(1, n - 1):
        if break_mask[i - 1] or break_mask[i] or break_mask[i + 1]:
            out[i] = vals[i]
            continue
        grad = max(abs(log_v[i] - log_v[i - 1]), abs(log_v[i + 1] - log_v[i]))
        if grad < log_slope_threshold:
            a, b, c = vals[i - 1], vals[i], vals[i + 1]
            out[i] = sorted((a, b, c))[1]
        else:
            out[i] = vals[i]
    return out

# ===========================================================================
# Population builder
# ===========================================================================
@dataclass
class SpeciesPopulations:
    name: str
    elem_label: str
    ion_stage: int
    nlev: int
    di: np.ndarray
    dc: np.ndarray
    t_write: np.ndarray
    target_pop: np.ndarray
    lte_pop: np.ndarray
    ion_frac_current: np.ndarray
    ion_frac_next: np.ndarray
    omitted_stage_frac: np.ndarray
    floored_mask: np.ndarray
    capped_mask: np.ndarray
    abundance_residual: np.ndarray
    smooth_diag: Dict = None

def _apply_floor(physical_frac, min_ion_frac):
    mask = physical_frac < min_ion_frac
    return np.where(mask, min_ion_frac, physical_frac), mask

def build_species(sp, ad, hydro, chi, gratio, gion, abundances, conv, rad=None, chi_to_next=None, smooth_di=False, elem_ion=None, physics_break_mask=None):
    label, mass, ion_stage = species_info(sp.name)
    nd = hydro.nd
    nlev = min(sp.n_f, ad.n_levels) if conv.write_level_count == "NF" else min(sp.n_s, sp.n_f)
    nlev = max(nlev, 1)
    
    t_gas = hydro.temperature
    t_ion = rad['t_ion'] if (conv.strategy == "GREY" and rad) else None
    w_dil = rad['w'] if (conv.strategy == "GREY" and rad) else None

    if elem_ion is not None:
        res_next = elem_ion.saha_result
        cur_idx_elem = elem_ion.stage_to_idx.get(ion_stage, 0)
        frac_cur_raw = elem_ion.fracs[cur_idx_elem].copy()
        next_stage = ion_stage + 1
        frac_next_raw = elem_ion.fracs[elem_ion.stage_to_idx[next_stage]].copy() if next_stage in elem_ion.stage_to_idx else elem_ion.omitted_next_frac.copy()
        di_next_raw = elem_ion.n_elem * frac_next_raw
    else:
        raise ValueError("build_species requires an element-level ionisation solution.")

    A = ELEMENT_INFO.get(label, (0, 1.0))[1]
    X = abundances.get(label, np.zeros(nd))
    n_elem = hydro.density * X / (A * M_H)
    
    smooth_source = elem_ion.fracs
    n_smooth = smooth_source.shape[0]
    cur_sidx = elem_ion.stage_to_idx.get(ion_stage, 0)
    next_sidx = elem_ion.stage_to_idx.get(ion_stage + 1, -1)

    if smooth_di and n_smooth > 1 and nd > 2:
        all_fracs_raw = smooth_source.copy()
        all_fracs_smoothed = np.zeros_like(all_fracs_raw)
        use_front = conv.smooth_mode == "front_aware" or (physics_break_mask is not None and np.any(physics_break_mask))
        for s in range(n_smooth):
            all_fracs_smoothed[s] = _median_filter_front_aware(all_fracs_raw[s], conv.front_log_slope_threshold, physics_break_mask) if use_front else _median_filter(all_fracs_raw[s])
        col_sum = np.maximum(np.sum(all_fracs_smoothed, axis=0), FLOOR_FRAC)
        all_fracs_smoothed /= col_sum[None, :]
        frac_cur = all_fracs_smoothed[cur_sidx]
        frac_next = all_fracs_smoothed[next_sidx] if next_sidx >= 0 else frac_next_raw
        smooth_diag = {"mode": conv.smooth_mode, "use_front_aware": bool(use_front)}
    else:
        frac_cur, frac_next = frac_cur_raw, frac_next_raw
        smooth_diag = {"mode": conv.smooth_mode, "use_front_aware": False}

    di_frac_floored, floored = _apply_floor(frac_next, conv.min_ion_frac)
    if np.any(floored): frac_next = di_frac_floored
    di = n_elem * di_frac_floored

    t_for_lte = t_ion if (conv.strategy == "GREY" and conv.write_t_ion_to_in and rad) else t_gas
    gion_next = gion.get(label, {}).get(ion_stage, float(sp.gion) if sp.gion > 0 else FLOOR_GION)

    if conv.use_w_in_ltepop and rad and conv.strategy == "GREY":
        w_occ = np.broadcast_to(rad['w'][None, :], (nlev, nd)).copy()
        lte_pop = ltepop_wld_v2(ad, nlev, nd, t_for_lte, hydro.ed, di, gion_next, w_occ)
    else:
        lte_pop = ltepop_wld_v2(ad, nlev, nd, t_for_lte, hydro.ed, di, gion_next)

    if conv.strategy == "GREY" and rad:
        if conv.t_exc_source == "t_ion": t_exc = rad['t_ion']
        elif conv.t_exc_source == "rad": t_exc = rad.get('t_rad', t_gas)
        elif conv.t_exc_source == "hybrid":
            tau_safe = np.maximum(rad.get('tau', np.zeros(nd)), 0.0)
            tau_scale = max(conv.tau_trans * conv.t_hybrid_tau_scale, 1e-30)
            weight = np.exp(-tau_safe / tau_scale)
            t_exc = t_gas + (rad.get('t_rad', t_gas) - t_gas) * weight
        else: t_exc = t_gas
    else: t_exc = t_gas
    t_exc = np.maximum(t_exc, conv.t_exc_min)

    n_current = n_elem * frac_cur
    dE = (ad.energy[:nlev] - ad.energy[0]) * H_PLANCK * SPEED_LIGHT
    arg = -dE[:, None] / (K_B * t_exc[None, :] * 1.0e4)
    arg = np.clip(arg, -LOG_CLIP, LOG_CLIP)
    Z = np.maximum(np.sum(ad.stat_wt[:nlev, None] * np.exp(arg), axis=0), ad.stat_wt[0] if ad.stat_wt.size else 1.0)
    boltzmann_pop = n_current[None, :] * ad.stat_wt[:nlev, None] * np.exp(arg) / Z[None, :]

    if conv.strategy == "LTE":
        target_pop = lte_pop.copy()
        dc = np.ones((nlev, nd))
        capped = np.zeros((nlev, nd), dtype=bool)
    else:
        target_pop = boltzmann_pop
        dc, capped, _ = dc_from_target(target_pop, lte_pop)

    if conv.write_level_count == "NF" and ad.f_to_s is not None and nlev == sp.n_f:
        dc = _replicate_via_f_to_s(dc, ad.f_to_s, sp.n_s)

    t_write = t_ion if (conv.strategy == "GREY" and conv.write_t_ion_to_in and rad) else t_gas
    sum_target = np.sum(target_pop, axis=0)
    abund_res = np.log10(np.maximum(sum_target, FLOOR_POP) / np.maximum(n_elem, FLOOR_POP))

    return SpeciesPopulations(sp.name, label, ion_stage, nlev, di, dc, t_write, target_pop, lte_pop, frac_cur, frac_next, res_next.omitted_stage_frac, floored, capped, abund_res, smooth_diag)

def _replicate_via_f_to_s(dc, f_to_s, n_s):
    nlev, nd = dc.shape
    out = dc.copy()
    for i in range(nd):
        ta = np.zeros(n_s + 1)
        for j in range(nlev):
            s = max(int(f_to_s[j]) if j < len(f_to_s) else 1, 1)
            if ta[s] == 0.0: ta[s] = dc[j, i]
        for j in range(nlev):
            s = max(int(f_to_s[j]) if j < len(f_to_s) else 1, 1)
            out[j, i] = ta[s] if ta[s] != 0.0 else dc[j, i]
    return out

def enforce_closure(pop, hydro, abundances):
    A = ELEMENT_INFO.get(pop.elem_label, (0, 1.0))[1]
    X = abundances.get(pop.elem_label, np.zeros(hydro.nd))
    n_elem = hydro.density * X / (A * M_H)
    n_current = n_elem * np.maximum(pop.ion_frac_current, FLOOR_FRAC)
    s = np.sum(pop.target_pop, axis=0)
    factor = np.where(s > 0, n_current / np.maximum(s, 1e-300), 1.0)
    pop.target_pop = pop.target_pop * factor[None, :]
    return factor

# ===========================================================================
# CMFGEN _IN round-trip
# ===========================================================================
@dataclass
class InFileData:
    format_date: str; r_nd: float; rlum: float; nold: int; ndold: int
    radius: np.ndarray; di: np.ndarray; ed: np.ndarray; t: np.ndarray; dc: np.ndarray

def read_in_file(filepath: str) -> Optional[InFileData]:
    if not os.path.exists(filepath): return None
    with open(filepath) as f: lines = f.readlines()
    format_date = ""
    i = 0
    while i < min(len(lines), 10):
        if '!Format date' in lines[i]:
            format_date = lines[i].split('!')[0].strip()
            i += 1; break
        i += 1
    else: return None
    while i < len(lines) and not lines[i].strip(): i += 1
    if i >= len(lines): return None
    hdr = lines[i].split()
    if len(hdr) < 4: return None
    try: r_nd, rlum, nold, ndold = float(hdr[0]), float(hdr[1]), int(hdr[2]), int(hdr[3])
    except ValueError: return None
    i += 1
    radius, di, ed, t = np.zeros(ndold), np.zeros(ndold), np.zeros(ndold), np.zeros(ndold)
    dc = np.zeros((nold, ndold))
    for k in range(ndold):
        while i < len(lines) and not lines[i].strip(): i += 1
        if i >= len(lines): return None
        toks = lines[i].split()
        if len(toks) < 4: return None
        try: radius[k], di[k], ed[k], t[k] = float(toks[0]), float(toks[1]), float(toks[2]), float(toks[3])
        except ValueError: return None
        i += 1
        dc_read = 0
        while dc_read < nold and i < len(lines):
            while i < len(lines) and not lines[i].strip(): i += 1
            if i >= len(lines): break
            for tok in lines[i].split():
                try: dc[dc_read, k] = float(tok); dc_read += 1
                except ValueError: break
                if dc_read >= nold: break
            i += 1
    return InFileData(format_date, r_nd, rlum, nold, ndold, radius, di, ed, t, dc)

def regrid_log_dc(old_r, new_r, old_dc, old_di, old_t, old_ed):
    ndold, nd, nold = old_r.size, new_r.size, old_dc.shape[0]
    if abs(old_r[-1] / new_r[-1] - 1.0) > 1e-4: old_r = old_r * (new_r[-1] / old_r[-1])
    log_old_r, log_new_r = np.log(np.maximum(old_r, 1e-30)), np.log(np.maximum(new_r, 1e-30))
    x_old, x_new = -log_old_r, -log_new_r
    log_old_dc = np.log(np.maximum(np.abs(old_dc), 1e-300))
    log_new_dc = np.zeros((nold, nd))
    for j in range(nold): log_new_dc[j] = np.interp(x_new, x_old, log_old_dc[j], left=log_old_dc[j, 0], right=log_old_dc[j, -1])
    new_dc = np.exp(log_new_dc)
    log_old_di = np.log(np.maximum(old_di, 1e-300))
    new_di = np.exp(np.interp(x_new, x_old, log_old_di, left=log_old_di[0], right=log_old_di[-1]))
    new_t = np.interp(x_new, x_old, old_t, left=old_t[0], right=old_t[-1])
    new_ed = np.interp(x_new, x_old, old_ed, left=old_ed[0], right=old_ed[-1])
    return new_dc, new_di, new_t, new_ed

def build_continuation_species(sp, ad, hydro, conv, prev_in, chi_to_next, gion, abundances):
    label, mass, ion_stage = species_info(sp.name)
    nd = hydro.nd
    nlev = min(sp.n_f, ad.n_levels) if conv.write_level_count == "NF" else min(sp.n_s, sp.n_f)
    nlev = max(nlev, 1)
    new_dc, new_di, _, _ = regrid_log_dc(prev_in.radius, hydro.radius, prev_in.dc, prev_in.di, prev_in.t, prev_in.ed)
    nlev = min(nlev, prev_in.nold)
    dc, di = new_dc[:nlev, :], new_di
    A = ELEMENT_INFO.get(label, (0, 1.0))[1]
    X = abundances.get(label, np.zeros(nd))
    n_elem = hydro.density * X / (A * M_H)
    di_frac_floored, floored = _apply_floor(di / np.maximum(n_elem, 1e-300), conv.min_ion_frac)
    di = n_elem * di_frac_floored
    gion_next = gion.get(label, {}).get(ion_stage, float(sp.gion) if sp.gion > 0 else 1.0)
    lte_pop = ltepop_wld_v2(ad, nlev, nd, hydro.temperature, hydro.ed, di, gion_next)
    target_pop = dc * lte_pop
    chi_list = [chi_to_next[label][ion_stage][0]] if label in chi_to_next and ion_stage in chi_to_next[label] else []
    omitted = np.ones(nd)
    if chi_list:
        t_K = np.maximum(hydro.temperature, conv.t_exc_min) * 1e4
        ne = np.maximum(hydro.ed, 1.0)
        log_r = _saha_log_ratios(chi_list, [1.0], t_K, ne, None)
        omitted = 1.0 / (1.0 + np.exp(-np.clip(log_r[0], -700, 700)))
    return SpeciesPopulations(sp.name, label, ion_stage, nlev, dc, di, hydro.temperature.copy(), lte_pop, lte_pop, target_pop, np.ones(nd), omitted, omitted, floored, np.zeros((nlev, nd), dtype=bool), np.zeros(nd))

def write_in_file(filepath, pop, hydro, conv):
    nd = hydro.nd
    r_nd = float(hydro.radius[nd - 1])
    os.makedirs(os.path.dirname(filepath) or '.', exist_ok=True)
    tmp = filepath + ".tmp"
    with open(tmp, 'w') as f:
        f.write("\n")
        f.write(f" {datetime.now().strftime('%d-%b-%Y').upper()}                           !Format date\n")
        f.write("\n")
        f.write(f" {r_nd:14.8E}     {0.0:12.4E}    {pop.nlev:4d}    {nd:4d}\n")
        for i in range(nd):
            f.write("\n")
            f.write(f" {hydro.radius[i]:17.10E}{pop.di[i]:17.8E}{hydro.ed[i]:17.8E}{pop.t_write[i]:17.8E}{0.0:17.8E}{hydro.velocity[i]:17.8E}{1.0:17.8E}  {i+1:4d}\n")
            for j0 in range(0, pop.nlev, 5):
                line = " "
                for j in range(j0, min(j0 + 5, pop.nlev)): line += f"{pop.dc[j, i]:17.7E}"
                f.write(line + "\n")
    os.replace(tmp, filepath)

def emulate_regrid_ltepop(pop, ad, hydro, gion_next):
    return pop.dc * ltepop_wld_v2(ad, pop.nlev, hydro.nd, pop.t_write, hydro.ed, pop.di, gion_next)

def roundtrip_error(pop, ad, hydro, gion_next):
    recon = emulate_regrid_ltepop(pop, ad, hydro, gion_next)
    target = pop.target_pop
    significant = np.isfinite(recon) & np.isfinite(target) & (target > 1e-100) & (recon > 1e-300)
    safe_target, safe_recon = np.maximum(np.abs(target), 1e-300), np.maximum(np.abs(recon), 1e-300)
    log10_ratio = np.zeros_like(target)
    log10_ratio[significant] = np.abs(np.log10(safe_recon[significant] / safe_target[significant]))
    rel = np.abs(recon - target) / safe_target
    dc_abs = np.abs(pop.dc)
    dc_finite = np.isfinite(dc_abs)
    return {
        'max_log10_ratio': float(np.max(log10_ratio[significant])) if np.any(significant) else 0.0,
        'p95_log10_ratio': float(np.percentile(log10_ratio[significant], 95)) if np.any(significant) else 0.0,
        'max_rel_error': float(np.max(rel[significant])) if np.any(significant) else 0.0,
        'n_nonfinite': int(np.sum(~np.isfinite(recon))),
        'n_failed': int(np.sum(rel[significant] > 1e-6)) if np.any(significant) else 0,
        'dc_min': float(np.min(dc_abs[dc_finite])) if np.any(dc_finite) else 0.0,
        'dc_max': float(np.max(dc_abs[dc_finite])) if np.any(dc_finite) else 0.0,
        'n_dc_floored': int(np.sum(pop.floored_mask)),
        'n_dc_capped': int(np.sum(pop.capped_mask)),
    }

# ===========================================================================
# Validation & Reporting
# ===========================================================================
def validate_run(pops, hydro, abundances, conv, elem_solutions=None):
    nd = hydro.nd
    issues = []
    for p in pops:
        if not np.all(np.isfinite(p.di)): issues.append(f"{p.name}: DI has non-finite values")
        if np.any(p.di < 0): issues.append(f"{p.name}: DI has negative values")
        if not np.all(np.isfinite(p.dc)): issues.append(f"{p.name}: DC has non-finite values")
        if np.any(p.dc < 0): issues.append(f"{p.name}: DC has negative values")
    ne_recomputed = np.zeros(nd)
    if elem_solutions:
        for sol in elem_solutions.values():
            for s in range(sol.saha_result.fracs.shape[0]): ne_recomputed += sol.n_elem * sol.saha_result.fracs[s] * s
    else:
        for p in pops:
            A = ELEMENT_INFO.get(p.elem_label, (0, 1.0))[1]
            n_elem = hydro.density * abundances.get(p.elem_label, np.zeros(nd)) / (A * M_H)
            ne_recomputed += n_elem * (p.ion_frac_current * (p.ion_stage - 1) + p.ion_frac_next * p.ion_stage)
    ne_ratio = ne_recomputed / np.maximum(hydro.ed, 1.0)
    return {'issues': issues, 'n_issues': len(issues), 'ne_ratio_max': float(np.max(np.abs(ne_ratio - 1.0)))}

def write_reports(pops, hydro, conv, abundances, vadat_abundances, output_dir, validation, roundtrip, strategy_log):
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'proxy_validation.txt'), 'w') as f:
        f.write(f"Strategy: {conv.strategy}\n")
        f.write(f"Saha Inversion: {conv.saha_inversion_for_ne}\n")
        f.write(f"Non-Local Reversals: {conv.non_local_reversals}\n")
        f.write(f"Validation issues: {validation['n_issues']}\n")
        for iss in validation['issues']: f.write(f"  - {iss}\n")
        f.write(f"Ne ratio max deviation: {validation['ne_ratio_max']:.3e}\n")

# ===========================================================================
# ADVANCED DIAGNOSTIC PLOTTING
# ===========================================================================
def _import_matplotlib():
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        return plt, True
    except ImportError:
        return None, False

def plot_velocity_and_reversals(hydro: HydroData, vel_diag: Optional[VelocityFieldDiagnostics], output_dir: str, plt) -> None:
    """Plot 1: Velocity field with shaded non-monotonic reversal zones."""
    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    r_cm = hydro.radius * 1e10
    v_km_s = hydro.velocity
    
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(r_cm, v_km_s, 'k-', linewidth=1.5, label='Velocity')
    ax1.set_xlabel('Radius (cm)', fontsize=12)
    ax1.set_ylabel('Velocity (km/s)', fontsize=12, color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.set_xscale('log')
    ax1.grid(True, alpha=0.3)
    
    if vel_diag is not None and np.any(vel_diag.extrema_mask | (vel_diag.dvds_signed < 0)):
        reversal_mask = vel_diag.extrema_mask | (vel_diag.dvds_signed < 0)
        ax1.fill_between(r_cm, 0, np.max(v_km_s)*1.1, where=reversal_mask, 
                         color='red', alpha=0.2, label='Non-monotonic / Reversal Zones')
        
        ax2 = ax1.twinx()
        ax2.plot(r_cm, vel_diag.dvds_signed, 'r--', alpha=0.6, linewidth=1, label='dv/dr (signed)')
        ax2.axhline(0, color='r', linestyle=':', alpha=0.5)
        ax2.set_ylabel('dv/dr (cm/s/cm)', fontsize=12, color='r')
        ax2.tick_params(axis='y', labelcolor='r')
        
    ax1.set_title('Velocity Field & Non-Monotonic Zones', fontsize=14)
    lines, labels = ax1.get_legend_handles_labels()
    if vel_diag is not None and np.any(vel_diag.extrema_mask | (vel_diag.dvds_signed < 0)):
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines + lines2, labels + labels2, loc='best')
    else:
        ax1.legend(loc='best')
        
    plt.tight_layout()
    path = os.path.join(plots_dir, 'POP_velocity_and_reversals.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)

def plot_temperature_decoupling(hydro: HydroData, rad: Optional[Dict], tau_cont: Optional[np.ndarray], output_dir: str, plt) -> None:
    """Plot 2: Temperature stratification showing gas/rad/ion decoupling."""
    plots_dir = os.path.join(output_dir, 'plots')
    r_cm = hydro.radius * 1e10
    tau = rad['tau'] if rad and 'tau' in rad else np.ones_like(r_cm)
    tau_plot = np.maximum(tau, 1e-10)  # floor for log-scale

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(tau_plot, hydro.temperature, 'k-', linewidth=1.5, label=r'$T_{gas}$')
    if rad:
        if 't_rad' in rad and rad['t_rad'] is not None:
            ax.plot(tau_plot, rad['t_rad'], 'r--', linewidth=1.5, label=r'$T_{rad}$')
        if 't_ion' in rad and rad['t_ion'] is not None:
            ax.plot(tau_plot, rad['t_ion'], 'b-.', linewidth=1.5, label=r'$T_{eff,ion}$')

    ax.set_xscale('log')
    ax.set_xlabel('Rosseland Optical Depth ($\tau$)', fontsize=12)
    ax.set_ylabel('Temperature ($10^4$ K)', fontsize=12)
    ax.set_title('Temperature Stratification & Decoupling', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(tau_plot.max(), tau_plot.min())  # Reverse: inner (high tau) on right
    
    plt.tight_layout()
    path = os.path.join(plots_dir, 'POP_temperature_decoupling.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)

def plot_electron_density_consistency(hydro: HydroData, elem_solutions: Dict[str, ElementIonization], output_dir: str, plt) -> None:
    """Plot 3: Verifies Saha inversion matched the STELLA n_e target."""
    plots_dir = os.path.join(output_dir, 'plots')
    r_cm = hydro.radius * 1e10
    ne_target = hydro.ed
    
    ne_saha = compute_ne_from_ionization(elem_solutions)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    ax1.plot(r_cm, ne_target, 'k-', linewidth=1.5, label='Target $n_e$ (STELLA/Hydro)')
    ax1.plot(r_cm, ne_saha, 'r--', linewidth=1.5, label='Saha Recomputed $n_e$')
    ax1.set_ylabel('Electron Density (cm$^{-3}$)', fontsize=12)
    ax1.set_title('Electron Density Consistency', fontsize=14)
    ax1.set_yscale('log')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3, which='both')
    
    ratio = ne_saha / np.maximum(ne_target, 1.0)
    ax2.plot(r_cm, ratio, 'b-', linewidth=1.5)
    ax2.axhline(1.0, color='k', linestyle=':', alpha=0.5)
    ax2.set_xlabel('Radius (cm)', fontsize=12)
    ax2.set_ylabel('Ratio (Saha / Target)', fontsize=12)
    ax2.set_xscale('log')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    path = os.path.join(plots_dir, 'POP_electron_density_consistency.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)

def plot_ionization_fractions(hydro: HydroData, elem_solutions: Dict[str, ElementIonization], output_dir: str, plt) -> None:
    """Plot 4: Ionization fractions vs velocity (top) and radius (bottom) for each element."""
    plots_dir = os.path.join(output_dir, 'plots')
    v_km_s = hydro.velocity
    r_cm = hydro.radius * 1e10

    # Element display names
    ELEM_NAMES = {'HYD': 'H', 'HE': 'He', 'OXY': 'O', 'SIL': 'Si', 'CAL': 'Ca', 'IRON': 'Fe',
                  'CARB': 'C', 'NIT': 'N', 'NEON': 'Ne', 'MG': 'Mg', 'S': 'S', 'NI': 'Ni'}

    # Distinct line styles to differentiate overlapping stages
    LINE_STYLES = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]

    for elem_label, sol in elem_solutions.items():
        fracs = np.asarray(sol.fracs)
        if fracs.ndim > 2:
            fracs = fracs.reshape(-1, fracs.shape[-1])
        stages = sol.stages
        n_stages = len(stages)
        elem_name = ELEM_NAMES.get(elem_label, elem_label)

        # Sort by velocity for the top subplot
        sort_idx = np.argsort(v_km_s)
        v_sorted = v_km_s[sort_idx]

        # Y-axis: span from the smallest positive fraction to 1
        # Use all finite positive values to set the lower bound
        positive_fracs = fracs[np.isfinite(fracs) & (fracs > 0)]
        if positive_fracs.size == 0:
            y_lo = 1e-30
        else:
            y_lo = max(positive_fracs.min() * 0.1, 1e-300)
        y_hi = 2.0

        colors = plt.cm.tab10(np.linspace(0, 0.9, max(n_stages, 2)))

        fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(12, 10), sharex=False)

        # --- Top: Ionization vs velocity ---
        for i, stage in enumerate(stages):
            label = f"{elem_name} {ROMAN.get(stage, stage)}"
            ls = LINE_STYLES[i % len(LINE_STYLES)]
            ax_top.plot(v_sorted, fracs[i][sort_idx], color=colors[i],
                        linewidth=1.8, linestyle=ls, label=label)

        ax_top.set_xlabel('Velocity (km/s)', fontsize=12)
        ax_top.set_ylabel('Ion Fraction', fontsize=12)
        ax_top.set_title(f'{elem_name} Ionization Balance vs Velocity', fontsize=14)
        ax_top.set_yscale('log')
        ax_top.set_ylim(y_lo, y_hi)
        ax_top.legend(fontsize=9, ncol=min(n_stages, 3), loc='best')
        ax_top.grid(True, alpha=0.3, which='both')

        # --- Bottom: Ionization vs radius ---
        for i, stage in enumerate(stages):
            label = f"{elem_name} {ROMAN.get(stage, stage)}"
            ls = LINE_STYLES[i % len(LINE_STYLES)]
            ax_bot.plot(r_cm, fracs[i], color=colors[i],
                        linewidth=1.8, linestyle=ls, label=label)

        ax_bot.set_xlabel('Radius (cm)', fontsize=12)
        ax_bot.set_ylabel('Ion Fraction', fontsize=12)
        ax_bot.set_title(f'{elem_name} Ionization Balance vs Radius', fontsize=14)
        ax_bot.set_xscale('log')
        ax_bot.set_yscale('log')
        ax_bot.set_ylim(y_lo, y_hi)
        ax_bot.legend(fontsize=9, ncol=min(n_stages, 3), loc='best')
        ax_bot.grid(True, alpha=0.3, which='both')

        plt.tight_layout()
        path = os.path.join(plots_dir, f'POP_ionization_{elem_label}.png')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)

def plot_dc_heatmap(pops: List[SpeciesPopulations], hydro: HydroData, output_dir: str, plt) -> None:
    """Plot 5: Heatmap of Departure Coefficients highlighting capped corrections."""
    plots_dir = os.path.join(output_dir, 'plots')
    
    max_caps = 0
    target_pop = None
    for p in pops:
        n_caps = np.sum(p.capped_mask)
        if n_caps > max_caps:
            max_caps = n_caps
            target_pop = p
            
    if target_pop is None or max_caps == 0:
        for p in pops:
            if p.nlev > 1:
                target_pop = p
                break
                
    if target_pop is None: return
        
    dc = target_pop.dc
    log_dc = np.log10(np.maximum(np.abs(dc), 1e-300))
    capped = target_pop.capped_mask
    
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(log_dc, aspect='auto', cmap='viridis', origin='lower', 
                   extent=[0, hydro.nd, 0, target_pop.nlev])
    
    if np.any(capped):
        cap_y, cap_x = np.where(capped)
        ax.scatter(cap_x, cap_y, color='red', marker='x', s=10, alpha=0.8, label='Capped (>100% corr)')
        
    ax.set_xlabel('Depth Index (Outer -> Inner)', fontsize=12)
    ax.set_ylabel('Level Index', fontsize=12)
    ax.set_title(f'Departure Coefficients ($\log_{{10}} b_i$) for {target_pop.name}', fontsize=14)
    fig.colorbar(im, ax=ax, label=r'$\log_{{10}}(DC)$')
    if np.any(capped):
        ax.legend(loc='upper right')
        
    plt.tight_layout()
    path = os.path.join(plots_dir, f'POP_dc_heatmap_{target_pop.name}.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)

def plot_diagnostics(
    pops: List[SpeciesPopulations],
    hydro: HydroData,
    elem_solutions: Dict[str, ElementIonization],
    rad: Optional[Dict],
    tau_cont: Optional[np.ndarray],
    vel_diag: Optional[VelocityFieldDiagnostics],
    output_dir: str,
    make_plots: bool = True,
) -> bool:
    if not make_plots: return False
    if not pops:
        print("WARNING: no populations to plot")
        return False
    plt, ok = _import_matplotlib()
    if not ok:
        print("WARNING: matplotlib not available, skipping plots")
        return False

    os.makedirs(os.path.join(output_dir, 'plots'), exist_ok=True)
    
    plot_velocity_and_reversals(hydro, vel_diag, output_dir, plt)
    plot_temperature_decoupling(hydro, rad, tau_cont, output_dir, plt)
    plot_electron_density_consistency(hydro, elem_solutions, output_dir, plt)
    plot_ionization_fractions(hydro, elem_solutions, output_dir, plt)
    plot_dc_heatmap(pops, hydro, output_dir, plt)
    
    print(f"Diagnostic plots written to {os.path.join(output_dir, 'plots')}/")
    return True

# ===========================================================================
# Top-level generation
# ===========================================================================
def generate(inputs, conv, output_dir, prev_model_dir=None, make_plots=False):
    conv.validate()
    os.makedirs(output_dir, exist_ok=True)
    
    # Abundances
    if inputs.vadat_abundances:
        abundances = {e: np.full(inputs.hydro.nd, f) for e, f in inputs.vadat_abundances.items()}
        for e, fracs in inputs.hydro.mass_fracs.items():
            if e not in abundances and np.any(fracs > 0): abundances[e] = fracs
    else:
        abundances = {e: fracs for e, fracs in inputs.hydro.mass_fracs.items()}

    strategy_log = {"requested": conv.strategy, "steps": []}
    chi, gratio, gion, chi_to_next = build_saha_chains(inputs.species, inputs.atom_data)

    if inputs.gammas_in_ed is not None and inputs.gammas_in_ed.size == inputs.hydro.nd:
        inputs.hydro.ed = inputs.gammas_in_ed.copy()

    rad = None
    vel_diag, temp_diag, tau_cont = None, None, None

    if conv.strategy == "GREY":
        tau_eff, tau_cont, vel_diag, tau_diag = compute_tau_enhanced(inputs.hydro, conv=conv)
        r_ph, _ = find_photosphere(inputs.hydro, tau_eff)
        w = grey_w_eff(inputs.hydro.radius, r_ph, tau_eff, conv.tau_trans)
        
        t_rad = grey_t_rad(inputs.hydro, tau_eff, r_ph, conv.tau_trans)
        if conv.use_t_ion or conv.t_ion_source in ("stella", "stella-direct", "grey"):
            t_ion = make_t_ion_from_radiation(inputs.hydro.temperature, t_rad, tau_eff, conv.tau_trans, "gas_floor", 1e-3)
        else:
            t_ion = inputs.hydro.temperature.copy()
            
        rad = {"tau": tau_eff, "tau_cont": tau_cont, "r_ph": r_ph, "w": w, "t_rad": t_rad, "t_ion": t_ion}

        # NEW OPTIMIZATION: Non-Local Radiation in Reversals
        if conv.non_local_reversals and vel_diag is not None:
            print("  Applying non-local radiation averaging in velocity reversals...")
            t_rad_nl, w_nl = apply_non_local_radiation_in_reversals(t_rad, w, vel_diag)
            t_rad, w = t_rad_nl, w_nl
            if conv.use_t_ion or conv.t_ion_source in ("stella", "stella-direct", "grey"):
                t_ion = make_t_ion_from_radiation(inputs.hydro.temperature, t_rad, tau_eff, conv.tau_trans, "gas_floor", 1e-3)
            rad = {"tau": tau_eff, "tau_cont": tau_cont, "r_ph": r_ph, "w": w, "t_rad": t_rad, "t_ion": t_ion}
            strategy_log["non_local_reversals"] = True

        # NEW OPTIMIZATION: Saha Inversion for n_e
        if conv.saha_inversion_for_ne:
            print("  Computing T_eff for Saha inversion to match STELLA n_e...")
            t_eff_ion = compute_t_eff_for_ne_target(inputs.hydro.ed, abundances, inputs.hydro, chi, gratio)
            rad['t_ion'] = t_eff_ion  # Override T_ion for the ionization solver
            strategy_log["saha_inversion"] = {"t_eff_min": float(np.min(t_eff_ion)), "t_eff_max": float(np.max(t_eff_ion))}

        ne_final, ne_log = inputs.hydro.ed, {"converged": True, "n_iter": 0}
        inputs.hydro.ed = ne_final
        strategy_log["ne_iteration"] = ne_log
    else:
        if conv.velocity_aware_smoothing or conv.tau_mode != "plane":
            vel_diag = compute_velocity_diagnostics(inputs.hydro, conv.velocity_dvdr_floor, conv.velocity_vth_cm_s)

    t_ion_for_elem = rad['t_ion'] if (conv.strategy == "GREY" and rad) else None
    w_dil_for_elem = rad['w'] if (conv.strategy == "GREY" and rad) else None

    elem_solutions = compute_all_element_ionizations(
        inputs.species, inputs.atom_data, inputs.hydro, chi, gratio, gion, abundances,
        t_override=t_ion_for_elem, w_dilute=w_dil_for_elem, t_floor_min=conv.t_exc_min,
        chi_to_next=chi_to_next, w_apply_mode=conv.w_apply_mode, floor_ne=conv.floor_ne
    )

    pops = []
    roundtrip = {}
    physics_break_mask = np.zeros(inputs.hydro.nd, dtype=bool)
    if conv.velocity_aware_smoothing and vel_diag: physics_break_mask |= vel_diag.extrema_mask | vel_diag.small_gradient_mask

    # Compute self-consistent Ne from proxy ionization before writing _IN files
    if elem_solutions:
        ne_self_consistent = compute_ne_from_ionization(elem_solutions)
        if ne_self_consistent.size == inputs.hydro.nd:
            inputs.hydro.ed = ne_self_consistent.copy()
            gammas_path = os.path.join(output_dir, 'GAMMAS_IN')
            update_gammas_in_ne(gammas_path, ne_self_consistent, inputs.hydro.nd)

    for sp in inputs.species:
        ad = inputs.atom_data.get(sp.name)
        if ad is None: continue
        
        pop = build_species(sp, ad, inputs.hydro, chi, gratio, gion, abundances, conv, rad, chi_to_next, conv.smooth, elem_solutions.get(species_info(sp.name)[0]), physics_break_mask)
        if conv.strategy != "LTE":
            enforce_closure(pop, inputs.hydro, abundances)
            pop.dc, pop.capped_mask, _ = dc_from_target(pop.target_pop, pop.lte_pop)
            if conv.write_level_count == "NF" and ad.f_to_s is not None and pop.nlev == sp.n_f:
                pop.dc = _replicate_via_f_to_s(pop.dc, ad.f_to_s, sp.n_s)
        
        gion_next = gion.get(pop.elem_label, {}).get(pop.ion_stage, float(sp.gion) if sp.gion > 0 else FLOOR_GION)
        rt = roundtrip_error(pop, ad, inputs.hydro, gion_next)
        roundtrip[sp.name] = rt
        pops.append(pop)
        write_in_file(os.path.join(output_dir, f"{sp.name}_IN"), pop, inputs.hydro, conv)

    validation = validate_run(pops, inputs.hydro, abundances, conv, elem_solutions)
    write_reports(pops, inputs.hydro, conv, abundances, inputs.vadat_abundances, output_dir, validation, roundtrip, strategy_log)
    
    if make_plots:
        plot_diagnostics(pops, inputs.hydro, elem_solutions, rad, tau_cont, vel_diag, output_dir, make_plots=True)
    
    return {"validation": validation, "roundtrip": roundtrip, "strategy_log": strategy_log, "n_species": len(pops)}

# ===========================================================================
# CLI
# ===========================================================================
def main(argv=None):
    p = argparse.ArgumentParser(description="Generate proxy CMFGEN _IN population files (v2 with non-monotonic optimizations & diagnostics).")
    p.add_argument("model_dir", help="CMFGEN model directory")
    p.add_argument("--hydro-data", help="Path to SN_HYDRO_DATA")
    p.add_argument("--output-dir", help="Output directory")
    p.add_argument("--strategy", choices=["LTE", "GREY", "CONTINUATION"], default="GREY")
    p.add_argument("--prev-model", help="Previous CMFGEN model dir (CONTINUATION)")
    p.add_argument("--write-t-ion", action="store_true")
    p.add_argument("--smooth", action="store_true")
    p.add_argument("--no-smooth", action="store_true")
    p.add_argument("--write-levels", choices=["NF", "NS"], default="NF")
    p.add_argument("--min-ion-frac", type=float, default=1e-100)
    p.add_argument("--t-exc-min", type=float, default=0.3)
    p.add_argument("--tau-trans", type=float, default=1.0)
    p.add_argument("--population-tolerance", type=float, default=1e-6)
    p.add_argument("--require-vadat-for-nlte", type=str, default="true", choices=["true", "false"])
    p.add_argument("--ne-max-iter", type=int, default=50)
    p.add_argument("--ne-tol", type=float, default=1e-3)
    p.add_argument("--use-w-in-ltepop", action="store_true")
    p.add_argument("--use-t-ion", action="store_true")
    p.add_argument("--w-apply-mode", choices=["first", "all", "sqrt", "none"], default="first")
    p.add_argument("--complex-model", action="store_true")
    p.add_argument("--kappa-fallback", choices=["thomson", "kramers"], default="thomson")
    p.add_argument("--kappa-bf-scale", type=float, default=1.0e-2)
    p.add_argument("--tau-mode", choices=["plane", "sobolev", "nonmonotonic"], default="nonmonotonic")
    p.add_argument("--velocity-vth", type=float, default=1.0e6)
    p.add_argument("--velocity-dvdr-floor", type=float, default=1.0e-12)
    p.add_argument("--velocity-trapping-factor", type=float, default=5.0)
    p.add_argument("--velocity-min-factor", type=float, default=0.05)
    p.add_argument("--velocity-aware-smoothing", action="store_true")
    p.add_argument("--temperature-aware-smoothing", action="store_true")
    p.add_argument("--temperature-log-slope-threshold", type=float, default=0.5)
    p.add_argument("--t-exc-source", choices=["gas", "t_ion", "rad", "hybrid"], default="hybrid")
    p.add_argument("--use-t-ion-for-levels", action="store_true")
    p.add_argument("--t-hybrid-tau-scale", type=float, default=1.0)
    p.add_argument("--smooth-mode", choices=["uniform", "front_aware"], default="front_aware")
    p.add_argument("--front-log-slope-threshold", type=float, default=2.0)
    p.add_argument("--ne-damping-mode", choices=["fixed", "adaptive"], default="fixed")
    p.add_argument("--ne-max-log-step", type=float, default=0.0)
    p.add_argument("--floor-ne", type=float, default=1.0)
    p.add_argument("--t-ion-source", choices=["gas", "grey", "stella", "stella-direct"], default="grey")
    p.add_argument("--stella-file", default=None)
    p.add_argument("--stella-radius-column", default=None)
    p.add_argument("--stella-velocity-column", default=None)
    p.add_argument("--stella-tau-column", default=None)
    p.add_argument("--stella-trad-column", default=None)
    p.add_argument("--stella-tion-column", default=None)
    p.add_argument("--stella-luminosity-column", default=None)
    p.add_argument("--stella-energy-density-column", default=None)
    p.add_argument("--stella-radius-unit", choices=["cm", "1e10cm", "rsun", "km"], default="cm")
    p.add_argument("--stella-temperature-unit", choices=["K", "1e4K", "eV"], default="K")
    p.add_argument("--stella-luminosity-unit", choices=["erg/s", "lsun"], default="erg/s")
    p.add_argument("--stella-interpolation-coordinate", choices=["radius", "velocity"], default="radius")
    p.add_argument("--stella-fallback", choices=["grey", "gas", "error"], default="grey")
    p.add_argument("--stella-t-ion-blend-mode", choices=["linear", "gas_floor"], default="linear")
    p.add_argument("--stella-t-ion-floor", type=float, default=1.0e-3)
    p.add_argument("--stella-t-ion-ceiling", type=float, default=0.0)
    p.add_argument("--mass-frac-floor", type=float, default=1.0e-8)
    p.add_argument("--force-overwrite-source", action="store_true")
    p.add_argument("--plots", action="store_true")
    
    # NEW OPTIMIZATION FLAGS
    p.add_argument("--saha-inversion-for-ne", action="store_true", help="Invert Saha equation to find T_eff that matches STELLA n_e")
    p.add_argument("--non-local-reversals", action="store_true", help="Replace T_rad and W in velocity reversals with nearest monotonic shells")

    args = p.parse_args(argv)
    if args.use_t_ion_for_levels and args.t_exc_source == p.get_default("t_exc_source"): args.t_exc_source = "t_ion"

    conv = CmfgenConventions(
        strategy=args.strategy, write_t_ion_to_in=args.write_t_ion, smooth=args.smooth,
        min_ion_frac=args.min_ion_frac, t_exc_min=args.t_exc_min, tau_trans=args.tau_trans,
        population_tolerance=args.population_tolerance, require_vadat_for_nlte=(args.require_vadat_for_nlte == "true"),
        ne_max_iter=args.ne_max_iter, ne_tol=args.ne_tol, use_w_in_ltepop=args.use_w_in_ltepop,
        use_t_ion=args.use_t_ion, w_apply_mode=args.w_apply_mode, write_level_count=args.write_levels,
        kappa_fallback=args.kappa_fallback, kappa_bf_scale=args.kappa_bf_scale, tau_mode=args.tau_mode,
        velocity_vth_cm_s=args.velocity_vth, velocity_dvdr_floor=args.velocity_dvdr_floor,
        velocity_trapping_factor=args.velocity_trapping_factor, velocity_min_factor=args.velocity_min_factor,
        velocity_aware_smoothing=args.velocity_aware_smoothing, temperature_aware_smoothing=args.temperature_aware_smoothing,
        temperature_log_slope_threshold=args.temperature_log_slope_threshold, t_exc_source=args.t_exc_source,
        t_hybrid_tau_scale=args.t_hybrid_tau_scale, smooth_mode=args.smooth_mode,
        front_log_slope_threshold=args.front_log_slope_threshold, ne_damping_mode=args.ne_damping_mode,
        ne_max_log_step=args.ne_max_log_step, floor_ne=args.floor_ne, t_ion_source=args.t_ion_source,
        stella_file=args.stella_file, stella_radius_column=args.stella_radius_column,
        stella_velocity_column=args.stella_velocity_column, stella_tau_column=args.stella_tau_column,
        stella_t_rad_column=args.stella_trad_column, stella_t_ion_column=args.stella_tion_column,
        stella_luminosity_column=args.stella_luminosity_column, stella_energy_density_column=args.stella_energy_density_column,
        stella_radius_unit=args.stella_radius_unit, stella_temperature_unit=args.stella_temperature_unit,
        stella_luminosity_unit=args.stella_luminosity_unit, stella_interpolation_coordinate=args.stella_interpolation_coordinate,
        stella_fallback=args.stella_fallback, stella_t_ion_blend_mode=args.stella_t_ion_blend_mode,
        stella_t_ion_floor=args.stella_t_ion_floor, stella_t_ion_ceiling=args.stella_t_ion_ceiling,
        saha_inversion_for_ne=args.saha_inversion_for_ne, non_local_reversals=args.non_local_reversals
    )

    if args.complex_model:
        d = p.get_default
        if args.kappa_fallback == d("kappa_fallback"): conv.kappa_fallback = "kramers"
        if args.tau_mode == d("tau_mode"): conv.tau_mode = "nonmonotonic"
        if not args.no_smooth:
            conv.smooth = True
            if args.smooth_mode == d("smooth_mode"): conv.smooth_mode = "front_aware"
        if not args.velocity_aware_smoothing: conv.velocity_aware_smoothing = True
        if not args.temperature_aware_smoothing: conv.temperature_aware_smoothing = True
        if args.ne_damping_mode == d("ne_damping_mode"): conv.ne_damping_mode = "adaptive"
        if args.ne_max_log_step == d("ne_max_log_step"): conv.ne_max_log_step = 1.0
        if args.floor_ne == d("floor_ne"): conv.floor_ne = 10.0
        if args.t_exc_source == d("t_exc_source"): conv.t_exc_source = "t_ion"
        if not args.use_t_ion: conv.use_t_ion = True
        # Enable new optimizations in complex model preset
        conv.saha_inversion_for_ne = True
        conv.non_local_reversals = True

    if args.no_smooth: conv.smooth = False

    output_dir = args.output_dir or (args.model_dir if args.force_overwrite_source else os.path.join(args.model_dir, "proxy_v2_out"))
    inputs = load_model(args.model_dir, hydro_path=args.hydro_data, conv=conv)
    result = generate(inputs, conv, output_dir, prev_model_dir=args.prev_model, make_plots=args.plots)
    
    print(f"Wrote {result['n_species']} _IN files to {output_dir}")
    v = result['validation']
    print(f"Validation issues: {v['n_issues']}")
    for iss in v['issues']: print(f"  - {iss}")
    print(f"Ne ratio max deviation: {v['ne_ratio_max']:.3e}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
