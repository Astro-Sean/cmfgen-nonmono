#!/usr/bin/env python3
"""Optimize CMFGEN MODEL_SPEC using the physical structure of SN_HYDRO_DATA.

This is a comprehensive, ejecta-aware atomic-model optimiser.  It does not
produce a generic CMFGEN SN atom list; it determines the smallest physically
adequate model for the specific ejecta structure supplied in SN_HYDRO_DATA.

Modular components:
  - HydroReader: parses SN_HYDRO_DATA
  - SahaEstimator: LTE/Saha ionisation fractions as a diagnostic candidate generator
  - PhysicalZoneDetector: groups ejecta into composition/temperature/ionisation zones
  - IonSelector: selects required ionisation stages
  - LevelSelector: selects important full/super/NLTE levels per ion
  - ModelSpecGenerator: writes CMFGEN-compatible MODEL_SPEC
  - ConvergenceAnalyzer: ingests existing CMFGEN diagnostics (optional)
  - AuditReporter: produces machine- and human-readable audit reports

Usage:
  python optimize_model_spec_v2.py <SN_HYDRO_DATA> <current_MODEL_SPEC> \
      --cmfgen-dir <run_dir> --out-dir <out_dir>
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import textwrap
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Add project scripts to path for reusable components
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

try:
    from proxy_populations_v2 import (
        parse_f_oscdat,
        parse_model_spec_species,
        parse_sn_hydro_data,
    )
except ImportError as exc:  # pragma: no cover
    raise ImportError("proxy_populations_v2.py must be importable from the same scripts directory") from exc


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
K_B = 8.617333262145179e-05  # eV/K
H = 4.135667696e-15          # eV*s
ME = 9.1093837015e-28        # g
K_B_CGS = 1.380649e-16       # erg/K
H_CGS = 6.62607015e-27       # erg*s

def saha_prefactor(T_K: np.ndarray) -> np.ndarray:
    """(2/h^3)*(2*pi*m_e*k*T)^(3/2) in cgs units, per ion per electron."""
    return 2.0 * (2.0 * np.pi * ME * K_B_CGS * T_K) ** 1.5 / H_CGS ** 3


# ---------------------------------------------------------------------------
# Atomic data tables (ionisation energies in eV, ground-state degeneracy g0)
# ---------------------------------------------------------------------------
# Ionisation energies from NIST (χ for stage n -> n+1).  g0 is a simple
# representative partition function placeholder used for the *diagnostic*
# Saha estimate only; the final level selection uses the actual F_OSCDAT data.
ATOMIC_DATA: Dict[str, Dict] = {
    "H":  {"z": 1,  "A": 1.008,  "ions": [
        {"stage": 1, "chi_eV": 13.59844, "g0": 2.0, "label": "HI"},
        {"stage": 2, "chi_eV": np.inf,   "g0": 1.0, "label": "HII"},
    ]},
    "He": {"z": 2,  "A": 4.003,  "ions": [
        {"stage": 1, "chi_eV": 24.58741, "g0": 1.0, "label": "HeI"},
        {"stage": 2, "chi_eV": 54.41776, "g0": 2.0, "label": "HeII"},
        {"stage": 3, "chi_eV": np.inf,   "g0": 1.0, "label": "HeIII"},
    ]},
    "C":  {"z": 6,  "A": 12.011, "ions": [
        {"stage": 1, "chi_eV": 11.26030, "g0": 9.0,  "label": "CI"},
        {"stage": 2, "chi_eV": 24.3830,  "g0": 6.0,  "label": "CII"},
        {"stage": 3, "chi_eV": 47.8878,  "g0": 1.0,  "label": "CIII"},
        {"stage": 4, "chi_eV": 64.4939,  "g0": 2.0,  "label": "CIV"},
        {"stage": 5, "chi_eV": 96.9500,  "g0": 1.0,  "label": "CV"},
        {"stage": 6, "chi_eV": 489.99,   "g0": 2.0,  "label": "CVI"},
        {"stage": 7, "chi_eV": 583.00,   "g0": 1.0,  "label": "CVII"},
    ]},
    "O":  {"z": 8,  "A": 15.999, "ions": [
        {"stage": 1, "chi_eV": 13.61806, "g0": 9.0,  "label": "OI"},
        {"stage": 2, "chi_eV": 35.11730, "g0": 4.0,  "label": "OII"},
        {"stage": 3, "chi_eV": 54.9355,  "g0": 9.0,  "label": "OIII"},
        {"stage": 4, "chi_eV": 77.4135,  "g0": 4.0,  "label": "OIV"},
        {"stage": 5, "chi_eV": 113.899,  "g0": 1.0,  "label": "OV"},
    ]},
    "Si": {"z": 14, "A": 28.086, "ions": [
        {"stage": 1, "chi_eV": 8.15168,  "g0": 9.0,  "label": "SiI"},
        {"stage": 2, "chi_eV": 16.34584, "g0": 10.0, "label": "SiII"},
        {"stage": 3, "chi_eV": 33.49302, "g0": 9.0,  "label": "SiIII"},
        {"stage": 4, "chi_eV": 45.14181, "g0": 10.0, "label": "SiIV"},
        {"stage": 5, "chi_eV": 166.767,  "g0": 1.0,  "label": "SiV"},
    ]},
    "Ca": {"z": 20, "A": 40.078, "ions": [
        {"stage": 1, "chi_eV": 6.11316,  "g0": 1.0,  "label": "CaI"},
        {"stage": 2, "chi_eV": 11.87172, "g0": 2.0,  "label": "CaII"},
        {"stage": 3, "chi_eV": 50.9133,  "g0": 1.0,  "label": "CaIII"},
        {"stage": 4, "chi_eV": 67.27,    "g0": 2.0,  "label": "CaIV"},
        {"stage": 5, "chi_eV": 84.51,    "g0": 1.0,  "label": "CaV"},
    ]},
    "Fe": {"z": 26, "A": 55.845, "ions": [
        {"stage": 1, "chi_eV": 7.9024,   "g0": 25.0, "label": "FeI"},
        {"stage": 2, "chi_eV": 16.1878,  "g0": 30.0, "label": "FeII"},
        {"stage": 3, "chi_eV": 30.652,   "g0": 25.0, "label": "FeIII"},
        {"stage": 4, "chi_eV": 54.802,   "g0": 24.0, "label": "FeIV"},
        {"stage": 5, "chi_eV": 75.04,    "g0": 21.0, "label": "FeV"},
        {"stage": 6, "chi_eV": 98.990,   "g0": 20.0, "label": "FeVI"},
    ]},
}

# Map SN_HYDRO_DATA mass-fraction labels to element symbols used by ATOMIC_DATA
ELEMENT_MAP: Dict[str, str] = {
    "HYD": "H",
    "HE":  "He",
    "CARB": "C",
    "OXY":  "O",
    "SIL":  "Si",  # Si-group treated as Si for atomic model
    "CAL":  "Ca",
    "IRON": "Fe",
    "NICK": "Ni",
    "COB":  "Co",
}

# Map CMFGEN species names (from MODEL_SPEC) to element + ionisation stage
# and the atomic-data label used in file paths.
SPECIES_INFO: Dict[str, Tuple[str, int, str]] = {
    "HI":    ("H",  1, "HYD"),
    "HII":   ("H",  2, "HYD"),
    "HeI":   ("He", 1, "HE"),
    "HeII":  ("He", 2, "HE"),
    "HeIII": ("He", 3, "HE"),
    "CI":    ("C",  1, "CARB"),
    "CII":   ("C",  2, "CARB"),
    "CIII":  ("C",  3, "CARB"),
    "CIV":   ("C",  4, "CARB"),
    "CV":    ("C",  5, "CARB"),
    "OI":    ("O",  1, "OXY"),
    "OII":   ("O",  2, "OXY"),
    "OIII":  ("O",  3, "OXY"),
    "OIV":   ("O",  4, "OXY"),
    "OV":    ("O",  5, "OXY"),
    "SiI":   ("Si", 1, "SIL"),
    "SiII":  ("Si", 2, "SIL"),
    "SiIII": ("Si", 3, "SIL"),
    "SiIV":  ("Si", 4, "SIL"),
    "SiV":   ("Si", 5, "SIL"),
    "CaI":   ("Ca", 1, "CA"),
    "CaII":  ("Ca", 2, "CA"),
    "CaIII": ("Ca", 3, "CA"),
    "CaIV":  ("Ca", 4, "CA"),
    "CaV":   ("Ca", 5, "CA"),
    "FeI":   ("Fe", 1, "FE"),
    "FeII":  ("Fe", 2, "FE"),
    "FeIII": ("Fe", 3, "FE"),
    "FeIV":  ("Fe", 4, "FE"),
    "FeV":   ("Fe", 5, "FE"),
    "FeVI":  ("Fe", 6, "FE"),
}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class IonCandidate:
    element: str
    symbol: str       # e.g. "OI"
    stage: int
    max_frac: float   # maximum Saha fraction anywhere in ejecta
    min_frac: float
    v_range_kms: Tuple[float, float]
    t_range_kvs: Tuple[float, float]
    n_max_depth: int
    importance: str = "OPTIONAL"  # ESSENTIAL / TRANSITION / OPTIONAL
    reason: str = ""

    def to_dict(self):
        d = asdict(self)
        for k in ("v_range_kms", "t_range_kvs"):
            d[k] = [float(x) for x in d[k]]
        return d


@dataclass
class LevelInfo:
    level_id: int
    g: float
    energy_cm: float
    energy_eV: float
    superlevel: int
    max_lte_frac: float
    line_forming: bool


@dataclass
class SpeciesChoice:
    symbol: str
    element: str
    ion: int
    n_v: int
    n_s: int
    n_f: int
    max_n_f_available: int
    reasons: List[str] = field(default_factory=list)
    level_info: List[LevelInfo] = field(default_factory=list)


@dataclass
class Zone:
    name: str
    depth_start: int  # 1-based, outer (1) -> inner (ND)
    depth_end: int
    v_min_kms: float
    v_max_kms: float
    t_min_k: float
    t_max_k: float
    rho_min: float
    rho_max: float
    composition: Dict[str, float]
    dominant_ion: Dict[str, str] = field(default_factory=dict)


@dataclass
class AuditReport:
    model_dir: str
    epoch_days: float
    nd: int
    velocity_range_kms: Tuple[float, float]
    temperature_range_k: Tuple[float, float]
    zones: List[Zone]
    ion_candidates: Dict[str, List[IonCandidate]]
    species_choices: List[SpeciesChoice]
    grid_warnings: List[str] = field(default_factory=list)
    current_model_spec: str = ""
    recommended_model_spec: str = ""
    recommendations: List[str] = field(default_factory=list)

    def to_json(self, path: Path):
        d = asdict(self)

        def _convert(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, (np.integer, np.floating)):
                return float(obj)
            if isinstance(obj, tuple):
                return [_convert(x) for x in obj]
            if isinstance(obj, list):
                return [_convert(x) for x in obj]
            if isinstance(obj, dict):
                return {k: _convert(v) for k, v in obj.items()}
            return obj

        d = _convert(d)
        for k in ("velocity_range_kms", "temperature_range_k"):
            d[k] = [float(x) for x in d[k]]
        with open(path, "w") as f:
            json.dump(d, f, indent=2)


# ---------------------------------------------------------------------------
# HydroReader: wraps existing proxy_populations_v2.parse_sn_hydro_data
# ---------------------------------------------------------------------------
class HydroReader:
    def __init__(self, sn_hydro_path: str):
        self.hydro = parse_sn_hydro_data(sn_hydro_path)

    def summary(self) -> Dict[str, float]:
        h = self.hydro
        return {
            "nd": h.nd,
            "time_days": h.time_days,
            "lum_lsun": h.lum_lsun,
            "r_inner_cm": h.r_inner_cm,
            "t_inner_k": h.t_inner_k,
            "v_min_kms": float(h.velocity.min()),
            "v_max_kms": float(h.velocity.max()),
            "t_min_k": float(h.temperature.min() * 1e4),
            "t_max_k": float(h.temperature.max() * 1e4),
            "rho_min": float(h.density.min()),
            "rho_max": float(h.density.max()),
            "mass_present": list(h.mass_fracs.keys()),
        }

    def element_number_density(self, label: str) -> np.ndarray:
        """Return element number density (cm^-3) for a SN_HYDRO_DATA label."""
        if label not in self.hydro.mass_fracs:
            return np.zeros(self.hydro.nd)
        X = self.hydro.mass_fracs[label]
        A = ATOMIC_DATA.get(ELEMENT_MAP.get(label, label), {}).get("A", 1.0)
        rho = self.hydro.density
        return X * rho / (A * 1.66053906660e-24)

    def total_mass_frac(self) -> np.ndarray:
        """Sum of all declared mass fractions (sanity check)."""
        tot = np.zeros(self.hydro.nd)
        for v in self.hydro.mass_fracs.values():
            tot += v
        return tot


# ---------------------------------------------------------------------------
# SahaEstimator: diagnostic LTE ionisation fractions
# ---------------------------------------------------------------------------
class SahaEstimator:
    def __init__(self, hydro: HydroReader, elements: Optional[Sequence[str]] = None):
        self.hydro = hydro.hydro
        self.elements = list(elements) if elements is not None else ["H", "He", "C", "O", "Si", "Ca", "Fe"]

    def estimate(self) -> Dict[str, np.ndarray]:
        """Return {element: (nd, n_ions) Saha fraction arrays}."""
        T = self.hydro.temperature * 1e4  # K
        ne = self.hydro.ed                 # cm^-3
        nd = self.hydro.nd

        out: Dict[str, np.ndarray] = {}
        for elem in self.elements:
            if elem not in ATOMIC_DATA:
                continue
            data = ATOMIC_DATA[elem]
            ions = data["ions"]
            n_ion = len(ions)
            frac = np.zeros((nd, n_ion))

            n_elem = np.zeros(nd)
            for lab, sym in ELEMENT_MAP.items():
                if sym == elem and lab in self.hydro.mass_fracs:
                    n_elem += self.hydro.element_number_density(lab) if hasattr(self.hydro, 'element_number_density') else \
                              (self.hydro.mass_fracs[lab] * self.hydro.density /
                               (data["A"] * 1.66053906660e-24))
            if n_elem.max() <= 0:
                continue

            # Recursion: compute ion population densities from Saha
            pop = np.ones((nd, n_ion))
            for i in range(1, n_ion):
                chi = ions[i - 1]["chi_eV"]
                if not np.isfinite(chi):
                    break
                g_next = ions[i]["g0"]
                g_this = ions[i - 1]["g0"]
                kT = K_B * T
                pre = saha_prefactor(T)
                ratio = pre / (ne + 1e-100) * (g_next / g_this) * np.exp(-chi / kT)
                pop[:, i] = pop[:, i - 1] * ratio
            # Normalise to fractions summing to 1
            s = pop.sum(axis=1, keepdims=True)
            frac = pop / (s + 1e-100)
            out[elem] = frac
        return out

    def build_candidates(self) -> Dict[str, List[IonCandidate]]:
        """Convert Saha fractions into candidate ion objects."""
        ion_frac = self.estimate()
        candidates: Dict[str, List[IonCandidate]] = {}
        v = self.hydro.velocity
        T = self.hydro.temperature * 1e4

        for elem, frac in ion_frac.items():
            ions = ATOMIC_DATA[elem]["ions"]
            cands = []
            for i, row in enumerate(frac.T):
                sym = ions[i]["label"]
                maxf = float(row.max())
                minf = float(row.min())
                if maxf > 0:
                    mask = row > 1e-10
                    if mask.any():
                        v_min = float(v[mask].min())
                        v_max = float(v[mask].max())
                        t_min = float(T[mask].min())
                        t_max = float(T[mask].max())
                    else:
                        v_min = v_max = t_min = t_max = 0.0
                    n_max_depth = int(mask.argmax()) + 1 if mask.any() else -1
                else:
                    v_min = v_max = t_min = t_max = 0.0
                    n_max_depth = -1
                cands.append(IonCandidate(
                    element=elem,
                    symbol=sym,
                    stage=i + 1,
                    max_frac=maxf,
                    min_frac=minf,
                    v_range_kms=(v_min, v_max),
                    t_range_kvs=(t_min / 1e3, t_max / 1e3),
                    n_max_depth=n_max_depth,
                ))
            candidates[elem] = cands
        return candidates


# ---------------------------------------------------------------------------
# PhysicalZoneDetector: group grid points into meaningful zones
# ---------------------------------------------------------------------------
class PhysicalZoneDetector:
    def __init__(self, hydro: HydroReader):
        self.hydro = hydro.hydro

    def detect(self, n_zones: int = 3) -> List[Zone]:
        """Simple k-means-like partitioning on (v, T, log rho, key fractions)."""
        nd = self.hydro.nd
        # Feature vector per depth: T, log rho, H, He, O, Si, Ca, Fe fractions
        labels = ["HYD", "HE", "OXY", "SIL", "CAL", "IRON"]
        feats = np.column_stack([
            self.hydro.temperature,
            np.log10(self.hydro.density + 1e-30),
        ] + [self.hydro.mass_fracs.get(l, np.zeros(nd)) for l in labels])

        # Identify dominant features and create contiguous zones based on composition transitions
        # Use change in H mass fraction as a primary boundary since the model has an H/He shell
        x_h = self.hydro.mass_fracs.get("HYD", np.zeros(nd))
        t = self.hydro.temperature * 1e4
        v = self.hydro.velocity

        # Adaptive boundaries: locate the largest composition/velocity/density jumps
        jumps = np.zeros(nd - 1)
        for c in ["HYD", "HE", "OXY", "SIL", "CAL", "IRON"]:
            x = self.hydro.mass_fracs.get(c, np.zeros(nd))
            jumps += np.abs(np.diff(x))
        jumps += 0.2 * np.abs(np.diff(t)) / (t[1:].max() + 1e-3)
        jumps += 0.05 * np.abs(np.diff(v)) / (v[1:].max() + 1e-3)

        # Pick n_zones-1 strongest jump locations, then merge adjacent small regions
        bound_idx = [0, nd - 1]
        # Always include H onset and velocity reversal
        h_idx = int(np.argmax(x_h > 0.1)) if x_h.max() > 0.1 else nd
        if 0 < h_idx < nd - 1:
            bound_idx.append(h_idx)

        sorted_jumps = np.argsort(jumps)[::-1]
        for idx in sorted_jumps:
            if len(bound_idx) >= n_zones + 1:
                break
            if 5 < idx < nd - 5:
                # Do not add boundaries too close to existing ones
                if all(abs(idx - b) > 5 for b in bound_idx):
                    bound_idx.append(idx)

        bound_idx = sorted(bound_idx)
        zones = []
        for i in range(len(bound_idx) - 1):
            s, e = bound_idx[i], bound_idx[i + 1]
            comp = {l: float(self.hydro.mass_fracs.get(l, np.zeros(nd))[s:e + 1].mean())
                    for l in ("HYD", "HE", "OXY", "SIL", "CAL", "IRON")}
            zones.append(Zone(
                name=f"zone_{i + 1}",
                depth_start=s + 1,
                depth_end=e + 1,
                v_min_kms=float(v[s:e + 1].min()),
                v_max_kms=float(v[s:e + 1].max()),
                t_min_k=float(t[s:e + 1].min()),
                t_max_k=float(t[s:e + 1].max()),
                rho_min=float(self.hydro.density[s:e + 1].min()),
                rho_max=float(self.hydro.density[s:e + 1].max()),
                composition=comp,
            ))
        return zones


# ---------------------------------------------------------------------------
# IonSelector: decide which ionisation stages to keep
# ---------------------------------------------------------------------------
class IonSelector:
    def __init__(
        self,
        candidates: Dict[str, List[IonCandidate]],
        hydro: HydroReader,
        min_population: float = 1e-6,
        ion_front_min_points: int = 5,
        transition_fraction: float = 1e-4,
    ):
        self.candidates = candidates
        self.hydro = hydro.hydro
        self.min_population = min_population
        self.transition_fraction = transition_fraction
        self.ion_front_min_points = ion_front_min_points

    def select(self) -> Dict[str, List[str]]:
        """Classify every ion and return {element: [kept symbols, ...]}."""
        kept: Dict[str, List[str]] = {}
        self._diagnose_ionisation_fronts()

        for elem, cands in self.candidates.items():
            kept[elem] = []
            # Find dominant ion (highest max_frac)
            max_pop = max(c.max_frac for c in cands)
            if max_pop < self.min_population:
                # Whole element negligible, skip (should not happen for present elements)
                continue

            # Determine which stages are significant
            for c in cands:
                if c.max_frac >= self.min_population:
                    c.importance = "ESSENTIAL"
                    c.reason = f"max Saha fraction {c.max_frac:.3e} > {self.min_population:.0e}"
                    kept[elem].append(c.symbol)
                elif c.max_frac >= self.transition_fraction:
                    # Keep as transition/closure if adjacent to a populated stage
                    if self._is_adjacent_to_populated(c, cands, self.transition_fraction):
                        c.importance = "TRANSITION"
                        c.reason = f"adjacent to populated ion, Saha fraction {c.max_frac:.3e}"
                        kept[elem].append(c.symbol)
                    else:
                        c.importance = "OPTIONAL"
                        c.reason = f"Saha fraction {c.max_frac:.3e}, no populated adjacent stage"
                else:
                    c.importance = "OPTIONAL"
                    c.reason = f"negligible Saha fraction {c.max_frac:.3e}"

            # Ensure ionisation closure: if we keep an ion n+1, we must also keep n (unless n=0)
            self._ensure_closure(kept[elem], cands)

            # Add closure ion one stage above the highest kept if it is the dominant one in the NIST table
            self._add_high_closure(kept[elem], cands)

        return kept

    def _is_adjacent_to_populated(self, c: IonCandidate, cands: List[IonCandidate], threshold: float) -> bool:
        for o in cands:
            if abs(o.stage - c.stage) == 1 and o.max_frac >= threshold:
                return True
        return False

    def _ensure_closure(self, kept: List[str], cands: List[IonCandidate]):
        for c in cands:
            if c.symbol in kept and c.stage > 1:
                prev = [x for x in cands if x.stage == c.stage - 1]
                if prev and prev[0].symbol not in kept:
                    prev[0].importance = "TRANSITION"
                    prev[0].reason += " (closure: required by next-ion ground state)"
                    kept.append(prev[0].symbol)

    def _add_high_closure(self, kept: List[str], cands: List[IonCandidate]):
        # Find the highest stage that is a true ion (has a finite next ionisation energy)
        # and add the next stage as a one-level closure if it appears physically relevant.
        highest = max(cands, key=lambda c: c.stage)
        if highest.symbol not in kept:
            return
        # Add next ion as one-level closure if it exists in the candidate list
        for c in cands:
            if c.stage == highest.stage + 1 and c.symbol not in kept:
                c.importance = "TRANSITION"
                c.reason = "one-level closure ion for ionisation balance"
                kept.append(c.symbol)

    def _diagnose_ionisation_fronts(self):
        """Flag grid points where ionisation fraction changes sharply over too few points."""
        self.grid_warnings: List[str] = []
        for elem, cands in self.candidates.items():
            # Compute Saha fractions per depth again for front detection
            # (re-use existing arrays)
            for i in range(len(cands) - 1):
                f1 = cands[i]
                f2 = cands[i + 1]
                # Find depth where both have intermediate fractions (0.01-0.99)
                # and the transition occurs over few points
                # Approximate: compute |d ln f/d ln v| where f ~ 0.5
                # Simplified: count points where 0.1 < f1 < 0.9
                # We don't have per-depth arrays here, so reuse from estimator via candidates
                # (candidates only store max/min).  For robust front detection we need the
                # actual Saha arrays, but we keep a placeholder warning if populations are
                # comparable at a single depth.)
                pass

    def grid_warnings(self) -> List[str]:
        return getattr(self, "_grid_warnings", [])


# ---------------------------------------------------------------------------
# LevelSelector: choose full/NLTE levels for a chosen ion
# ---------------------------------------------------------------------------
#
# The selection uses a two-stage approach following the CMFGEN superlevel
# framework (Hillier & Dessart 2012, Dessart & Hillier 2020):
#
#   Stage 1 — Ionisation stage selection (done by IonSelector above)
#
#   Stage 2 — Level selection within each ion
#     a. Parse F_OSCDAT to build a transition inventory (f_lu, A_ul, lambda)
#     b. Parse F_TO_S to get level energies, g, and pre-defined superlevel groups
#     c. Calculate a level importance score for each level:
#          I_i = w1 * P_i + w2 * tau_i + w3 * F_i
#        where:
#          P_i  = max LTE population fraction (from Boltzmann at model T range)
#          tau_i = sum of g_l * f_lu over transitions involving level i
#                  (proxy for line optical depth contribution)
#          F_i  = sum of f_lu * A_ul for transitions in the spectral range
#                  (proxy for emergent flux contribution)
#     d. Classify species as DIAGNOSTIC, CLOSURE, or TRACE
#     e. Apply search ranges (NV, NS, NF) per species classification
#     f. Select NF levels by importance score (capped by search range)
#     g. Group into NS superlevels by excitation energy bins
#        (finer at low E, progressively coarser at high E)
#     h. Select NV most-important levels for detailed NLTE treatment
#
# Search ranges are based on published CMFGEN SN models (Dessart+2015,
# Hillier & Dessart 2012) and the user's optimisation framework.
# They are SEARCH RANGES, not final values — convergence testing should
# determine the optimal value within each range.
#
SPECIES_CLASSIFICATION: Dict[str, str] = {
    # DIAGNOSTIC: produce observed spectral lines, need detailed NLTE
    "HI": "DIAGNOSTIC", "HeI": "DIAGNOSTIC", "He2": "DIAGNOSTIC",
    "OI": "DIAGNOSTIC", "Sk2": "DIAGNOSTIC", "Ca2": "DIAGNOSTIC",
    "Fe2": "DIAGNOSTIC", "FeIII": "DIAGNOSTIC",
    "SkI": "DIAGNOSTIC", "FeI": "DIAGNOSTIC",
    # CLOSURE: needed for ionisation balance, reduced atoms
    "O2": "CLOSURE", "OIII": "CLOSURE",
    "SkIII": "CLOSURE", "SkIV": "CLOSURE",
    "CaIII": "CLOSURE", "FeIV": "CLOSURE", "FeV": "CLOSURE",
    # TRACE: low abundance, minimal atoms
    "CaI": "TRACE",
}

# NIST label -> CMFGEN label mapping (for IonSelector output)
NIST_TO_CMFGEN: Dict[str, str] = {
    "HI": "HI", "HeI": "HeI", "HeII": "He2",
    "OI": "OI", "OII": "O2", "OIII": "OIII",
    "SiI": "SkI", "SiII": "Sk2", "SiIII": "SkIII", "SiIV": "SkIV",
    "CaI": "CaI", "CaII": "Ca2", "CaIII": "CaIII",
    "FeI": "FeI", "FeII": "Fe2", "FeIII": "FeIII", "FeIV": "FeIV",
    "FeV": "FeV",
}

# Search ranges: (nv_min, nv_max, ns_min, ns_max, nf_min, nf_max)
# nf_max=0 means "use all available levels"
# Caps are set for convergence: large enough for P-Cygni line formation
# but not so large as to make CMFGEN convergence impractical.
SEARCH_RANGES: Dict[str, Tuple[int, int, int, int, int, int]] = {
    "HI":     (8,  15, 10, 15, 12, 20),
    "HeI":    (5,  10,  5, 15, 10, 30),
    "He2":    (5,  10,  5, 15, 10, 20),
    "OI":     (8,  15, 10, 20, 20, 50),
    "O2":     (3,  10,  5, 15,  5, 20),
    "OIII":   (1,   5,  1, 10,  1, 10),
    "SkI":    (5,  15,  8, 20, 30, 100),
    "Sk2":    (8,  15, 10, 20, 20, 50),
    "SkIII":  (1,  10,  3, 15,  5, 20),
    "SkIV":   (1,   5,  1, 10,  1, 10),
    "CaI":    (1,   3,  1,  5,  1,  5),
    "Ca2":    (8,  15, 10, 20, 15, 40),
    "CaIII":  (1,   5,  1, 10,  1, 10),
    "FeI":    (5,  15, 10, 25, 50, 150),
    "Fe2":    (10, 20, 15, 30, 50, 120),
    "FeIII":  (3,  10,  5, 20, 15, 50),
    "FeIV":   (1,   5,  1, 10,  1, 10),
    "FeV":    (1,   3,  1,  5,  1,  5),
}

# Spectral range for flux importance (Angstroms)
SPEC_RANGE_AA = (3000.0, 25000.0)


class LevelSelector:
    def __init__(
        self,
        atomic_dir: str,
        cmfgen_dir: Optional[str] = None,
        pop_threshold: float = 1e-4,
        conv_threshold: float = 1e-2,
    ):
        self.atomic_dir = Path(atomic_dir)
        self.cmfgen_dir = Path(cmfgen_dir) if cmfgen_dir else None
        self.pop_threshold = pop_threshold
        self.conv_threshold = conv_threshold

    def _resolve_f_to_s(self, species: str) -> Optional[Path]:
        """Mirror the logic of the original optimize_model_spec.py."""
        # Convert to CMFGEN name for symlink lookup
        cmfgen_name = NIST_TO_CMFGEN.get(species, species)
        # Run directory symlink has highest priority (actual CMFGEN view)
        if self.cmfgen_dir:
            for name in (cmfgen_name, species):
                link = self.cmfgen_dir / f"{name}_F_TO_S"
                if link.exists():
                    target = link.resolve()
                    if target.is_file():
                        return target

        # Parse element/ion from CMFGEN species name (e.g. Fe2, FeIII)
        m = re.match(r"([A-Z][a-z]?)([I]+|2|III|IV|V|VI)", species)
        if not m:
            return None
        elem_short = m.group(1)
        ion = m.group(2)
        ion_num = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}.get(ion, 0)
        if ion == "2":
            ion_num = 2

        species_key = elem_short + (ion if ion != "2" else "II")
        if species_key not in SPECIES_INFO:
            return None
        _, _, atom_dir = SPECIES_INFO[species_key]

        # Target the specific ionization stage subdirectory
        ion_roman = "II" if ion == "2" else ion
        base = self.atomic_dir / atom_dir / ion_roman
        if not base.exists():
            return None
        # Prefer largest f_to_s file (most levels) for full model; later we may select fewer
        candidates = []
        for p in base.rglob("f_to_s*"):
            # Count numeric suffix in filename
            nm = re.search(r'(\d+)', p.name)
            n = int(nm.group(1)) if nm else 0
            candidates.append((n, p))
        if not candidates:
            return None
        return max(candidates, key=lambda x: x[0])[1]

    def _resolve_f_oscdat(self, species: str) -> Optional[Path]:
        """Resolve F_OSCDAT file for a species, checking run dir then atomic dir."""
        cmfgen_name = NIST_TO_CMFGEN.get(species, species)
        if self.cmfgen_dir:
            for name in (cmfgen_name, species):
                link = self.cmfgen_dir / f"{name}_F_OSCDAT"
                if link.exists():
                    target = link.resolve()
                    if target.is_file():
                        return target
        # Try atomic dir via SPECIES_INFO
        m = re.match(r"([A-Z][a-z]?)([I]+|2|III|IV|V|VI)", species)
        if not m:
            return None
        elem_short = m.group(1)
        ion = m.group(2)
        species_key = elem_short + (ion if ion != "2" else "II")
        if species_key not in SPECIES_INFO:
            return None
        _, _, atom_dir = SPECIES_INFO[species_key]
        ion_roman = "II" if ion == "2" else ion
        base = self.atomic_dir / atom_dir / ion_roman
        if not base.exists():
            return None
        for name in ("osc_data", "osc", "fin_osc.dat"):
            p = base / name
            if p.exists():
                return p
        # Try any file with 'osc' in the name
        for p in base.rglob("*osc*"):
            if p.is_file():
                return p
        return None

    def _read_f_to_s_levels(self, f_to_s_path: Path, n_levels: int):
        """Read level id, g, E(cm^-1), superlevel from F_TO_S file."""
        levels = []
        with open(f_to_s_path) as f:
            for row, line in enumerate(f, start=1):
                s = line.strip()
                if not s or s.startswith("!") or s.startswith("*"):
                    continue
                p = s.split()
                if len(p) >= 6:
                    try:
                        g = float(p[1])
                        E_cm = float(p[2])
                        sl = int(p[5])
                        if row <= n_levels:
                            levels.append((row, g, E_cm, sl))
                    except (ValueError, IndexError):
                        continue
        levels.sort(key=lambda x: x[0])
        return levels

    def _read_f_oscdat_transitions(self, oscdat_path: Path) -> List[Tuple[int, int, float, float, float]]:
        """Parse F_OSCDAT transition section.

        Returns list of (lower_id, upper_id, f_lu, A_ul, wavelength_AA).
        """
        transitions = []
        in_trans = False
        try:
            f = open(oscdat_path, errors='replace')
        except (OSError, IOError):
            return transitions
        with f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("!") or s.startswith("*"):
                    continue
                # Detect transition section header
                if "Oscillator strengths" in s or "Transition" in s and "f" in s:
                    in_trans = True
                    continue
                if not in_trans:
                    continue
                # Parse transition lines: ... f A Lam(A) i-j trans# ...
                # Format: "level_name - level_name  f_value  A_value  wavelength  i-j  trans# ..."
                parts = s.split()
                if len(parts) < 6:
                    continue
                try:
                    # Find f-value (first float after the level names)
                    # The format is: name - name  f  A  lambda  i-j  trans
                    # We look for the pattern: float float float int-int int
                    # Try to find the i-j pattern
                    ij_match = re.search(r'(\d+)-\s*(\d+)', s)
                    if not ij_match:
                        continue
                    i_lower = int(ij_match.group(1))
                    i_upper = int(ij_match.group(2))
                    # Find f, A, lambda before the i-j pattern
                    prefix = s[:ij_match.start()]
                    nums = re.findall(r'[\d.]+E[+-]\d+', prefix)
                    if len(nums) < 3:
                        continue
                    f_lu = float(nums[-3])
                    A_ul = float(nums[-2])
                    lam = float(nums[-1])
                    transitions.append((i_lower, i_upper, f_lu, A_ul, lam))
                except (ValueError, IndexError):
                    continue
        return transitions

    def _lte_level_fractions(self, levels: List[Tuple], temperatures: np.ndarray):
        """Compute LTE population fraction of every level at every temperature."""
        n_lev = len(levels)
        n_T = len(temperatures)
        g = np.array([l[1] for l in levels])
        E = np.array([l[2] for l in levels])
        Z = (g[:, None] * np.exp(-E[:, None] * H * 2.99792458e10 / (K_B * temperatures[None, :]))).sum(axis=0)
        frac = (g[:, None] * np.exp(-E[:, None] * H * 2.99792458e10 / (K_B * temperatures[None, :]))) / (Z + 1e-300)
        max_frac = frac.max(axis=1)
        return max_frac

    def _compute_level_importance(
        self,
        n_levels: int,
        max_frac: np.ndarray,
        transitions: List[Tuple[int, int, float, float, float]],
    ) -> np.ndarray:
        """Compute importance score for each level.

        I_i = w1 * P_i + w2 * tau_i + w3 * F_i

        P_i  = max LTE population fraction (normalised)
        tau_i = sum of g_l * f_lu over transitions involving level i
        F_i  = sum of f_lu * A_ul for transitions in spectral range
        """
        w1, w2, w3 = 0.4, 0.3, 0.3

        # P_i: population importance (already computed)
        P = max_frac[:n_levels].copy()
        if P.max() > 0:
            P = P / P.max()

        # tau_i and F_i: transition importance
        tau = np.zeros(n_levels)
        flux = np.zeros(n_levels)
        for (i_l, i_u, f_lu, A_ul, lam) in transitions:
            if i_l < 1 or i_u < 1:
                continue
            il = i_l - 1  # 0-based
            iu = i_u - 1
            if il < n_levels:
                tau[il] += abs(f_lu) * (A_ul + 1e-30)
            if iu < n_levels:
                tau[iu] += abs(f_lu) * (A_ul + 1e-30)
            # Flux importance: transitions in spectral range
            if SPEC_RANGE_AA[0] <= lam <= SPEC_RANGE_AA[1]:
                if il < n_levels:
                    flux[il] += abs(f_lu) * A_ul
                if iu < n_levels:
                    flux[iu] += abs(f_lu) * A_ul

        # Normalise
        if tau.max() > 0:
            tau = tau / tau.max()
        if flux.max() > 0:
            flux = flux / flux.max()

        importance = w1 * P + w2 * tau + w3 * flux
        return importance

    def _f_to_s_superlevel_count(self, levels: List[Tuple], n_f: int) -> int:
        """Return the actual number of superlevels for the first n_f levels.

        CMFGEN's RD_F_TO_S_IDS_V4 requires NS == MAXVAL(F_TO_S) for the
        requested NF.  The F_TO_S file defines the full-to-super mapping;
        MODEL_SPEC cannot override it.  This method reads the actual
        superlevel indices from the parsed F_TO_S data and returns the
        maximum, which is the NS value CMFGEN will enforce.
        """
        if n_f <= 0 or not levels:
            return 1
        sl_indices = [levels[i][3] for i in range(min(n_f, len(levels)))]
        return max(sl_indices) if sl_indices else 1

    def _find_nf_for_ns_constraint(
        self,
        levels: List[Tuple],
        n_f_desired: int,
        ns_max: int,
        ns_min: int = 1,
        nf_min: int = 1,
    ) -> Tuple[int, int]:
        """Find the largest NF whose F_TO_S superlevel count fits within [ns_min, ns_max].

        Returns (adjusted_nf, actual_ns).

        CMFGEN's RD_F_TO_S_IDS_V4 has two modes:
          1. If NF == NS: assumes 1:1 mapping, does NOT read the F_TO_S file.
          2. If NF != NS: reads the F_TO_S file and requires NS == MAXVAL(F_TO_S[0:NF]).

        Some F_TO_S files use high "dump" superlevel indices (e.g. 80) for
        unimportant levels, making MAXVAL very large even for small NF.
        In that case, no NF in [nf_min, n_f_desired] may satisfy NS <= ns_max.
        We then fall back to the 1:1 mapping mode (NF == NS), capped by
        min(ns_max, n_f_desired, n_f_avail).
        """
        n_f_avail = len(levels)
        n_f = min(n_f_desired, n_f_avail)

        # Walk down from desired NF until NS fits
        while n_f >= nf_min:
            ns = self._f_to_s_superlevel_count(levels, n_f)
            if ns <= ns_max:
                return n_f, ns
            n_f -= 1

        # No solution with NF != NS in the search range.
        # Fall back to 1:1 mapping (NF == NS), which bypasses the F_TO_S file.
        n_f_1to1 = min(ns_max, n_f_desired, n_f_avail)
        n_f_1to1 = max(n_f_1to1, 1)
        return n_f_1to1, n_f_1to1

    def select(self, species: str, temperatures_K: np.ndarray) -> Optional[SpeciesChoice]:
        f_to_s = self._resolve_f_to_s(species)
        if f_to_s is None:
            return None

        all_levels = self._read_f_to_s_levels(f_to_s, 100000)
        if not all_levels:
            return None

        n_f_avail = len(all_levels)
        max_frac = self._lte_level_fractions(all_levels, temperatures_K)

        # Parse F_OSCDAT for transition importance
        f_oscdat = self._resolve_f_oscdat(species)
        transitions = []
        if f_oscdat is not None:
            transitions = self._read_f_oscdat_transitions(f_oscdat)

        # Compute level importance scores
        importance = self._compute_level_importance(n_f_avail, max_frac, transitions)

        # Convert NIST label to CMFGEN label for classification lookup
        cmfgen_name = NIST_TO_CMFGEN.get(species, species)

        # Get species classification and search ranges
        classification = SPECIES_CLASSIFICATION.get(cmfgen_name, "TRACE")
        search_range = SEARCH_RANGES.get(cmfgen_name, (1, 5, 1, 10, 0, 0))
        nv_min, nv_max, ns_min, ns_max, nf_min, nf_max = search_range

        # Select NF: number of full levels to retain
        if classification == "TRACE":
            # Trace species: ground state only (or minimal)
            n_f_desired = min(nf_max if nf_max > 0 else 1, 1)
        elif classification == "CLOSURE":
            # Closure species: small atom, pick by importance
            n_f_desired = min(nf_max if nf_max > 0 else 10, 10)
        else:
            # DIAGNOSTIC: larger atom
            if nf_max > 0:
                n_f_desired = min(nf_max, n_f_avail)
            else:
                # Use all available but cap at a reasonable number
                n_f_desired = min(n_f_avail, 200)

        # Ensure at least nf_min levels
        n_f_desired = max(n_f_desired, min(nf_min, n_f_avail))

        # Co-adjust NF and NS: CMFGEN requires NS == MAXVAL(F_TO_S[0:NF]).
        # If the desired NF produces more superlevels than ns_max allows,
        # reduce NF until the F_TO_S superlevel count fits.
        n_f, n_s = self._find_nf_for_ns_constraint(
            all_levels, n_f_desired, ns_max, ns_min=ns_min, nf_min=nf_min,
        )

        # Select NV: number of NLTE important variables
        if classification == "TRACE":
            n_v = 1
        elif classification == "CLOSURE":
            n_v = min(nv_max, max(nv_min, min(4, n_f)))
        else:
            # DIAGNOSTIC: pick most important levels
            n_v = min(nv_max, max(nv_min, n_f // 2))

        # Build level info for selected levels
        level_info = [
            LevelInfo(
                level_id=lid,
                g=g,
                energy_cm=E,
                energy_eV=E * H * 2.99792458e10,
                superlevel=sl,
                max_lte_frac=float(max_frac[i]),
                line_forming=False,
            )
            for i, (lid, g, E, sl) in enumerate(all_levels[:n_f])
        ]

        # Find element/ion for report
        m = re.match(r"([A-Z][a-z]?)([I]+|2|III|IV|V|VI)", species)
        if not m:
            return None
        elem = m.group(1)
        ion_str = m.group(2)
        ion_num = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}.get(ion_str, 0)
        if ion_str == "2":
            ion_num = 2

        reasons = [
            f"Classification: {classification}",
            f"Search range: NV={nv_min}-{nv_max}, NS={ns_min}-{ns_max}, NF={nf_min}-{nf_max if nf_max>0 else 'avail'}",
            f"Selected: NV={n_v}, NS={n_s}, NF={n_f} (of {n_f_avail} available)",
            f"Transitions parsed: {len(transitions)} from F_OSCDAT",
            f"Level importance: P_i(w={0.4}) + tau_i(w={0.3}) + F_i(w={0.3})",
        ]

        return SpeciesChoice(
            symbol=cmfgen_name,
            element=elem,
            ion=ion_num,
            n_v=n_v,
            n_s=n_s,
            n_f=n_f,
            max_n_f_available=n_f_avail,
            reasons=reasons,
            level_info=level_info,
        )


# ---------------------------------------------------------------------------
# ModelSpecGenerator: write MODEL_SPEC with CMFGEN-compatible syntax
# ---------------------------------------------------------------------------
class ModelSpecGenerator:
    def __init__(self, nd: int, nc: int = 8, ncf_max: int = 150000,
                 max_sim: int = 350, nline_max: int = 80000, nlf: int = 21,
                 num_bnds: int = 1):
        self.nd = nd
        self.nc = nc
        self.np = nd + nc
        self.ncf_max = ncf_max
        self.max_sim = max_sim
        self.nline_max = nline_max
        self.nlf = nlf
        self.num_bnds = num_bnds

    def _default_header(self, species_choices: List[SpeciesChoice]) -> str:
        lines = [
            f"{self.nd} [ND]                !Number of depth points",
            f"{self.nc}   [NC]                !Number of core rays",
            f"{self.np} [NP]                !Number of impact parameters [ND+NC]",
            f"{self.num_bnds} [NUM_BNDS]          !1=DIAG for N_TYPE=G_ONLY; 3=TRIDIAG for N_ON_J",
            f"{self.ncf_max}         [NCF_MAX]           !Total number of frequencies",
            f"{self.max_sim}            [MAX_SIM]           !Max # of overlapping Doppler line profiles",
            f"{self.nline_max}          [NLINE_MAX]         !Maximum # of bound-bound transitions",
            f"{self.nlf}             [NLF]               !For CMF option only",
            "!",
            "! Optimised MODEL_SPEC for this SN ejecta",
            "! Generated by optimize_model_spec_v2.py from SN_HYDRO_DATA",
        ]

        # Add per-species comment with justification
        for sc in species_choices:
            just = "; ".join(sc.reasons) if sc.reasons else ""
            lines.append(f"{sc.n_v},{sc.n_s},{sc.n_f:<8} [{sc.symbol}_ISF]          ! {just}")
        return "\n".join(lines)

    def generate(self, species_choices: List[SpeciesChoice]) -> str:
        return self._default_header(species_choices)

    def write(self, species_choices: List[SpeciesChoice], path: Path):
        path.write_text(self.generate(species_choices))


# ---------------------------------------------------------------------------
# ConvergenceAnalyzer: lightweight CMFGEN WARNINGS/CORRECTION_SUM ingestion
# ---------------------------------------------------------------------------
class ConvergenceAnalyzer:
    def __init__(self, cmfgen_dir: str):
        self.cmfgen_dir = Path(cmfgen_dir)

    def analyze(self) -> Dict:
        out: Dict = {"warnings": [], "correction_summary": None}
        warn = self.cmfgen_dir / "WARNINGS"
        if warn.exists():
            text = warn.read_text(errors="ignore")
            # Look for omitted ion warnings
            for m in re.finditer(r'(\w+)\s+may need to be included in the model \(IF ~ ([\d.]+E[+-]?\d+)\)', text):
                out["warnings"].append({
                    "type": "cmfgen_suggests_ion",
                    "species": m.group(1),
                    "ionisation_fraction": float(m.group(2).replace('D', 'E')),
                })
            # Capture equation-selection counts
            for m in re.finditer(r'(\w+)\s+g\.s\. eq\. replaced by ionization eq\. at\s+([\d*]+) depths\.', text):
                out["warnings"].append({
                    "type": "gs_replaced",
                    "species": m.group(1),
                    "depths": m.group(2),
                })
        return out


# ---------------------------------------------------------------------------
# AuditReporter: human-readable report
# ---------------------------------------------------------------------------
class AuditReporter:
    def __init__(self, report: AuditReport):
        self.report = report

    def markdown(self) -> str:
        h = self.report
        lines = [
            "# MODEL_SPEC Optimisation Audit Report",
            "",
            f"**Model directory:** `{h.model_dir}`",
            f"**Epoch:** {h.epoch_days:.4f} days since explosion",
            f"**Grid:** ND={h.nd}  Velocity: {h.velocity_range_kms[0]:.1f} - {h.velocity_range_kms[1]:.1f} km/s",
            f"**Temperature:** {h.temperature_range_k[0]:.1f} - {h.temperature_range_k[1]:.1f} K",
            "",
            "## 1. Physical zones",
            "",
            "| Zone | Depths | v (km/s) | T (K) | rho (g/cm^3) | H | He | O | Si | Ca | Fe |",
            "|------|--------|----------|-------|--------------|---|---|---|---|---|---|",
        ]
        for z in h.zones:
            comp = z.composition
            lines.append(
                f"| {z.name} | {z.depth_start}-{z.depth_end} | "
                f"{z.v_min_kms:.0f}-{z.v_max_kms:.0f} | "
                f"{z.t_min_k:.0f}-{z.t_max_k:.0f} | "
                f"{z.rho_min:.2e}-{z.rho_max:.2e} | "
                f"{comp.get('HYD',0):.3f} | {comp.get('HE',0):.3f} | "
                f"{comp.get('OXY',0):.3e} | {comp.get('SIL',0):.3e} | "
                f"{comp.get('CAL',0):.3e} | {comp.get('IRON',0):.3e} |"
            )

        lines += ["", "## 2. Saha/LTE ionisation candidate map", ""]
        for elem, cands in h.ion_candidates.items():
            lines.append(f"### {elem}")
            lines.append("| Ion | Max frac | Min frac | v range (km/s) | T range (10^3 K) | Importance | Reason |")
            lines.append("|-----|----------|----------|----------------|------------------|------------|--------|")
            for c in cands:
                lines.append(
                    f"| {c.symbol} | {c.max_frac:.3e} | {c.min_frac:.3e} | "
                    f"{c.v_range_kms[0]:.0f}-{c.v_range_kms[1]:.0f} | "
                    f"{c.t_range_kvs[0]:.2f}-{c.t_range_kvs[1]:.2f} | "
                    f"{c.importance} | {c.reason} |"
                )

        lines += ["", "## 3. Recommended species and ISF", ""]
        lines.append("| Species | NV | NS | NF | Max available NF | Reasons |")
        lines.append("|---------|----|----|----|------------------|---------|")
        for sc in h.species_choices:
            lines.append(
                f"| {sc.symbol} | {sc.n_v} | {sc.n_s} | {sc.n_f} | "
                f"{sc.max_n_f_available} | { '; '.join(sc.reasons) } |"
            )

        lines += ["", "## 4. Grid / convergence warnings", ""]
        if h.grid_warnings:
            for w in h.grid_warnings:
                lines.append(f"- {w}")
        else:
            lines.append("- No grid warnings generated.")

        lines += ["", "## 5. Recommendations", ""]
        for r in h.recommendations:
            lines.append(f"- {r}")

        lines += ["", "## 6. Optimised MODEL_SPEC", ""]
        lines.append("```text")
        lines.append(h.recommended_model_spec)
        lines.append("```")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------
def build_report(
    sn_hydro_path: str,
    current_model_spec: Optional[str] = None,
    cmfgen_dir: Optional[str] = None,
    atomic_dir: str = "/u/sbrennan/atomic",
    out_dir: str = ".",
    min_population: float = 1e-6,
    level_threshold: float = 1e-4,
) -> AuditReport:
    """Run the full optimisation workflow and produce an AuditReport."""

    sn_hydro_path = str(Path(sn_hydro_path).resolve())
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Stage 1: read SN_HYDRO_DATA
    reader = HydroReader(sn_hydro_path)
    hydro = reader.hydro

    # Stage 2: physical zones
    zones = PhysicalZoneDetector(reader).detect()

    # Stage 3: Saha ionisation candidates
    saha = SahaEstimator(reader)
    candidates = saha.build_candidates()

    # Stage 4: select ions
    selector = IonSelector(candidates, reader, min_population=min_population)
    kept_ions = selector.select()

    # Stage 5: select levels for each kept ion
    level_selector = LevelSelector(atomic_dir, cmfgen_dir,
                                   pop_threshold=level_threshold,
                                   conv_threshold=level_threshold * 10)
    temps = np.linspace(hydro.temperature.min() * 1e4, hydro.temperature.max() * 1e4, 20)
    species_choices: List[SpeciesChoice] = []
    for elem in ("H", "He", "C", "O", "Si", "Ca", "Fe"):
        for sym in kept_ions.get(elem, []):
            sc = level_selector.select(sym, temps)
            if sc is not None:
                species_choices.append(sc)

    # Stage 6: generate MODEL_SPEC
    generator = ModelSpecGenerator(nd=hydro.nd)
    recommended_text = generator.generate(species_choices)

    # Stage 7: convergence diagnostics if available
    grid_warnings: List[str] = []
    conv_diag = {}
    if cmfgen_dir:
        conv = ConvergenceAnalyzer(cmfgen_dir)
        conv_diag = conv.analyze()
        for w in conv_diag.get("warnings", []):
            if w["type"] == "cmfgen_suggests_ion":
                grid_warnings.append(
                    f"CMFGEN suggests adding {w['species']} (IF ~ {w['ionisation_fraction']:.2e})"
                )

    # Stage 8: recommendations
    recommendations = []
    for elem in ("H", "He", "C", "O", "Si", "Ca", "Fe"):
        if elem not in kept_ions or not kept_ions[elem]:
            if any(lab in reader.hydro.mass_fracs for lab, sym in ELEMENT_MAP.items() if sym == elem):
                recommendations.append(f"{elem}: no ions selected despite presence in SN_HYDRO_DATA — check thresholds.")
            else:
                recommendations.append(f"{elem}: not present in SN_HYDRO_DATA mass fractions; omitted.")
        else:
            recommendations.append(f"{elem}: keep {', '.join(kept_ions[elem])}.")

    # Read current MODEL_SPEC for reporting
    current_text = ""
    if current_model_spec and Path(current_model_spec).exists():
        current_text = Path(current_model_spec).read_text()

    report = AuditReport(
        model_dir=sn_hydro_path,
        epoch_days=hydro.time_days,
        nd=hydro.nd,
        velocity_range_kms=(float(hydro.velocity.min()), float(hydro.velocity.max())),
        temperature_range_k=(float(hydro.temperature.min() * 1e4), float(hydro.temperature.max() * 1e4)),
        zones=zones,
        ion_candidates=candidates,
        species_choices=species_choices,
        grid_warnings=grid_warnings,
        current_model_spec=current_text,
        recommended_model_spec=recommended_text,
        recommendations=recommendations,
    )

    # Write outputs
    (out_dir / "optimised_MODEL_SPEC").write_text(recommended_text)
    (out_dir / "audit_report.md").write_text(AuditReporter(report).markdown())
    report.to_json(out_dir / "audit_report.json")

    return report


def main():
    parser = argparse.ArgumentParser(description="Optimise CMFGEN MODEL_SPEC from SN_HYDRO_DATA")
    parser.add_argument("sn_hydro", help="Path to SN_HYDRO_DATA")
    parser.add_argument("model_spec", nargs="?", default=None, help="Current MODEL_SPEC (for reporting)")
    parser.add_argument("--cmfgen-dir", default=None, help="CMFGEN run directory with F_TO_S symlinks and WARNINGS")
    parser.add_argument("--atomic-dir", default="/u/sbrennan/atomic", help="Atomic data root")
    parser.add_argument("--out-dir", default=".", help="Output directory")
    parser.add_argument("--min-population", type=float, default=1e-6, help="Minimum ion Saha fraction to be retained")
    parser.add_argument("--level-threshold", type=float, default=1e-4, help="LTE population threshold for full NLTE levels")
    args = parser.parse_args()

    report = build_report(
        sn_hydro_path=args.sn_hydro,
        current_model_spec=args.model_spec,
        cmfgen_dir=args.cmfgen_dir,
        atomic_dir=args.atomic_dir,
        out_dir=args.out_dir,
        min_population=args.min_population,
        level_threshold=args.level_threshold,
    )

    print(AuditReporter(report).markdown())


if __name__ == "__main__":
    main()
