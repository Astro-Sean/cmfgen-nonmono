#!/usr/bin/env python3
"""Robust CMFGEN convergence diagnostics.

CMFGEN convergence is NOT necessarily monotonic iteration-by-iteration.
A monotonic solver can show temporary increases or oscillations while
still converging. This module distinguishes:

  CONVERGING          — corrections generally decrease, occasional increases
  DAMPED_OSCILLATION  — oscillates but amplitude decreases
  PERSISTENT_OSCILLATION — oscillates without amplitude decrease
  STAGNATING          — corrections plateau above tolerance
  DIVERGING           — corrections increase over the window
  CONVERGED           — corrections below tolerance and stable

The diagnostic uses a configurable rolling window (default 10 iterations)
and examines:
  - absolute magnitude of corrections
  - median correction
  - maximum correction
  - recent minimum / maximum
  - overall slope/trend (linear fit in log space)
  - oscillation amplitude (peak-to-peak)
  - whether approaching a stable floor

CMFGEN reports (per iteration):
  Maximum % increase at depth D is V (LAMBDA/BA) --- iteration N
  Maximum % decrease at depth D is V (LAMBDA/BA) --- iteration N
  Maximm changes as returned by SOLVEBA_V13 is V
  Luminosity of star (d=1,ND)(iteration N) is: L_inner L_outer
  DEC_VEC / INC_VEC (top-10 corrections per depth bucket)
  CORRECTION_SUM (per-depth histogram of correction magnitudes)

Key CMFGEN semantics (solveba_v13.f):
  - STEQ is +ve for population decreases, -ve for increases
  - "Maximum % increase" is the unscaled correction residual (NOT the applied update)
  - LAM_SCALE_OPT=LIMIT clips decreases >1.1 to 0.999
  - MAX_LAM (CHANGE_LIM) limits applied increases to ~factor of CHANGE_LIM
  - SCALE_OPT=MAJOR excludes trace populations from the scale factor
  - The reported "Maximm changes" is max(|increase|, converted_decrease)
  - DEC_VEC/INC_VEC are the top-NV largest corrections across all depths

A low-abundance species may have huge fractional corrections because its
absolute population is negligible. The diagnostic must not let numerical
noise from trace species dominate the global convergence metric.
"""

from __future__ import annotations

import re
import math
import json
import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict, Any

import numpy as np


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

ITER_RE = re.compile(r'Current great iteration count is\s+(\d+)')
LUM_RE = re.compile(
    r'Luminosity of star \(d=1,ND\)\(iteration\s+(\d+)\)\s+is:\s+'
    r'([\d.eE+\-]+)\s+([\d.eE+\-]+)'
)
MAX_INC_RE = re.compile(
    r'Maximum % increase at depth\s+(\d+)\s+is\s+([\d.eE+\-]+)'
)
MAX_DEC_RE = re.compile(
    r'Maximum % decrease at depth\s+(\d+)\s+is\s+([\d.eE+\-]+)'
)
MAXCH_RE = re.compile(
    r'Maximm changes as returned by SOLVEBA_V13 is\s+([\d.eE+\-]+)'
)
SCALE_RE = re.compile(
    r'minimum value of scale for Major species is:\s+([\d.eE+\-]+)'
)
NG_RE = re.compile(r'Ng acceleration (used|requested|performed)', re.IGNORECASE)
CRASH_RE = re.compile(
    r'forrtl: severe|segmentation fault|SIGSEGV|NaN|NaN detected|'
    r'Error in MON_INTERP|Mean intensity blowing up|'
    r'bad initial population guesses|Predicted changes are too large|'
    r'Error in SETPOP_WLD_V2|zero RGU'
)
COMPLETE_RE = re.compile(
    r'CMFGEN completed|Normal termination|END OF PROGRAM'
)
# DEC_VEC / INC_VEC lines (top-10 per depth bucket)
DEC_VEC_RE = re.compile(r'DEC_VEC:\s+(.+)')
INC_VEC_RE = re.compile(r'INC_VEC:\s+(.+)')
# CORRECTION_SUM header
NT_RE = re.compile(r'NT=\s+(\d+)')
# BA_RES_ACC rejection markers.  A masked/vetoed BA step zeroes every
# depth correction, so the reported "Maximum % increase/decrease" is
# exactly 0.0 for that iteration — a rejected step, NOT convergence.
# ("BA_RES_ACC PARTIAL REJ it  N" header, and "ITER_SUM: it=N ... ba=MASKED".)
BA_REJ_RE = re.compile(r'BA_RES_ACC\s+\S+\s+REJ\s+it\s+(\d+)')
ITERSUM_RE = re.compile(r'ITER_SUM:\s+it=\s*(\d+)\s+.*?ba=(\S+)')
# Reported corrections at/below this floor mean the solver applied zero
# update (a vetoed/masked step); they must not count as convergence.
ZERO_CORR_TOL = 1.0e-12


@dataclass
class IterationRecord:
    """One CMFGEN iteration's diagnostics."""
    iteration: int
    max_inc: float
    max_inc_depth: int
    max_dec: float
    max_dec_depth: int
    maxch: float
    lum_inner: float
    lum_outer: float
    major_scale: Optional[float] = None
    ng_used: bool = False
    rejected: bool = False  # BA_RES_ACC masked/vetoed step (ba=MASKED)


@dataclass
class ConvergenceHistory:
    """Full parsed convergence history from OUTGEN."""
    iterations: List[IterationRecord] = field(default_factory=list)
    crashed: bool = False
    completed: bool = False
    nt: Optional[int] = None  # total levels (from CORRECTION_SUM)

    @property
    def n_iterations(self) -> int:
        return len(self.iterations)

    @property
    def max_corrections(self) -> np.ndarray:
        """Array of MAXCH (maximum change) per iteration."""
        return np.array([it.maxch for it in self.iterations])

    @property
    def max_increases(self) -> np.ndarray:
        return np.array([it.max_inc for it in self.iterations])

    @property
    def max_decreases(self) -> np.ndarray:
        return np.array([it.max_dec for it in self.iterations])

    @property
    def lum_inner(self) -> np.ndarray:
        return np.array([it.lum_inner for it in self.iterations])

    @property
    def lum_outer(self) -> np.ndarray:
        return np.array([it.lum_outer for it in self.iterations])

    @property
    def major_scales(self) -> np.ndarray:
        vals = [it.major_scale for it in self.iterations if it.major_scale is not None]
        return np.array(vals) if vals else np.array([])


@dataclass
class CorrectionSumRecord:
    """Per-depth correction histogram from CMFGEN's CORRECTION_SUM file.

    Each depth has cumulative counts of variables (levels/populations)
    exceeding each correction threshold.  A depth with many variables
    at the 10% or 100% level is a problematic depth that may need
    additional grid points for the transfer solution to converge.
    """
    depth: int                          # CMFGEN depth index (1=outer, ND=inner)
    n_vars_100pct: int                  # count with >100% correction
    n_vars_10pct: int                   # count with >10% correction
    n_vars_1pct: int                    # count with >1% correction
    n_vars_0_1pct: int                  # count with >0.1% correction
    n_vars_0_01pct: int                 # count with >0.01% correction
    n_vars_0_001pct: int                # count with >0.001% correction
    n_vars_0_0001pct: int               # count with >0.0001% correction
    nt: int = 0                         # total variables

    @property
    def frac_10pct(self) -> float:
        """Fraction of variables with >10% correction."""
        return self.n_vars_10pct / max(self.nt, 1)

    @property
    def frac_1pct(self) -> float:
        """Fraction of variables with >1% correction."""
        return self.n_vars_1pct / max(self.nt, 1)

    @property
    def severity(self) -> float:
        """Severity score: weighted sum of correction fractions.

        Higher = more problematic.  Weights emphasise large corrections
        (>10%, >1%) over small ones (>0.1%, >0.01%).
        """
        return (10.0 * self.frac_10pct +
                1.0 * self.frac_1pct +
                0.1 * (self.n_vars_0_1pct / max(self.nt, 1)))


def parse_correction_sum(correction_sum_path: Path) -> List[CorrectionSumRecord]:
    """Parse a CMFGEN CORRECTION_SUM file into per-depth records.

    The CORRECTION_SUM file has the format:

        Summary of changes at each depth
        NT=   115

          Depth  100.0%   10.0%    1.0%    0.1%   0.01%  0.001% 0.0001%
              1      16      72      85     109     113     113     113
              2      12      47      86     108     112     113     113
              ...

    The columns are cumulative counts of variables with corrections
    exceeding each threshold.  NT is the total number of variables.

    Parameters
    ----------
    correction_sum_path : Path
        Path to the CORRECTION_SUM file.

    Returns
    -------
    List[CorrectionSumRecord], one per depth (outer->inner, depth 1=outer).
    """
    if not correction_sum_path.exists():
        return []

    lines = correction_sum_path.read_text(errors='replace').splitlines()
    records = []
    nt = 0

    # Find NT and the header line
    in_data = False
    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Parse NT
        nt_m = NT_RE.search(line)
        if nt_m and not in_data:
            nt = int(nt_m.group(1))
            continue

        # Look for the header line with "Depth" and percentage columns
        if 'Depth' in line and '100.0%' in line:
            in_data = True
            continue

        if not in_data:
            continue

        # Parse data lines: depth, then 7 cumulative counts
        parts = line.split()
        if len(parts) >= 8:
            try:
                depth = int(parts[0])
                counts = [int(x) for x in parts[1:8]]
            except ValueError:
                continue
            records.append(CorrectionSumRecord(
                depth=depth,
                n_vars_100pct=counts[0],
                n_vars_10pct=counts[1],
                n_vars_1pct=counts[2],
                n_vars_0_1pct=counts[3],
                n_vars_0_01pct=counts[4],
                n_vars_0_001pct=counts[5],
                n_vars_0_0001pct=counts[6],
                nt=nt,
            ))

    return records


def identify_problematic_depths(records: List[CorrectionSumRecord],
                                frac_10pct_threshold: float = 0.1,
                                frac_1pct_threshold: float = 0.5,
                                min_severity: float = 0.5) -> List[dict]:
    """Identify depths with poor convergence from CORRECTION_SUM records.

    A depth is flagged as problematic if:
    - >10% of variables have corrections >10% (frac_10pct_threshold)
    - OR >50% of variables have corrections >1% (frac_1pct_threshold)
    - OR severity score > 0.5 (min_severity)

    Parameters
    ----------
    records : List[CorrectionSumRecord]
        Parsed CORRECTION_SUM records.
    frac_10pct_threshold : float
        Flag depths where >this fraction of variables have >10% corrections.
    frac_1pct_threshold : float
        Flag depths where >this fraction of variables have >1% corrections.
    min_severity : float
        Flag depths with severity above this value.

    Returns
    -------
    List of dicts with keys: 'depth', 'severity', 'frac_10pct', 'frac_1pct',
    'reason'.  Sorted by severity (most problematic first).
    """
    problematic = []
    for rec in records:
        reasons = []
        if rec.frac_10pct > frac_10pct_threshold:
            reasons.append(f"frac_10pct={rec.frac_10pct:.2f}")
        if rec.frac_1pct > frac_1pct_threshold:
            reasons.append(f"frac_1pct={rec.frac_1pct:.2f}")
        if rec.severity > min_severity:
            reasons.append(f"severity={rec.severity:.2f}")

        if reasons:
            problematic.append({
                'depth': rec.depth,
                'severity': rec.severity,
                'frac_10pct': rec.frac_10pct,
                'frac_1pct': rec.frac_1pct,
                'n_vars_10pct': rec.n_vars_10pct,
                'n_vars_1pct': rec.n_vars_1pct,
                'nt': rec.nt,
                'reason': '; '.join(reasons),
            })

    problematic.sort(key=lambda x: x['severity'], reverse=True)
    return problematic


def map_depths_to_radii(problematic_depths: List[dict],
                         rdinr_path: Path) -> List[dict]:
    """Map CMFGEN depth indices to physical radii using the RDINR file.

    CMFGEN depth 1 = outermost, depth ND = innermost.
    RDINR format: header line, then R(1) R(ND) ND, then per-depth lines
    with R (1e10 cm), V (km/s), and other columns.

    Parameters
    ----------
    problematic_depths : List[dict]
        Output from ``identify_problematic_depths``.
    rdinr_path : Path
        Path to the RDINR file.

    Returns
    -------
    The same list with added keys 'r_cm' and 'v_kms' for each entry.
    """
    if not rdinr_path.exists():
        return problematic_depths

    lines = rdinr_path.read_text(errors='replace').splitlines()
    # Parse radii: skip header (1 line), then R(1) R(ND) ND line, then per-depth
    radii = []
    velocities = []
    for i, line in enumerate(lines):
        if i < 2:
            continue
        parts = line.split()
        if len(parts) >= 6:
            try:
                r = float(parts[0])  # 1e10 cm
                v = float(parts[5])  # km/s
                radii.append(r * 1e10)  # convert to cm
                velocities.append(v)
            except (ValueError, IndexError):
                continue

    if not radii:
        return problematic_depths

    # RDINR is written inner->outer (R increasing), but CMFGEN depths are
    # outer->inner (depth 1 = outer).  Check ordering.
    if radii[0] < radii[-1]:
        # inner->outer: reverse to get outer->inner
        radii = radii[::-1]
        velocities = velocities[::-1]

    for p in problematic_depths:
        depth_idx = p['depth'] - 1  # 0-based
        if 0 <= depth_idx < len(radii):
            p['r_cm'] = radii[depth_idx]
            p['v_kms'] = velocities[depth_idx]
        else:
            p['r_cm'] = None
            p['v_kms'] = None

    return problematic_depths


def group_problematic_regions(problematic_depths: List[dict],
                              max_gap: int = 3) -> List[dict]:
    """Group contiguous problematic depths into refinement regions.

    Parameters
    ----------
    problematic_depths : List[dict]
        Output from ``identify_problematic_depths`` (with radii mapped).
    max_gap : int
        Maximum gap (in depth indices) to consider two problematic depths
        as part of the same region.

    Returns
    -------
    List of dicts with keys: 'depth_start', 'depth_end', 'r_outer', 'r_inner',
    'n_depths', 'max_severity', 'mean_severity'.
    """
    if not problematic_depths:
        return []

    depths = sorted([p['depth'] for p in problematic_depths])
    regions = []
    start = depths[0]
    prev = depths[0]
    for d in depths[1:]:
        if d - prev > max_gap:
            regions.append((start, prev))
            start = d
        prev = d
    regions.append((start, prev))

    # Build region dicts with radii
    depth_map = {p['depth']: p for p in problematic_depths}
    result = []
    for r_start, r_end in regions:
        depths_in_region = [depth_map[d] for d in range(r_start, r_end + 1)
                           if d in depth_map]
        if not depths_in_region:
            continue
        severities = [p['severity'] for p in depths_in_region]
        r_outer = min((p.get('r_cm', float('inf')) for p in depths_in_region),
                      default=None)
        r_inner = max((p.get('r_cm', 0) for p in depths_in_region),
                      default=None)
        result.append({
            'depth_start': r_start,
            'depth_end': r_end,
            'n_depths': r_end - r_start + 1,
            'r_outer_cm': r_outer,
            'r_inner_cm': r_inner,
            'max_severity': max(severities),
            'mean_severity': sum(severities) / len(severities),
        })

    return result


def parse_outgen(outgen_path: Path) -> ConvergenceHistory:
    """Parse an OUTGEN file into a ConvergenceHistory.

    Reads the full file (not just the tail) so the complete iteration
    history is available for trend analysis.
    """
    history = ConvergenceHistory()
    if not outgen_path.exists():
        return history

    text = outgen_path.read_text(errors='replace')

    # Check for crash / completion
    history.crashed = bool(CRASH_RE.search(text))
    history.completed = bool(COMPLETE_RE.search(text))

    # NT from CORRECTION_SUM
    nt_m = NT_RE.search(text)
    if nt_m:
        history.nt = int(nt_m.group(1))

    # Parse iteration records
    # We match Luminosity lines (which have iteration numbers) and then
    # find the corresponding max inc/dec/maxch for that iteration.
    lum_matches = list(LUM_RE.finditer(text))
    inc_matches = list(MAX_INC_RE.finditer(text))
    dec_matches = list(MAX_DEC_RE.finditer(text))
    maxch_matches = list(MAXCH_RE.finditer(text))
    scale_matches = list(SCALE_RE.finditer(text))
    ng_matches = list(NG_RE.finditer(text))
    rej_matches = list(BA_REJ_RE.finditer(text))
    itersum_matches = list(ITERSUM_RE.finditer(text))

    # Build iteration -> data mapping
    # Luminosity lines have the iteration number explicitly.
    # Max inc/dec/maxch lines also have iteration numbers at the end.
    # We parse by iteration number.
    iter_data: Dict[int, Dict[str, Any]] = {}

    for m in lum_matches:
        it = int(m.group(1))
        if it not in iter_data:
            iter_data[it] = {}
        iter_data[it]['lum_inner'] = float(m.group(2))
        iter_data[it]['lum_outer'] = float(m.group(3))

    for m in inc_matches:
        # The line ends with "--- iteration N"
        # We need to find the iteration number from context
        # The regex doesn't capture it, so we search backwards for the
        # iteration number on the same line
        line_end = text.find('--- iteration', m.start())
        if line_end > 0:
            it_match = re.search(r'iteration\s+(\d+)', text[line_end:line_end+30])
            if it_match:
                it = int(it_match.group(1))
                if it not in iter_data:
                    iter_data[it] = {}
                iter_data[it]['max_inc'] = float(m.group(2))
                iter_data[it]['max_inc_depth'] = int(m.group(1))

    for m in dec_matches:
        line_end = text.find('--- iteration', m.start())
        if line_end > 0:
            it_match = re.search(r'iteration\s+(\d+)', text[line_end:line_end+30])
            if it_match:
                it = int(it_match.group(1))
                if it not in iter_data:
                    iter_data[it] = {}
                iter_data[it]['max_dec'] = float(m.group(2))
                iter_data[it]['max_dec_depth'] = int(m.group(1))

    for m in maxch_matches:
        # Maximm changes line doesn't have an explicit iteration number,
        # but it appears right after the max inc/dec lines for that iteration.
        # We find the nearest preceding "iteration N" reference.
        search_start = max(0, m.start() - 500)
        it_match = None
        for im in re.finditer(r'iteration\s+(\d+)', text[search_start:m.start()]):
            it_match = im
        if it_match:
            it = int(it_match.group(1))
            if it not in iter_data:
                iter_data[it] = {}
            iter_data[it]['maxch'] = float(m.group(1))

    for m in scale_matches:
        search_start = max(0, m.start() - 500)
        it_match = None
        for im in re.finditer(r'iteration\s+(\d+)', text[search_start:m.start()]):
            it_match = im
        if it_match:
            it = int(it_match.group(1))
            if it not in iter_data:
                iter_data[it] = {}
            iter_data[it]['major_scale'] = float(m.group(1))

    for m in ng_matches:
        search_start = max(0, m.start() - 200)
        it_match = None
        for im in re.finditer(r'iteration\s+(\d+)', text[search_start:m.start()]):
            it_match = im
        if it_match:
            it = int(it_match.group(1))
            if it not in iter_data:
                iter_data[it] = {}
            iter_data[it]['ng_used'] = True

    # Rejection markers carry an explicit iteration number.
    for m in rej_matches:
        it = int(m.group(1))
        if it not in iter_data:
            iter_data[it] = {}
        iter_data[it]['rejected'] = True

    for m in itersum_matches:
        if m.group(2).upper() == 'MASKED':
            it = int(m.group(1))
            if it not in iter_data:
                iter_data[it] = {}
            iter_data[it]['rejected'] = True

    # Build sorted iteration list
    for it_num in sorted(iter_data.keys()):
        d = iter_data[it_num]
        rec = IterationRecord(
            iteration=it_num,
            max_inc=d.get('max_inc', 0.0),
            max_inc_depth=d.get('max_inc_depth', 0),
            max_dec=d.get('max_dec', 0.0),
            max_dec_depth=d.get('max_dec_depth', 0),
            maxch=d.get('maxch', max(d.get('max_inc', 0.0), d.get('max_dec', 0.0))),
            lum_inner=d.get('lum_inner', 0.0),
            lum_outer=d.get('lum_outer', 0.0),
            major_scale=d.get('major_scale'),
            ng_used=d.get('ng_used', False),
            rejected=d.get('rejected', False),
        )
        history.iterations.append(rec)

    return history


# ---------------------------------------------------------------------------
# Convergence classification
# ---------------------------------------------------------------------------

@dataclass
class ConvergenceAssessment:
    """Multi-iteration convergence assessment."""
    classification: str  # CONVERGING, DAMPED_OSCILLATION, etc.
    window_size: int
    n_total_iterations: int

    # Correction metrics (in the rolling window)
    max_correction: float
    median_correction: float
    min_correction: float
    recent_min: float
    recent_max: float

    # Trend
    slope_log10: float  # slope of log10(correction) vs iteration
    trend: str  # "decreasing", "increasing", "flat"

    # Oscillation
    oscillation_amplitude: float  # log10(max/min) in window
    oscillation_indicator: str  # "none", "damped", "persistent"

    # Stagnation
    stagnation_indicator: str  # "none", "stagnating"

    # Luminosity
    lum_inner: float
    lum_outer: float
    lum_outer_negative: bool
    lum_trend: str  # "stable", "improving", "worsening", "unstable"

    # Worst depth/species
    worst_depth: int
    worst_species: str  # "unknown" (CMFGEN doesn't report species in the summary)

    # Physical validity
    negative_luminosity: bool
    crashed: bool
    completed: bool

    # Detailed message
    message: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            'classification': self.classification,
            'window_size': self.window_size,
            'n_total_iterations': self.n_total_iterations,
            'max_correction': self.max_correction,
            'median_correction': self.median_correction,
            'min_correction': self.min_correction,
            'recent_min': self.recent_min,
            'recent_max': self.recent_max,
            'slope_log10': self.slope_log10,
            'trend': self.trend,
            'oscillation_amplitude': self.oscillation_amplitude,
            'oscillation_indicator': self.oscillation_indicator,
            'stagnation_indicator': self.stagnation_indicator,
            'lum_inner': self.lum_inner,
            'lum_outer': self.lum_outer,
            'lum_outer_negative': self.lum_outer_negative,
            'lum_trend': self.lum_trend,
            'worst_depth': self.worst_depth,
            'worst_species': self.worst_species,
            'negative_luminosity': self.negative_luminosity,
            'crashed': self.crashed,
            'completed': self.completed,
            'message': self.message,
        }


def assess_convergence(
    history: ConvergenceHistory,
    window: int = 10,
    converged_tol: float = 10.0,
    stagnation_tol: float = 0.5,
    divergence_slope: float = 0.15,
    oscillation_ratio: float = 2.0,
) -> ConvergenceAssessment:
    """Assess convergence using a rolling window.

    Parameters
    ----------
    history : ConvergenceHistory
        Parsed OUTGEN history.
    window : int
        Rolling window size (number of recent iterations to analyze).
    converged_tol : float
        Correction below this (in %) is considered converged.
    stagnation_tol : float
        Relative change below this (in log10 space) is considered stagnant.
    divergence_slope : float
        Slope of log10(correction) vs iteration above which is diverging.
    oscillation_ratio : float
        max/min ratio (in linear space) above which oscillation is significant.

    Returns
    -------
    ConvergenceAssessment
    """
    n = history.n_iterations
    corrections = history.max_corrections

    # Default values
    cls = "STARTING"
    msg = ""
    slope = 0.0
    trend = "unknown"
    osc_amp = 0.0
    osc_ind = "none"
    stag_ind = "none"
    lum_trend = "unknown"

    worst_depth = 0
    if history.iterations:
        worst_depth = history.iterations[-1].max_inc_depth

    lum_inner = history.iterations[-1].lum_inner if history.iterations else 0.0
    lum_outer = history.iterations[-1].lum_outer if history.iterations else 0.0
    lum_outer_negative = lum_outer < 0

    # Crash and completion take priority over iteration count
    if history.crashed:
        cls = "CRASHED"
        msg = "Fatal error detected in OUTGEN (NaN, segfault, etc.)."
        return ConvergenceAssessment(
            classification=cls, window_size=0, n_total_iterations=n,
            max_correction=0.0, median_correction=0.0, min_correction=0.0,
            recent_min=0.0, recent_max=0.0, slope_log10=0.0, trend="unknown",
            oscillation_amplitude=0.0, oscillation_indicator="none",
            stagnation_indicator="none", lum_inner=lum_inner,
            lum_outer=lum_outer, lum_outer_negative=lum_outer_negative,
            lum_trend="unknown", worst_depth=worst_depth,
            worst_species="unknown", negative_luminosity=lum_outer_negative,
            crashed=True, completed=history.completed, message=msg,
        )

    if n == 0:
        cls = "STARTING"
        msg = "No iterations parsed yet."
    elif n < 3:
        cls = "STARTING"
        msg = f"Only {n} iterations — insufficient data for trend analysis."
    else:
        # Use the rolling window
        w = min(window, n)
        corr_w = corrections[-w:]
        log_corr = np.log10(np.maximum(corr_w, 1e-30))

        max_corr = float(np.max(corr_w))
        med_corr = float(np.median(corr_w))
        min_corr = float(np.min(corr_w))
        recent_min = float(np.min(corr_w[-3:])) if w >= 3 else min_corr
        recent_max = float(np.max(corr_w[-3:])) if w >= 3 else max_corr

        # Trend: linear fit of log10(correction) vs iteration index
        x = np.arange(w, dtype=float)
        if w > 1:
            slope = float(np.polyfit(x, log_corr, 1)[0])
        else:
            slope = 0.0

        if slope < -divergence_slope:
            trend = "decreasing"
        elif slope > divergence_slope:
            trend = "increasing"
        else:
            trend = "flat"

        # Also check first-half vs second-half median for a more robust
        # trend direction that doesn't depend on the slope threshold.
        first_half = corr_w[:w // 2]
        second_half = corr_w[w // 2:]
        first_half_med = float(np.median(first_half))
        second_half_med = float(np.median(second_half))
        # "Generally decreasing" means the median of the second half is
        # smaller than the first half by at least a factor of 1.5.
        generally_decreasing = (
            second_half_med < first_half_med / 1.5
        )
        generally_increasing = (
            second_half_med > first_half_med * 1.5
        )

        # Oscillation: peak-to-peak in the window
        osc_amp = float(log_corr.max() - log_corr.min())

        # Check for oscillation pattern (alternating up/down)
        # Only classify as oscillation if the amplitude is significant
        # relative to the overall magnitude AND there are sign changes
        osc_log_threshold = math.log10(oscillation_ratio)
        if w >= 4 and osc_amp > osc_log_threshold:
            diffs = np.diff(log_corr)
            n_sign_changes = int(np.sum(np.diff(np.sign(diffs)) != 0))
            sign_change_ratio = n_sign_changes / max(len(diffs) - 1, 1)

            if sign_change_ratio > 0.3:
                # Compare first-half vs second-half amplitude AND median
                first_half_amp = float(
                    log_corr[:w // 2].max() - log_corr[:w // 2].min()
                )
                second_half_amp = float(
                    log_corr[w // 2:].max() - log_corr[w // 2:].min()
                )

                # Damped: overall level is decreasing (median drops) AND
                # amplitude is not growing.
                if generally_decreasing and second_half_amp <= first_half_amp * 1.2:
                    osc_ind = "damped"
                # Persistent: overall level is NOT decreasing AND
                # amplitude is NOT decreasing.
                elif not generally_decreasing and second_half_amp >= first_half_amp * 0.8:
                    osc_ind = "persistent"
                # If level is decreasing but amplitude is growing, still damped
                # (the overall trend matters more than local amplitude).
                elif generally_decreasing:
                    osc_ind = "damped"
                else:
                    osc_ind = "persistent"
        elif osc_amp > osc_log_threshold:
            osc_ind = "present"

        # Stagnation: corrections plateau above tolerance, no oscillation,
        # no downward trend.
        if max_corr > converged_tol and osc_ind == "none":
            if osc_amp < stagnation_tol and abs(slope) < stagnation_tol * 0.5:
                stag_ind = "stagnating"

        # Luminosity trend
        if n >= 3:
            lum_w = history.lum_outer[-min(window, n):]
            if np.all(np.abs(lum_w) < 1e-30):
                lum_trend = "zero"
            else:
                lum_rel = np.abs(np.diff(lum_w)) / np.maximum(np.abs(lum_w[:-1]), 1e-30)
                lum_rel_recent = float(np.median(lum_rel[-3:])) if len(lum_rel) >= 3 else float(np.median(lum_rel))
                if lum_rel_recent < 0.01:
                    lum_trend = "stable"
                elif lum_rel_recent < 0.1:
                    lum_trend = "improving"
                else:
                    lum_trend = "unstable"

        # Classification (priority order):
        # 1. Completed + recent convergence → CONVERGED
        # 2. Completed + not converged → STAGNATING
        # 3. Below tolerance (recent) → CONVERGED
        # 4. Diverging (increasing slope) → DIVERGING
        # 5. Persistent oscillation (not damped, not decreasing) → PERSISTENT_OSCILLATION
        # 6. Stagnation → STAGNATING
        # 7. Damped oscillation → DAMPED_OSCILLATION
        # 8. Decreasing trend → CONVERGING
        # 9. Otherwise → CONVERGING (insufficient trend)

        # For "completed" and "below tolerance" checks, use the recent
        # minimum (last few iterations), not the window max. A model that
        # started with 1e20 corrections but ended at 2% has converged.
        #
        # False-convergence veto (m15_2p82d_jfnk_hyb, job 12073582): a
        # reported max-% change of ~0 does NOT mean converged — it means
        # the solver applied zero correction because the step was vetoed
        # (BA_RES_ACC rejected the step and masked every depth
        # correction).  A genuinely converged iteration still applies a
        # nonzero correction (or the run converged earlier).  Only
        # nonzero, non-rejected iterations count as below-tolerance
        # evidence.
        recent_recs = history.iterations[-3:]
        real_recent = [r.maxch for r in recent_recs
                       if not r.rejected and r.maxch > ZERO_CORR_TOL]
        recent_below_tol = (bool(real_recent)
                            and min(real_recent) < converged_tol)

        if history.completed and recent_below_tol:
            cls = "CONVERGED"
            msg = (f"CMFGEN completed with recent max correction {recent_max:.2f}% "
                   f"(< {converged_tol}%).")
        elif history.completed and not recent_below_tol:
            cls = "STAGNATING"
            msg = (f"CMFGEN completed but recent max correction {recent_max:.2f}% "
                   f"exceeds tolerance {converged_tol}%.")
        elif recent_below_tol:
            cls = "CONVERGED"
            msg = f"Max correction {max_corr:.2f}% below tolerance {converged_tol}%."
        elif trend == "increasing" and slope > divergence_slope:
            cls = "DIVERGING"
            msg = (f"Corrections increasing (slope={slope:.3f} dex/iter, "
                   f"max={max_corr:.2e}%). Model is diverging.")
        elif osc_ind == "persistent" and not generally_decreasing:
            cls = "PERSISTENT_OSCILLATION"
            msg = (f"Corrections oscillating without damping "
                   f"(amplitude={osc_amp:.2f} dex, max={max_corr:.2e}%). "
                   f"Investigate species/depth {worst_depth}.")
        elif stag_ind == "stagnating":
            cls = "STAGNATING"
            msg = (f"Corrections plateaued at {max_corr:.2e}% "
                   f"(slope={slope:.3f}, amplitude={osc_amp:.2f} dex). "
                   f"Model is stalled above tolerance.")
        elif osc_ind == "damped":
            cls = "DAMPED_OSCILLATION"
            msg = (f"Corrections oscillating with decreasing amplitude "
                   f"(amp={osc_amp:.2f} dex, slope={slope:.3f}, max={max_corr:.2e}%). "
                   f"Convergence is likely but slow.")
        elif trend == "decreasing" or generally_decreasing:
            cls = "CONVERGING"
            msg = (f"Corrections decreasing (slope={slope:.3f} dex/iter, "
                   f"max={max_corr:.2e}%, median={med_corr:.2e}%). "
                   f"Occasional increases are normal.")
        else:
            cls = "CONVERGING"
            msg = (f"Max correction {max_corr:.2e}%, slope={slope:.3f}. "
                   f"Insufficient trend to classify definitively.")

        # Physical validity warnings
        n_vetoed = sum(1 for r in recent_recs
                       if r.rejected or r.maxch <= ZERO_CORR_TOL)
        if n_vetoed:
            msg += (f" NOTE: {n_vetoed} of the last {len(recent_recs)} "
                    f"iterations applied zero correction (vetoed/masked "
                    f"BA step) — not counted as convergence evidence.")
        if lum_outer_negative and abs(lum_outer) > 1e-3 * max(abs(lum_inner), 1.0):
            msg += f" WARNING: outer luminosity is negative ({lum_outer:.2e})."
        if lum_trend == "unstable":
            msg += f" WARNING: luminosity is unstable."

        return ConvergenceAssessment(
            classification=cls,
            window_size=w,
            n_total_iterations=n,
            max_correction=max_corr,
            median_correction=med_corr,
            min_correction=min_corr,
            recent_min=recent_min,
            recent_max=recent_max,
            slope_log10=slope,
            trend=trend,
            oscillation_amplitude=osc_amp,
            oscillation_indicator=osc_ind,
            stagnation_indicator=stag_ind,
            lum_inner=lum_inner,
            lum_outer=lum_outer,
            lum_outer_negative=lum_outer_negative,
            lum_trend=lum_trend,
            worst_depth=worst_depth,
            worst_species="unknown",
            negative_luminosity=lum_outer_negative,
            crashed=history.crashed,
            completed=history.completed,
            message=msg,
        )

    # Fallback for early iterations
    return ConvergenceAssessment(
        classification=cls,
        window_size=min(window, n) if n > 0 else 0,
        n_total_iterations=n,
        max_correction=float(corrections[-1]) if n > 0 else 0.0,
        median_correction=float(corrections[-1]) if n > 0 else 0.0,
        min_correction=float(corrections[-1]) if n > 0 else 0.0,
        recent_min=float(corrections[-1]) if n > 0 else 0.0,
        recent_max=float(corrections[-1]) if n > 0 else 0.0,
        slope_log10=slope,
        trend=trend,
        oscillation_amplitude=osc_amp,
        oscillation_indicator=osc_ind,
        stagnation_indicator=stag_ind,
        lum_inner=lum_inner,
        lum_outer=lum_outer,
        lum_outer_negative=lum_outer_negative,
        lum_trend=lum_trend,
        worst_depth=worst_depth,
        worst_species="unknown",
        negative_luminosity=lum_outer_negative,
        crashed=history.crashed,
        completed=history.completed,
        message=msg,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def format_assessment(assessment: ConvergenceAssessment) -> str:
    """Format a ConvergenceAssessment as a human-readable report."""
    lines = [
        "CMFGEN CONVERGENCE STATUS",
        "=" * 60,
        f"Total iterations:    {assessment.n_total_iterations}",
        f"Rolling window:      {assessment.window_size}",
        "",
        f"Maximum correction:  {assessment.max_correction:.4e}%",
        f"Median correction:   {assessment.median_correction:.4e}%",
        f"Recent min:          {assessment.recent_min:.4e}%",
        f"Recent max:          {assessment.recent_max:.4e}%",
        f"Trend (slope):       {assessment.slope_log10:.4f} dex/iter ({assessment.trend})",
        f"Oscillation:         amplitude={assessment.oscillation_amplitude:.2f} dex "
        f"({assessment.oscillation_indicator})",
        f"Stagnation:          {assessment.stagnation_indicator}",
        "",
        f"Luminosity (inner):  {assessment.lum_inner:.4e}",
        f"Luminosity (outer):  {assessment.lum_outer:.4e} "
        f"({'NEGATIVE' if assessment.lum_outer_negative else 'positive'})",
        f"Luminosity trend:    {assessment.lum_trend}",
        "",
        f"Worst depth:         {assessment.worst_depth}",
        f"Worst species:       {assessment.worst_species}",
        "",
        f"Crashed:             {assessment.crashed}",
        f"Completed:           {assessment.completed}",
        "",
        f"Classification:      {assessment.classification}",
        f"Message:             {assessment.message}",
        "=" * 60,
    ]
    return "\n".join(lines)


def assess_outgen(outgen_path: Path, window: int = 10) -> ConvergenceAssessment:
    """Convenience: parse OUTGEN and assess convergence."""
    history = parse_outgen(outgen_path)
    return assess_convergence(history, window=window)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Assess CMFGEN convergence from OUTGEN."
    )
    parser.add_argument('outgen', type=Path, help='Path to OUTGEN file.')
    parser.add_argument('--window', type=int, default=10,
                        help='Rolling window size (default: 10).')
    parser.add_argument('--json', action='store_true',
                        help='Output JSON instead of text.')
    args = parser.parse_args()

    assessment = assess_outgen(args.outgen, window=args.window)
    if args.json:
        print(json.dumps(assessment.to_dict(), indent=2))
    else:
        print(format_assessment(assessment))


if __name__ == '__main__':
    main()
