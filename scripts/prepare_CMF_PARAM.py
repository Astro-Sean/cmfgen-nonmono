#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HPC-side script: enforce correct CMF_FLUX_PARAM settings for the OBS
formal solution (cmf_flux.exe).

This is the OBS analogue of prepare_VADAT.py.  Where prepare_VADAT.py
dynamically enforces non-monotonic solver settings in the main CMFGEN
VADAT, this script enforces the *same* non-monotonic settings in the
OBS CMF_FLUX_PARAM file that cmf_flux.exe reads.

The static template (CMF_FLUX_PARAM_INIT) can drift out of sync with
the main VADAT — exactly the bug that caused the OBS run to use
G_ONLY / SOB / FG_J_CMF_V12 while the main run used N_ON_J / BLANK /
FG_J_CMF_NONMONO_V1.  This script eliminates that drift by enforcing
the correct values at runtime, just as prepare_VADAT.py does.

Non-monotonic consistency (Knop et al. 2009):
  - N_TYPE = N_ON_J          stable at velocity reversals where H=0
                             (G_ONLY uses N/H and crashes when H=0)
  - GLOBAL_LINE = BLANK      blanketing mode; lines handled in continuum
                             frequency set, not via separate Sobolev pass
  - DO_SOB_LINES = F         no Sobolev equivalent-widths
  - SOB_FREQ_IN_OBS = F      no Sobolev frequencies in observed spectrum
  - LAM_SET = F              do not switch to SOB for long wavelengths
  - USE_LAM_ES = T            formal-solution J for e.s. (always T in OBS:
                             the OBS binary has no MOM_J_CMF_NONMONO_V1 —
                             a main run that converged with USE_LAM_ES=F
                             still gets its OBS J from the formal solve)
  - NON_MONO_SOLVER = T       activate FG_J_CMF_NONMONO_V1
  - INC_GRID = F              required by non-monotonic solver
  - FG_OPT = INT/INS          required by non-monotonic solver (not DIFF)
  - NEG_OPAC_OPT = SRCE_CHK   negative opacity handling
  - All TRANS_* = BLANK       blanketing mode for every species
  - USE_FIXED_J = F           incompatible with NON_MONO_SOLVER
  - USE_FRM_REL = F           incompatible with NON_MONO_SOLVER
  - USE_J_REL = F             incompatible with NON_MONO_SOLVER
  - USE_DJDT_RTE = F          incompatible with NON_MONO_SOLVER
  - PP_MOD = F                incompatible with NON_MONO_SOLVER
  - PP_NOV = F                incompatible with NON_MONO_SOLVER
  - TRAP_J = T                trapezoidal weights for J (always)

These forced settings are a subset of prepare_VADAT.py:NON_MONO_FORCED_VADAT,
restricted to keys that cmf_flux.exe actually reads (rd_cmf_flux_controls.f).
Main-run-only keys (VEL_LAW, RD_IN_R_GRID, PURE_HUB, AUTO_SMOOTH,
MAX_CHNG, LAM_SCALE_OPT, etc.) are not applicable to the OBS stage.

Usage:
  # Default: read CMF_FLUX_PARAM_INIT in cwd, enforce, overwrite in-place
  python prepare_CMF_PARAM.py

  # Explicit template and output
  python prepare_CMF_PARAM.py --template /path/to/CMF_FLUX_PARAM_INIT \
                               --output  /path/to/CMF_FLUX_PARAM_INIT

  # Read YAML config for optional OBS overrides (cmf_flux_param_overrides)
  python prepare_CMF_PARAM.py --vadat-overrides config.yaml

The script is called from run_cmfgen_OBS.slurm AFTER the template is
copied from obs_template/ and BEFORE batobs.sh runs.  batobs.sh then
copies CMF_FLUX_PARAM_INIT -> CMF_FLUX_PARAM (with TRAP_J sed, now a
no-op since TRAP_J is already T).

Non-coherent electron scattering: when the converged model ran COH_ES=F
(read from ../VADAT, falling back to the YAML vadat_overrides), the OBS
step automatically gets NUM_ES>=2 (iteration 1 computes J coherently and
writes ES_J_CONV, iteration 2 applies the convolved J) and the ES wing
extents propagated from VADAT.  Explicit cmf_flux_param_overrides still
win.

Frequency range: MIN_CF/MAX_CF are propagated from ../VADAT (expand-only)
because SET_CONT_FREQ_V4 STOPs when the OBS grid does not cover the
model's NU range from RVTJ.
"""

import argparse
import re
import sys
from pathlib import Path

# ==============================================================================
# Forced settings — must match prepare_VADAT.py:NON_MONO_FORCED_VADAT
# for all keys that cmf_flux.exe reads (rd_cmf_flux_controls.f).
# ==============================================================================
NON_MONO_FORCED_CMF_PARAM = {
    # ── Non-monotonic solver core ──────────────────────────────────
    'NON_MONO_SOLVER': 'T',
    'USE_LAM_ES': 'T',
    'N_TYPE': 'N_ON_J',
    'FG_OPT': 'INT/INS',
    'INC_GRID': 'F',

    # ── Line treatment (BLANK = blanketing mode) ───────────────────
    'GLOBAL_LINE': 'BLANK',
    'DO_SOB_LINES': 'F',
    'SOB_FREQ_IN_OBS': 'F',
    'LAM_SET': 'F',

    # ── Incompatible with NON_MONO_SOLVER (all must be F) ───────────
    'USE_FIXED_J': 'F',
    'USE_FRM_REL': 'F',
    'USE_J_REL': 'F',
    'USE_DJDT_RTE': 'F',
    'PP_MOD': 'F',
    'PP_NOV': 'F',

    # ── Numerical / opacity ────────────────────────────────────────
    'NEG_OPAC_OPT': 'SRCE_CHK',
    'CHK_L_POS': 'T',

    # ── Always-on ──────────────────────────────────────────────────
    'TRAP_J': 'T',

    # ── Per-species line treatment: all BLANK ──────────────────────
    'TRANS_HI': 'BLANK',
    'TRANS_HeI': 'BLANK',
    'TRANS_He2': 'BLANK',
    'TRANS_CI': 'BLANK',
    'TRANS_C2': 'BLANK',
    'TRANS_CIV': 'BLANK',
    'TRANS_NI': 'BLANK',
    'TRANS_N2': 'BLANK',
    'TRANS_NIII': 'BLANK',
    'TRANS_NIV': 'BLANK',
    'TRANS_NV': 'BLANK',
    'TRANS_OI': 'BLANK',
    'TRANS_O2': 'BLANK',
    'TRANS_OIII': 'BLANK',
    'TRANS_OIV': 'BLANK',
    'TRANS_OV': 'BLANK',
    'TRANS_OSIX': 'BLANK',
    'TRANS_NeI': 'BLANK',
    'TRANS_Ne2': 'BLANK',
    'TRANS_NeIII': 'BLANK',
    'TRANS_NeIV': 'BLANK',
    'TRANS_NeV': 'BLANK',
    'TRANS_NeSIX': 'BLANK',
    'TRANS_NeSEV': 'BLANK',
    'TRANS_SkI': 'BLANK',
    'TRANS_Sk2': 'BLANK',
    'TRANS_SkIII': 'BLANK',
    'TRANS_SkIV': 'BLANK',
    'TRANS_NaI': 'BLANK',
    'TRANS_MgI': 'BLANK',
    'TRANS_Mg2': 'BLANK',
    'TRANS_AlI': 'BLANK',
    'TRANS_Al2': 'BLANK',
    'TRANS_AlIII': 'BLANK',
    'TRANS_SI': 'BLANK',
    'TRANS_S2': 'BLANK',
    'TRANS_SIII': 'BLANK',
    'TRANS_SIV': 'BLANK',
    'TRANS_SV': 'BLANK',
    'TRANS_SSIX': 'BLANK',
    'TRANS_ArI': 'BLANK',
    'TRANS_Ar2': 'BLANK',
    'TRANS_ArIII': 'BLANK',
    'TRANS_ArIV': 'BLANK',
    'TRANS_ArV': 'BLANK',
    'TRANS_ArSIX': 'BLANK',
    'TRANS_ArSEV': 'BLANK',
    'TRANS_CaI': 'BLANK',
    'TRANS_Ca2': 'BLANK',
    'TRANS_CaIII': 'BLANK',
    'TRANS_CaIV': 'BLANK',
    'TRANS_CaV': 'BLANK',
    'TRANS_CaSIX': 'BLANK',
    'TRANS_CaSEV': 'BLANK',
    'TRANS_CaVIII': 'BLANK',
    'TRANS_FeI': 'BLANK',
    'TRANS_Fe2': 'BLANK',
    'TRANS_FeIII': 'BLANK',
    'TRANS_FeIV': 'BLANK',
    'TRANS_FeV': 'BLANK',
    'TRANS_FeSIX': 'BLANK',
    'TRANS_FeSEV': 'BLANK',
    'TRANS_FeVIII': 'BLANK',
    'TRANS_FeIX': 'BLANK',
    'TRANS_FeX': 'BLANK',
    'TRANS_NkI': 'BLANK',
    'TRANS_Nk2': 'BLANK',
    'TRANS_NkIII': 'BLANK',
    'TRANS_NkIV': 'BLANK',
    'TRANS_NkV': 'BLANK',
    'TRANS_NkSIX': 'BLANK',
    'TRANS_NkSEV': 'BLANK',
    'TRANS_NkVIII': 'BLANK',
    'TRANS_NkIX': 'BLANK',
    'TRANS_NkX': 'BLANK',
}

# Keys that are forced and cannot be overridden by YAML.
# OBS-specific settings (frequency range, velocity spacing, profiles,
# Eddington factors, etc.) remain in the template and are NOT forced.
NON_OVERRIDABLE = set(NON_MONO_FORCED_CMF_PARAM.keys())


def _apply_forced_keys(lines, forced, tag="NON_MONO"):
    """Replace or append CMF_FLUX_PARAM keys.

    Mirrors prepare_VADAT.py:_apply_forced_vadat_keys.  Any TRANS_* key
    not explicitly listed in *forced* is still forced to BLANK.
    """
    key_re = re.compile(r'^\s*([^\s]+)\s+\[([A-Za-z0-9_/]+)\]')
    found = set()
    new_lines = []
    force_trans_blank = any(k.startswith('TRANS_') for k in forced)
    for line in lines:
        m = key_re.match(line)
        if not m:
            new_lines.append(line)
            continue
        key = m.group(2)
        if key in forced:
            found.add(key)
            forced_val = forced[key]
            current = m.group(1).strip()
            if current.upper() != str(forced_val).upper():
                print(f"  [{tag}] {key} forced: {current} -> {forced_val}")
            new_lines.append(f"{forced_val} [{key}]\n")
        elif force_trans_blank and key.startswith('TRANS_'):
            found.add(key)
            current = m.group(1).strip()
            if current.upper() != 'BLANK':
                print(f"  [{tag}] {key} forced: {current} -> BLANK")
            new_lines.append(f"BLANK [{key}]\n")
        else:
            new_lines.append(line)
    # Append any forced keys that were not already present
    for key, val in forced.items():
        if key not in found:
            new_lines.append(f"{val} [{key}]\n")
            print(f"  [{tag}] {key} added: {val}")
    return new_lines


def sob_mode_active(config_path):
    """True when the config requests the dev7 Sobolev non-monotonic path.

    Mirrors prepare_VADAT.py mode resolution: cmfgen.non_mono_mode='SOB'
    (primary), legacy non_mono_line_treatment='SOB', or an explicit
    vadat_overrides DO_NONMONO_SOB=T.  The OBS flag only switches the
    Sobolev EW estimator to SIGMA_TRUE; the spectrum is unchanged.
    """
    try:
        import yaml
    except ImportError:
        return False

    config_path = Path(config_path) if config_path else None
    if config_path is None or not config_path.is_file():
        return False

    with open(config_path) as f:
        cfg = yaml.safe_load(f) or {}

    cmfgen_cfg = cfg.get('cmfgen', {}) or {}
    overrides = cmfgen_cfg.get('vadat_overrides', {}) or {}
    if str(overrides.get('DO_NONMONO_SOB', '')).strip().upper().startswith('T'):
        return True
    mode = str(cmfgen_cfg.get('non_mono_mode', '')).strip().upper()
    if mode:
        return mode == 'SOB'
    lt = str(cmfgen_cfg.get('non_mono_line_treatment', '')).strip().upper()
    return lt == 'SOB'


def knop_mode_active(config_path):
    """True when the config requests the dev7 Knop non-monotonic path.

    Mirrors sob_mode_active: cmfgen.non_mono_mode='KNOP' (primary), or
    an explicit vadat_overrides DO_NONMONO_KNOP/DO_NOMONO_KNOP=T.  The
    OBS flag switches the CMF formal J solve to FG_J_CMF_NONMONO_V1;
    required when sigma < -1, where FG_J_CMF_V12 aborts.
    """
    try:
        import yaml
    except ImportError:
        return False

    config_path = Path(config_path) if config_path else None
    if config_path is None or not config_path.is_file():
        return False

    with open(config_path) as f:
        cfg = yaml.safe_load(f) or {}

    cmfgen_cfg = cfg.get('cmfgen', {}) or {}
    overrides = cmfgen_cfg.get('vadat_overrides', {}) or {}
    for key in ('DO_NONMONO_KNOP', 'DO_NOMONO_KNOP'):
        if str(overrides.get(key, '')).strip().upper().startswith('T'):
            return True
    mode = str(cmfgen_cfg.get('non_mono_mode', '')).strip().upper()
    if mode:
        return mode == 'KNOP'
    lt = str(cmfgen_cfg.get('non_mono_line_treatment', '')).strip().upper()
    return lt == 'KNOP'


def _param_value(lines, key):
    """Return the token (value) for [key] in the param lines, or None."""
    key_re = re.compile(r'^\s*([^\s]+)\s+\[([A-Za-z0-9_/]+)\]')
    for line in lines:
        m = key_re.match(line)
        if m and m.group(2) == key:
            return m.group(1).strip()
    return None


def main_run_es_settings(template_path, config_path=None, keys=None):
    """Return settings of the converged main run.

    ../VADAT relative to the obs directory is the executed truth; the YAML
    cmfgen.vadat_overrides are the fallback when VADAT is unreadable.
    Returns a dict with the keys found (default: COH_ES, ES_WING_EXT,
    R_CMF_WING_EXT) plus '_source'; empty dict when undetermined.
    """
    if keys is None:
        keys = ('COH_ES', 'ES_WING_EXT', 'R_CMF_WING_EXT')
    key_re = re.compile(r'^\s*([^\s]+)\s+\[([A-Za-z0-9_/]+)\]')
    out = {}
    for cand in (Path(template_path).parent.parent / 'VADAT',
                 Path('..') / 'VADAT'):
        if cand.is_file():
            with open(cand) as f:
                for line in f:
                    m = key_re.match(line)
                    if m and m.group(2) in keys and m.group(2) not in out:
                        out[m.group(2)] = m.group(1).strip()
            if out:
                out['_source'] = str(cand)
                return out
    if config_path:
        try:
            import yaml
        except ImportError:
            return out
        config_path = Path(config_path)
        if config_path.is_file():
            with open(config_path) as f:
                cfg = yaml.safe_load(f) or {}
            ov = (cfg.get('cmfgen', {}) or {}).get('vadat_overrides', {}) or {}
            for k in keys:
                if k in ov:
                    out[k] = str(ov[k])
            if out:
                out['_source'] = f'{config_path} (vadat_overrides)'
    return out


def _num_or_none(tok):
    """Parse a Fortran-style numeric token (D exponent) or return None."""
    if tok is None:
        return None
    try:
        return float(str(tok).strip().upper().replace('D', 'E'))
    except ValueError:
        return None


def h_atom_split_n2(template_path):
    """True when the model's H I atom uses split n=2 levels.

    Reads level names from ../HI_F_OSCDAT (linked into the model dir by
    batch.sh assign): the row whose last field is the level index 2 is
    named '2___' for combined atoms and e.g. '2s_2Se' when n=2 is split.
    cmf_flux_sub_v5.f STOPs unless PROF_HI is DOPPLER for split atoms.
    Returns False when the atom file cannot be read or parsed.
    """
    for cand in (Path(template_path).parent.parent / 'HI_F_OSCDAT',
                 Path('..') / 'HI_F_OSCDAT'):
        if not cand.is_file():
            continue
        rows = cand.read_text(errors='replace').splitlines()
        idx = next((i for i, l in enumerate(rows)
                    if 'energy levels' in l.lower()), None)
        if idx is None:
            continue
        for l in rows[idx:idx + 60]:
            parts = l.split()
            if len(parts) >= 6 and parts[-1] == '2':
                return parts[0] != '2___'
    return False


def build_noncoherent_es_companions(lines, es_main):
    """OBS companions when the converged model ran COH_ES=F.

    cmf_flux_sub_v5.f: with COH_ES=F the formal solution iterates the ES
    convolution NUM_ES times (COH_ES=T forces NUM_ES_ITERATIONS=1
    internally).  ES iteration 1 computes J coherently and writes
    ES_J_CONV; iteration 2+ applies the convolved J, so NUM_ES>=2 is
    required for the spectrum to carry any non-coherent ES at all.  The
    wing extents are propagated from VADAT (floored at the values the IIn
    grid needs) so the observed-frame frequency grid covers the same
    scattering wings.  OLD_J is not touched: it reads a literal OLD_J_FILE
    that nothing in the pipeline produces.
    """
    companions = {'COH_ES': 'F'}
    for key, floor in (('ES_WING_EXT', 12000.0), ('R_CMF_WING_EXT', 4.0)):
        vals = [floor]
        for tok in (es_main.get(key), _param_value(lines, key)):
            v = _num_or_none(tok)
            if v is not None:
                vals.append(v)
        companions[key] = f"{max(vals):.1f}"
    cur_es = _num_or_none(_param_value(lines, 'NUM_ES'))
    if cur_es is None or cur_es < 2:
        companions['NUM_ES'] = '2'
    return companions


def load_yaml_overrides(config_path):
    """Load cmf_flux_param_overrides from YAML config (if present).

    Returns a dict of key->value strings.  Non-overridable keys
    (NON_MONO_FORCED_CMF_PARAM) are silently ignored.
    """
    try:
        import yaml
    except ImportError:
        return {}

    config_path = Path(config_path)
    if not config_path.is_file():
        return {}

    with open(config_path) as f:
        cfg = yaml.safe_load(f) or {}

    cmfgen_cfg = cfg.get('cmfgen', {})
    overrides = cmfgen_cfg.get('cmf_flux_param_overrides', {})
    if not isinstance(overrides, dict):
        return {}

    # Strip non-overridable keys
    clean = {}
    for k, v in overrides.items():
        if k in NON_OVERRIDABLE:
            print(f"  WARNING: '{k}' is non-overridable (forced to "
                  f"{NON_MONO_FORCED_CMF_PARAM[k]}); ignoring YAML value {v}")
            continue
        clean[k] = str(v)
    return clean


def apply_overrides(lines, overrides):
    """Apply user-configurable overrides (non-forced keys only)."""
    if not overrides:
        return lines
    key_re = re.compile(r'^\s*([^\s]+)\s+\[([A-Za-z0-9_/]+)\]')
    found = set()
    new_lines = []
    for line in lines:
        m = key_re.match(line)
        if not m:
            new_lines.append(line)
            continue
        key = m.group(2)
        if key in overrides:
            found.add(key)
            new_val = overrides[key]
            current = m.group(1).strip()
            if current != new_val:
                print(f"  [OBS_OVERRIDE] {key}: {current} -> {new_val}")
            new_lines.append(f"{new_val} [{key}]\n")
        else:
            new_lines.append(line)
    for key, val in overrides.items():
        if key not in found:
            new_lines.append(f"{val} [{key}]\n")
            print(f"  [OBS_OVERRIDE] {key} added: {val}")
    return new_lines


def verify_consistency(lines):
    """Verify that forced keys have correct values after enforcement.

    Returns list of error messages (empty if all correct).
    """
    key_re = re.compile(r'^\s*([^\s]+)\s+\[([A-Za-z0-9_/]+)\]')
    values = {}
    for line in lines:
        m = key_re.match(line)
        if m:
            values[m.group(2)] = m.group(1).strip()

    errors = []
    for key, expected in NON_MONO_FORCED_CMF_PARAM.items():
        actual = values.get(key, "")
        if actual.upper() != str(expected).upper():
            errors.append(f"  {key}: expected {expected}, got '{actual}'")

    # Cross-check: NON_MONO_SOLVER=T requires USE_LAM_ES=T in OBS even
    # when the main run used USE_LAM_ES=F.  MOM_J_CMF_NONMONO_V1 exists
    # only in the main binary; the OBS F path would fall back to the
    # clamped-sigma MOM_J_CMF_V8, inconsistent with the non-monotonic
    # model.  J=TC from the formal solution is the correct OBS field.
    if values.get('NON_MONO_SOLVER', '').upper() == 'T':
        if values.get('USE_LAM_ES', '').upper() != 'T':
            errors.append("  NON_MONO_SOLVER=T but USE_LAM_ES != T "
                          "(OBS requires formal-solution J)")

    # Cross-check: NON_MONO_SOLVER=T requires FG_OPT != DIFF
    if values.get('NON_MONO_SOLVER', '').upper() == 'T':
        fg = values.get('FG_OPT', '').upper()
        if fg.startswith('DIFF'):
            errors.append(f"  NON_MONO_SOLVER=T but FG_OPT={fg} (must be INT/INS)")

    return errors


def main():
    parser = argparse.ArgumentParser(
        description="Enforce non-monotonic OBS settings in CMF_FLUX_PARAM")
    parser.add_argument(
        '--template', default='CMF_FLUX_PARAM_INIT',
        help="Input template file (default: CMF_FLUX_PARAM_INIT in cwd)")
    parser.add_argument(
        '--output', default=None,
        help="Output file (default: overwrite --template in-place)")
    parser.add_argument(
        '--vadat-overrides', default=None,
        help="YAML config file for cmf_flux_param_overrides (optional)")
    args = parser.parse_args()

    template_path = Path(args.template)
    if not template_path.is_file():
        print(f"ERROR: template not found: {template_path}", file=sys.stderr)
        return 1

    output_path = Path(args.output) if args.output else template_path

    print(f"prepare_CMF_PARAM: enforcing non-monotonic OBS settings")
    print(f"  template: {template_path}")
    print(f"  output:   {output_path}")

    with open(template_path) as f:
        lines = f.readlines()

    # 1. Apply forced non-monotonic settings (non-overridable)
    lines = _apply_forced_keys(lines, NON_MONO_FORCED_CMF_PARAM)

    # 1b. dev7 Sobolev flag: T only when the main run used DO_NONMONO_SOB.
    # The flag must NOT be stamped unconditionally — flag-off runs must
    # remain bit-identical to the monotonic path.
    if sob_mode_active(args.vadat_overrides):
        lines = _apply_forced_keys(
            lines, {'DO_NONMONO_SOB': 'T'}, tag="NON_MONO_SOB")
        print("  non-mono SOB mode: DO_NONMONO_SOB=T stamped "
              "(EW estimators use SIGMA_TRUE)")

    # dev7 Knop flag: T only when the main run used FG_J_CMF_NONMONO_V1.
    # Must not be stamped unconditionally — monotonic OBS runs must keep
    # the standard FG_J_CMF_V12 path.
    if knop_mode_active(args.vadat_overrides):
        lines = _apply_forced_keys(
            lines, {'DO_NONMONO_KNOP': 'T'}, tag="NON_MONO_KNOP")
        print("  non-mono KNOP mode: DO_NONMONO_KNOP=T stamped "
              "(formal J uses FG_J_CMF_NONMONO_V1)")

    # 1c. Non-coherent ES companions: when the converged model ran
    # COH_ES=F (cmfgen_sub.f -- ES source term uses J convolved with the
    # e.s. kernel), the formal solution needs matching physics.  Applied
    # before user overrides so explicit cmf_flux_param_overrides still win.
    es_main = main_run_es_settings(template_path, args.vadat_overrides)
    if str(es_main.get('COH_ES', '')).strip().upper() == 'F':
        companions = build_noncoherent_es_companions(lines, es_main)
        lines = _apply_forced_keys(lines, companions, tag="ES_NONCOH")
        print(f"  non-coherent ES model ({es_main.get('_source', '?')}): "
              f"OBS companions applied: {companions}")

    # 1d. The OBS continuum grid must contain the model's frequency range —
    # SET_CONT_FREQ_V4 STOPs with "MAX_FREQ too small" when MAX_CF is below
    # the highest NU read from RVTJ.  Expand (never shrink) the template
    # MIN_CF/MAX_CF to the converged run's VADAT values so obs_fin_10 spans
    # the full model range.  Applied before user overrides so explicit
    # cmf_flux_param_overrides still win.
    freq_main = main_run_es_settings(
        template_path, args.vadat_overrides, keys=('MIN_CF', 'MAX_CF'))
    freq_fix = {}
    v = _num_or_none(freq_main.get('MIN_CF'))
    if v is not None:
        cur = _num_or_none(_param_value(lines, 'MIN_CF'))
        if cur is None or v < cur:
            freq_fix['MIN_CF'] = freq_main['MIN_CF']
    v = _num_or_none(freq_main.get('MAX_CF'))
    if v is not None:
        cur = _num_or_none(_param_value(lines, 'MAX_CF'))
        if cur is None or v > cur:
            freq_fix['MAX_CF'] = freq_main['MAX_CF']
    if freq_fix:
        lines = _apply_forced_keys(lines, freq_fix, tag="FREQ_RANGE")
        print(f"  OBS freq range expanded to cover model grid "
              f"({freq_main.get('_source', '?')}): {freq_fix}")

    # 2. Apply user-configurable overrides from YAML (if any)
    if args.vadat_overrides:
        overrides = load_yaml_overrides(args.vadat_overrides)
        if overrides:
            print(f"  Applying {len(overrides)} user overrides from YAML")
            lines = apply_overrides(lines, overrides)

    # 2b. cmf_flux implements only Doppler profiles for split n=2 H atoms
    # (hard STOP in cmf_flux_sub_v5.f), so this wins over user overrides.
    # GLOBAL_PROF=LIST stamps PROF_TYPE='LIST' on every line, which also
    # fails the 'DOP' check — GLOBAL_PROF itself must be DOPPLER.
    if h_atom_split_n2(template_path):
        split_fix = {}
        for key in ('PROF_HI', 'GLOBAL_PROF'):
            cur = _param_value(lines, key)
            if cur is not None and 'DOP' not in cur.strip().upper():
                split_fix[key] = 'DOPPLER'
        if split_fix:
            lines = _apply_forced_keys(lines, split_fix, tag="H_SPLIT")
            print(f"  H atom has split n=2 levels: {split_fix}")

    # 3. Verify consistency
    errors = verify_consistency(lines)
    if errors:
        print("ERROR: consistency check failed:", file=sys.stderr)
        for e in errors:
            print(e, file=sys.stderr)
        return 2

    # 4. Write output
    with open(output_path, 'w') as f:
        f.writelines(lines)

    print(f"  CMF_FLUX_PARAM written: {output_path}")
    print(f"  Non-monotonic OBS settings verified: "
          f"N_ON_J, BLANK, FG_J_CMF_NONMONO_V1")
    return 0


if __name__ == '__main__':
    sys.exit(main())
