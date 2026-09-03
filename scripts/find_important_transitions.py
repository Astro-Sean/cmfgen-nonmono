#!/usr/bin/env python3
"""Find the 100 most important bound-bound transitions per species and suggest/update MODEL_SPEC.

Reads F_OSCDAT files for each species in the CMFGEN model directory, ranks transitions by
oscillator strength (f-value), and determines the smallest ISF (NV,NS,NF) settings that include
all levels involved in the top 100 transitions while remaining consistent with the F_TO_S file.

Usage:
    python find_important_transitions.py <model_dir> [--update] [--n-transitions 100]
"""
import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple


SPECIES_TO_OSCDAT = {
    "HI": "HI_F_OSCDAT",
    "HeI": "HeI_F_OSCDAT",
    "OI": "OI_F_OSCDAT",
    "SkI": "SkI_F_OSCDAT",
    "Sk2": "Sk2_F_OSCDAT",
    "CaI": "CaI_F_OSCDAT",
    "Ca2": "Ca2_F_OSCDAT",
    "FeI": "FeI_F_OSCDAT",
    "Fe2": "Fe2_F_OSCDAT",
    "FeIII": "FeIII_F_OSCDAT",
}

SPECIES_TO_FTOS = {
    "HI": "HI_F_TO_S",
    "HeI": "HeI_F_TO_S",
    "OI": "OI_F_TO_S",
    "SkI": "SkI_F_TO_S",
    "Sk2": "Sk2_F_TO_S",
    "CaI": "CaI_F_TO_S",
    "Ca2": "Ca2_F_TO_S",
    "FeI": "FeI_F_TO_S",
    "Fe2": "Fe2_F_TO_S",
    "FeIII": "FeIII_F_TO_S",
}


def parse_f_oscdat(filepath: Path) -> Tuple[int, List[Tuple[int, int, float, float, float]]]:
    """Parse CMFGEN F_OSCDAT file. Returns (n_levels, transitions) where each transition is
    (lower_index, upper_index, f_value, a_value, wavelength_A)."""
    with open(filepath) as f:
        lines = f.readlines()

    # Skip until the !Date line (level section)
    i = 0
    while i < len(lines) and '!Date' not in lines[i]:
        i += 1
    i += 1  # date line

    n_levels = int(lines[i].split()[0])
    i += 1  # ionization energy
    i += 1  # screened nuclear charge
    i += 1  # number of transitions
    n_transitions = int(lines[i].split()[0])
    i += 1

    # Skip blank/comment lines to level data
    while i < len(lines) and ('!' in lines[i] or lines[i].strip() == ''):
        i += 1

    # Skip level data
    i += n_levels

    # After the level data there is a brief header; the actual transition
    # data follows. We scan the rest of the file and collect the first
    # n_transitions valid transition rows. Some files include an extra
    # 'Trans. #' token; others end with '| Lam(obs) | % Acc |'. The only
    # required columns are the transition label, f, A, wavelength, and the
    # i-j pair, so we locate the i-j from the right end of the token list.

    transitions = []
    for j in range(i, len(lines)):
        if len(transitions) >= n_transitions:
            break
        parts = lines[j].split()
        if len(parts) < 6:
            continue

        found = False
        # Scan from the right for an i-j pair: either a single token 'l-u'
        # (merged) or two tokens 'l-' and 'u' (split). Anything to the right
        # of that (Trans. #, |, Lam(obs), % Acc) is ignored.
        for k in range(len(parts) - 1, 3, -1):
            # Try merged i-j in parts[k]
            m_merged = re.match(r'^(\d+)-\s*(\d+)$', parts[k])
            if m_merged:
                lower = int(m_merged.group(1))
                upper = int(m_merged.group(2))
                f = float(parts[k - 3])
                a = float(parts[k - 2])
                lam = float(parts[k - 1])
                found = True
                break

            # Try split i-j: parts[k] is the upper index (pure int),
            # parts[k - 1] is the lower index ending in '-' (e.g. '1-').
            m_lower = re.match(r'^(\d+)-$', parts[k - 1])
            if m_lower is None:
                continue
            try:
                lower = int(m_lower.group(1))
                upper = int(parts[k])
                a = float(parts[k - 3])
                lam = float(parts[k - 2])
            except ValueError:
                continue
            # f is the next token to the left of the wavelength
            f = float(parts[k - 4])
            found = True
            break

        if not found or lower < 1 or upper < 1:
            continue
        transitions.append((lower, upper, f, a, lam))

    return n_levels, transitions


def parse_f_to_s(filepath: Path) -> Dict[int, int]:
    """Parse F_TO_S file and return mapping {full_level: superlevel}.
    The data lines mirror the F_OSCDAT level format with a trailing 0 and the
    full level index at the end: name ... superlevel 0 full_level."""
    mapping = {}
    with open(filepath) as f:
        for line in f:
            if '!' in line or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                sl = int(parts[-3])
                fl = int(parts[-1])
            except ValueError:
                continue
            mapping[fl] = sl
    return mapping


def find_top_transitions(
    model_dir: Path,
    species: str,
    n: int = 100,
    lam_min: float = 1000.0,
    lam_max: float = 20000.0,
    max_nf: int = 0,
) -> Optional[Dict]:
    """Return analysis for top n transitions of a species.
    If max_nf > 0, only transitions with lower/upper levels <= max_nf are
    considered, so the resulting ISF stays within a practical atom size."""
    oscdat_path = model_dir / SPECIES_TO_OSCDAT.get(species, f"{species}_F_OSCDAT")
    if not oscdat_path.is_file():
        return None

    n_levels, transitions = parse_f_oscdat(oscdat_path)
    if not transitions:
        return None

    # Drop physically suspect / continuum-like transitions and restrict to the
    # requested wavelength window (default: 1000-20000 A, optical/NIR).
    physical = []
    for lower, upper, f, a, lam in transitions:
        if not (0.0 < f < 10.0):
            continue
        if not (0.0 < lam < 1.0e6) or not (lam < float('inf')):
            continue
        if not (lam_min <= lam <= lam_max):
            continue
        if max_nf > 0 and (lower > max_nf or upper > max_nf):
            continue
        physical.append((lower, upper, f, a, lam))

    if not physical:
        return None

    # Rank by f-value (descending)
    ranked = sorted(physical, key=lambda x: x[2], reverse=True)[:n]

    levels_involved = set()
    for lower, upper, f, a, lam in ranked:
        levels_involved.add(lower)
        levels_involved.add(upper)

    max_full_level = max(levels_involved)
    n_unique_levels = len(levels_involved)

    # Determine superlevel mapping from F_TO_S
    ftos_path = model_dir / SPECIES_TO_FTOS.get(species, f"{species}_F_TO_S")
    sl_values = []
    bypass = False
    if ftos_path.is_file():
        mapping = parse_f_to_s(ftos_path)
        for fl in range(1, max_full_level + 1):
            sl = mapping.get(fl, fl)
            sl_values.append(sl)
        if any(sl == 80 for sl in sl_values):
            # LTE marker in range; safer to use 1-to-1 bypass
            bypass = True

    if bypass:
        # Use 1-to-1 mapping for first max_full_level levels
        nv = max_full_level
        ns = max_full_level
        nf = max_full_level
        strategy = "1-to-1 (SL=80 marker in range)"
    else:
        # Use the F_TO_S grouping up to the top levels
        sls_in_range = [mapping.get(i, i) for i in range(1, max_full_level + 1) if i in levels_involved]
        unique_sls = len(set(sls_in_range))
        max_sl = max(mapping.get(i, i) for i in range(1, max_full_level + 1))
        nv = unique_sls
        ns = max_sl
        nf = max_full_level
        strategy = f"F_TO_S grouping (maxSL={max_sl})"

    return {
        "species": species,
        "n_levels": n_levels,
        "n_transitions": len(transitions),
        "n_unique_levels": n_unique_levels,
        "max_full_level": max_full_level,
        "nv": nv,
        "ns": ns,
        "nf": nf,
        "strategy": strategy,
        "top": ranked[:10],
    }


def parse_model_spec(filepath: Path) -> List[str]:
    with open(filepath) as f:
        return f.read().splitlines()


def update_model_spec(lines: List[str], changes: Dict[str, str]) -> List[str]:
    out = []
    for line in lines:
        m = re.match(r'^(\s*)(\d+,\d+,\d+)(\s+\[([A-Za-z0-9_]+)_ISF\])', line)
        if m:
            sp = m.group(4)
            if sp in changes:
                new_isf = changes[sp]
                out.append(f"{m.group(1)}{new_isf}{m.group(3)}{line[m.end():]}")
                continue
        out.append(line)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('model_dir', type=Path, help='CMFGEN model directory')
    p.add_argument('--n-transitions', type=int, default=100)
    p.add_argument('--lam-min', type=float, default=1000.0, help='Minimum wavelength (Angstroms)')
    p.add_argument('--lam-max', type=float, default=20000.0, help='Maximum wavelength (Angstroms)')
    p.add_argument('--max-nf', type=int, default=200, help='Maximum full levels to consider (0=no cap)')
    p.add_argument('--model-spec', type=Path, default=None, help='MODEL_SPEC to read/update')
    p.add_argument('--update', action='store_true', help='Apply updates to MODEL_SPEC')
    p.add_argument('--species', default='', help='Comma-separated species (default: all in MODEL_SPEC)')
    args = p.parse_args()

    model_dir = args.model_dir
    model_spec_path = args.model_spec or (model_dir / 'MODEL_SPEC')

    species_list = [s.strip() for s in args.species.split(',') if s.strip()]
    if not species_list and model_spec_path.is_file():
        text = model_spec_path.read_text()
        species_list = re.findall(r'\[([A-Z][a-zA-Z0-9_]*)_ISF\]', text)

    print(f"Analyzing top {args.n_transitions} transitions in "
          f"{args.lam_min}-{args.lam_max} A for: {', '.join(species_list)}\n")

    all_changes = {}
    for species in species_list:
        result = find_top_transitions(
            model_dir, species, n=args.n_transitions,
            lam_min=args.lam_min, lam_max=args.lam_max,
            max_nf=args.max_nf,
        )
        if result is None:
            print(f"{species}: F_OSCDAT not found, skipping")
            continue

        top_desc = "\n  ".join(
            f"{lower:4d}-{upper:4d} f={f:10.3e} A={a:10.3e} lam={lam:10.3f} A"
            for lower, upper, f, a, lam in result["top"]
        )
        print(
            f"{result['species']}: n_levels={result['n_levels']}, "
            f"top covers {result['n_unique_levels']} unique levels, max_full_level={result['max_full_level']}\n"
            f"  Suggested ISF: {result['nv']},{result['ns']},{result['nf']}  ({result['strategy']})\n"
            f"  Top 10 transitions:\n  {top_desc}\n"
        )
        all_changes[species] = f"{result['nv']},{result['ns']},{result['nf']}"

    if args.update:
        if not model_spec_path.is_file():
            print(f"ERROR: {model_spec_path} not found")
            return 1
        backup = model_spec_path.with_suffix('.bak')
        backup.write_text(model_spec_path.read_text())
        lines = parse_model_spec(model_spec_path)
        new_lines = update_model_spec(lines, all_changes)
        model_spec_path.write_text('\n'.join(new_lines) + '\n')
        print(f"Updated {model_spec_path}; backup at {backup}")
    else:
        print("Dry run: pass --update to write changes to MODEL_SPEC")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
