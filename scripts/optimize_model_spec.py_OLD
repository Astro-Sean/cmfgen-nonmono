#!/usr/bin/env python3
"""Optimize MODEL_SPEC ISF values based on SN_HYDRO_DATA physical conditions.

Computes LTE level populations at the snapshot's temperature range and determines
which levels have significant population (> threshold). Levels with negligible
LTE population can safely be placed in superlevels without affecting the spectrum.

Two ISF sets are generated:
  - Stage 1 (convergence): Reduced NLTE levels for stable convergence
  - Stage 2 (spectrum): Full NLTE levels, restarted from Stage 1 populations

Usage:
  python optimize_model_spec.py <SN_HYDRO_DATA> <MODEL_SPEC> [--threshold 1e-4]
                              [--atomic-dir /u/sbrennan/atomic]
                              [--write-stage1] [--write-stage2]
                              [--stage1-suffix _conv] [--stage2-suffix _full]
"""
import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np

K_B = 8.617e-5  # eV/K
CM1_TO_EV = 1.2398e-4


def read_sn_hydro_temperatures(sn_hydro_path):
    """Read temperature range from SN_HYDRO_DATA."""
    with open(sn_hydro_path) as f:
        lines = f.readlines()

    nd = int(lines[0].split(":")[1].strip())
    t_inner = float(lines[6].split(":")[1].strip())

    # Parse all numeric data after header
    data = []
    for line in lines[7:]:
        parts = line.strip().split()
        if parts:
            try:
                data.extend([float(x) for x in parts])
            except ValueError:
                pass

    n_vars = len(data) // nd
    # T is in 10^4 K (var index 3: after R, V, sigma)
    t_col = 3
    if n_vars <= t_col:
        # Try var 2 (R, V, T without sigma)
        t_col = 2
    t_10k = np.array(data[t_col * nd:(t_col + 1) * nd])
    t_kelvin = t_10k * 1e4

    return t_kelvin.min(), t_kelvin.max(), nd


def read_f_to_s_levels(f_to_s_path, n_levels):
    """Read level energies and superlevel assignments from F_TO_S file.

    Returns list of (level_id, g, E_cm, sl_group) tuples.
    """
    with open(f_to_s_path) as f:
        lines = f.readlines()

    levels = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("!") or stripped.startswith("*"):
            continue
        parts = stripped.split()
        if len(parts) >= 8:
            try:
                # Format: name  g  E(cm^-1)  freq  eV  SL  ?  ID
                g = float(parts[1])
                E_cm = float(parts[2])
                sl = int(parts[5])  # superlevel assignment (column 6, 0-indexed=5)
                level_id = int(parts[7])  # level ID (column 8, 0-indexed=7)
                if level_id <= n_levels:
                    levels.append((level_id, g, E_cm, sl))
            except (ValueError, IndexError):
                pass

    levels.sort(key=lambda x: x[0])
    return levels


def compute_lte_populations(levels, temperatures):
    """Compute LTE population fraction for each level at each temperature.

    Returns: (fractions, max_fraction) where fractions is [n_levels x n_T]
    """
    n_lev = len(levels)
    n_T = len(temperatures)

    fractions = np.zeros((n_lev, n_T))
    for j, T in enumerate(temperatures):
        kT = K_B * T
        # Partition function
        Z = sum(g * np.exp(-E_cm * CM1_TO_EV / kT) for _, g, E_cm, _ in levels)
        for i, (lid, g, E_cm, _) in enumerate(levels):
            fractions[i, j] = g * np.exp(-E_cm * CM1_TO_EV / kT) / Z

    max_frac = fractions.max(axis=1)
    return fractions, max_frac


def determine_optimal_isf(levels, max_frac, threshold=1e-4):
    """Determine optimal ISF (NLTE, NS, NT) based on LTE population significance.

    ISF format: NLTE, NS, NT where:
      NLTE = number of levels with full NLTE treatment (the first NLTE levels)
      NS   = number of superlevel groups (MUST match F_TO_S max SL link for NT levels)
      NT   = total levels in the model atom

    NS is fixed by the F_TO_S file — it must equal the max superlevel link
    for the first NT levels. We only optimize NLTE based on the threshold.

    Returns: (nlte, ns, nt)
    """
    n_levels = len(levels)

    # NS: max superlevel link from F_TO_S for the first NT levels only (fixed by atomic data)
    sl_groups = set()
    for i in range(min(n_levels, len(levels))):
        sl_groups.add(levels[i][3])
    ns = max(sl_groups) if sl_groups else 1

    # NLTE: levels with significant population get full NLTE treatment
    significant = max_frac > threshold
    nlte = 0
    for i in range(n_levels):
        if significant[i]:
            nlte = i + 1

    # Ensure at least 1 NLTE level
    nlte = max(nlte, 1)

    # NLTE cannot exceed NT
    nlte = min(nlte, n_levels)

    return nlte, ns, n_levels


def determine_convergence_isf(levels, max_frac, threshold=1e-2):
    """Determine reduced ISF for Stage 1 (convergence).

    Uses a higher threshold (1%) so only the most populated levels get NLTE.
    This ensures stable convergence. The full NLTE is recovered in Stage 2.
    """
    return determine_optimal_isf(levels, max_frac, threshold=threshold)


def parse_model_spec(model_spec_path):
    """Parse MODEL_SPEC file and return list of (species, isf, comment, line)."""
    species_list = []
    with open(model_spec_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("!"):
                continue
            # Format: NF,NS,NT  [SPECIES_ISF]  !comment
            match = re.match(r"(\d+),(\d+),(\d+)\s+\[(\w+)\]\s*(.*)", stripped)
            if match:
                nf, ns, nt = int(match.group(1)), int(match.group(2)), int(match.group(3))
                species = match.group(4)
                comment = match.group(5)
                species_list.append({
                    "species": species,
                    "nf": nf, "ns": ns, "nt": nt,
                    "comment": comment,
                    "original_line": line.rstrip(),
                })
    return species_list


def resolve_f_to_s_path(species, atomic_dir, cmfgen_dir=None):
    """Resolve the F_TO_S file path for a given species.

    Priority:
    1. If cmfgen_dir is given, check for the actual symlink in the run directory
       (e.g. Fe2_F_TO_S -> /u/sbrennan/atomic/FE/II/19apr23/f_to_s_135)
    2. Fall back to searching the atomic data directory
    """
    sp = species.replace("_ISF", "")

    # Check for actual symlink in CMFGEN run directory
    if cmfgen_dir:
        link_name = Path(cmfgen_dir) / f"{sp}_F_TO_S"
        if link_name.exists():
            target = link_name.resolve()
            if target.exists():
                return str(target)

    # Parse species name: e.g. Fe2 -> FE/II, OI -> OXY/I, HI -> HYD/I
    # Map species to element/ion
    ion_map = {"I": 1, "2": 2, "III": 3, "IV": 4, "V": 5}
    element_map = {
        "H": "HYD", "He": "HE", "C": "CARB", "O": "OXY", "Si": "SIL",
        "Ca": "CA", "Fe": "FE", "Ni": "NICK", "Co": "COB",
    }

    # Parse: element name + ion stage
    match = re.match(r"([A-Z][a-z]?)(I+|2|III|IV|V)", sp)
    if not match:
        return None

    elem = match.group(1)
    ion = match.group(2)
    ion_num = ion_map.get(ion, 1)
    ion_roman = {1: "I", 2: "II", 3: "III", 4: "IV", 5: "V"}[ion_num]
    elem_dir = element_map.get(elem, elem.upper())

    # Try to find any f_to_s file
    base = Path(atomic_dir) / elem_dir / ion_roman
    if base.exists():
        for p in sorted(base.rglob("f_to_s*")):
            return str(p)

    return None


def main():
    parser = argparse.ArgumentParser(description="Optimize MODEL_SPEC ISF based on SN_HYDRO_DATA")
    parser.add_argument("sn_hydro", help="Path to SN_HYDRO_DATA file")
    parser.add_argument("model_spec", help="Path to MODEL_SPEC file")
    parser.add_argument("--threshold", type=float, default=1e-4,
                        help="Population threshold for Stage 2 (full NLTE): 1e-4 = 0.01%%")
    parser.add_argument("--conv-threshold", type=float, default=1e-2,
                        help="Population threshold for Stage 1 (convergence): 1e-2 = 1%%")
    parser.add_argument("--atomic-dir", default="/u/sbrennan/atomic",
                        help="Atomic data directory")
    parser.add_argument("--cmfgen-dir", default=None,
                        help="CMFGEN run directory (to read actual F_TO_S symlinks)")
    parser.add_argument("--write-stage1", action="store_true",
                        help="Write Stage 1 MODEL_SPEC (reduced NLTE for convergence)")
    parser.add_argument("--write-stage2", action="store_true",
                        help="Write Stage 2 MODEL_SPEC (full NLTE for spectrum)")
    parser.add_argument("--stage1-suffix", default="_conv",
                        help="Suffix for Stage 1 MODEL_SPEC file")
    parser.add_argument("--stage2-suffix", default="_full",
                        help="Suffix for Stage 2 MODEL_SPEC file")
    parser.add_argument("--n-temps", type=int, default=10,
                        help="Number of temperature samples")
    args = parser.parse_args()

    # Read temperature range from SN_HYDRO_DATA
    t_min, t_max, nd = read_sn_hydro_temperatures(args.sn_hydro)
    temperatures = np.linspace(t_min, t_max, args.n_temps)
    print(f"SN_HYDRO_DATA: ND={nd}, T={t_min:.0f}-{t_max:.0f} K")
    print(f"Stage 1 threshold: {args.conv_threshold:.0e} (convergence)")
    print(f"Stage 2 threshold: {args.threshold:.0e} (full NLTE)")
    print()

    # Parse MODEL_SPEC
    species_list = parse_model_spec(args.model_spec)
    if not species_list:
        print("ERROR: No species found in MODEL_SPEC")
        sys.exit(1)

    # Analyze each species
    stage1_isf = {}
    stage2_isf = {}
    print(f"{'Species':<12} {'Current':<12} {'Stage1':<12} {'Stage2':<12} {'Status':<20}")
    print("-" * 70)

    for sp in species_list:
        species = sp["species"]
        nt = sp["nt"]
        current = f"{sp['nf']},{sp['ns']},{sp['nt']}"

        # Find F_TO_S file
        f_to_s_path = resolve_f_to_s_path(species, args.atomic_dir, args.cmfgen_dir)
        if f_to_s_path is None:
            print(f"{species:<12} {current:<12} {'?':<12} {'?':<12} {'F_TO_S not found':<20}")
            stage1_isf[species] = (sp["nf"], sp["ns"], sp["nt"])
            stage2_isf[species] = (sp["nf"], sp["ns"], sp["nt"])
            continue

        # Read levels
        levels = read_f_to_s_levels(f_to_s_path, nt)
        if not levels:
            print(f"{species:<12} {current:<12} {'?':<12} {'?':<12} {'No levels read':<20}")
            stage1_isf[species] = (sp["nf"], sp["ns"], sp["nt"])
            stage2_isf[species] = (sp["nf"], sp["ns"], sp["nt"])
            continue

        # Compute LTE populations
        _, max_frac = compute_lte_populations(levels, temperatures)

        # Determine optimal ISF
        nf1, ns1, _ = determine_convergence_isf(levels, max_frac, args.conv_threshold)
        nf2, ns2, _ = determine_optimal_isf(levels, max_frac, args.threshold)

        stage1_isf[species] = (nf1, ns1, nt)
        stage2_isf[species] = (nf2, ns2, nt)

        s1 = f"{nf1},{ns1},{nt}"
        s2 = f"{nf2},{ns2},{nt}"

        if s1 == current and s2 == current:
            status = "optimal"
        elif s2 == current:
            status = "already full NLTE"
        elif s1 == current:
            status = "needs Stage 2"
        else:
            status = "needs both stages"

        print(f"{species:<12} {current:<12} {s1:<12} {s2:<12} {status:<20}")

    # Write output files
    if args.write_stage1 or args.write_stage2:
        spec_path = Path(args.model_spec)
        stem = spec_path.stem
        parent = spec_path.parent

        if args.write_stage1:
            out1 = parent / f"{stem}{args.stage1_suffix}{spec_path.suffix}"
            with open(args.model_spec) as f:
                lines = f.readlines()
            with open(out1, "w") as f:
                for line in lines:
                    stripped = line.strip()
                    match = re.match(r"(\d+),(\d+),(\d+)\s+\[(\w+)\]", stripped)
                    if match and not stripped.startswith("!"):
                        species = match.group(4)
                        nf, ns, nt = stage1_isf.get(species, (int(match.group(1)), int(match.group(2)), int(match.group(3))))
                        comment = line.split("!", 1)[1] if "!" in line else ""
                        f.write(f"{nf},{ns},{nt:>5}        [{species}]{'':>4}! Stage 1 (conv): {nf} NLTE, {ns} SL, {nt} total{comment.rstrip()}\n")
                    else:
                        f.write(line)
            print(f"\nWrote Stage 1 MODEL_SPEC: {out1}")

        if args.write_stage2:
            out2 = parent / f"{stem}{args.stage2_suffix}{spec_path.suffix}"
            with open(args.model_spec) as f:
                lines = f.readlines()
            with open(out2, "w") as f:
                for line in lines:
                    stripped = line.strip()
                    match = re.match(r"(\d+),(\d+),(\d+)\s+\[(\w+)\]", stripped)
                    if match and not stripped.startswith("!"):
                        species = match.group(4)
                        nf, ns, nt = stage2_isf.get(species, (int(match.group(1)), int(match.group(2)), int(match.group(3))))
                        comment = line.split("!", 1)[1] if "!" in line else ""
                        f.write(f"{nf},{ns},{nt:>5}        [{species}]{'':>4}! Stage 2 (full): {nf} NLTE, {ns} SL, {nt} total{comment.rstrip()}\n")
                    else:
                        f.write(line)
            print(f"\nWrote Stage 2 MODEL_SPEC: {out2}")

    print("\nTwo-stage workflow:")
    print("  1. Run CMFGEN with Stage 1 MODEL_SPEC (reduced NLTE) → converges")
    print("  2. Rename *OUT to *_IN (departure coefficients from Stage 1)")
    print("  3. Run CMFGEN with Stage 2 MODEL_SPEC (full NLTE) + LIN_INT=T")
    print("  4. Stage 2 converges quickly from Stage 1 populations")


if __name__ == "__main__":
    main()
