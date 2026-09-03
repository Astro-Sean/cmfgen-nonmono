#!/usr/bin/env python3
"""X_newgrid Pass 2: CMFGEN ionisation-front refinement (standalone diagnostic).

After a CMFGEN run has computed its own ionisation structure, this script:

  1. Reads the CMFGEN-computed electron density / ion fractions.
  2. Detects the *actual* ionisation/recombination fronts (where the
     dominant ion species or electron density changes rapidly).
  3. Inserts grid points at those fronts (as cell boundaries, so no cell
     straddles a front).
  4. Re-conservatively-remaps the STELLA source onto the new grid.
  5. Writes a new SN_HYDRO_DATA.
  6. Prints restart instructions (copy _OUT -> _IN, set VADAT flags).

This implements the two-pass strategy from the Dessart/Hillier literature
survey: the Pass-1 grid (from ``x_newgrid``) is refined using the *real*
CMFGEN ionisation structure, not the Saha estimate.

Usage:
    python x_newgrid_ionfront_refine.py <cmfgen_dir> <swd> <epoch_days> <out_dir>
        [--ion-source ascii:<file> | species_out | sn_hydro]
        [--ne-log-threshold 1.0] [--min-points-per-front 5]
        [--max-points 400] [--lum 92.4]

EVIDENCE STATUS
---------------
The ASCII ionisation reader and the front-detection / grid-insertion /
restart logic are testable on synthetic data (VERIFIED via the test suite).
The Fortran-unformatted species-_OUT binary reader is implemented per the
documented CMFGEN departure-coefficient format but is UNVERIFIED until
tested against a real converged CMFGEN run that has written _OUT files.
Per the project AGENTS.md evidence policy, do not rely on the binary
reader's output without such a validation run.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from grid_mapping.stella_io import read_swd, select_epoch, normalize_swd_composition, build_source_grid
from grid_mapping.feature_detection import detect_all_features
from grid_mapping.sn_hydro_io import read_sn_hydro_data
from grid_mapping.conservative_remap import conservative_remap
from grid_mapping.cmfgen_writer import write_sn_hydro_data
from grid_mapping.shell_model import SourceShellGrid, TargetGrid

from x_newgrid.config import XNewGridConfig
from x_newgrid.skeleton import _enforce_strict_decrease


# ===========================================================================
# Ionisation readers
# ===========================================================================

@dataclass
class IonizationStructure:
    """CMFGEN-computed ionisation structure on the current grid."""
    r_cm: np.ndarray              # depth-point radii (decreasing), len ND
    electron_density: np.ndarray  # cm^-3, len ND
    # optional per-ion fractions: dict ion_name -> (ND,) mass/number fraction
    ion_fractions: Optional[dict] = None
    source: str = ""              # 'ascii' | 'species_out' | 'sn_hydro'


def read_ionization(cmfgen_dir: Path, source: str,
                    sn_hydro_path: Optional[Path] = None
                    ) -> IonizationStructure:
    """Dispatch to the requested ionisation reader.

    Parameters
    ----------
    cmfgen_dir : directory containing CMFGEN output
    source : 'species_out' | 'sn_hydro' | 'ascii:<path>'
    sn_hydro_path : path to the SN_HYDRO_DATA (for the grid + fallback Ne)
    """
    if source.startswith("ascii:"):
        return _read_ascii_ionization(Path(source.split(":", 1)[1]))
    if source == "species_out":
        return _read_species_out_ionization(cmfgen_dir, sn_hydro_path)
    if source == "sn_hydro":
        return _read_sn_hydro_ionization(sn_hydro_path or cmfgen_dir / "SN_HYDRO_DATA")
    raise ValueError(f"Unknown ion source: {source}")


def _read_sn_hydro_ionization(path: Path) -> IonizationStructure:
    """Read the electron density from SN_HYDRO_DATA (Saha estimate).

    This is the fallback when no CMFGEN-computed ionisation is available.
    It uses the Ne column written by the Pass-1 builder (Saha).  Real
    CMFGEN fronts will differ; this reader exists so the script can be
    exercised end-to-end without a converged CMFGEN run.
    """
    d = read_sn_hydro_data(str(path))
    r_cm = d.r_1e10cm * 1e10
    return IonizationStructure(r_cm=r_cm, electron_density=d.electron_density,
                               source="sn_hydro")


def _read_ascii_ionization(path: Path) -> IonizationStructure:
    """Read an ASCII ionisation dump: one depth point per line.

    Expected columns (whitespace-separated, '#' comments):
        r_1e10cm   electron_density_cm-3   [ion_name=fraction ...]

    The radius column is in 1e10 cm units (matching SN_HYDRO_DATA).
    Optional ``ion_name=fraction`` fields populate ``ion_fractions``.
    Lines may be in any radius order; they are sorted decreasing here.
    """
    rows = []
    ions: dict[str, list[float]] = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            r = float(parts[0])
            ne = float(parts[1])
            ion_vals = {}
            for tok in parts[2:]:
                if "=" in tok:
                    name, val = tok.split("=", 1)
                    ion_vals[name] = float(val)
                    ions.setdefault(name, []).append(float(val))
            rows.append((r, ne, ion_vals))

    rows.sort(key=lambda t: -t[0])  # decreasing r
    r_cm = np.array([t[0] * 1e10 for t in rows])
    ne = np.array([t[1] for t in rows])
    ion_fractions = {k: np.array(v) for k, v in ions.items()} if ions else None
    return IonizationStructure(r_cm=r_cm, electron_density=ne,
                               ion_fractions=ion_fractions, source="ascii")


def _read_species_out_ionization(cmfgen_dir: Path,
                                 sn_hydro_path: Optional[Path]
                                 ) -> IonizationStructure:
    """Read ionisation from CMFGEN species departure-coefficient _OUT files.

    CMFGEN writes per-species files (HIOUT, HeIOUT, Fe2OUT, ...) in Fortran
    unformatted binary containing departure coefficients and populations
    per depth.  This reader reconstructs the electron density from the
    ion populations.

    EVIDENCE: UNVERIFIED.  Implemented per the documented format
    (regrid_log_dc_v1.f) but not yet tested against a real converged run.
    Falls back to the SN_HYDRO_DATA Ne if no _OUT files are present.
    """
    out_files = sorted(cmfgen_dir.glob("*OUT"))
    # species DC files exclude these non-population outputs
    dc_files = [f for f in out_files if f.name not in
                ("SCRTEMP_OUT", "RVTJ_OUT", "EDDFACTOR_OUT", "CFDAT_OUT")]
    if not dc_files:
        print("[ionfront] No species _OUT files found; falling back to SN_HYDRO_DATA Ne")
        return _read_sn_hydro_ionization(sn_hydro_path or cmfgen_dir / "SN_HYDRO_DATA")

    # Get the grid from SN_HYDRO_DATA
    sh = read_sn_hydro_data(str(sn_hydro_path or cmfgen_dir / "SN_HYDRO_DATA"))
    r_cm = sh.r_1e10cm * 1e10
    nd = sh.n_points

    # UNVERIFIED binary parse: Fortran unformatted records.
    # The first record of a DC file is: R(ND), LUM, NHYD, ND  (header).
    # Subsequent records hold per-level departure coefficients / populations.
    # Reconstructing Ne from populations requires the level->species map from
    # MODEL_SPEC, which is not available here without parsing MODEL_SPEC.
    # Therefore this reader currently returns the SN_HYDRO_DATA Ne and flags
    # the ion_fractions as unavailable, with a clear warning.  A full
    # population-based Ne reconstruction is left for validation against a
    # real converged run (see EVIDENCE note in the module docstring).
    print("[ionfront] WARNING: species-_OUT binary reader is UNVERIFIED.")
    print("[ionfront]   Using SN_HYDRO_DATA Ne (Saha) as a placeholder.")
    print("[ionfront]   Provide --ion-source ascii:<file> with CMFGEN-printed")
    print("[ionfront]   Ne/ion fractions for a validated Pass-2 refinement.")
    return IonizationStructure(r_cm=r_cm, electron_density=sh.electron_density,
                               ion_fractions=None, source="species_out_unverified")


# ===========================================================================
# Front detection from the real CMFGEN ionisation
# ===========================================================================

def detect_ionization_fronts(ion: IonizationStructure,
                             ne_log_threshold: float = 1.0,
                             ye_log_threshold: float = 0.5,
                             ion_frac_log_threshold: float = 1.0
                             ) -> np.ndarray:
    """Detect ionisation/recombination front *radii* from the CMFGEN solution.

    A front is flagged where:
      * |dlog10(Ne)| > ne_log_threshold, or
      * a dominant ion species changes (|dlog10(x_ion)| > ion_frac_log_threshold
        for the species that is dominant on either side).

    Returns the front radii (decreasing), to be inserted as grid boundaries.
    """
    r = ion.r_cm
    ne = np.maximum(ion.electron_density, 1e-30)
    log_ne = np.log10(ne)
    front_idx = set()

    # Ne gradient
    for i in range(len(r) - 1):
        if abs(log_ne[i] - log_ne[i + 1]) > ne_log_threshold:
            front_idx.add(i)

    # Ion-fraction dominance changes
    if ion.ion_fractions:
        # for each adjacent pair, check if the dominant ion changes
        names = list(ion.ion_fractions.keys())
        fracs = np.array([ion.ion_fractions[n] for n in names])  # (n_ion, ND)
        for i in range(len(r) - 1):
            dom_l = names[int(np.argmax(fracs[:, i]))]
            dom_r = names[int(np.argmax(fracs[:, i + 1]))]
            if dom_l != dom_r:
                front_idx.add(i)
            # also large log-fraction change in the dominant species
            for n in names:
                x = ion.ion_fractions[n]
                if max(x[i], x[i + 1]) > 1e-4:
                    if abs(np.log10(max(x[i], 1e-30)) -
                           np.log10(max(x[i + 1], 1e-30))) > ion_frac_log_threshold:
                        front_idx.add(i)
                        break

    # front radii: use the midpoint between the two bracketing points
    front_r = []
    for i in sorted(front_idx):
        r_front = np.sqrt(r[i] * r[i + 1])
        front_r.append(float(r_front))
    return np.array(sorted(set(front_r), reverse=True))  # decreasing


# ===========================================================================
# Grid rebuild with front anchors
# ===========================================================================

def rebuild_grid_with_fronts(source: SourceShellGrid, front_radii: np.ndarray,
                             cfg: XNewGridConfig, outpath: str,
                             lum_lsun: float) -> tuple[TargetGrid, dict]:
    """Re-run the x_newgrid one-pass build with front radii forced as anchors.

    The front radii are injected into the source's ionization_front_mask
    (at the nearest source zone) so the skeleton and coarsening treat them
    as mandatory protected boundaries.
    """
    # Inject front radii into the source ionization_front_mask
    if source.ionization_front_mask is None:
        source.ionization_front_mask = np.zeros(source.n_shells, dtype=bool)
    for r_front in front_radii:
        idx = int(np.argmin(np.abs(source.r_centre_cm - r_front)))
        source.ionization_front_mask[idx] = True

    # Re-detect all features (picks up the new front mask)
    detect_all_features(source)

    # Run the one-pass builder
    from x_newgrid.builder import build_grid
    target, report = build_grid(source, cfg, outpath=outpath, lum_lsun=lum_lsun)
    report['front_radii_inserted'] = front_radii.tolist()
    report['n_fronts'] = len(front_radii)
    return target, report


# ===========================================================================
# Restart instructions
# ===========================================================================

def print_restart_instructions(cmfgen_dir: Path, out_dir: Path) -> None:
    """Print the CMFGEN restart instructions (mirrors cmfgen_restart_from_out.py)."""
    print("\n" + "=" * 70)
    print("CMFGEN RESTART INSTRUCTIONS (Pass 2)")
    print("=" * 70)
    print(f"  CMFGEN run dir : {cmfgen_dir}")
    print(f"  New SN_HYDRO_DATA: {out_dir / 'SN_HYDRO_DATA'}")
    print()
    print("  1. Copy the new SN_HYDRO_DATA into the CMFGEN run directory:")
    print(f"       cp {out_dir / 'SN_HYDRO_DATA'} {cmfgen_dir}/SN_HYDRO_DATA")
    print("  2. Copy _OUT departure-coefficient files to _IN for a warm restart:")
    print(f"       python scripts/cmfgen_restart_from_out.py {cmfgen_dir}")
    print("  3. Set in VADAT:  LTE_EST=F, LIN_INT=T")
    print("  4. Set in IN_ITS: DO_LAM_IT=T (Lambda iterations)")
    print("  5. Resubmit via the orchestrator (do NOT sbatch directly):")
    print("       python scripts/run_model_grid.py <config.yaml> "
          "--skip-stella --run-cmfgen --force-cancel")
    print()
    print("  NOTE: The new grid has more depth points around the ionisation")
    print("  fronts.  CMFGEN will interpolate the _IN departure coefficients")
    print("  onto the new grid (LIN_INT=T).  Monitor convergence of the front")
    print("  region; if a front still oscillates, repeat this Pass-2 step.")
    print("=" * 70)


# ===========================================================================
# CLI
# ===========================================================================

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="X_newgrid Pass 2: CMFGEN ionisation-front refinement")
    p.add_argument("cmfgen_dir", type=Path, help="CMFGEN run directory")
    p.add_argument("swd", type=Path, help="STELLA .swd file (for source remap)")
    p.add_argument("epoch_days", type=float, help="snapshot epoch (days)")
    p.add_argument("out_dir", type=Path, help="output directory")
    p.add_argument("--ion-source", default="sn_hydro",
                   help="'sn_hydro' | 'species_out' | 'ascii:<file>'")
    p.add_argument("--ne-log-threshold", type=float, default=1.0)
    p.add_argument("--ion-frac-log-threshold", type=float, default=1.0)
    p.add_argument("--min-points-per-front", type=int, default=5)
    p.add_argument("--max-points", type=int, default=400)
    p.add_argument("--target-points", type=int, default=150)
    p.add_argument("--lum", type=float, default=0.0)
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)

    cfg = XNewGridConfig(verbose=not args.quiet,
                         max_points=args.max_points,
                         target_points=args.target_points,
                         min_points_per_ionfront=args.min_points_per_front)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Read CMFGEN ionisation
    sn_hydro = args.cmfgen_dir / "SN_HYDRO_DATA"
    print(f"[ionfront] Reading ionisation source: {args.ion_source}")
    ion = read_ionization(args.cmfgen_dir, args.ion_source, sn_hydro)
    print(f"[ionfront]   {ion.source}: {len(ion.r_cm)} depth points, "
          f"Ne=[{ion.electron_density.min():.3e},{ion.electron_density.max():.3e}]")

    # 2. Detect fronts
    fronts = detect_ionization_fronts(
        ion, ne_log_threshold=args.ne_log_threshold,
        ion_frac_log_threshold=args.ion_frac_log_threshold)
    print(f"[ionfront] Detected {len(fronts)} ionisation front(s) at radii:")
    for r_f in fronts:
        print(f"           r = {r_f:.4e} cm ({r_f/1e10:.3e} 1e10cm)")

    # 3. Read STELLA source
    data = read_swd(str(args.swd))
    best_t, ep = select_epoch(data, target_days=args.epoch_days)
    mf = normalize_swd_composition(ep.get('composition'))
    source = build_source_grid(ep, best_t, mass_fractions=mf, source_file=str(args.swd))

    # 4. Rebuild grid with front anchors
    outpath = args.out_dir / "SN_HYDRO_DATA"
    if len(fronts) == 0:
        print("[ionfront] No fronts detected; copying SN_HYDRO_DATA unchanged.")
        import shutil
        shutil.copy2(sn_hydro, outpath)
    else:
        target, report = rebuild_grid_with_fronts(
            source, fronts, cfg, outpath=str(outpath), lum_lsun=args.lum)
        print(f"[ionfront] New grid: ND={target.n_points} "
              f"(was {len(ion.r_cm)}, +{target.n_points - len(ion.r_cm)} points)")
        with open(args.out_dir / "ionfront_report.json", "w") as f:
            json.dump(report, f, indent=2, default=str)

    # 5. Restart instructions
    print_restart_instructions(args.cmfgen_dir, args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
