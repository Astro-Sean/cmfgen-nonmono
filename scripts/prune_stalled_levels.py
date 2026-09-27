#!/usr/bin/env python3
"""Find and remove superlevels that stall CMFGEN convergence.

CMFGEN's MAXCH is an unweighted max-norm over all statistical-equilibrium
variables. A superlevel with negligible population can pin MAXCH while the
spectrum-forming levels are already converged: such levels typically want
~100% decreases every iteration but are throttled by the per-depth SCALE
cap, so they collapse geometrically over hundreds of iterations.

Two modes:

  * Default (post-run): parses the last K STEQ SOLUTION ARRAY blocks in
    STEQ_VALS to find which SE variables (species, superlevel, depth)
    still carry the largest fractional corrections.

  * --prior (pre-run): no STEQ_VALS needed. A superlevel is a prune
    candidate when it is Boltzmann-insignificant at every hydro
    temperature, carries no transition inside the diagnostic spectral
    range (F_OSCDAT), and sits more than --edge-guard eV below the
    highest level in its F_TO_S file (edge levels feed recombination
    cascades and are kept even when dead).

Shared steps:

  1. Maps SE rows to (species, superlevel) using the NS entries of
     MODEL_SPEC. Each element block contributes sum(NS) rows plus one
     row for the next ionisation stage's ground state; ED and T follow.
  2. Reads each species' F_TO_S file (run-directory symlink preferred,
     else the atomic archive via --atomic-dir) to get the full levels
     inside each superlevel and their energies.
  3. Reads the temperature profile from SN_HYDRO_DATA and computes the
     maximum LTE (Boltzmann) fraction of every superlevel.
  4. Proposes an ISF trim NV,NS,NF only when the doomed superlevels
     form a contiguous block at the top of the atom. RD_F_TO_S_IDS_V4
     requires MAXVAL(F_TO_S[1:NF]) == NS with every SL occupied, so only
     tail truncation of NF is legal; interior superlevels cannot be
     removed through MODEL_SPEC alone.

Usage:
    prune_stalled_levels.py RUN_DIR [options]

Options:
    --prior             a-priori mode: prune Boltzmann-dead, non-line-
                        forming tail SLs without needing STEQ_VALS
    --atomic-dir PATH   atomic archive root for F_TO_S/F_OSCDAT fallback
                        (default: $ATOMIC or /u/sbrennan/atomic)
    --extra-data-dir P  additional data dir (project atomic overrides)
    --spec-range LO,HI  diagnostic wavelength window in A (3000,25000)
    --edge-guard EV     keep SLs within EV of the atom's top level (2.0)
    --corr-min FRAC     flag SE rows with |correction| >= FRAC (default 0.02)
    --boltzmann-max F   max LTE fraction for a SL to be prunable (1e-4)
    --last K            number of STEQ blocks to analyse (default 8)
    --top N             rows shown in the stall table (default 15)
    --trim SP:NS        force species SP to keep NS superlevels (repeatable)
    --write             rewrite MODEL_SPEC in place (keeps .preprune backup)
    --output PATH       write pruned MODEL_SPEC here instead of in place
"""
import argparse
import fnmatch
import os
import re
import sys
from pathlib import Path

K_B_EV = 8.617333e-5  # eV/K
CM1_TO_EV = 1.2398e-4
DEFAULT_SPEC_RANGE = (3000.0, 25000.0)  # diagnostic window, Angstrom

ISF_RE = re.compile(r"^(\d+),(\d+),(\d+)\s+\[(\w+)\]")
ROMAN_RE = re.compile(r"^([A-Za-z]+?)(I{1,3}|IV|V|VI{0,3}|[2-9])$")

ROMAN_STR = {1: "I", 2: "II", 3: "III", 4: "IV", 5: "V", 6: "VI", 7: "VII"}
ROMAN_NUM = {v: k for k, v in ROMAN_STR.items()}

# Element label -> subdirectory name in the atomic archive. "Sk" is the
# pipeline's alias for silicon (avoids a SiI/Si filename collision).
# Verified against /u/sbrennan/atomic; species outside this map are
# skipped rather than guessed.
ELEM_DIRS = {
    "H": "HYD",
    "He": "HE",
    "C": "CARB",
    "N": "NIT",
    "O": "OXY",
    "F": "FLU",
    "Ne": "NEON",
    "Na": "NA",
    "Mg": "MG",
    "Al": "AL",
    "Sk": "SIL",
    "Si": "SIL",
    "P": "PHOS",
    "S": "SUL",
    "Cl": "CHL",
    "Ar": "ARG",
    "K": "POT",
    "Ca": "CA",
    "Ti": "TIT",
    "V": "VAN",
    "Cr": "CHRO",
    "Mn": "MAN",
    "Fe": "FE",
    "Co": "COB",
    "Ni": "NICK",
    "Ba": "BAR",
}


def element_of(key):
    """Species key to element label: O2_ISF -> O, FeIII_ISF -> Fe."""
    name = key.removesuffix("_ISF")
    m = ROMAN_RE.match(name)
    return m.group(1) if m else name


def ion_stage_of(key):
    """Species key to ionisation stage number: O2_ISF -> 2, FeIII -> 3."""
    name = key.removesuffix("_ISF")
    m = ROMAN_RE.match(name)
    if not m:
        return None
    suffix = m.group(2)
    if suffix.isdigit():
        return int(suffix)
    return ROMAN_NUM.get(suffix)


def read_batch_links(batch_sh, atomic_dir, project_dir):
    """Parse batch.sh `ln -sf TARGET LINK` lines into {link_name: Path}.

    batch.sh is the authoritative record of which atomic dataset a run
    links; the archive holds several variants per ion (f_to_s vs
    f_to_s_181), so a filename glob cannot reproduce the choice.
    """
    links = {}
    if batch_sh is None or not Path(batch_sh).exists():
        return links
    for line in Path(batch_sh).read_text(errors="replace").splitlines():
        s = line.strip()
        if s.startswith("#") or "ln -sf" not in s:
            continue
        toks = s.split()
        try:
            i = toks.index("ln")
        except ValueError:
            continue
        toks = [t for t in toks[i + 1 :] if not t.startswith("-")]
        if len(toks) < 2:
            continue
        target, link = toks[-2], toks[-1]
        target = target.replace("$ATOMIC", str(atomic_dir)).replace(
            "${ATOMIC}", str(atomic_dir)
        )
        if project_dir:
            target = target.replace("$PROJECT_DIR", str(project_dir)).replace(
                "${PROJECT_DIR}", str(project_dir)
            )
        links[link] = Path(target)
    return links


def resolve_ion_file(name, kind, run_dir, atomic_dir, extra_dirs, links=None):
    """Locate a species' F_TO_S or F_OSCDAT data file.

    Search order: the run directory symlink (the file CMFGEN actually
    uses), the batch.sh link table, project data dirs that override
    archive defaults, then dated datasets under ATOMIC/<ELEM>/<ION>/.
    Returns None when unresolved.
    """
    local = run_dir / f"{name}_{kind}"
    if local.exists():
        return local.resolve()
    if links and f"{name}_{kind}" in links:
        linked = links[f"{name}_{kind}"]
        if linked.exists():
            return linked
    # Project data dirs hold per-ion files, so the species name must
    # match; a dated archive dir holds one ion only, so a generic
    # pattern suffices there (e.g. OXY/I/19apr23/f_to_s_58).
    proj_pattern = {
        "F_TO_S": f"*{name.lower()}*f_to_s*",
        "F_OSCDAT": f"*{name.lower()}*osc*",
    }[kind]
    arch_pattern = {"F_TO_S": "*f_to_s*", "F_OSCDAT": "osc*"}[kind]
    candidates = []
    for base in (Path(d) for d in extra_dirs):
        if not base.is_dir():
            continue
        candidates += [
            f
            for f in base.iterdir()
            if f.is_file() and fnmatch.fnmatch(f.name.lower(), proj_pattern)
        ]
    elem_dir = ELEM_DIRS.get(element_of(name))
    roman = ROMAN_STR.get(ion_stage_of(name) or -1)
    if atomic_dir and elem_dir and roman:
        ion_base = Path(atomic_dir) / elem_dir / roman
        if ion_base.is_dir():
            subs = sorted(
                (d for d in ion_base.iterdir() if d.is_dir()),
                key=lambda d: (d.name != "19apr23", -d.stat().st_mtime),
            )
            for sub in subs:
                candidates += [
                    f
                    for f in sub.iterdir()
                    if f.is_file() and fnmatch.fnmatch(f.name.lower(), arch_pattern)
                ]
    if not candidates:
        return None
    preferred = f"{name.lower()}_f_to_s.dat" if kind == "F_TO_S" else "osc_data"
    candidates.sort(key=lambda f: (f.name.lower() != preferred, -f.stat().st_mtime))
    return candidates[0]


def read_line_levels(oscdat_path, lam_lo, lam_hi):
    """Level ids joined by any transition inside [lam_lo, lam_hi] A.

    A superlevel whose member levels appear here can form or feed
    diagnostic lines, so it stays even when its LTE fraction is tiny.
    """
    in_range = set()
    in_trans = False
    saw_trans = False
    try:
        lines = Path(oscdat_path).read_text(errors="replace").splitlines()
    except OSError:
        return None
    for line in lines:
        s = line.strip()
        if not s or s.startswith("!") or s.startswith("*"):
            continue
        if "Oscillator strengths" in s:
            in_trans = True
            saw_trans = True
            continue
        if not in_trans:
            continue
        # The i-j index pair is the last "N-M" field on the line; f, A
        # and wavelength are the three tokens before it. Wavelengths are
        # plain decimals (e.g. 3728.815), not E-notation.
        matches = list(re.finditer(r"(\d+)\s*-\s*(\d+)", s))
        if not matches:
            continue
        ij = matches[-1]
        fields = s[: ij.start()].split()
        if len(fields) < 3:
            continue
        try:
            lam = float(fields[-1])
        except ValueError:
            continue
        if lam_lo <= lam <= lam_hi:
            in_range.add(int(ij.group(1)))
            in_range.add(int(ij.group(2)))
    # A file without the expected header is an unknown format: return
    # None so the caller treats line-forming status as unverifiable
    # rather than absent.
    return in_range if saw_trans else None


def read_all_f_to_s_levels(path):
    """Return [(level_id, g, E_cm1, sl, name)] for every level in F_TO_S.

    Some datasets omit the trailing level-id column (e.g. the HI 15-level
    file); there the row order is the level id.
    """
    levels = []
    seq = 0
    for line in Path(path).read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("!") or s.startswith("*"):
            continue
        p = s.split()
        if len(p) < 7:
            continue
        try:
            g, e_cm, sl = float(p[1]), float(p[2]), int(p[5])
            lid = int(p[7]) if len(p) >= 8 else seq + 1
        except ValueError:
            continue
        seq += 1
        levels.append((lid, g, e_cm, sl, p[0]))
    return levels


def read_model_spec(path):
    """Return list of ISF entries plus the raw file lines."""
    entries = []
    lines = Path(path).read_text(errors="replace").splitlines()
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("!"):
            continue
        m = ISF_RE.match(stripped)
        if m:
            entries.append(
                {
                    "key": m.group(4),
                    "nv": int(m.group(1)),
                    "ns": int(m.group(2)),
                    "nf": int(m.group(3)),
                    "line_idx": i,
                }
            )
    return entries, lines


def build_row_map(entries):
    """Map SE row index -> (species key, SL index, is_ion_ground).

    Ordering follows MODEL_SPEC: each element block is its ions' NS rows
    in listed order, then one row for the next ion stage's ground state.
    ED and T occupy the final two rows.
    """
    row_map = []
    i = 0
    while i < len(entries):
        elem = element_of(entries[i]["key"])
        j = i
        while j < len(entries) and element_of(entries[j]["key"]) == elem:
            e = entries[j]
            for sl in range(1, e["ns"] + 1):
                row_map.append((e["key"], sl, False))
            j += 1
        row_map.append((f"{elem}ion", 0, True))
        i = j
    row_map.append(("ED", 0, False))
    row_map.append(("T", 0, False))
    return row_map


def read_steq_solutions(path, n_blocks):
    """Return list of {row: {depth: corr}} for the last n_blocks arrays."""
    if not Path(path).exists():
        return []
    lines = Path(path).read_text(errors="replace").splitlines()
    starts = [i for i, line in enumerate(lines) if "STEQ SOLUTION ARRAY" in line]
    if not starts:
        return []
    blocks = []
    row_re = re.compile(r"^\s*(\d+)\(\s*(\d+)\)#(.*)")
    for start in starts[-n_blocks:]:
        rows = {}
        for line in lines[start + 1 : start + 5000]:
            m = row_re.match(line)
            if m:
                row, d0 = int(m.group(1)), int(m.group(2))
                try:
                    vals = [float(x) for x in m.group(3).split()]
                except ValueError:
                    continue
                rec = rows.setdefault(row, {})
                for k, v in enumerate(vals):
                    rec[d0 + k] = v
            elif "STEQ" in line and rows:
                break
        if rows:
            blocks.append(rows)
    return blocks


def read_f_to_s(path, nf):
    """Return {sl: [(level_id, g, E_cm1, name)]} for the first nf levels."""
    groups = {}
    seq = 0
    for line in Path(path).read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("!") or s.startswith("*"):
            continue
        p = s.split()
        if len(p) < 7:
            continue
        try:
            g, e_cm, sl = float(p[1]), float(p[2]), int(p[5])
            lid = int(p[7]) if len(p) >= 8 else seq + 1
        except ValueError:
            continue
        seq += 1
        if 1 <= lid <= nf and sl > 0:
            groups.setdefault(sl, []).append((lid, g, e_cm, p[0]))
    for sl in groups:
        groups[sl].sort()
    return groups


def read_hydro_temperatures(path):
    """Return (T_min, T_max) in K from SN_HYDRO_DATA's temperature block."""
    lines = Path(path).read_text(errors="replace").splitlines()
    vals = []
    in_block = False
    for line in lines:
        if line.startswith("Temperature"):
            in_block = True
            continue
        if in_block:
            if re.match(r"^[A-Za-z]", line):
                break
            vals.extend(float(x) for x in line.split())
    if not vals:
        raise ValueError(f"no Temperature block in {path}")
    temps = [v * 1.0e4 for v in vals]
    return min(temps), max(temps)


def lte_fractions(all_levels, t_kelvin):
    """LTE population fraction of each level at t_kelvin.

    Boltzmann factor normalized by the partition sum over the levels in
    the model atom. NLTE departure coefficients can move this by orders
    of magnitude, so the result is a significance guide, not a guarantee.
    """
    k_t = K_B_EV * t_kelvin
    weights = [
        g * pow(2.718281828459045, -e * CM1_TO_EV / k_t) for _, g, e, _ in all_levels
    ]
    z = sum(weights)
    return [w / z for w in weights]


def sl_boltzmann_fraction(members, all_levels, t_min, t_max, n_t=8):
    """Max LTE fraction of a superlevel over the hydro T range.

    Returns the largest fraction the superlevel reaches at any sampled
    temperature. High levels are most populated at T_max, but sampling
    the range guards against odd partition-sum behaviour at low T.
    """
    best = 0.0
    for i in range(n_t):
        t = t_min + (t_max - t_min) * i / max(n_t - 1, 1)
        frac = lte_fractions(all_levels, t)
        total = sum(frac[lid - 1] for lid, _, _, _ in members)
        best = max(best, total)
    return best


def parse_args():
    p = argparse.ArgumentParser(
        description="Prune stalling superlevels from a CMFGEN MODEL_SPEC."
    )
    p.add_argument("run_dir", help="CMFGEN run directory (MODEL_SPEC, STEQ_VALS)")
    p.add_argument(
        "--spec", default=None, help="MODEL_SPEC path (default: RUN_DIR/MODEL_SPEC)"
    )
    p.add_argument(
        "--prior",
        action="store_true",
        help="a-priori mode: no STEQ_VALS needed; prune dead tail SLs",
    )
    p.add_argument(
        "--atomic-dir",
        default=os.environ.get("ATOMIC", "/u/sbrennan/atomic"),
        help="atomic archive root for F_TO_S/F_OSCDAT fallback",
    )
    p.add_argument(
        "--extra-data-dir",
        action="append",
        default=[],
        help="project data dir searched before the archive (repeatable)",
    )
    p.add_argument(
        "--batch-sh",
        default=None,
        help="batch.sh to read F_TO_S/F_OSCDAT link targets from "
        "(default: RUN_DIR/batch.sh)",
    )
    p.add_argument(
        "--project-dir",
        default=None,
        help="project root for $PROJECT_DIR expansion in batch.sh links",
    )
    p.add_argument(
        "--spec-range",
        default=f"{DEFAULT_SPEC_RANGE[0]:.0f},{DEFAULT_SPEC_RANGE[1]:.0f}",
        help="diagnostic wavelength window LO,HI in A",
    )
    p.add_argument(
        "--edge-guard",
        type=float,
        default=2.0,
        help="keep SLs within EV of the atom's top level (default 2.0)",
    )
    p.add_argument(
        "--corr-min", type=float, default=0.02, help="|corr| stall flag (default 0.02)"
    )
    p.add_argument(
        "--boltzmann-max",
        type=float,
        default=1e-4,
        help="max LTE fraction for a prunable SL (default 1e-4)",
    )
    p.add_argument(
        "--last", type=int, default=8, help="STEQ blocks analysed (default 8)"
    )
    p.add_argument(
        "--top", type=int, default=15, help="rows in stall table (default 15)"
    )
    p.add_argument(
        "--trim",
        action="append",
        default=[],
        metavar="SP:NS",
        help="force species SP to keep NS superlevels (repeatable)",
    )
    p.add_argument("--write", action="store_true", help="rewrite MODEL_SPEC in place")
    p.add_argument("--output", default=None, help="write pruned MODEL_SPEC to PATH")
    return p.parse_args()


def main():
    args = parse_args()
    run_dir = Path(args.run_dir)
    spec_path = Path(args.spec) if args.spec else run_dir / "MODEL_SPEC"

    entries, spec_lines = read_model_spec(spec_path)
    if not entries:
        sys.exit(f"ERROR: no ISF entries in {spec_path}")
    row_map = build_row_map(entries)
    print(f"MODEL_SPEC: {len(entries)} ions, {len(row_map)} SE rows expected")

    blocks = read_steq_solutions(run_dir / "STEQ_VALS", args.last)
    if not blocks and not args.prior:
        sys.exit(f"ERROR: no STEQ SOLUTION ARRAY in {run_dir / 'STEQ_VALS'}")
    stall_rows = {}
    if blocks:
        nt = max(max(b) for b in blocks)
        print(f"STEQ_VALS: {len(blocks)} block(s), NT={nt}")
        if nt != len(row_map):
            print(
                f"WARNING: NT={nt} != predicted rows={len(row_map)}; "
                "row map may be off"
            )

        # Worst correction per row over the last blocks.
        worst = {}
        for rows in blocks:
            for row, dv in rows.items():
                for d, v in dv.items():
                    if row not in worst or abs(v) > abs(worst[row][0]):
                        worst[row] = (v, d)

        ranked = sorted(worst.items(), key=lambda kv: -abs(kv[1][0]))
        print(f"\nTop {args.top} SE rows by |correction| (last block max):")
        print(f"  {'row':>4} {'species':<14} {'SL':>3} {'depth':>5} {'corr':>9}")
        for idx, (row, (v, d)) in enumerate(ranked):
            if row - 1 < len(row_map):
                key, sl, ground = row_map[row - 1]
            else:
                key, sl, ground = ("?", 0, False)
            if abs(v) >= args.corr_min and not ground and sl > 0:
                stall_rows.setdefault(key, set()).add(sl)
            if idx < args.top:
                tag = "gnd" if ground else f"SL{sl}"
                print(f"  {row:>4} {key:<14} {tag:>3} {d:>5} {v:>+9.4f}")
    else:
        print("STEQ_VALS: not used (a-priori mode)")

    # Temperature range for the Boltzmann test.
    t_min, t_max = read_hydro_temperatures(run_dir / "SN_HYDRO_DATA")
    print(f"\nSN_HYDRO_DATA: T = {t_min:.0f} - {t_max:.0f} K")

    forced = {}
    for item in args.trim:
        key, keep = item.rsplit(":", 1)
        key = key.removesuffix("_ISF")
        forced[key] = int(keep)

    batch_sh = Path(args.batch_sh) if args.batch_sh else run_dir / "batch.sh"
    links = read_batch_links(batch_sh, args.atomic_dir, args.project_dir)
    lam_lo, lam_hi = (float(x) for x in args.spec_range.split(","))
    print(
        f"\nSuperlevel audit (stall flag = |corr|>={args.corr_min:.0%}; "
        f"Boltzmann bound at T_max={t_max:.0f} K; "
        f"line window {lam_lo:.0f}-{lam_hi:.0f} A; "
        f"edge guard {args.edge_guard:.1f} eV):"
    )
    proposals = {}
    for e in entries:
        key = e["key"]
        name = key.removesuffix("_ISF")
        f_to_s = resolve_ion_file(
            name, "F_TO_S", run_dir, args.atomic_dir, args.extra_data_dir, links
        )
        if f_to_s is None:
            print(
                f"  {key:<10} {e['nv']},{e['ns']},{e['nf']:<4} "
                "F_TO_S not found - skipped"
            )
            continue
        groups = read_f_to_s(f_to_s, e["nf"])
        occupied = sorted(groups)
        contiguous = occupied == list(range(1, len(occupied) + 1))
        one_to_one = e["nf"] == e["ns"]
        all_levels = [m for g in groups.values() for m in g]
        file_levels = read_all_f_to_s_levels(f_to_s)
        e_edge = max(m[2] for m in file_levels) * CM1_TO_EV if file_levels else None
        oscdat = resolve_ion_file(
            name, "F_OSCDAT", run_dir, args.atomic_dir, args.extra_data_dir, links
        )
        line_levels = (
            read_line_levels(oscdat, lam_lo, lam_hi) if oscdat is not None else None
        )

        print(
            f"  {key:<10} current {e['nv']},{e['ns']},{e['nf']}"
            + ("  [1:1 mapping - F_TO_S bypassed]" if one_to_one else "")
        )
        prunable = {}
        for sl in occupied:
            members = groups[sl]
            e_min = min(m[2] for m in members) * CM1_TO_EV
            e_max = max(m[2] for m in members) * CM1_TO_EV
            boltz = sl_boltzmann_fraction(members, all_levels, t_min, t_max)
            dead = boltz < args.boltzmann_max
            has_lines = (
                any(m[0] in line_levels for m in members)
                if line_levels is not None
                else None
            )
            edge_ok = e_edge is not None and (e_edge - e_max) > args.edge_guard
            st = " STALL" if sl in stall_rows.get(key, set()) else ""
            sig = " negligible" if dead else ""
            lin = " lines" if has_lines else (" lines?" if has_lines is None else "")
            edg = "" if edge_ok else " EDGE"
            print(
                f"    SL{sl:>3}: {len(members):>3} lvls {members[0][0]:>3}-{members[-1][0]:>3}"
                f"  E {e_min:7.2f}-{e_max:7.2f} eV  gB {boltz:9.2e}{st}{sig}{lin}{edg}"
                f"    e.g. {members[0][3]}"
            )
            stalled = sl in stall_rows.get(key, set())
            if args.prior:
                prunable[sl] = dead and has_lines is False and edge_ok
            else:
                prunable[sl] = stalled and dead and has_lines is not True and edge_ok

        # A trim is legal only if the pruned SLs form a contiguous block at
        # the top of the occupied SL list and the F_TO_S check passes.
        keep = forced.get(name)
        if keep is None:
            if one_to_one or not contiguous:
                continue
            keep = len(occupied)
            while keep > 1 and prunable.get(occupied[keep - 1], False):
                keep -= 1
            if keep == len(occupied):
                continue
        if keep < 1 or keep > len(occupied):
            print(f"    -> requested keep={keep} invalid; skipped")
            continue
        new_ns = occupied[keep - 1] if contiguous else keep
        retained_max_level = max(m[0] for s in occupied[:keep] for m in groups[s])
        new_nv = min(e["nv"], new_ns)
        proposals[key] = (new_nv, new_ns, retained_max_level)
        print(
            f"    -> PROPOSE {key} {e['nv']},{e['ns']},{e['nf']} -> "
            f"{new_nv},{new_ns},{retained_max_level}"
        )

    if not proposals:
        print("\nNo legal trims found.")
        return

    if not args.write and not args.output:
        print("\nDry run - rerun with --write or --output PATH to apply.")
        return

    new_lines = list(spec_lines)
    for e in entries:
        if e["key"] not in proposals:
            continue
        nv, ns, nf = proposals[e["key"]]
        old = spec_lines[e["line_idx"]]
        comment = old.split("!", 1)[1].rstrip() if "!" in old else ""
        new_lines[e["line_idx"]] = (
            f"{nv},{ns},{nf:<5}       [{e['key']}]"
            f"          !{comment} | pruned {e['nv']},{e['ns']},{e['nf']} "
            f"by prune_stalled_levels.py"
        )

    out = Path(args.output) if args.output else spec_path
    if args.write and args.output is None:
        backup = spec_path.with_suffix(spec_path.suffix + ".preprune")
        backup.write_text(spec_path.read_text(errors="replace"))
        print(f"Backup written: {backup}")
    out.write_text("\n".join(new_lines) + "\n")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
