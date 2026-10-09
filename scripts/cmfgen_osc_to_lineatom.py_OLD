#!/usr/bin/env python3
"""
Convert CMFGEN atomic data (osc_data / *_osc files in $ATOMIC) into a STELLA
``lineatom.dat`` line list for the vladsf expansion-opacity machinery
(``lineexpop`` / ``xronfict.exe``).

STELLA record format (rlnexp.f, format 30):

    format(f10.4,1x,a4,1x,f7.3,f4.1,f11.3,f4.1,1x,f11.3,f7.2,a8,2x,a8)

    wl(nm) ref loggf J_l E_l(cm^-1) J_u E_u(cm^-1) Z.(stage-1) namelw nameup

with an ``i10`` line-count header as the first record.  ``ionstage`` is
decoded as ``int(100*(code - Z) + 0.1) + 1`` (stage 1 = neutral).  The J
values are read into scalars that are never used by the opacity
calculation; only wl, log(gf), E_l, ion stage and Z enter the Sobolev
expansion opacity (gf carries the lower-level statistical weight).

CMFGEN ``osc_data`` format (Hillier atomic tree):

    header ... !Number of energy levels / !Ionization energy /
               !Screened nuclear charge / !Number of transitions
    level lines:   name  g  E(cm^-1)  [various]  ID
    transition lines (after a '  Transition' header):
                   lowname - upname  f  A  Lam(A)  i-j  Trans.# ...

gf for STELLA = g_lower * f_absorption (standard absorption convention).

Modes:
  * ``--osc Z:STAGE:PATH`` repeated: convert those files only.
  * ``--from-batch-sh PATH``: extract every ``*_F_OSCDAT`` link target from
    a CMFGEN batch.sh and derive (Z, stage) from the $ATOMIC directory
    layout ``<EL>/<ROMAN>/<date>/<file>``.
  * ``--merge-with PATH``: keep lines from an existing lineatom.dat for
    (Z, stage) NOT covered by the converted set (merge-mode 'replace',
    default), so e.g. Ni/Co/Ti blanketing is preserved.  'only' writes a
    pure CMFGEN-converted list.

A JSON provenance sidecar (<output>.provenance.json) records sources and
line counts.
"""

import argparse
import json
import re
import sys
from pathlib import Path

# Hillier $ATOMIC element directory names -> nuclear charge Z.
ELEMENT_DIR_Z = {
    'HYD': 1, 'HE': 2, 'LIT': 3, 'BER': 4, 'BOR': 5, 'CARB': 6, 'NIT': 7,
    'OXY': 8, 'FLU': 9, 'NEON': 10, 'NA': 11, 'MG': 12, 'AL': 13, 'SIL': 14,
    'PHOS': 15, 'SUL': 16, 'CHL': 17, 'ARG': 18, 'POT': 19, 'CA': 20,
    'SCAN': 21, 'TIT': 22, 'VAN': 23, 'CHRO': 24, 'MAN': 25, 'FE': 26,
    'COB': 27, 'NICK': 28, 'COP': 29, 'ZINC': 30,
}

ROMAN = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7,
         'VIII': 8, 'IX': 9, 'X': 10, 'XI': 11, 'XII': 12}

# STELLA RT grid is ~1 A - 5e4 A (ronfitab.trf: wlmin=1, wlmax=5e4 AA).
# Keep a generous margin; the f10.4 field caps at 99999.9999 nm anyway.
WL_MIN_NM = 1.0e-4
WL_MAX_NM = 99999.0

REF_TAG = 'CMFG'

_IJ_RE = re.compile(r'^(\d+)-\s*(\d+)$')


def roman_to_stage(name):
    """'IV' -> 4 etc."""
    return ROMAN.get(name.upper())


def ion_from_atomic_path(path):
    """Derive (Z, stage) from $ATOMIC/<EL>/<ROMAN>/<date>/<file> layout.

    Returns (Z, stage) or None if the path does not match the layout.
    """
    parts = Path(path).parts
    for i, p in enumerate(parts):
        z = ELEMENT_DIR_Z.get(p.upper())
        if z is not None and i + 1 < len(parts):
            stage = roman_to_stage(parts[i + 1])
            if stage is not None:
                return z, stage
    return None


def parse_osc_data(path, z, stage):
    """Parse a CMFGEN osc_data file.

    Returns (transitions, stats) where transitions is a list of
    (wl_nm, loggf, j_l, e_l, j_u, e_u, name_l, name_u) and stats is a
    dict of counters for reporting.
    """
    levels = {}          # id -> (g, E_cm, name)
    n_levels = None
    n_trans = None
    stats = {'levels': 0, 'transitions_listed': 0, 'kept': 0,
             'dropped_dup': 0, 'dropped_f': 0, 'dropped_wl': 0,
             'dropped_level': 0, 'dropped_parse': 0}

    lines = Path(path).read_text(errors='replace').splitlines()

    # --- locate section sizes ---
    for idx, ln in enumerate(lines):
        if '!Number of energy levels' in ln:
            n_levels = int(ln.split('!')[0].strip())
        elif '!Number of transitions' in ln:
            n_trans = int(ln.split('!')[0].strip())
    if n_levels is None or n_trans is None:
        raise ValueError(f'{path}: could not find level/transition counts')

    # --- parse level block: first n_levels non-blank lines after the
    # '!Number of transitions' marker (level table follows the header
    # numerics in all observed variants) ---
    trans_marker = next(i for i, ln in enumerate(lines)
                        if '!Number of transitions' in ln)
    got = 0
    scanned = 0
    i = trans_marker + 1
    while i < len(lines) and got < n_levels and scanned < 2 * n_levels + 50:
        ln = lines[i].strip()
        i += 1
        scanned += 1
        if not ln:
            continue
        if ln.startswith('Transition') or ln.startswith('*'):
            # reached the transition table (or a separator) before the
            # declared level count — stop rather than misparse transitions
            break
        tok = ln.split()
        # transition rows with merged names ('low-up f A lam i-j ...')
        # must not be taken for level rows
        if re.search(r'\S-\S', tok[0]):
            continue
        # level line: name g E(cm^-1) ... ID ... [optional ARAD C4 C6]
        # ID is the first bare-integer token after tok[2] (may be
        # negative); trailing columns (ARAD, C4, C6) are floats.
        try:
            g = float(tok[1])
            e = float(tok[2])
        except (ValueError, IndexError):
            continue
        lev_id = None
        for t in tok[3:]:
            if re.fullmatch(r'-?\d+', t):
                lev_id = int(t)
                break
        if lev_id is None:
            continue
        # Transitions index levels by row number (1..n_levels); the ID
        # column's sign is metadata (negative marks e.g. collapsed
        # levels), so key by row index, not the ID value.
        levels[got + 1] = (g, e, tok[0])
        got += 1
    stats['levels'] = len(levels)

    # --- parse transitions: every line containing a ' i- j ' anchor ---
    # The i-j field is followed by a Trans.# column, a '|' column, or EOL
    # depending on the file variant.
    trans_re = re.compile(r'(-?\d+)-\s*(-?\d+)(?=\s*(?:\||\d|$))')
    seen = set()
    transitions = []
    for ln in lines:
        if '-' not in ln:
            continue
        m = trans_re.search(ln)
        if not m:
            continue
        i_lo, i_hi = int(m.group(1)), int(m.group(2))
        # fields immediately before the i-j anchor: f, A, Lam(A)
        pre = ln[:m.start()].split()
        if len(pre) < 3:
            stats['dropped_parse'] += 1
            continue
        try:
            lam_a = float(pre[-1])
            f_val = float(pre[-3])
        except ValueError:
            stats['dropped_parse'] += 1
            continue
        stats['transitions_listed'] += 1

        key = (i_lo, i_hi)
        if key in seen:
            stats['dropped_dup'] += 1
            continue
        seen.add(key)

        if i_lo not in levels or i_hi not in levels:
            stats['dropped_level'] += 1
            continue
        g_l, e_l, name_l = levels[i_lo]
        g_u, e_u, name_u = levels[i_hi]
        if e_l > e_u:  # ensure lower/upper ordering
            g_l, g_u = g_u, g_l
            e_l, e_u = e_u, e_l
            name_l, name_u = name_u, name_l

        if f_val <= 0.0:
            stats['dropped_f'] += 1
            continue
        wl_nm = lam_a / 10.0
        if not (WL_MIN_NM < wl_nm < WL_MAX_NM):
            stats['dropped_wl'] += 1
            continue

        loggf = __import__('math').log10(f_val * g_l)
        # J is unused by lineexpop; keep it valid for the f4.1 field.
        j_l = min((g_l - 1.0) / 2.0, 99.9)
        j_u = min((g_u - 1.0) / 2.0, 99.9)
        transitions.append((wl_nm, loggf, j_l, abs(e_l), j_u, abs(e_u),
                            name_l[:8], name_u[:8]))
        stats['kept'] += 1

    return transitions, stats


def format_record(wl, loggf, j_l, e_l, j_u, e_u, code, name_l, name_u,
                  ref=REF_TAG):
    """One STELLA lineatom record per rlnexp.f format 30."""
    return (f'{wl:10.4f} {ref:<4s} {loggf:7.3f}{j_l:4.1f}{e_l:11.3f}'
            f'{j_u:4.1f} {e_u:11.3f}{code:7.2f}{name_l:<8s}  {name_u:<8s}')


_CODE_RE = re.compile(r'\b(\d{1,2}\.\d{2})\b')

# Ion code field per rlnexp.f format 30 (0-based slice, cols 55-61).
_CODE_SLICE = slice(54, 61)


def _code_of(line):
    """Ion code (Z.(stage-1)) from a lineatom record, or None.

    Primary: the f7.2 field at fixed columns.  Fallback: the last
    Z.nn-shaped token on the line (some legacy rows deviate from the
    fixed layout).
    """
    if len(line) >= 61:
        try:
            return float(line[_CODE_SLICE])
        except ValueError:
            pass
    matches = _CODE_RE.findall(line)
    if matches:
        return float(matches[-1])
    return None


def iter_original_ions(path):
    """Yield (line_text, z, stage) for each data line of a lineatom.dat,
    skipping the i10 header.  (z, stage) may be None if the ion code
    cannot be located on the line."""
    with open(path, errors='replace') as fh:
        header = fh.readline()
        for ln in fh:
            ln = ln.rstrip('\n')
            if not ln.strip():
                continue
            code = _code_of(ln)
            z = stage = None
            if code is not None:
                z = int(code + 0.001)
                stage = int(100.0 * (code - z) + 0.1) + 1
            yield ln, z, stage


def extract_osc_paths_from_batch_sh(batch_sh, atomic_dir):
    """Return ordered [(z, stage, path)] for every *_F_OSCDAT link in
    a CMFGEN batch.sh, resolved against atomic_dir."""
    out = []
    txt = Path(batch_sh).read_text(errors='replace')
    # Any token that looks like $ATOMIC/... or an absolute path
    # under the atomic tree, on a line containing '_F_OSCDAT'.
    for ln in txt.splitlines():
        if '_F_OSCDAT' not in ln or 'ln' not in ln:
            continue
        if ln.lstrip().startswith('#'):
            continue  # commented-out link
        for tok in ln.split():
            tok = tok.strip('"\'')
            if tok.startswith('$ATOMIC/') or tok.startswith('${ATOMIC}/'):
                rel = tok.split('/', 1)[1]
                cand = Path(atomic_dir) / rel
            elif tok.startswith('/') and 'OSCDAT' not in tok:
                cand = Path(tok)
            else:
                continue
            ion = ion_from_atomic_path(tok)
            if ion is None:
                print(f'  WARNING: cannot derive ion from {tok}',
                      file=sys.stderr)
                continue
            out.append((ion[0], ion[1], cand))
    # dedup preserving order
    seen, res = set(), []
    for z, s, p in out:
        if (z, s) not in seen:
            seen.add((z, s))
            res.append((z, s, p))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--osc', action='append', default=[],
                    metavar='Z:STAGE:PATH',
                    help='CMFGEN osc_data file with explicit ion tag '
                         '(repeatable)')
    ap.add_argument('--from-batch-sh', metavar='BATCH_SH',
                    help='Extract osc_data paths from *_F_OSCDAT links in '
                         'a CMFGEN batch.sh')
    ap.add_argument('--atomic-dir', default=None,
                    help='$ATOMIC tree root (required with --from-batch-sh)')
    ap.add_argument('--merge-with', metavar='LINEATOM',
                    help='Original lineatom.dat; keep its lines for '
                         '(Z,stage) not covered by the converted set')
    ap.add_argument('--merge-mode', choices=['replace', 'only'],
                    default='replace',
                    help="replace (default): converted ions replace the "
                         "original list's same-ion lines; other ions kept. "
                         "only: pure CMFGEN-converted list.")
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    sources = []  # (z, stage, path)
    for spec in args.osc:
        try:
            z_s, stage_s, p = spec.split(':', 2)
            sources.append((int(z_s), int(stage_s), Path(p)))
        except ValueError:
            ap.error(f'bad --osc spec: {spec!r} (want Z:STAGE:PATH)')
    if args.from_batch_sh:
        if not args.atomic_dir:
            ap.error('--from-batch-sh requires --atomic-dir')
        sources.extend(extract_osc_paths_from_batch_sh(args.from_batch_sh,
                                                       args.atomic_dir))
    if not sources:
        ap.error('no input files: give --osc or --from-batch-sh')

    covered = set()
    converted = []
    provenance = {'sources': [], 'merge_with': args.merge_with,
                  'merge_mode': args.merge_mode, 'ref_tag': REF_TAG}
    for z, stage, path in sources:
        if not path.exists():
            print(f'  WARNING: osc_data not found: {path} '
                  f'(Z={z} stage={stage}) — skipped', file=sys.stderr)
            continue
        trans, stats = parse_osc_data(path, z, stage)
        code = z + (stage - 1) / 100.0
        for t in trans:
            converted.append((*t[:6], code, t[6], t[7]))
        covered.add((z, stage))
        provenance['sources'].append({
            'file': str(path), 'z': z, 'stage': stage, 'code': code,
            **stats})
        print(f'  Z={z:2d} stage={stage:2d}  {path.name:24s} '
              f'levels={stats["levels"]:5d} trans={stats["transitions_listed"]:7d} '
              f'kept={stats["kept"]:7d} dup={stats["dropped_dup"]:5d} '
              f'bad_f={stats["dropped_f"]:4d} wl={stats["dropped_wl"]:4d} '
              f'nolvl={stats["dropped_level"]:4d} parse={stats["dropped_parse"]:4d}')

    out_lines = [format_record(wl, loggf, j_l, e_l, j_u, e_u, code, nl, nu)
                 for (wl, loggf, j_l, e_l, j_u, e_u, code, nl, nu)
                 in converted]

    n_kept_orig = 0
    if args.merge_with and args.merge_mode == 'replace':
        kept_orig = []
        n_unparsed = 0
        for ln, z, stage in iter_original_ions(args.merge_with):
            if z is None:
                n_unparsed += 1
                kept_orig.append(ln)   # conservative: keep unclassifiable
            elif (z, stage) not in covered:
                kept_orig.append(ln)
            n_kept_orig += 1
        out_lines.extend(kept_orig)
        provenance['original_lines_kept'] = len(kept_orig)
        provenance['original_lines_unparsed'] = n_unparsed
        print(f'  merge: kept {len(kept_orig)} original lines '
              f'({n_unparsed} unclassifiable kept conservatively)')

    # sort by wavelength (first f10.4 field), matching the original file's
    # ascending-wavelength order that lineexpop's incremental log assumes
    def _wl_of(line):
        try:
            return float(line[:10])
        except ValueError:
            return float('inf')
    out_lines.sort(key=_wl_of)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w') as fh:
        fh.write(f'{len(out_lines):10d}\n')
        fh.write('\n'.join(out_lines) + '\n')
    prov_path = out.with_suffix(out.suffix + '.provenance.json')
    prov_path.write_text(json.dumps(provenance, indent=2))
    print(f'  wrote {out}: {len(out_lines)} lines '
          f'({len(converted)} converted), provenance -> {prov_path.name}')


if __name__ == '__main__':
    main()
