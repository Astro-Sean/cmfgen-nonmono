#!/usr/bin/env python3
"""Convert a converged CMFGEN run's POP* files into per-ion <ION>_IN files.

The _IN files let a cold start seed populations through departure
coefficients (LTE_EST=F, LIN_INT=T) while T and ED come from T_IN /
GAMMAS_IN.  This is how to restart at the hydro temperature when the
checkpoint carries a previously fitted T: the POP files store XzV_F
(linear departure coefficients) which are temperature independent; CMFGEN
rebuilds the level populations as DC x LTE(T_IN) at start-up.

Format written here matches prepare_PROXY_POPULATIONS/utils/in_files.py:
a '!Format date' header, one 'RP_OLD 0 NOLD ND' line, then per depth
(outer to inner) a structure row 'R DI ED T 0 V CLUMP idx' followed by
NOLD departure coefficients, five per line.
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

# DI=0 kills LTEPOP_WLD_V2; DC=0 is fine to write (the reader takes LOG
# only when the file predates the log-format date, and ours is new), but
# floor both anyway so downstream LOG interpolations stay finite.
DI_FLOOR = 1e-30
DC_FLOOR = 1e-30
RVTJ_FLOATS_PER_LINE = 8


def read_nd_floats(lines, start, count):
    """Read `count` floats starting at lines[start], returning (arr, next_line)."""
    vals = []
    i = start
    while len(vals) < count and i < len(lines):
        vals.extend(float(x) for x in lines[i].split())
        i += 1
    if len(vals) < count:
        raise ValueError(f"expected {count} floats at line {start}, got {len(vals)}")
    return np.array(vals[:count]), i


def parse_rvtj(path, nd):
    """Extract R(1e10 cm), V(km/s), ED(cm^-3), T(1e4 K) sections from RVTJ."""
    lines = path.read_text().splitlines()
    sections = {}
    labels = {
        "radius": "Radius (10^10 cm)",
        "velocity": "Velocity (km/s)",
        "ed": "Electron density",
        "t": "Temperature (10^4K)",
    }
    for idx, line in enumerate(lines):
        label = line.strip()
        for key, text in labels.items():
            if key not in sections and label == text:
                arr, _ = read_nd_floats(lines, idx + 1, nd)
                sections[key] = arr
    missing = [k for k in labels if k not in sections]
    if missing:
        raise ValueError(f"{path}: missing sections {missing}")
    return sections


def parse_pop_file(path, nd):
    """Yield (ion_name, dc[nlev, nd], di[nd]) blocks from a POP<elem> file."""
    lines = path.read_text().splitlines()
    blocks = []
    hdr_re = re.compile(r"Number of (\S+) levels:\s*(\d+)")
    i = 0
    while i < len(lines):
        m = hdr_re.search(lines[i])
        if not m:
            i += 1
            continue
        ion, nlev = m.group(1), int(m.group(2))
        # RITE_ASC writes only the header line when CI_PRES is false.
        if nlev == 0:
            i += 1
            continue
        # Next line is the oscillator date; the DC block follows it.
        flat, i = read_nd_floats(lines, i + 2, nlev * nd + nd)
        dc = flat[: nlev * nd].reshape(nd, nlev).T  # Fortran column-major
        di = flat[nlev * nd :]
        blocks.append((ion, dc, di))
    return blocks


def write_in_file(path, ion, dc, di, rvtj, nd):
    """Emit one <ION>_IN file in the format REGRID_LOG_DC_V1 consumes."""
    di = np.where(np.isfinite(di) & (di > 0), di, DI_FLOOR)
    dc = np.where(np.isfinite(dc), np.maximum(dc, DC_FLOOR), 1.0)
    with open(path, "w") as f:
        f.write("\n")
        f.write(
            f" {datetime.now().strftime('%d-%b-%Y').upper()}"
            "                           !Format date\n"
        )
        f.write("\n")
        nlev = dc.shape[0]
        f.write(f" {rvtj['radius'][nd - 1]:14.8E}     {0.0:12.4E}    {nlev:4d}    {nd:4d}\n")
        for i in range(nd):
            f.write("\n")
            f.write(
                f" {rvtj['radius'][i]:17.10E}{di[i]:17.8E}{rvtj['ed'][i]:17.8E}"
                f"{rvtj['t'][i]:17.8E}{0.0:17.8E}{rvtj['velocity'][i]:17.8E}"
                f"{1.0:17.8E}  {i + 1:4d}\n"
            )
            for j0 in range(0, nlev, 5):
                line = " "
                for j in range(j0, min(j0 + 5, nlev)):
                    line += f"{dc[j, i]:17.7E}"
                f.write(line + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir", required=True, help="dir containing POP* and RVTJ")
    ap.add_argument("--out-dir", required=True, help="where to write <ION>_IN files")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rvtj_path = run_dir / "RVTJ"
    m = re.search(r"ND:\s*(\d+)", rvtj_path.read_text()[:2000])
    if not m:
        sys.exit("RVTJ: could not find ND")
    nd = int(m.group(1))
    rvtj = parse_rvtj(rvtj_path, nd)

    n_ions = 0
    for pop in sorted(run_dir.glob("POP*")):
        if pop.suffix or not pop.is_file():
            continue
        for ion, dc, di in parse_pop_file(pop, nd):
            write_in_file(out_dir / f"{ion}_IN", ion, dc, di, rvtj, nd)
            n_ions += 1
            print(f"  {ion}_IN: nlev={dc.shape[0]} "
                  f"dc[{np.min(dc):.2e},{np.max(dc):.2e}] "
                  f"di[{np.min(di):.2e},{np.max(di):.2e}]")
    print(f"Wrote {n_ions} _IN files to {out_dir} (ND={nd})")


if __name__ == "__main__":
    main()
