#!/usr/bin/env python3
"""Score SN_B double-peak model lightcurves against SN 2024pba targets.

Targets (from observed SN 2024pba):
  1. Early plateau / first peak ~ -15.5 to -16 mag (r-band pre-main-peak)
  2. Main ATLAS o-band peak -18.54 mag with double-peak morphology
  3. Slow post-peak tail (~+2 mag/100d decline)

For each model .tt file reports:
  - r-band first-peak magnitude and time (local max before main peak)
  - o-band global peak magnitude and time (brightest point)
  - dip depth between first and second peak (mag)
  - post-peak decline rate (mag / 100 d) from o-band peak to end

.tt column layout (24-col runs): cols 13-15 SDSS g/r/i, cols 20-21 ATLAS c/o.
Older 14-col runs: cols ~10-12 Johnson V/R/I (use MR col as r proxy).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

OBS_O_PEAK = -18.54
OBS_FIRST_PEAK = -15.6   # r-band plateau level before main rise
OBS_DECAY = 1.97         # mag/100d


def load_tt(path):
    """Return (t_days, mag_columns dict) from a STELLA .tt file."""
    rows = []
    for ln in open(path):
        if not ln.strip():
            continue
        parts = ln.split()
        try:
            t = float(parts[0])
            vals = [float(x) for x in parts[1:]]
        except (ValueError, IndexError):
            continue
        if len(vals) >= 12:
            rows.append((t, vals))
    if not rows:
        return None, {}
    t = np.array([r[0] for r in rows])
    n = len(rows[0][1])
    return t, {i: np.array([r[1][i] for r in rows]) for i in range(n)}


def band_cols(ncol):
    """Pick o-band and r-band column indices (0-based within mag block)."""
    # .tt cols (1-indexed incl. time): g=13,r=14,i=15 -> mag idx 11,12,13
    # ATLAS c=20,o=21 -> mag idx 18,19. Johnson V=10,R=11 -> mag idx 9,10.
    if ncol >= 21:
        return 19, 12      # o-band, r-band
    if ncol >= 13:
        return 13, 12      # i-band as o proxy? use r
    return ncol - 1, ncol - 2


def local_extrema(t, mag, t_lo, t_hi):
    """Brightest (most negative) mag in window -> (t*, mag*)."""
    m = (t >= t_lo) & (t <= t_hi)
    if not m.any():
        return None, None
    i = np.argmin(mag[m])
    return t[m][i], mag[m][i]


def score_model(path):
    t, cols = load_tt(path)
    if t is None:
        return None
    nmag = len(cols)
    o_idx, r_idx = band_cols(nmag)
    o = cols.get(o_idx)
    r = cols.get(r_idx)
    if o is None or r is None:
        return None

    # main (brightest) peak
    i_pk = int(np.argmin(o))
    t_pk, m_pk = float(t[i_pk]), float(o[i_pk])

    # first peak: local minimum in mag before main peak (search t < t_pk)
    # restrict to t in [0, t_pk]
    pre = t < t_pk
    m_first = t_first = None
    if pre.any():
        # first local brightening: take brightest r-band point before t_pk*0.7
        early = (t < t_pk * 0.7)
        if early.any():
            i_f = np.argmin(r[early])
            t_first, m_first = float(t[early][i_f]), float(r[early][i_f])

    # dip between first peak and second peak (faintest mag between them)
    dip = None
    if t_first is not None:
        mid = (t > t_first) & (t < t_pk)
        if mid.any():
            i_d = np.argmax(o[mid])
            dip = float(m_first - o[mid][i_d]) if m_first is not None else None
            # dip depth = second_peak - dip_level
            dip = float(m_pk - o[mid][i_d])

    # post-peak decline: slope from peak to last point (mag/100d)
    decay = None
    if t[-1] > t_pk + 5:
        decay = float((o[-1] - m_pk) / (t[-1] - t_pk) * 100.0)

    return {
        't_first': t_first, 'm_first': m_first,
        't_peak': t_pk, 'm_peak': m_pk,
        'dip_mag': dip, 'decay': decay,
        't_end': float(t[-1]),
        'err_peak': m_pk - OBS_O_PEAK,
        'err_first': (m_first - OBS_FIRST_PEAK) if m_first else None,
        'err_decay': (decay - OBS_DECAY) if decay else None,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dirs', nargs='*')
    args = ap.parse_args()
    print(f"{'model':44s} {'m1st':>7s} {'dip':>6s} {'mPk':>7s} {'dPk':>6s} {'dec':>6s}")
    print('-' * 90)
    for d in args.dirs:
        d = Path(d)
        tts = list(d.glob('stella/*.tt')) or list(d.glob('**/*.tt'))
        if not tts:
            print(f"{d.name:44s}  no .tt")
            continue
        s = score_model(tts[0])
        if s is None:
            print(f"{d.name:44s}  parse fail")
            continue
        mf = f"{s['m_first']:.2f}" if s['m_first'] is not None else "  --"
        dp = f"{s['dip_mag']:.2f}" if s['dip_mag'] is not None else "  --"
        dc = f"{s['decay']:.2f}" if s['decay'] is not None else "  --"
        print(f"{d.name:44s} {mf:>7s} {dp:>6s} {s['m_peak']:7.2f} {s['err_peak']:+6.2f} {dc:>6s}")


if __name__ == '__main__':
    main()
