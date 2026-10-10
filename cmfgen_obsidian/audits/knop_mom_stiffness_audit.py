#!/usr/bin/env python3
"""Reproduce the stiffness diagnosis for the dev7 MOM_KNOP_REFRESH_V1 run.

Inputs:
  --hydro   path to SN_HYDRO_DATA (70-point source profile)
  --errors  path to MOM_J_NONMONO_ERRORS (per-cell diagnostic dump)

Outputs (stdout):
  - standard-grid and expanded-grid sigma ranges and locations
  - mapping of the worst-conditioned expanded depths to R, V, parent
  - evidence from the error dump that the ill-conditioned cells sit at
    the zero-GAM resonance 1 + sigma*K_on_J ~= 0.
"""
import argparse
import re
from pathlib import Path
import numpy as np


def _read_block(lines, i):
    vals = []
    for line in lines[i + 1:]:
        s = line.strip()
        if not s or s[0].isalpha():
            break
        vals.extend(map(float, s.split()))
    return np.array(vals)


def read_hydro(path):
    text = Path(path).read_text().splitlines()
    data = {}
    for i, line in enumerate(text):
        if 'Radius grid' in line:
            data['R'] = _read_block(text, i)
        elif line.strip().startswith('Velocity (km/s)'):
            data['V'] = _read_block(text, i)
        elif line.strip().startswith('Sigma (dlnV/dlnr-1)'):
            data['sigma_file'] = _read_block(text, i)
    return data


def sigma_velocity(r, v):
    """VEL_LAW=12 finite-difference sigma used by CMFGEN."""
    n = len(r)
    s = np.zeros(n)
    # boundaries one-sided
    s[0] = (r[0] / v[0]) * (v[0] - v[1]) / (r[0] - r[1]) - 1.0
    s[-1] = (r[-1] / v[-1]) * (v[-2] - v[-1]) / (r[-2] - r[-1]) - 1.0
    for i in range(1, n - 1):
        s[i] = (r[i] / v[i]) * (v[i - 1] - v[i + 1]) / (r[i - 1] - r[i + 1]) - 1.0
    return s


def build_expanded_grid(r, v, vdop_frac=1.0, t2_override=None, max_ins=256):
    """Replicate MOM_J_CMF_NONMONO_V1 grid expansion."""
    # VDOP_VEC is not in SN_HYDRO_DATA; use a typical T2 if not given.
    if t2_override is None:
        t2 = vdop_frac * np.min(np.abs(v))
    else:
        t2 = t2_override
    inserts = []
    for i in range(len(r) - 1):
        it = int(abs(v[i] - v[i + 1]) / t2)
        if it > max_ins:
            it = max_ins
        if it % 2 != 0:
            it += 1
        inserts.append(it)
    r_nm = [r[0]]
    pnt = [1]
    for i in range(len(r) - 1):
        it = inserts[i]
        if it > 0:
            dr = (r[i + 1] - r[i]) / (it + 1)
            for _ in range(it):
                r_nm.append(r_nm[-1] + dr)
                pnt.append(i + 1)
        r_nm.append(r[i + 1])
        pnt.append(i + 1)
    return np.array(r_nm), np.array(pnt, dtype=int), inserts, t2


def steffen_coeffs(r, y):
    """MON_INT_FUNS_V2: Steffen monotonic cubic coefficients."""
    n = len(r)
    h = np.diff(r)
    s = np.diff(y) / h
    d = np.empty(n)
    d[0] = s[0] + (s[0] - s[1]) * h[0] / (h[0] + h[1])
    for i in range(1, n - 1):
        d[i] = (s[i - 1] * h[i] + s[i] * h[i - 1]) / (h[i - 1] + h[i])
    d[-1] = s[-1] + (s[-1] - s[-2]) * h[-1] / (h[-2] + h[-1])
    d[0] = (np.sign(s[0]) + np.sign(d[0])) * min(abs(s[0]), 0.5 * abs(d[0]))
    for i in range(1, n - 1):
        d[i] = (np.sign(s[i - 1]) + np.sign(s[i])) * min(
            abs(s[i - 1]), abs(s[i]), 0.5 * abs(d[i])
        )
    d[-1] = (np.sign(s[-1]) + np.sign(d[-1])) * min(
        abs(s[-1]), 0.5 * abs(d[-1])
    )
    coef = np.zeros((n, 4))
    for i in range(n - 1):
        coef[i, 0] = (d[i] + d[i + 1] - 2.0 * s[i]) / (h[i] * h[i])
        coef[i, 1] = (3.0 * s[i] - 2.0 * d[i] - d[i + 1]) / h[i]
        coef[i, 2] = d[i]
        coef[i, 3] = y[i]
    coef[-1, 3] = y[-1]
    coef[-1, 2] = d[-1]
    return coef


def eval_steffen(r, coef, x):
    if x >= r[0]:
        i = 0
    elif x <= r[-1]:
        i = len(r) - 2
    else:
        i = np.searchsorted(-r, -x, side='right') - 1
    t = x - r[i]
    c = coef[i]
    return ((c[0] * t + c[1]) * t + c[2]) * t + c[3]


def parse_error_dump(path):
    """Read MOM_J_NONMONO_ERRORS and return list of dicts for dumped cells."""
    rows = []
    if not Path(path).exists():
        return rows
    text = Path(path).read_text()
    for m in re.finditer(
        r'^\s*(\d+)\s+([\d\.E\+-]+)\s+([\d\.E\+-]+)\s+([\d\.E\+-]+)\s+'
        r'([\d\.E\+-]+)\s+([\d\.E\+-]+)\s+([\d\.E\+-]+)\s+([\d\.E\+-]+)\s+'
        r'([\d\.E\+-]+)',
        text, re.MULTILINE
    ):
        groups = m.groups()
        rows.append({
            'depth': int(groups[0]),
            'freq': float(groups[1]),
            'xm': float(groups[2]),
            'eta': float(groups[3]),
            'chi': float(groups[4]),
            'esec': float(groups[5]),
            'k_on_j': float(groups[6]),
            'nmid_on_j': float(groups[7]),
            'nmid_on_hmid': float(groups[8]),
        })
    # Append the sigma_true line that follows each row.
    sigma_lines = re.findall(
        r'sigma_true,gam,chi,chi_hat=\s+([\d\.E\+-]+)\s+([\d\.E\+-]+)\s+'
        r'([\d\.E\+-]+)\s+([\d\.E\+-]+)',
        text
    )
    for row, sg in zip(rows, sigma_lines):
        row['sigma_true'] = float(sg[0])
        row['gam'] = float(sg[1])
        row['chi_hat'] = float(sg[3])
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hydro', required=True)
    parser.add_argument('--errors', required=True)
    parser.add_argument('--t2', type=float, default=50.0,
                        help='Grid-expansion velocity scale (km/s); '
                             'OUTGEN reports 5.0E+01 for the M15 run.')
    args = parser.parse_args()

    h = read_hydro(args.hydro)
    R, V = h['R'], h['V']
    SIG = sigma_velocity(R, V)
    print('Standard grid (hydro profile)')
    print(f'  ND = {len(R)}')
    print(f'  R range  = {R.min():.3f} -- {R.max():.3f}  (1e10 cm)')
    print(f'  V range  = {V.min():.3f} -- {V.max():.3f}  (km/s)')
    print(f'  sigma min/max = {SIG.min():.3f} / {SIG.max():.3f}')
    imin = int(SIG.argmin()) + 1
    print(f'  sigma min at standard depth {imin}, R={R[imin-1]:.3f}, V={V[imin-1]:.3f}')

    R_nm, Pnt, inserts, t2 = build_expanded_grid(R, V, t2_override=args.t2)
    print(f'\nExpanded grid: {len(R)} -> {len(R_nm)}  (T2={t2:.2f} km/s)')
    print(f'  Insertions per interval (first 40): {inserts[:40]}')

    c_v = steffen_coeffs(R, V)
    c_s = steffen_coeffs(R, SIG)
    V_nm = np.array([eval_steffen(R, c_v, x) for x in R_nm])
    SIG_nm = np.array([eval_steffen(R, c_s, x) for x in R_nm])
    print(f'  Expanded V range   = {V_nm.min():.3f} -- {V_nm.max():.3f}')
    print(f'  Expanded sigma min/max = {SIG_nm.min():.3f} / {SIG_nm.max():.3f}')
    imin_nm = int(SIG_nm.argmin()) + 1
    print(f'  Expanded sigma min at depth {imin_nm}, R={R_nm[imin_nm-1]:.3f}')

    print('\nMapping of OUTGEN worst-conditioned depths (22, 33) and the')
    print('shock-transition band (103-116):')
    for d in [22, 33, 100, 103, 110, 115, 116, 180, 190, 200, 210, 212]:
        if d <= len(R_nm):
            print(f'  depth {d:3d}: R={R_nm[d-1]:10.3f}  V={V_nm[d-1]:9.3f}  '
                  f'sigma={SIG_nm[d-1]:8.3f}  parent={Pnt[d-1]:2d}')

    rows = parse_error_dump(args.errors)
    if rows:
        print(f'\nParsed {len(rows)} dumped negative-J cells from {args.errors}')
        print('Cells in the shock-transition band (depth 100-120) closest to')
        print('the zero-GAM resonance 1 + sigma*K_on_J = 0:')
        band = [r for r in rows if 100 <= r['depth'] <= 120]
        band.sort(key=lambda r: abs(1.0 + r.get('sigma_true', np.nan) * r['k_on_j']))
        for row in band[:10]:
            s = row.get('sigma_true', np.nan)
            k = row['k_on_j']
            gam_fac = 1.0 + s * k
            print(f"  d={row['depth']:3d} f={row['freq']:.4e}  "
                  f"K_on_J={k:.4f}  sigma={s:7.3f}  "
                  f"1+sigma*K={gam_fac:9.4f}  chi/chi_hat="
                  f"{row['chi']/row.get('chi_hat',np.nan):.3e}")
        print('\nAll dumped cells with |1+sigma*K_on_J| < 0.1:')
        near = [r for r in rows
                if abs(1.0 + r.get('sigma_true', np.nan) * r['k_on_j']) < 0.1]
        near.sort(key=lambda r: abs(1.0 + r.get('sigma_true', np.nan) * r['k_on_j']))
        for row in near[:10]:
            s = row.get('sigma_true', np.nan)
            k = row['k_on_j']
            gam_fac = 1.0 + s * k
            print(f"  d={row['depth']:3d} f={row['freq']:.4e}  "
                  f"K_on_J={k:.4f}  sigma={s:7.3f}  "
                  f"1+sigma*K={gam_fac:9.4f}")
    else:
        print(f'\nNo error dump rows parsed from {args.errors}')


if __name__ == '__main__':
    main()
