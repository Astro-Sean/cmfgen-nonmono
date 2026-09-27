#!/usr/bin/env python3
"""Compare SN94W_IIn CMFGEN spectrum to Dessart+ (2009) day-49.5 reference.

Reference: observational_data/sn94w_reference/e1_v2_s2_7_Bs3_numes40_vt10.fl
(+ _cont.fl continuum). Luc's model is parameterized (rho~r^-10, V~830,
T_phot~7500, L~4.7e9). This is a morphology benchmark, not an exact match.

Usage:
    python scripts/compare_sn94w.py <model_spectrum> [--ref REF.fl] \
        [--ref-cont REF_cont.fl] [--out plot.png] [--ha] 
"""
import argparse
import sys

import numpy as np


C_KMS = 2.99792458e5
HA = 6562.8


def load_spec(path):
    """Load a 2-col wavelength(Angstrom)/flux spectrum."""
    d = np.loadtxt(path, comments='#')
    if d.ndim != 2 or d.shape[1] < 2:
        raise ValueError(f"{path}: expected >=2 columns")
    return d[:, 0], d[:, 1]


def profile_metrics(lam, fl, lam_c=None, fl_c=None):
    """Narrow-core FWHM and broad-wing extent around Halpha (km/s)."""
    v = (lam - HA) / HA * C_KMS
    m = np.abs(v) < 25000.0
    vv, ff = v[m], fl[m]
    if vv.size < 10:
        return {}
    pk = ff.max()
    out = {'peak_flux': float(pk)}
    # excess above continuum (if provided) else above local minimum
    if lam_c is not None and fl_c is not None:
        cc = np.interp(lam, lam_c, fl_c)[m]
    else:
        cc = np.full_like(ff, np.percentile(ff, 20))
    ex = ff - cc
    # narrow core FWHM on the excess
    half = ex.max() * 0.5
    iw = np.where(ex > half)[0]
    if iw.size:
        out['core_fwhm_kms'] = float(vv[iw].max() - vv[iw].min())
    # broad-wing extent: velocity where excess > 1% and 5% of peak
    for frac, key in [(0.05, 'wing_5pct_kms'), (0.01, 'wing_1pct_kms')]:
        iw = np.where(ex > ex.max() * frac)[0]
        if iw.size:
            out[key] = float(vv[iw].max() - vv[iw].min())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('model', help='model spectrum (wavelength Angstrom, flux)')
    ap.add_argument('--ref', default='observational_data/sn94w_reference/'
                    'e1_v2_s2_7_Bs3_numes40_vt10.fl')
    ap.add_argument('--ref-cont', default='observational_data/sn94w_reference/'
                    'e1_v2_s2_7_Bs3_numes40_vt10_cont.fl')
    ap.add_argument('--model-cont', default=None,
                    help='model continuum spectrum (optional)')
    ap.add_argument('--out', default=None, help='output plot path')
    ap.add_argument('--ha', action='store_true',
                    help='also zoom the Halpha profile')
    args = ap.parse_args()

    rl, rf = load_spec(args.ref)
    rcl, rcf = load_spec(args.ref_cont)
    ml, mf = load_spec(args.model)
    mcl = mcf = None
    if args.model_cont:
        mcl, mcf = load_spec(args.model_cont)

    # --- metrics ---
    print('=== Halpha morphology (velocity frame, km/s) ===')
    rm = profile_metrics(rl, rf, rcl, rcf)
    mm = profile_metrics(ml, mf, mcl, mcf)
    print(f"{'metric':>16}  {'Luc(d49.5)':>12}  {'model':>12}")
    for k in ['core_fwhm_kms', 'wing_5pct_kms', 'wing_1pct_kms', 'peak_flux']:
        rv = f"{rm.get(k, float('nan')):.3g}" if k in rm else 'n/a'
        mv = f"{mm.get(k, float('nan')):.3g}" if k in mm else 'n/a'
        print(f"{k:>16}  {rv:>12}  {mv:>12}")

    # --- continuum level comparison near Ha and bolometric scaling ---
    rcont = np.interp(6563, rcl, rcf)
    print(f"\nLuc continuum F(6563) = {rcont:.3e}")
    print(f"model peak F(6563)    = {np.interp(6563, ml, mf):.3e}")

    # --- plot ---
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f'(matplotlib unavailable: {e})')
        return
    fig, ax = plt.subplots(1, 2 if args.ha else 1, figsize=(13 if args.ha else 8, 5))
    if not args.ha:
        ax = [ax]
    a = ax[0]
    a.plot(rl, rf, lw=0.7, color='k', label='Luc d49.5 (e1_v2_s2_7)')
    a.plot(rcl, rcf, lw=0.7, color='gray', ls='--', label='Luc cont')
    a.plot(ml, mf, lw=0.8, color='tab:red', label='SN94W_IIn model')
    a.set_xlim(3000, 10000)
    a.set_xlabel('wavelength [A]')
    a.set_ylabel('flux')
    a.set_title('SN1994W day-49.5 benchmark')
    a.legend(fontsize=7)
    if args.ha:
        b = ax[1]
        vr = (rl - HA) / HA * C_KMS
        vm = (ml - HA) / HA * C_KMS
        b.plot(vr, rf, lw=0.8, color='k', label='Luc d49.5')
        b.plot(vr, rcf * 0 + np.interp(rl, rcl, rcf), lw=0.7, color='gray',
               ls='--', label='Luc cont')
        b.plot(vm, mf, lw=0.9, color='tab:red', label='model')
        b.set_xlim(-20000, 20000)
        m = (vr > -20000) & (vr < 20000)
        b.set_ylim(0, rf[m].max() * 1.15)
        b.set_xlabel('velocity [km/s]')
        b.set_title('Halpha profile: narrow core + broad es wings')
        b.legend(fontsize=7)
    fig.tight_layout()
    out = args.out or 'sn94w_compare.png'
    fig.savefig(out, dpi=150)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
