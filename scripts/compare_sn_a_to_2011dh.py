#!/usr/bin/env python3
"""Compare SN_A grid model lightcurves with actual SN 2011dh observations
(Ergon et al. 2015, J/A+A/580/A142).

Uses the real photometric (tablec4.dat: UBVRI) and bolometric (tablec9.dat)
data rather than hardcoded approximations.  Computes a weighted chi-squared
for V-band, B-V colour, U-V colour and bolometric magnitude, then ranks all
completed SN_A models.

Usage:
    python scripts/compare_sn_a_to_2011dh.py [--plot] [--top N]
"""
import os, sys, glob, argparse
import numpy as np
import pandas as pd

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
OBS_DIR = f"{BASE}/observational_data/SN2011dh"

# SN 2011dh parameters (Ergon et al. 2014/2015, Bersten et al. 2012)
DISTANCE_MPC = 8.4          # Mpc
MU = 5.0 * np.log10(DISTANCE_MPC * 1e6 / 10.0)   # distance modulus ~29.62
EBV = 0.04                   # Milky Way reddening
RV = 3.1
AV = EBV * RV                # ~0.124
# Band extinction ratios (Cardelli et al. 1989, RV=3.1)
A_BAND = {
    'U': EBV * 4.1,
    'B': EBV * 4.1,
    'V': AV,
    'R': EBV * 2.5,
    'I': EBV * 1.8,
}
# Bolometric: M_bol_sun
M_BOL_SUN = 4.74
L_SUN = 3.828e33  # erg/s


def load_2011dh_ubvri():
    """Parse tablec4.dat (fixed-width UBVRI photometry)."""
    path = f"{OBS_DIR}/tablec4.dat"
    records = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if len(line) < 16:
                continue
            phase = float(line[10:16].strip())
            vals = {}
            for band, s0, s1 in [('U', 18, 32), ('B', 34, 48),
                                 ('V', 50, 64), ('R', 66, 80),
                                 ('I', 82, 96)]:
                mag_str = line[s0:s0+6].strip()
                err_str = line[s0+6:s1].strip()
                if mag_str:
                    vals[band] = float(mag_str)
                    vals[f'e_{band}'] = float(err_str) if err_str else 0.02
            if 'V' in vals:
                rec = {'phase': phase}
                rec.update(vals)
                records.append(rec)
    df = pd.DataFrame(records)
    # Convert apparent -> absolute (dereddened)
    for band in ['U', 'B', 'V', 'R', 'I']:
        if band in df.columns:
            df[f'M_{band}'] = df[band] - MU - A_BAND[band]
    if 'B' in df.columns and 'V' in df.columns:
        df['MB_MV'] = df['M_B'] - df['M_V']
    if 'U' in df.columns and 'V' in df.columns:
        df['MU_MV'] = df['M_U'] - df['M_V']
    return df


def load_2011dh_bol():
    """Parse tablec9.dat (bolometric luminosity in 10^41 erg/s)."""
    path = f"{OBS_DIR}/tablec9.dat"
    records = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if len(line) < 16:
                continue
            phase = float(line[10:16].strip())
            lbol_str = line[18:24].strip()
            elbol_str = line[26:32].strip()
            if lbol_str:
                lbol = float(lbol_str) * 1e41  # erg/s
                elbol = float(elbol_str) * 1e41 if elbol_str else 0.05e41
                m_bol = -2.5 * np.log10(lbol / L_SUN) + M_BOL_SUN
                e_mbol = 2.5 * elbol / (lbol * np.log(10))
                records.append({'phase': phase, 'lbol': lbol,
                                'M_bol': m_bol, 'e_M_bol': e_mbol})
    return pd.DataFrame(records)


def load_model(model_dir):
    """Load model lightcurve CSV."""
    name = os.path.basename(model_dir)
    csv = f"{model_dir}/output/{name}_lightcurve_complete.csv"
    if not os.path.exists(csv):
        return None
    df = pd.read_csv(csv)
    df = df[df['bol_mag'] > -30].copy()
    df = df[df['time'] >= 0].copy()
    return df


def interp_model(df, t_obs, col):
    """Interpolate model column at observed phases."""
    t = df['time'].values
    v = df[col].values
    return np.interp(t_obs, t, v)


def compute_chisq(model_df, obs_ubvri, obs_bol,
                  t_min=7.0, t_max=45.0):
    """Compute weighted RMS deviation comparing model to observations.

    Uses RMS (root-mean-square) in magnitudes rather than formal chi-squared,
    because the tiny observational errors (~0.01 mag) would make chi-squared
    dominated by systematic model-observation offsets rather than the
    scientifically meaningful shape/peak comparison.  The SBO flash (t<7 d)
    is excluded because SN 2011dh was not observed at shock breakout.
    """
    # V-band
    obs_v = obs_ubvri[(obs_ubvri['phase'] >= t_min) &
                      (obs_ubvri['phase'] <= t_max)].copy()
    mv = interp_model(model_df, obs_v['phase'].values, 'V')
    dv = mv - obs_v['M_V'].values
    rms_v = np.sqrt(np.mean(dv ** 2)) if len(dv) > 0 else np.nan
    ndof_v = len(obs_v)

    # B-V colour
    has_bv = 'B' in obs_ubvri.columns and obs_ubvri['B'].notna().any()
    if has_bv:
        obs_bv = obs_ubvri[(obs_ubvri['phase'] >= t_min) &
                           (obs_ubvri['phase'] <= t_max)].copy()
        obs_bv = obs_bv.dropna(subset=['MB_MV'])
        if len(obs_bv) > 0:
            mbv_model = (interp_model(model_df, obs_bv['phase'].values, 'B') -
                         interp_model(model_df, obs_bv['phase'].values, 'V'))
            dbv = mbv_model - obs_bv['MB_MV'].values
            rms_bv = np.sqrt(np.mean(dbv ** 2))
            ndof_bv = len(obs_bv)
        else:
            rms_bv, ndof_bv = np.nan, 0
    else:
        rms_bv, ndof_bv = np.nan, 0

    # U-V colour
    has_uv = 'U' in obs_ubvri.columns and obs_ubvri['U'].notna().any()
    if has_uv:
        obs_uv = obs_ubvri[(obs_ubvri['phase'] >= t_min) &
                           (obs_ubvri['phase'] <= t_max)].copy()
        obs_uv = obs_uv.dropna(subset=['MU_MV'])
        if len(obs_uv) > 0:
            muv_model = (interp_model(model_df, obs_uv['phase'].values, 'U') -
                         interp_model(model_df, obs_uv['phase'].values, 'V'))
            duv = muv_model - obs_uv['MU_MV'].values
            rms_uv = np.sqrt(np.mean(duv ** 2))
            ndof_uv = len(obs_uv)
        else:
            rms_uv, ndof_uv = np.nan, 0
    else:
        rms_uv, ndof_uv = np.nan, 0

    # Bolometric
    obs_b = obs_bol[(obs_bol['phase'] >= t_min) &
                    (obs_bol['phase'] <= t_max)].copy()
    mbol_model = interp_model(model_df, obs_b['phase'].values, 'bol_mag')
    dbol = mbol_model - obs_b['M_bol'].values
    rms_bol = np.sqrt(np.mean(dbol ** 2)) if len(dbol) > 0 else np.nan
    ndof_bol = len(obs_b)

    total_ndof = ndof_v + ndof_bv + ndof_uv + ndof_bol
    if total_ndof == 0:
        return np.nan, {}, 0

    # Weighted combination of RMS deviations (in mag):
    # Bolometric and V get higher weight (most diagnostic for IIb)
    w_v, w_bv, w_uv, w_bol = 2.0, 1.5, 1.0, 3.0
    parts = [(rms_v, w_v), (rms_bv, w_bv),
             (rms_uv, w_uv), (rms_bol, w_bol)]
    w_sum = sum(w for r, w in parts if not np.isnan(r))
    combined = sum(r * w for r, w in parts if not np.isnan(r)) / w_sum

    details = {
        'rcs_v': rms_v, 'rcs_bv': rms_bv, 'rcs_uv': rms_uv,
        'rcs_bol': rms_bol,
        'ndof_v': ndof_v, 'ndof_bv': ndof_bv,
        'ndof_uv': ndof_uv, 'ndof_bol': ndof_bol,
    }
    return combined, details, total_ndof


def get_model_params(model_dir):
    """Extract key parameters from config for reporting."""
    cfg = f"{model_dir}/config.yaml"
    if not os.path.exists(cfg):
        return {}
    params = {}
    with open(cfg) as f:
        text = f.read()
    import re
    for key, label in [('explosion_energy', 'E_burst'),
                       ('h_ext_mass_msun', 'M_env'),
                       ('h_ext_r_outer_rsun', 'R_env'),
                       ('h_ext_T_base_K', 'T_base'),
                       ('h_ext_n_rho', 'n_rho'),
                       ('h_ext_fraction', 'X_H')]:
        m = re.search(rf'{key}:\s*([\d.]+)', text)
        if m:
            params[label] = float(m.group(1))
    # Ni56 mass
    m = re.search(r'nickel56:.*?total_mass_msun:\s*([\d.]+)', text, re.DOTALL)
    if m:
        params['M_Ni'] = float(m.group(1))
    return params


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot', action='store_true',
                        help='Generate comparison plot')
    parser.add_argument('--top', type=int, default=5,
                        help='Number of top models to show in detail')
    parser.add_argument('--t-max', type=float, default=45.0,
                        help='Max phase (days) for comparison')
    parser.add_argument('--t-min', type=float, default=7.0,
                        help='Min phase (days) for comparison (excludes SBO)')
    args = parser.parse_args()

    print("Loading SN 2011dh observations (Ergon et al. 2015)...")
    print(f"  Distance: {DISTANCE_MPC} Mpc (mu={MU:.2f}), "
          f"E(B-V)={EBV}, A_V={AV:.3f}")
    obs_ubvri = load_2011dh_ubvri()
    obs_bol = load_2011dh_bol()
    print(f"  UBVRI points: {len(obs_ubvri)}, "
          f"Bolometric points: {len(obs_bol)}")
    print(f"  V range: {obs_ubvri['phase'].min():.1f} - "
          f"{obs_ubvri['phase'].max():.1f} d")
    print(f"  Bol range: {obs_bol['phase'].min():.1f} - "
          f"{obs_bol['phase'].max():.1f} d")

    # Find all completed SN_A grid models
    csvs = sorted(glob.glob(f"{BASE}/SN_A_g_*/output/*_lightcurve_complete.csv"))
    print(f"\nFound {len(csvs)} completed SN_A grid models\n")

    results = []
    for csv_path in csvs:
        model_name = os.path.basename(
            os.path.dirname(os.path.dirname(csv_path)))
        model_dir = f"{BASE}/{model_name}"
        df = load_model(model_dir)
        if df is None or len(df) < 10:
            continue

        score, details, ndof = compute_chisq(
            df, obs_ubvri, obs_bol, t_min=args.t_min, t_max=args.t_max)

        # Peak properties
        peak_idx = df['bol_mag'].idxmin()
        v_peak_idx = df['V'].idxmin()
        early = df[df['time'] < 5]
        sbo_peak = early['bol_mag'].min() if len(early) > 0 else 0
        sbo_time = (early.loc[early['bol_mag'].idxmin(), 'time']
                    if len(early) > 0 else -99)

        params = get_model_params(model_dir)

        results.append({
            'name': model_name,
            'score': score,
            'bol_peak_t': df.loc[peak_idx, 'time'],
            'bol_peak': df.loc[peak_idx, 'bol_mag'],
            'v_peak_t': df.loc[v_peak_idx, 'time'],
            'v_peak': df.loc[v_peak_idx, 'V'],
            'sbo_peak_t': sbo_time,
            'sbo_peak': sbo_peak,
            'ndof': ndof,
            'params': params,
            **details,
        })

    # Sort by score (lower = better)
    results.sort(key=lambda x: x['score'] if not np.isnan(x['score']) else 9e9)

    # Print ranking table
    print(f"{'Model':<22} {'Score':>7} {'rms_V':>7} {'rms_BV':>7} "
          f"{'rms_UV':>7} {'rms_Bol':>7} {'BolPk_t':>8} {'BolPk':>7} "
          f"{'VPk_t':>7} {'VPk':>7}")
    print("-" * 100)
    for r in results:
        print(f"{r['name']:<22} {r['score']:7.2f} {r['rcs_v']:7.2f} "
              f"{r['rcs_bv']:7.2f} {r['rcs_uv']:7.2f} {r['rcs_bol']:7.2f} "
              f"{r['bol_peak_t']:8.1f} {r['bol_peak']:7.2f} "
              f"{r['v_peak_t']:7.1f} {r['v_peak']:7.2f}")

    print(f"\nSN 2011dh targets: BolPk ~20d @ M_bol~-16.9, "
          f"VPk ~18-20d @ M_V~-17.2")

    # Print parameters for top models
    print(f"\n=== Top {args.top} models — parameters ===")
    for r in results[:args.top]:
        p = r['params']
        pstr = ", ".join(f"{k}={v}" for k, v in sorted(p.items()))
        print(f"  {r['name']:<22} score={r['score']:.2f}  {pstr}")

    # Detailed epoch comparison for top model
    print(f"\n=== Top model detailed epoch comparison: {results[0]['name']} ===")
    best_dir = f"{BASE}/{results[0]['name']}"
    best_df = load_model(best_dir)
    if best_df is not None:
        print(f"{'t':>5} | {'Model V':>8} {'2011dh V':>8} {'dV':>6} | "
              f"{'Model B-V':>9} {'2011dh B-V':>10} {'dB-V':>6} | "
              f"{'Model Bol':>9} {'2011dh Bol':>10} {'dBol':>6}")
        for t in [3, 5, 7, 10, 15, 18, 20, 25, 30, 35, 40]:
            mv = np.interp(t, best_df['time'], best_df['V'])
            mb = np.interp(t, best_df['time'], best_df['B'])
            mbol = np.interp(t, best_df['time'], best_df['bol_mag'])
            # Observed
            ov = obs_ubvri.iloc[(obs_ubvri['phase'] - t).abs().idxmin()]
            ob = obs_bol.iloc[(obs_bol['phase'] - t).abs().idxmin()]
            dv = mv - ov['M_V']
            dbv = (mb - mv) - ov.get('MB_MV', np.nan)
            dbol = mbol - ob['M_bol']
            print(f"{t:5d} | {mv:8.2f} {ov['M_V']:8.2f} {dv:6.2f} | "
                  f"{mb-mv:9.2f} {ov.get('MB_MV', float('nan')):10.2f} "
                  f"{dbv:6.2f} | "
                  f"{mbol:9.2f} {ob['M_bol']:10.2f} {dbol:6.2f}")

    # Plot
    if args.plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not available, skipping plot")
            return

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        # V-band
        ax = axes[0, 0]
        ax.scatter(obs_ubvri['phase'], obs_ubvri['M_V'], c='k',
                   s=10, zorder=5, label='SN 2011dh')
        colors = plt.cm.tab10(np.linspace(0, 1, min(len(results), args.top)))
        for i, r in enumerate(results[:args.top]):
            df = load_model(f"{BASE}/{r['name']}")
            if df is not None:
                ax.plot(df['time'], df['V'], alpha=0.8, color=colors[i],
                        label=f"{r['name']} ({r['score']:.2f})")
        ax.set_xlabel('Phase (days)')
        ax.set_ylabel('Absolute V magnitude')
        ax.set_title('V-band')
        ax.invert_yaxis()
        ax.legend(fontsize=7)
        ax.set_xlim(0, args.t_max)

        # B-V
        ax = axes[0, 1]
        obs_bv = obs_ubvri.dropna(subset=['MB_MV'])
        ax.scatter(obs_bv['phase'], obs_bv['MB_MV'], c='k', s=10, zorder=5)
        for i, r in enumerate(results[:args.top]):
            df = load_model(f"{BASE}/{r['name']}")
            if df is not None:
                bv = df['B'] - df['V']
                ax.plot(df['time'], bv, alpha=0.8, color=colors[i])
        ax.set_xlabel('Phase (days)')
        ax.set_ylabel('B-V colour')
        ax.set_title('B-V colour')
        ax.set_xlim(0, args.t_max)

        # Bolometric
        ax = axes[1, 0]
        ax.scatter(obs_bol['phase'], obs_bol['M_bol'], c='k', s=10, zorder=5)
        for i, r in enumerate(results[:args.top]):
            df = load_model(f"{BASE}/{r['name']}")
            if df is not None:
                ax.plot(df['time'], df['bol_mag'], alpha=0.8, color=colors[i])
        ax.set_xlabel('Phase (days)')
        ax.set_ylabel('Absolute bolometric magnitude')
        ax.set_title('Bolometric')
        ax.invert_yaxis()
        ax.set_xlim(0, args.t_max)

        # U-V
        ax = axes[1, 1]
        obs_uv = obs_ubvri.dropna(subset=['MU_MV'])
        ax.scatter(obs_uv['phase'], obs_uv['MU_MV'], c='k', s=10, zorder=5)
        for i, r in enumerate(results[:args.top]):
            df = load_model(f"{BASE}/{r['name']}")
            if df is not None:
                uv = df['U'] - df['V']
                ax.plot(df['time'], uv, alpha=0.8, color=colors[i])
        ax.set_xlabel('Phase (days)')
        ax.set_ylabel('U-V colour')
        ax.set_title('U-V colour')
        ax.set_xlim(0, args.t_max)

        fig.suptitle('SN_A grid models vs SN 2011dh (Ergon et al. 2015)',
                     fontsize=14)
        plt.tight_layout()
        out = f"{BASE}/outputs/sn_a_vs_2011dh_comparison.png"
        plt.savefig(out, dpi=150)
        print(f"\nPlot saved to {out}")


if __name__ == "__main__":
    main()
