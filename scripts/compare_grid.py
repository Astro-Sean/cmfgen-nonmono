#!/usr/bin/env python3
"""Compare grid model lightcurves with SN 2011dh observations."""
import os, sys, glob
import numpy as np
import pandas as pd

BASE = "/viper/ptmp1/sbrennan/SN2024pba"

# SN 2011dh observed properties (approximate, from literature)
obs_2011dh = {
    1:  {"B-V": -0.10, "U-V": -0.50, "V": -14.5, "bol": -15.5},
    3:  {"B-V":  0.00, "U-V": -0.30, "V": -15.3, "bol": -16.0},
    5:  {"B-V":  0.15, "U-V":  0.00, "V": -16.0, "bol": -16.5},
    10: {"B-V":  0.35, "U-V":  0.50, "V": -16.8, "bol": -17.2},
    18: {"B-V":  0.55, "U-V":  1.30, "V": -17.2, "bol": -17.5},
    20: {"B-V":  0.60, "U-V":  1.50, "V": -17.2, "bol": -17.4},
    30: {"B-V":  0.85, "U-V":  2.00, "V": -16.8, "bol": -17.0},
    40: {"B-V":  0.90, "U-V":  2.20, "V": -16.3, "bol": -16.5},
}

def load_model(model_dir):
    """Load model lightcurve CSV."""
    csv = f"{model_dir}/output/{os.path.basename(model_dir)}_lightcurve_complete.csv"
    if not os.path.exists(csv):
        return None
    df = pd.read_csv(csv)
    df = df[df['bol_mag'] > -30]  # remove artifacts
    return df

def get_epochs(df, epochs):
    """Get model values at specific epochs."""
    results = {}
    for t in epochs:
        idx = (df['time'] - t).abs().idxmin()
        row = df.loc[idx]
        results[t] = {
            'time': row['time'],
            'U': row['U'], 'B': row['B'], 'V': row['V'],
            'R': row['R'], 'I': row['I'],
            'bol': row['bol_mag'],
            'B-V': row['B'] - row['V'],
            'U-V': row['U'] - row['V'],
        }
    return results

def score_model(df):
    """Score how well a model matches SN 2011dh (lower = better)."""
    epochs_data = get_epochs(df, list(obs_2011dh.keys()))
    score = 0
    n = 0
    for t, obs in obs_2011dh.items():
        if t not in epochs_data:
            continue
        m = epochs_data[t]
        # Weight bolometric and colors
        score += (m['bol'] - obs['bol'])**2 * 1.0
        score += (m['B-V'] - obs['B-V'])**2 * 2.0
        score += (m['U-V'] - obs['U-V'])**2 * 2.0
        score += (m['V'] - obs['V'])**2 * 0.5
        n += 1
    return np.sqrt(score / n) if n > 0 else 999

def main():
    # Find all grid models
    models = sorted(glob.glob(f"{BASE}/SN_A_g_*/output/*_lightcurve_complete.csv"))
    if not models:
        print("No grid model outputs found yet.")
        return

    print(f"Found {len(models)} completed models\n")

    # Load all and score
    results = []
    for csv_path in models:
        model_name = os.path.basename(os.path.dirname(os.path.dirname(csv_path)))
        df = pd.read_csv(csv_path)
        df = df[df['bol_mag'] > -30]
        if len(df) < 10:
            continue

        peak_idx = df['bol_mag'].idxmin()
        v_peak_idx = df['V'].idxmin()

        # Early SBO peak
        early = df[df['time'] < 5]
        sbo_peak = early['bol_mag'].min() if len(early) > 0 else 0
        sbo_time = early.loc[early['bol_mag'].idxmin(), 'time'] if len(early) > 0 else -99

        score = score_model(df)

        results.append({
            'name': model_name,
            'score': score,
            'bol_peak_t': df.loc[peak_idx, 'time'],
            'bol_peak': df.loc[peak_idx, 'bol_mag'],
            'v_peak_t': df.loc[v_peak_idx, 'time'],
            'v_peak': df.loc[v_peak_idx, 'V'],
            'sbo_peak_t': sbo_time,
            'sbo_peak': sbo_peak,
        })

    results.sort(key=lambda x: x['score'])

    print(f"{'Model':<22} {'Score':>7} {'BolPk_t':>8} {'BolPk':>7} {'VPk_t':>7} {'VPk':>7} {'SBO_t':>7} {'SBO':>7}")
    print("-" * 80)
    for r in results:
        print(f"{r['name']:<22} {r['score']:7.2f} {r['bol_peak_t']:8.1f} {r['bol_peak']:7.2f} "
              f"{r['v_peak_t']:7.1f} {r['v_peak']:7.2f} {r['sbo_peak_t']:7.1f} {r['sbo_peak']:7.2f}")

    print(f"\nTarget (SN 2011dh): BolPk~18-20d @ -17.2, VPk~18-20d @ -17.2, SBO< -16")

    # Detailed comparison for top 3
    print("\n\n=== Top 3 models detailed comparison ===")
    for r in results[:3]:
        model_dir = f"{BASE}/{r['name']}"
        df = load_model(model_dir)
        if df is None:
            continue
        epochs_data = get_epochs(df, list(obs_2011dh.keys()))

        print(f"\n--- {r['name']} (score={r['score']:.2f}) ---")
        print(f"{'t':>5} | {'Model':>35} | {'2011dh':>20} | {'Diff':>14}")
        print(f"{'':>5} | {'B-V':>7} {'U-V':>7} {'V':>7} {'bol':>7} | {'B-V':>7} {'U-V':>7} {'bol':>7} | {'dB-V':>5} {'dU-V':>5}")
        for t in [1, 3, 5, 10, 18, 20, 30, 40]:
            if t not in epochs_data:
                continue
            m = epochs_data[t]
            obs = obs_2011dh[t]
            print(f"{t:5d} | {m['B-V']:7.2f} {m['U-V']:7.2f} {m['V']:7.2f} {m['bol']:7.2f} | "
                  f"{obs['B-V']:7.2f} {obs['U-V']:7.2f} {obs['bol']:7.2f} | "
                  f"{m['B-V']-obs['B-V']:5.2f} {m['U-V']-obs['U-V']:5.2f}")

if __name__ == "__main__":
    main()
