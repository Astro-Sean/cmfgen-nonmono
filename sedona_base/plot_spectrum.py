import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
import sys
import os
import glob
import re

c = 2.998e10

# Distance constants
PC_CM = 3.085677581491367e18

# Distance to SN2024pba — model spectra are placed at this distance
# so they are directly comparable to the observed spectrum (no scaling of observed data)
D_MPC = 83.58
D_CM = D_MPC * 1e6 * PC_CM
FOUR_PI_D2 = 4.0 * np.pi * D_CM**2  # 4*pi*d^2 in cm^2

# Observed spectrum file
OBS_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN2024pba_60875.145_PHASE25.1_Nordic-Optical-Telescope_ALFOSC.txt'


def load_and_convert(fname):
    """Load SEDONA spectrum -> (wavelength_A, flux_at_10pc_per_A).

    SEDONA output: freq (Hz), luminosity (erg/s/Hz), count.
    Convert to flux at 10 pc: F = L / (4*pi*d^2), then /d_lambda to get per-Angstrom.
    """
    freq, lum_per_hz, _ = np.loadtxt(fname, skiprows=1, unpack=True)
    lam = c / freq * 1e8  # Angstrom
    idx = np.argsort(lam)
    lam = lam[idx]
    lum_per_hz = lum_per_hz[idx]
    # L_per_A = L_per_Hz * d_nu/d_lambda = L_per_Hz * c / lambda^2
    lum_per_a = lum_per_hz * (c * 1e8) / (lam ** 2)
    # Flux at SN2024pba distance (83.58 Mpc)
    flux_per_a = lum_per_a / FOUR_PI_D2
    return lam, flux_per_a


def load_observed(fname):
    """Load observed spectrum: wavelength(A), flux(erg/s/cm2/A), error.

    Observed data is already at 83.58 Mpc — no scaling needed.
    """
    lam, flux, err = np.loadtxt(fname, unpack=True)
    idx = np.argsort(lam)
    return lam[idx], flux[idx], err[idx]


def get_model_label():
    """Extract model name and snapshot time for the plot legend."""
    cwd = os.getcwd()
    parts = cwd.replace('\\', '/').split('/')
    model_name = ''
    if 'sedona' in parts:
        idx = parts.index('sedona')
        if idx > 0:
            model_name = parts[idx - 1]
    if not model_name:
        mod_files = glob.glob('*.mod')
        if mod_files:
            model_name = os.path.splitext(mod_files[0])[0]

    snapshot_days = None
    if os.path.isfile('param.lua'):
        with open('param.lua') as f:
            text = f.read()
        m = re.search(r'tstep_time_start\s*=\s*([^\n]+)', text)
        if m:
            expr = m.group(1).strip()
            expr = expr.split('--')[0].strip()
            expr = expr.rstrip(',;').strip()
            try:
                t_start_sec = float(eval(expr, {'__builtins__': {}}, {}))
                snapshot_days = t_start_sec / (3600.0 * 24.0)
            except Exception as e:
                print(f"Warning: could not evaluate tstep_time_start='{expr}': {e}")

    if model_name and snapshot_days is not None:
        if abs(snapshot_days - round(snapshot_days)) < 0.5:
            day_str = f"{int(round(snapshot_days))} days"
        else:
            day_str = f"{snapshot_days:.1f} days"
        return f"{model_name} {day_str}"
    elif model_name:
        return model_name
    else:
        return "SEDONA"


def get_plots_dir():
    """Find the model's plots/ directory from the sedona run directory."""
    cwd = os.getcwd()
    parts = cwd.replace('\\', '/').split('/')
    if 'sedona' in parts:
        idx = parts.index('sedona')
        if idx > 0:
            model_dir = '/'.join(parts[:idx])
            plots_dir = os.path.join(model_dir, 'plots')
            os.makedirs(plots_dir, exist_ok=True)
            return plots_dir
    return '.'


# Rest-wavelength annotations for key lines (Å)
LINE_MARKS = [
    (3890, 'H$\\zeta$', 'C0'),
    (3970, 'H$\\epsilon$', 'C0'),
    (4102, 'H$\\delta$', 'C0'),
    (4340, 'H$\\gamma$', 'C0'),
    (4861, 'H$\\beta$', 'C0'),
    (6563, 'H$\\alpha$', 'C0'),
    (4924, 'Fe II', 'C1'),
    (6300, '[O I]', 'C2'),
    (6347, 'Si II', 'C3'),
    (7231, 'C II', 'C4'),
    (7774, 'O I', 'C5'),
    (8498, 'Ca II', 'C6'),
]

LABEL_POSITIONS = {
    'H$\\alpha$': (6563, 'left'),
    'H$\\beta$': (4861, 'left'),
    'H$\\gamma$': (4340, 'left'),
    'H$\\delta$': (4102, 'left'),
    'H$\\epsilon$': (3970, 'left'),
    'H$\\zeta$': (3890, 'left'),
    'Fe II': (4924, 'left'),
    '[O I]': (6300, 'right'),
    'Si II': (6347, 'right'),
    'C II': (7231, 'left'),
    'O I': (7774, 'left'),
    'Ca II': (8498, 'right'),
}


def add_line_annotations(ax):
    """Add vertical dashed lines and labels for key spectral features."""
    for wl, label, color in LINE_MARKS:
        if 3000 <= wl <= 10000:
            ax.axvline(wl, color=color, linestyle='--', linewidth=0.7, alpha=0.6, zorder=1)
            if label in LABEL_POSITIONS:
                ref_wl, ha_align = LABEL_POSITIONS[label]
                if wl == ref_wl:
                    y_pos = ax.get_ylim()[1] * 0.95
                    if label in ['Fe II', '[O I]', 'Si II', 'Ca II']:
                        y_pos = ax.get_ylim()[1] * 0.85
                    ax.text(wl, y_pos, label, rotation=90, va='top', ha=ha_align,
                            fontsize=6, color=color, alpha=0.8)


# ── Main ──────────────────────────────────────────────────────────────

# Find all spectrum files and select the final one
spec_files = sorted(glob.glob('spectrum_*.dat'), key=lambda f: int(f.split('_')[1].split('.')[0]))
if not spec_files:
    print("No spectrum_*.dat files found in the current directory.")
    sys.exit(1)

final_file = spec_files[-1]
print(f"Plotting final spectrum: {final_file} (of {len(spec_files)} iterations)")

model_label = get_model_label()
print(f"Model label for legend: {model_label}")

plots_dir = get_plots_dir()
print(f"Plots directory: {plots_dir}")

# Extract epoch label from model_label for filenames (e.g. "X_1_fix 30d" -> "30d")
epoch_tag = ''
parts = model_label.split()
for p in parts:
    if p.endswith('d') and p[:-1].replace('.', '').isdigit():
        epoch_tag = p
        break
if not epoch_tag:
    # Fallback: try to read from param.lua
    if os.path.isfile('param.lua'):
        with open('param.lua') as f:
            text = f.read()
        m = re.search(r'tstep_time_start\s*=\s*([^\n]+)', text)
        if m:
            expr = m.group(1).strip().split('--')[0].strip().rstrip(',;').strip()
            try:
                t_days = float(eval(expr, {'__builtins__': {}}, {})) / 86400.0
                epoch_tag = f"{t_days:g}d"
            except Exception:
                pass
if not epoch_tag:
    epoch_tag = 'snapshot'
print(f"Epoch tag for filenames: {epoch_tag}")

# Load model spectrum
lam, flux = load_and_convert(final_file)
mask = (lam >= 3000) & (lam <= 10000)
lam = lam[mask]
flux = flux[mask]

# ── Plot 1: SEDONA spectrum only ──────────────────────────────────────

fig, ax = plt.subplots(figsize=(12, 6))
ax.step(lam, flux, where='mid', color='C0', linewidth=1.0,
        label=f'SEDONA {model_label}', zorder=5)
add_line_annotations(ax)
ax.set_xlabel(r'Wavelength ($\AA$)')
ax.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
ax.set_xlim(3000, 10000)
ax.set_title(f'SEDONA Spectrum (d={D_MPC} Mpc) - {model_label}')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper right', fontsize='small')
# Auto-scale y-lims to the visible data within xlims
pos_flux = flux[flux > 0]
if len(pos_flux) > 0:
    ax.set_ylim(pos_flux.min() * 0.3, pos_flux.max() * 3)
plt.tight_layout()
out1 = os.path.join(plots_dir, f'sedona_{epoch_tag}_spectrum.png')
plt.savefig(out1, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved spectrum plot: {out1}')

# ── Plot 2: SEDONA (smoothed) vs SN2024pba, dual y-axis ──────────────

fig, ax1 = plt.subplots(figsize=(12, 6))

# Smooth the model spectrum slightly for visibility (running mean over ~15 A)
smooth_window = max(3, int(round(15.0 / np.median(np.diff(lam)))))
flux_smooth = uniform_filter1d(flux, size=smooth_window)

# Model spectrum on left y-axis (at SN2024pba distance)
ax1.step(lam, flux_smooth, where='mid', color='C0', linewidth=1.2,
         label=f'SEDONA {model_label} (smoothed)', zorder=5)
ax1.set_xlabel(r'Wavelength ($\AA$)')
ax1.set_ylabel(r'Model flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)', color='C0')
ax1.tick_params(axis='y', labelcolor='C0')
ax1.set_xlim(3000, 10000)
# Auto-scale y-lims to the visible model data within xlims
pos_flux = flux_smooth[flux_smooth > 0]
if len(pos_flux) > 0:
    ax1.set_ylim(pos_flux.min() * 0.3, pos_flux.max() * 3)

# Observed spectrum on right y-axis (at observed distance, no scaling)
has_obs = False
if os.path.isfile(OBS_FILE):
    try:
        obs_lam, obs_flux, obs_err = load_observed(OBS_FILE)
        mask = (obs_lam >= 3000) & (obs_lam <= 10000)
        obs_lam = obs_lam[mask]
        obs_flux = obs_flux[mask]
        obs_err = obs_err[mask]

        ax2 = ax1.twinx()
        ax2.step(obs_lam, obs_flux, where='mid', color='black', linewidth=0.8,
                 label='SN2024pba +25.1d (NOT/ALFOSC)', zorder=4)
        ax2.fill_between(obs_lam, obs_flux - obs_err, obs_flux + obs_err,
                         color='black', alpha=0.15, step='mid', zorder=2)
        ax2.set_ylabel(r'Observed flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)',
                       color='black')
        ax2.tick_params(axis='y', labelcolor='black')
        has_obs = True
        print(f"Observed spectrum loaded: {len(obs_lam)} points, "
              f"{obs_lam.min():.0f}-{obs_lam.max():.0f} A")
    except Exception as e:
        print(f"Warning: could not read observed spectrum: {e}")
else:
    print(f"Warning: observed spectrum not found at {OBS_FILE}")

add_line_annotations(ax1)
ax1.set_title(f'SEDONA {model_label} vs SN2024pba +25.1d')
ax1.grid(True, alpha=0.3)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
if has_obs:
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right',
               fontsize='small')
else:
    ax1.legend(loc='upper right', fontsize='small')

plt.tight_layout()
out2 = os.path.join(plots_dir, f'sedona_{epoch_tag}_spectrum_comparison.png')
plt.savefig(out2, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved comparison plot: {out2}')

# Backward-compatible single-axis plot
fig, ax = plt.subplots(figsize=(12, 6))
ax.step(lam, flux_smooth, where='mid', color='C0', linewidth=1.2,
        label=f'SEDONA {model_label} (smoothed)', zorder=5)
if os.path.isfile(OBS_FILE):
    try:
        obs_lam, obs_flux, obs_err = load_observed(OBS_FILE)
        mask = (obs_lam >= 3000) & (obs_lam <= 10000)
        ax.step(obs_lam[mask], obs_flux[mask], where='mid', color='black',
                linewidth=0.8, label=f'SN2024pba +25.1d', zorder=4)
        ax.fill_between(obs_lam[mask], (obs_flux - obs_err)[mask],
                        (obs_flux + obs_err)[mask], color='black', alpha=0.15,
                        step='mid', zorder=2)
    except Exception:
        pass
add_line_annotations(ax)
ax.set_xlabel(r'Wavelength ($\AA$)')
ax.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
ax.set_xlim(3000, 10000)
ax.set_title(f'SEDONA {model_label} vs SN2024pba +25.1d (d={D_MPC} Mpc)')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper right', fontsize='small')
# Auto-scale y-lims to the visible data within xlims
pos_flux = flux_smooth[flux_smooth > 0]
if len(pos_flux) > 0:
    ax.set_ylim(pos_flux.min() * 0.3, pos_flux.max() * 3)
plt.tight_layout()
plt.savefig('spectrum_auto_plot.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved backward-compatible plot: spectrum_auto_plot.png')
