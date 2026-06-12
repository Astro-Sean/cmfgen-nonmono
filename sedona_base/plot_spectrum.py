import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import glob

c = 2.998e10

# Observed spectrum file - update path as needed
OBS_FILE = '/ptmp/sbrennan/SN2024pba/SN2024pba_60875.145_PHASE25.1_Nordic-Optical-Telescope_ALFOSC.txt'

# Distance to SN2024pba: scale observed flux -> intrinsic luminosity to match SEDONA output
D_MPC = 83.58
D_CM = D_MPC * 1e6 * 3.085677581e18  # cm
FLUX_TO_LUM = 4 * np.pi * D_CM**2      # multiply observed flux by this to get L_λ

def load_and_convert(fname):
    """Load a SEDONA spectrum and return wavelength (Angstrom) and flux per Angstrom."""
    freq, flux, _ = np.loadtxt(fname, skiprows=1, unpack=True)
    lam = c / freq * 1e8
    idx = np.argsort(lam)
    lam = lam[idx]
    flux = flux[idx]
    flux_per_a = flux * (c * 1e8) / (lam ** 2)
    return lam, flux_per_a

def load_observed(fname):
    """Load observed spectrum: wavelength(A), flux(erg/s/cm2/A), error."""
    lam, flux, err = np.loadtxt(fname, unpack=True)
    idx = np.argsort(lam)
    return lam[idx], flux[idx], err[idx]

# Find all spectrum files and sort by iteration number
spec_files = sorted(glob.glob('spectrum_*.dat'), key=lambda f: int(f.split('_')[1].split('.')[0]))
if not spec_files:
    print("No spectrum_*.dat files found in the current directory.")
    sys.exit(1)

print(f"Found {len(spec_files)} spectrum file(s): {spec_files}")

fig, ax = plt.subplots(figsize=(12, 6))
colors = plt.cm.viridis(np.linspace(0, 1, len(spec_files)))

for fname, color in zip(spec_files, colors):
    try:
        lam, flux_per_a = load_and_convert(fname)
        # Only plot 3000-10000 A for correct y-axis scaling
        mask = (lam >= 3000) & (lam <= 10000)
        lam = lam[mask]
        flux_per_a = flux_per_a[mask]
        label = fname.replace('.dat', '')
        ax.step(lam, flux_per_a, where='mid', color=color, linewidth=0.8, label=label, zorder=3)
    except Exception as e:
        print(f"Warning: could not read {fname}: {e}")

# Overlay observed spectrum
if os.path.isfile(OBS_FILE):
    try:
        obs_lam, obs_flux, obs_err = load_observed(OBS_FILE)
        # Only plot 3000-10000 A for correct y-axis scaling
        mask = (obs_lam >= 3000) & (obs_lam <= 10000)
        obs_lam = obs_lam[mask]
        obs_flux = obs_flux[mask]
        obs_err = obs_err[mask]
        obs_flux_scaled = obs_flux * FLUX_TO_LUM
        obs_err_scaled = obs_err * FLUX_TO_LUM
        ax.step(obs_lam, obs_flux_scaled, where='mid', color='black', linewidth=0.8,
                label=f'SN2024pba +25.1d (NOT/ALFOSC, scaled by 4πd², d={D_MPC} Mpc)', zorder=4)
        ax.fill_between(obs_lam, obs_flux_scaled - obs_err_scaled, obs_flux_scaled + obs_err_scaled,
                        color='black', alpha=0.2, step='mid', zorder=2)
        print(f"Observed spectrum loaded: {len(obs_lam)} points, "
              f"{obs_lam.min():.0f}–{obs_lam.max():.0f} Å")
    except Exception as e:
        print(f"Warning: could not read observed spectrum: {e}")
else:
    print(f"Warning: observed spectrum not found at {OBS_FILE}")

# Rest-wavelength annotations for key lines (Å)
# Balmer series: H-alpha 6563, H-beta 4861, H-gamma 4340, H-delta 4102, H-epsilon 3970, H-zeta 3890
# 10 strongest Fe II optical lines (approximate wavelengths)
LINE_MARKS = [
    # Balmer series
    (3890, 'H$\\zeta$', 'C0'),
    (3970, 'H$\\epsilon$', 'C0'),
    (4102, 'H$\\delta$', 'C0'),
    (4340, 'H$\\gamma$', 'C0'),
    (4861, 'H$\\beta$', 'C0'),
    (6563, 'H$\\alpha$', 'C0'),
    # Fe II - 10 strongest optical lines
    (4233, 'Fe II', 'C1'),
    (4385, 'Fe II', 'C1'),
    (4924, 'Fe II', 'C1'),
    (5018, 'Fe II', 'C1'),
    (5169, 'Fe II', 'C1'),
    (5276, 'Fe II', 'C1'),
    (5316, 'Fe II', 'C1'),
    (5425, 'Fe II', 'C1'),
    (5588, 'Fe II', 'C1'),
    (6149, 'Fe II', 'C1'),
    # Other lines
    (6300, '[O I]', 'C2'),
    (6347, 'Si II', 'C3'),
    (6371, 'Si II', 'C3'),
    (7231, 'C II', 'C4'),
    (7774, 'O I', 'C5'),
    (8498, 'Ca II', 'C6'),
    (8542, 'Ca II', 'C6'),
    (8662, 'Ca II', 'C6'),
]

# Strategic label positioning to avoid overlap
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

for wl, label, color in LINE_MARKS:
    if 3000 <= wl <= 10000:
        ax.axvline(wl, color=color, linestyle='--', linewidth=0.7, alpha=0.6, zorder=1)
        
        # Add labels only for designated representative lines
        if label in LABEL_POSITIONS:
            ref_wl, ha_align = LABEL_POSITIONS[label]
            if wl == ref_wl:
                # Alternate vertical positions to avoid overlap
                y_pos = ax.get_ylim()[1] * 0.95
                if label in ['Fe II', '[O I]', 'Si II', 'Ca II']:
                    y_pos = ax.get_ylim()[1] * 0.85
                
                ax.text(wl, y_pos, label, rotation=90, va='top', ha=ha_align,
                        fontsize=6, color=color, alpha=0.8)

ax.set_xlabel(r'Wavelength ($\AA$)')
ax.set_ylabel(r'Luminosity (erg s$^{-1}$ $\AA^{-1}$)')
ax.set_xlim(3000, 10000)
ax.set_title('SEDONA Spectrum vs SN2024pba +25.1d (NOT/ALFOSC)')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper right', fontsize='small')
plt.tight_layout()
plt.savefig('spectrum_auto_plot.png', dpi=150, bbox_inches='tight')
print('Spectrum plot saved to spectrum_auto_plot.png')
