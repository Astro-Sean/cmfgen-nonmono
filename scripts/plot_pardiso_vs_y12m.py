#!/usr/bin/env python3
"""Compare baseline Y12M vs PARDISO parallel STELLA at matched epochs."""
import importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import glob

PROJ = Path('/viper/ptmp1/sbrennan/SN2024pba')
spec = importlib.util.spec_from_file_location(
    'conv', PROJ/'scripts'/'convert_stella_swd_to_cmfgen.py')
conv = importlib.util.module_from_spec(spec); spec.loader.exec_module(conv)

def load(path):
    ep = conv.read_swd(path)
    return {t: ep[t] for t in sorted(ep) if ep[t] is not None}

bswd = glob.glob(str(PROJ/'parallel_test_baseline/stella/strad/*.swd'))[0]
pswd = str(PROJ/'pardiso_direct/strad/parallel_test_baseline-debug.swd')
B = load(bswd); P = load(pswd)
bt, pt = sorted(B), sorted(P)
print(f'baseline {len(bt)} ep  {bt[0]:.2f}..{bt[-1]:.2f}')
print(f'parallel {len(pt)} ep  {pt[0]:.2f}..{pt[-1]:.2f}')

# matched epochs
common = [t for t in bt if any(abs(t-p) < 0.05 for p in pt)]
print(f'{len(common)} common epochs (|dt|<0.05d)')

fig, axes = plt.subplots(2, 2, figsize=(15, 11))
picks = common[::max(1, len(common)//6)] or common

# density & velocity overlay at matched times
for tt in picks:
    pm = min(pt, key=lambda p: abs(p-tt))
    a = 0.3 + 0.6*(tt-bt[0])/(bt[-1]-bt[0])
    axes[0,0].plot(B[tt]['r']/1e15, B[tt]['rho'], 'k-', alpha=a, lw=1.6)
    axes[0,0].plot(P[pm]['r']/1e15, P[pm]['rho'], 'r--', alpha=a, lw=0.9)
    axes[0,1].plot(B[tt]['r']/1e15, B[tt]['v'], 'k-', alpha=a, lw=1.6)
    axes[0,1].plot(P[pm]['r']/1e15, P[pm]['v'], 'r--', alpha=a, lw=0.9)
axes[0,0].plot([],[],'k-',label='Y12M'); axes[0,0].plot([],[],'r--',label='PARDISO')
axes[0,0].set_yscale('log'); axes[0,0].set_xlabel('R (1e15 cm)'); axes[0,0].set_ylabel('rho'); axes[0,0].legend(); axes[0,0].set_title('Density evolution')
axes[0,1].set_xlabel('R (1e15 cm)'); axes[0,1].set_ylabel('v (km/s)'); axes[0,1].set_title('Velocity evolution')

# lightcurve (lum col index 11 -> 'flux' in 1e40 erg/s)
tb_ = [t for t in bt]; lb = [B[t]['flux'][0] for t in tb_]
tp_ = [t for t in pt]; lp = [P[t]['flux'][0] for t in tp_]
axes[1,0].plot(tb_, lb, 'k-', lw=1.5, label='Y12M')
axes[1,0].plot(tp_, lp, 'r--', lw=1.0, label='PARDISO')
axes[1,0].set_yscale('log'); axes[1,0].set_xlabel('t (d)'); axes[1,0].set_ylabel('L (1e40 erg/s)'); axes[1,0].legend(); axes[1,0].set_title('Bolometric lightcurve')

# relative differences at matched epochs vs time
rt, rr, rv = [], [], []
for tt in common:
    pm = min(pt, key=lambda p: abs(p-tt))
    nb = min(len(B[tt]['rho']), len(P[pm]['rho']))
    rt.append(tt)
    rr.append(np.abs(B[tt]['rho'][:nb]-P[pm]['rho'][:nb]).max()/(np.abs(B[tt]['rho'][:nb]).max()+1e-300))
    rv.append(np.abs(B[tt]['v'][:nb]-P[pm]['v'][:nb]).max()/(np.abs(B[tt]['v'][:nb]).max()+1e-300))
axes[1,1].semilogy(rt, rr, 'o-', ms=3, label='rho max rel diff')
axes[1,1].semilogy(rt, rv, 's-', ms=3, label='v max rel diff')
axes[1,1].set_xlabel('t (d)'); axes[1,1].set_ylabel('max relative diff'); axes[1,1].legend(); axes[1,1].set_title('Y12M vs PARDISO agreement'); axes[1,1].grid(alpha=0.3)

for ax in axes.ravel(): ax.grid(alpha=0.3)
fig.suptitle('STELLA: original Y12M vs parallel PARDISO (matched epochs)', fontsize=13)
fig.tight_layout()
out = PROJ/'parallel_test_plots'/'y12m_vs_pardiso.png'
out.parent.mkdir(exist_ok=True)
fig.savefig(out, dpi=150)
print('saved', out)
