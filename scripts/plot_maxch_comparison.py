#!/usr/bin/env python3
"""Plot the CMFGEN population-correction diagnostic (MAXCH, i.e. "Maximum %
increase" reported by SOLVEBA_V13 each LAMBDA iteration) across the
different initial-population methods tried for SN_12C_proxy / X_1 5-day
snapshots.

MAXCH is CMFGEN's own measure of how much the level populations changed
during that iteration's correction step — it is the most direct available
proxy for "population change" driven by the choice of initial condition
(proxy_populations variant vs. built-in LTE_EST) and grid/VADAT settings.

Reads real per-iteration MAXCH traces from three OUTGEN logs (full runs) and
overlays single-point markers for the proxy_populations variants that were tested
and diverged/crashed before producing a usable iteration series (their
first-iteration MAXCH is recorded directly from the batch logs discussed in
AGENTS.md).
"""
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

MAXCH_RE = re.compile(
    r'Maximum % increase at depth\s+(\d+)\s+is\s+([\d.]+E[+-]?\d+)\s+\(LAMBDA\)\s+---\s+iteration\s+(\d+)'
)


def parse_maxch(path):
    """Return (iterations, maxch_values, depths) from an OUTGEN file."""
    iters, vals, depths = [], [], []
    seen = set()
    with open(path) as f:
        for line in f:
            m = MAXCH_RE.search(line)
            if m:
                depth, val, it = int(m.group(1)), float(m.group(2)), int(m.group(3))
                if it not in seen:
                    seen.add(it)
                    iters.append(it)
                    vals.append(val)
                    depths.append(depth)
    return iters, vals, depths


# --- Full iteration-series runs (real OUTGEN traces) -----------------------
RUNS = [
    {
        'label': 'LTE_EST baseline (no AMR fix), MAX_LAM=1.5',
        'path': '/viper/ptmp1/sbrennan/SN2024pba/SN_12C_LTE_EST/cmfgen/5.0days/OUTGEN',
        'color': 'tab:red',
    },
    {
        'label': 'LTE_EST + AMR/Ne-priority fix, MAX_LAM=1.5',
        'path': '/viper/ptmp1/sbrennan/SN2024pba/SN_12C_proxy/cmfgen/5.0days/OUTGEN',
        'color': 'tab:orange',
    },
    {
        'label': 'LTE_EST + AMR fix + outer-boundary fix\n(X_1 day 5), MAX_LAM=1.05',
        'path': '/viper/ptmp1/sbrennan/SN2024pba/models/X_1/cmfgen/5.0days/OUTGEN',
        'color': 'tab:green',
    },
]

# --- Single-point references: proxy_populations variants that diverged/crashed ------
# First-iteration MAXCH values recorded directly from the batch logs while
# testing each approach (see AGENTS.md "proxy_populations vs LTE_EST" section).
PROXY_POINTS = [
    {
        'label': 'proxy_populations: original bug\n(Saha DI at T_local, wrong-sign DC)',
        'value': 6.75e23,
        'color': 'tab:purple',
        'marker': 's',
    },
    {
        'label': 'proxy_populations: dilute-Saha DI (T_ion)\n+ DC=1 nebular approx.',
        'value': 4.71e25,
        'color': 'tab:brown',
        'marker': '^',
    },
    {
        'label': 'proxy_populations: "exact" analytic DC\n(DI at T_ion, DC~1e13-1e26)',
        'value': 9.53e149,
        'color': 'black',
        'marker': 'x',
    },
]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={'width_ratios': [2.3, 1]})

max_iter_overall = 1
first_vals = {}  # label -> first-iteration MAXCH, for the bar chart panel
for run in RUNS:
    iters, vals, depths = parse_maxch(run['path'])
    if not iters:
        print(f"WARNING: no MAXCH data found in {run['path']}")
        continue
    max_iter_overall = max(max_iter_overall, max(iters))
    first_vals[run['label'].split('\n')[0]] = vals[0]
    ax.plot(iters, vals, color=run['color'], lw=1.8, label=run['label'])

    # Mark cascade point: value increases after a preceding decrease
    for i in range(2, len(vals)):
        if vals[i] > vals[i - 1] * 1.5 and vals[i - 1] < vals[i - 2]:
            ax.scatter([iters[i]], [vals[i]], color=run['color'], marker='*',
                       s=220, zorder=5, edgecolor='black', linewidth=0.8)
            ax.annotate(f'cascade, iter {iters[i]}\n(depth {depths[i-1]}\u2192{depths[i]})',
                        xy=(iters[i], vals[i]), xytext=(iters[i] + max_iter_overall * 0.015, vals[i] * 2.2),
                        fontsize=8, color=run['color'], ha='left')
            break

# Plot the two proxy_populations single-point references that are within the visible
# range of the zoomed panel (the "exact DC" point at ~1e150 is off-scale here
# and is shown only in the right-hand bar-chart panel instead).
for p in PROXY_POINTS:
    if p['value'] > 1e30:
        continue
    ax.scatter([1], [p['value']], color=p['color'], marker=p['marker'], s=140,
               zorder=6, edgecolor='black', linewidth=0.8, label=p['label'])

ax.set_yscale('log')
ax.set_ylim(1e9, 1e30)
ax.set_xlim(0, max_iter_overall * 1.02)
ax.set_xlabel('LAMBDA iteration', fontsize=12)
ax.set_ylabel('MAXCH: "Maximum % increase" in populations\n(CMFGEN SOLVEBA_V13)', fontsize=11)
ax.set_title('Convergence trajectory (zoomed)', fontsize=12)
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=7.5, loc='upper right', framealpha=0.9, ncol=1)

# --- Right panel: bar chart of first-iteration MAXCH across ALL methods,
# including the catastrophic "exact DC" case (off-scale in the left panel) ---
bar_items = list(first_vals.items()) + [(p['label'].split('\n')[0], p['value']) for p in PROXY_POINTS]
bar_colors = [r['color'] for r in RUNS] + [p['color'] for p in PROXY_POINTS]
names = [n for n, _ in bar_items]
values = [v for _, v in bar_items]
y_pos = range(len(names))
ax2.barh(y_pos, values, color=bar_colors, edgecolor='black', linewidth=0.6)
ax2.set_xscale('log')
ax2.set_yticks(list(y_pos))
ax2.set_yticklabels(names, fontsize=8)
ax2.set_xlabel('First-iteration MAXCH', fontsize=11)
ax2.set_title('Full-range comparison\n(all methods, iteration 1-2)', fontsize=12)
ax2.axvline(1.0, color='gray', ls=':', lw=1, alpha=0.6)
for y, v in zip(y_pos, values):
    ax2.text(v * 1.3, y, f'{v:.2e}', va='center', fontsize=7.5)
ax2.grid(True, which='major', axis='x', alpha=0.3)
ax2.set_xlim(1e19, 1e155)

fig.suptitle('SN_12C_proxy / X_1 5-day snapshot: population-correction diagnostic (MAXCH)\n'
             'across different initial-condition methods', fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.93])
out_path = '/viper/ptmp1/sbrennan/SN2024pba/plots/maxch_comparison_methods.png'
fig.savefig(out_path, dpi=150)
print(f"Saved: {out_path}")
