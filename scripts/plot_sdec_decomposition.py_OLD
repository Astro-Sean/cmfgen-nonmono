#!/usr/bin/env python3
"""
Build a pseudo-SDEC (spectral element decomposition) plot from CMFGEN
species-ablation OBS runs produced by submit_sdec_ablation.py.

Reads the baseline spectrum (SDEC/ALL/obs_fin_10 preferred, else
obs/obs_fin_10) plus each SDEC/<Elem>/obs_fin_10
ablated spectrum.  For every element the difference

    delta = F_all - F_noElem

is its net contribution: where delta > 0 the element adds flux (net
emission -- removing it lowers the spectrum), where delta < 0 it
removes flux (net absorption).  The plot stacks clip(delta,0) upward
from zero as emission and clip(-delta,0) upward from the model
spectrum as absorption -- the top of the absorption stack is the
approximate 'unabsorbed' flux, TARDIS SDEC style.

Caveat: the per-species deltas are NOT additive -- line overlap and
blanketing mean F_all - sum(F_all - F_noX) != sum of the components.
Populations are frozen at the converged state, same approximation as
TARDIS SDEC's last-interaction tag.  The default 'additive' attribution
normalizes the positive deltas at each wavelength so they partition
F_all exactly -- the colored stacks fill the area under the spectrum,
each element's share proportional to its delta; bins where no element
emits are left unpainted.  '--attribution exclusive'
assigns each bin to the single species with the largest delta;
'--attribution shared' stacks the raw overlapping deltas instead.

Baseline provenance: the deltas assume the baseline and ablated runs
restored the same FGINEXT_GS I_NEXT generation -- the emergent flux
contains an ESEC*J term evaluated on that field.  A mismatch injects a
species-independent offset into every delta, which shows up as the same
fake emission/absorption under several elements.  The script therefore
prefers SDEC/ALL/obs_fin_10 (same submission batch) over
obs/obs_fin_10 and warns when the restored tags differ.
"""

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_cmfgen_obs_spectrum import read_obs_fin_10  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402


I_NEXT_TAG_RE = re.compile(r'I_NEXT GENERATION \(tag=\s*(\d+)\)')


def named_spec_path(sdec_dir, label):
    """Descriptive copy written by newer submit_sdec_ablation jobs:
    SDEC/output_spectrum_<model>_<epoch>_<label>.dat.  The baseline
    has no suffix; the label is always the last _-separated token."""
    suffix_re = re.compile(
        r'^output_spectrum_.+_' + re.escape(label) + r'\.dat$')
    cands = [f for f in sorted(sdec_dir.glob('output_spectrum_*.dat'))
             if suffix_re.match(f.name)
             and f.is_file() and f.stat().st_size]
    return cands[0] if cands else None


def find_baseline(cmfgen_dir):
    # SDEC/ALL is preferred: it is submitted in the same batch as the
    # ablation jobs and therefore restores the same FGINEXT_GS I_NEXT
    # generation.  A pipeline obs/obs_fin_10 can carry a different
    # generation (e.g. if the OBS job read FGINEXT_GS while FINAL was
    # still writing it), which adds a species-independent offset to
    # every delta through the ESEC*J emissivity term.  The named
    # output_spectrum_* copy of ALL is equivalent provenance.
    sdec = cmfgen_dir / 'SDEC'
    cands = [f for f in sorted(sdec.glob('output_spectrum_*.dat'))
             if not re.search(r'_[A-Z][a-z]?\.dat$', f.name)
             and '_CONT.dat' not in f.name and '_ALL.dat' not in f.name]
    for cand in cands + [sdec / 'ALL' / 'obs_fin_10',
                         cmfgen_dir / 'obs' / 'obs_fin_10']:
        if cand.is_file() and cand.stat().st_size:
            return cand
    return None


def restored_tag(run_dir):
    """I_NEXT generation tag an OBS run restored from FGINEXT_GS, else
    None if OUT_FLUX is missing or the tag line is absent."""
    out_flux = Path(run_dir) / 'OUT_FLUX'
    if not out_flux.is_file():
        return None
    m = I_NEXT_TAG_RE.search(out_flux.read_text(errors='replace'))
    return int(m.group(1)) if m else None


ELEMENT_ORDER = (
    'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn '
    'Fe Co Ni Cu Zn Ba').split()
KNOWN_ELEMENTS = frozenset(ELEMENT_ORDER)
_Z_INDEX = {e: i for i, e in enumerate(ELEMENT_ORDER)}

# ColorBrewer Paired: qualitative palette, cycled if more species
# are plotted than it has colours.
SPECIES_PALETTE = list(plt.get_cmap('Paired').colors)


def load_ablated(cmfgen_dir, only=None):
    """Return ordered list of (elem, path).  The descriptive
    SDEC/output_spectrum_*_<elem>.dat copies are preferred over
    SDEC/<elem>/obs_fin_10 (identical contents, clearer provenance)."""
    sdec = cmfgen_dir / 'SDEC'
    out = []
    for d in sorted(sdec.iterdir()):
        if not d.is_dir() or d.name in ('ALL', 'CONT'):
            continue
        if d.name not in KNOWN_ELEMENTS:
            print(f'  skipping {d.name}: not an element symbol')
            continue
        f = named_spec_path(sdec, d.name) or d / 'obs_fin_10'
        if f.is_file() and f.stat().st_size:
            if only is None or d.name in only:
                out.append((d.name, f))
    out.sort(key=lambda t: _Z_INDEX[t[0]])
    return out


def draw_stacks(ax, wl_m, emis, abso, names, colors, abs_base=None,
                abs_colors=None):
    """Stack emission upward from zero and absorption upward from
    abs_base (the model spectrum) on ax; return the cumulative arrays
    so callers can set y limits.  abs_colors may differ from colors
    when the emission stack carries a non-species base layer."""
    cum_e = np.zeros_like(wl_m)
    for e, name, c in zip(emis, names, colors):
        ax.fill_between(wl_m, cum_e, cum_e + e, color=c, alpha=0.85,
                        label=name, linewidth=0)
        cum_e += e
    cum_a = np.zeros_like(wl_m) if abs_base is None else abs_base.copy()
    for a, c in zip(abso, abs_colors if abs_colors is not None else colors):
        ax.fill_between(wl_m, cum_a, cum_a + a, color=c, alpha=0.85,
                        linewidth=0)
        cum_a += a
    return cum_e, cum_a


def main():
    ap = argparse.ArgumentParser(
        description='Pseudo-SDEC stacked decomposition plot for CMFGEN')
    ap.add_argument('cmfgen_dir',
                    help='CMFGEN run dir containing obs/ and SDEC/')
    ap.add_argument('--baseline', default=None,
                    help='baseline obs_fin_10 (default: obs/obs_fin_10 '
                         'then SDEC/ALL/obs_fin_10)')
    ap.add_argument('--species', nargs='+', default=None,
                    help='only these element dirs')
    ap.add_argument('--lam-min', type=float, default=3500.0)
    ap.add_argument('--lam-max', type=float, default=9000.0)
    ap.add_argument('--continuum', action=argparse.BooleanOptionalAction,
                    default=True,
                    help='draw SDEC/CONT/obs_fin_10 (all-species '
                         'ablation, i.e. the OBS-computed continuum) '
                         'as a dashed line, when present '
                         '(default: on)')
    ap.add_argument('--attribution',
                    choices=('additive', 'exclusive', 'shared'),
                    default='additive',
                    help='additive: normalize positive deltas per bin so '
                         'they partition F_all exactly -- colored area '
                         'fills the spectrum (default); exclusive: each '
                         'wavelength bin is assigned to the single '
                         'species with the largest delta; shared: stack '
                         'all per-species deltas (overlapping, '
                         'non-additive)')
    ap.add_argument('--out', default=None,
                    help='output PNG (default <dir>/SDEC/sdec_decomposition.png)')
    ap.add_argument('--write-csv', action=argparse.BooleanOptionalAction,
                    default=True,
                    help='write sdec_components.csv next to the plot '
                         '(wavelength, model flux, continuum, '
                         'per-element contributions; default: on)')
    args = ap.parse_args()

    cmfgen_dir = Path(args.cmfgen_dir).resolve()
    if args.baseline:
        base_path = Path(args.baseline)
    else:
        base_path = find_baseline(cmfgen_dir)
    if base_path is None or not base_path.is_file():
        sys.exit('ERROR: no baseline obs_fin_10 (checked obs/ and SDEC/ALL/)')

    ablated = load_ablated(cmfgen_dir, args.species)
    if not ablated:
        sys.exit(f'ERROR: no SDEC/<elem>/obs_fin_10 under {cmfgen_dir}/SDEC')

    print(f'baseline: {base_path}')
    if base_path.parent.name == 'obs':
        print('WARNING: baseline is obs/obs_fin_10, not a same-batch '
              'SDEC/ALL run -- check I_NEXT tags below')
    # OUT_FLUX lives in the per-label run dir, not beside the named
    # output_spectrum_* copy at SDEC/ level.
    base_tag = restored_tag(cmfgen_dir / 'SDEC' / 'ALL') \
        if base_path.parent.name == 'SDEC' \
        else restored_tag(base_path.parent)
    wl, f_all = read_obs_fin_10(str(base_path))
    order = np.argsort(wl)
    wl, f_all = wl[order], f_all[order]

    mask = (wl >= args.lam_min) & (wl <= args.lam_max)
    wl_m, f_m = wl[mask], f_all[mask]
    f_pos = np.maximum(f_m, 0.0)

    # OBS-computed continuum: a same-batch CONT run has every species
    # SCL_*_ABUND-scaled, so only electron scattering remains.  Same
    # I_NEXT generation is required for consistency with the deltas.
    cont = None
    if args.continuum:
        sdec_dir = cmfgen_dir / 'SDEC'
        cont_path = (named_spec_path(sdec_dir, 'CONT')
                     or sdec_dir / 'CONT' / 'obs_fin_10')
        if cont_path.is_file() and cont_path.stat().st_size:
            ctag = restored_tag(sdec_dir / 'CONT')
            if base_tag is not None and ctag is not None and \
                    ctag != base_tag:
                print(f'WARNING: I_NEXT generation mismatch: baseline '
                      f'tag {base_tag}, CONT tag {ctag}')
            wl_c, f_c = read_obs_fin_10(str(cont_path))
            o = np.argsort(wl_c)
            cont = np.interp(wl_m, wl_c[o], f_c[o],
                             left=np.nan, right=np.nan)
            cont = np.nan_to_num(cont, nan=0.0)
            print(f'continuum: {cont_path} (I_NEXT tag {ctag})')

    emis, abso, names = [], [], []
    for elem, path in ablated:
        tag = restored_tag(cmfgen_dir / 'SDEC' / elem)
        if base_tag is not None and tag is not None and \
                tag != base_tag:
            print(f'WARNING: I_NEXT generation mismatch: baseline tag '
                  f'{base_tag}, {elem} tag {tag} -- delta_{elem} carries '
                  'a species-independent offset')
        wl_i, f_i = read_obs_fin_10(str(path))
        o = np.argsort(wl_i)
        f_i = np.interp(wl, wl_i[o], f_i[o], left=np.nan, right=np.nan)
        delta = f_all - f_i
        delta = np.nan_to_num(delta, nan=0.0, posinf=0.0, neginf=0.0)
        emis.append(np.clip(delta, 0.0, None)[mask])
        abso.append(np.clip(-delta, 0.0, None)[mask])
        names.append(elem)
        print(f'  {elem}: {path} (I_NEXT tag {tag})')

    if args.attribution == 'additive':
        # Partition F_all among the positive deltas: species k gets
        # fraction e_k/sum(e) of the flux at each wavelength, so the
        # colored stacks fill the area under the spectrum wherever an
        # element ablation changed the flux.  Bins with no delta stay
        # unpainted -- the dashed continuum line marks the CONT level.
        e_all = np.asarray(emis)
        wsum = e_all.sum(axis=0)
        frac = np.divide(e_all, wsum, out=np.zeros_like(e_all),
                         where=wsum > 0)
        emis_plot = [frac[k] * f_pos for k in range(len(names))]
        names_plot = list(names)
        abso_plot = abso
    elif args.attribution == 'exclusive':
        names_plot = list(names)
        # Each bin belongs to the species with the largest delta --
        # the deltas overlap, so stacking all of them assigns one
        # feature to several species.
        e_all = np.asarray(emis)
        a_all = np.asarray(abso)
        win_e = np.argmax(e_all, axis=0)
        win_a = np.argmax(a_all, axis=0)
        emis_plot = [np.where(win_e == k, e_all[k], 0.0)
                     for k in range(len(names))]
        abso_plot = [np.where(win_a == k, a_all[k], 0.0)
                     for k in range(len(names))]
    else:
        names_plot = list(names)
        emis_plot, abso_plot = emis, abso

    # Paired qualitative palette, cycled if there are more species
    # than it has colours.
    species_colors = [SPECIES_PALETTE[k % len(SPECIES_PALETTE)]
                      for k in range(len(names))]
    colors = species_colors
    flabel = r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ A$^{-1}$)'

    # Legend order: species in atomic order first, the continuum
    # line entry last.
    cont_label = r'Photosphere ($\tau_{es} \sim 2/3$)'
    tail_labels = (cont_label,)

    def legend_top(axis):
        handles, labels = axis.get_legend_handles_labels()
        order = ([i for i, l in enumerate(labels) if l not in tail_labels]
                 + [i for i, l in enumerate(labels) if l in tail_labels])
        axis.legend([handles[i] for i in order],
                    [labels[i] for i in order],
                    fontsize=9, ncol=max(len(labels), 1),
                    loc='lower center', bbox_to_anchor=(0.5, 1.0),
                    frameon=False)

    figsize = (16.0, 6.0)
    fig, ax = plt.subplots(figsize=figsize)
    cum_e, cum_a = draw_stacks(ax, wl_m, emis_plot, abso_plot,
                               names_plot, colors, abs_base=f_pos,
                               abs_colors=species_colors)
    ax.plot(wl_m, f_m, color='black', linewidth=0.6)
    if cont is not None:
        ax.plot(wl_m, cont, color='0.35', linewidth=0.9,
                linestyle='--', label=cont_label)
    ax.axhline(0.0, color='gray', linewidth=0.5)

    finite = f_m[np.isfinite(f_m)]
    if finite.size:
        ax.set_ylim(min(finite.min(), 0.0),
                    max(np.percentile(finite, 99) * 1.3,
                        cum_e.max() * 1.1, cum_a.max() * 1.1))
    ax.set_xlim(args.lam_min, args.lam_max)
    ax.set_xlabel('Wavelength (A)')
    ax.set_ylabel(flabel)
    ax.xaxis.set_minor_locator(MultipleLocator(500))
    ax.tick_params(axis='x', which='minor', length=4)
    legend_top(ax)
    ax.grid(True, alpha=0.2)

    out = Path(args.out) if args.out else \
        cmfgen_dir / 'SDEC' / 'sdec_decomposition.png'
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Plot saved to {out}')

    # Log version: single panel. Emission fills the area under the
    # spectrum, absorption stacks above the line toward the
    # 'unabsorbed' flux -- everything is positive so a plain log axis
    # suffices.
    fig2, ax_b = plt.subplots(figsize=figsize)

    cum_e, cum_a = draw_stacks(ax_b, wl_m, emis_plot, abso_plot,
                               names_plot, colors, abs_base=f_pos,
                               abs_colors=species_colors)
    ax_b.plot(wl_m, f_m, color='black', linewidth=0.6)
    if cont is not None:
        ax_b.plot(wl_m, cont, color='0.35', linewidth=0.9,
                  linestyle='--', label=cont_label)
    pos = f_m[np.isfinite(f_m) & (f_m > 0)]
    if pos.size:
        ax_b.set_yscale('log')
        ax_b.set_ylim(max(pos.min(), 1e-3),
                      max(cum_a.max(), cum_e.max(), pos.max()) * 1.3)
    ax_b.set_xlim(args.lam_min, args.lam_max)
    ax_b.set_xlabel('Wavelength (A)')
    ax_b.set_ylabel(flabel)
    ax_b.xaxis.set_minor_locator(MultipleLocator(500))
    ax_b.tick_params(axis='x', which='minor', length=4)
    legend_top(ax_b)
    ax_b.grid(True, alpha=0.2)

    log_out = out.with_name(out.stem + '_log' + out.suffix)
    plt.tight_layout()
    plt.savefig(log_out, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    print(f'Log plot saved to {log_out}')

    if args.write_csv:
        # The emis/abs columns are the plotted contributions (the
        # additive partition by default), so they reproduce the figure.
        # The raw ablation deltas are included separately since they
        # are the physically meaningful quantities -- the partition is
        # bookkeeping only.
        csv = out.with_name('sdec_components.csv')
        header = ['wavelength_A', 'F_all']
        cols = [wl_m, f_m]
        if cont is not None:
            header.append('F_cont')
            cols.append(cont)
        for n, e in zip(names_plot, emis_plot):
            header.append(f'emis_{n}')
            cols.append(e)
        for n, a in zip(names, abso_plot):
            header.append(f'abs_{n}')
            cols.append(a)
        for n, e, a in zip(names, emis, abso):
            header.append(f'delta_{n}')
            cols.append(e - a)
        np.savetxt(csv, np.column_stack(cols), delimiter=',',
                   header=','.join(header), comments='')
        print(f'Components written to {csv}')


if __name__ == '__main__':
    main()
