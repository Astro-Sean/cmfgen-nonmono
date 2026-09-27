#!/usr/bin/env python3
"""
Submit CMFGEN OBS (cmf_flux) species-ablation jobs for a pseudo-SDEC
spectral decomposition of a converged model.

For each requested element this creates <cmfgen_dir>/SDEC/<Elem>/ and
submits an OBS job whose CMF_FLUX_PARAM sets SCL_<species>_ABUND to a
tiny value, removing that element's bound-bound / bound-free / free-free
opacity from the formal solution while T, n_e and all other populations
stay frozen at the converged state.  The difference F_all - F_noX is the
element's net spectral contribution (emission where negative, absorption
where positive); plot_sdec_decomposition.py renders it TARDIS-SDEC style.

An unablated ALL job is always submitted in the same batch (unless
--no-submit-all) because the baseline must be computed with the same
FGINEXT_GS I_NEXT generation as the ablations: the emergent flux
contains an ESEC*J electron-scattering emissivity evaluated on the
restored radiation field.  A pipeline obs/obs_fin_10 can carry a
different generation (e.g. the OBS job read FGINEXT_GS while FINAL was
still writing it), which injects a species-independent offset into
every delta and produces spurious multi-species features in the
decomposition.  plot_sdec_decomposition.py therefore prefers
SDEC/ALL/obs_fin_10 over obs/obs_fin_10 and cross-checks the restored
I_NEXT tag of every run.

Why SCL_*_ABUND and not TRANS_*=NONE:
  - The non-monotonic OBS solver requires GLOBAL_LINE=BLANK and all
    TRANS_*=BLANK (prepare_CMF_PARAM.py enforces this).  Under
    GLOBAL_LINE=BLANK, per-ion TRANS settings are overwritten anyway
    (cmf_flux_sub_v5.f), so TRANS_*=NONE cannot ablate an ion.
  - SCL_*_ABUND scales the species populations inside cmf_flux only
    (cmf_flux_sub_v5.f / par_opac_sub.f); no new CMFGEN solve is done.
  - Electron-scattering opacity uses ED from RVTJ (ESEC=6.65e-15*ED), so
    SCL_HYD/SCL_HE remove H/He lines and edges without destroying the
    e.s. continuum -- exactly the "no <elem>" spectrum SDEC needs.

Directory layout (mirrors the pipeline obs/ convention; cmf_flux
resolves MODEL, MODEL_SPEC and POP<SPECIES> via DIR_NAME derived from
the ../RVTJ stdin path, and batobs.sh needs ../VADAT + ../batch.sh):

    <cmfgen_dir>/SDEC/
        VADAT RVTJ batch.sh MODEL_SPEC MODEL POP* FGINEXT_GS
            -> symlinks to ../<name>
        <Elem>/
            run_sdec_<Elem>.slurm
            (batobs.sh, CMF_FLUX_PARAM_INIT, SCRATCH/ made at runtime)

Submission: job i is submitted with --dependency=afterok:<job_{i-N}> so
at most --max-parallel (default 2) ablation jobs run concurrently, in
order.  A final job depends afterok on all ablation jobs and runs
plot_sdec_decomposition.py.

These are OBS-only diagnostics on an existing converged model (RVTJ +
POP* are reused, no populations are recomputed), not new CMFGEN models.
"""

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

# SPECIES_ABR -> element symbol (authoritative: cmfgen_source par_opac.f
# / cmf_flux_v5.f SPECIES_ABR table).
ABR_TO_ELEM = {
    'H': 'H', 'He': 'He', 'C': 'C', 'N': 'N', 'O': 'O', 'F': 'F',
    'Ne': 'Ne', 'Na': 'Na', 'Mg': 'Mg', 'Al': 'Al', 'Sk': 'Si',
    'P': 'P', 'S': 'S', 'Cl': 'Cl', 'Ar': 'Ar', 'K': 'K', 'Ca': 'Ca',
    'Sc': 'Sc', 'Tk': 'Ti', 'V': 'V', 'Cr': 'Cr', 'Mn': 'Mn',
    'Fe': 'Fe', 'Co': 'Co', 'Nk': 'Ni', 'Cu': 'Cu', 'Zn': 'Zn',
    'Ba': 'Ba',
}

# element symbol -> SPECIES name used in SCL_<SPECIES>_ABUND keys.
ELEM_TO_SPECIES = {
    'H': 'HYD', 'He': 'HE', 'C': 'CARB', 'N': 'NIT', 'O': 'OXY',
    'F': 'FLU', 'Ne': 'NEON', 'Na': 'SOD', 'Mg': 'MAG', 'Al': 'ALUM',
    'Si': 'SIL', 'P': 'PHOS', 'S': 'SUL', 'Cl': 'CHL', 'Ar': 'ARG',
    'K': 'POT', 'Ca': 'CAL', 'Sc': 'SCAN', 'Ti': 'TIT', 'V': 'VAN',
    'Cr': 'CHRO', 'Mn': 'MAN', 'Fe': 'IRON', 'Co': 'COB', 'Ni': 'NICK',
    'Cu': 'COP', 'Zn': 'ZINC', 'Ba': 'BAR',
}

ROMAN_SUFFIXES = ('VIII', 'SEV', 'SIX', 'III', 'IV', 'IX', 'II', '2',
                  'I', 'V', 'X')

# Files cmf_flux / batobs.sh resolve through ../ (DIR_NAME from ../RVTJ).
REQUIRED_LINKS = ['VADAT', 'RVTJ', 'batch.sh', 'MODEL_SPEC']
OPTIONAL_LINK_GLOBS = ['MODEL', 'FGINEXT_GS', 'POP*', 'MOD_SUM_DATA']


def ion_label_to_element(label):
    """Map a MODEL_SPEC ion label (e.g. Sk2, FeIII, HI) to an element
    symbol, or None if unrecognised."""
    for abr in sorted(ABR_TO_ELEM, key=len, reverse=True):
        if label.startswith(abr):
            rest = label[len(abr):]
            if rest in ROMAN_SUFFIXES:
                return ABR_TO_ELEM[abr]
    return None


def parse_model_spec(path):
    """Return ordered unique element symbols present in MODEL_SPEC."""
    elems = []
    for line in path.read_text().splitlines():
        m = re.search(r'\[([A-Za-z0-9]+)_ISF\]', line)
        if m:
            elem = ion_label_to_element(m.group(1))
            if elem and elem not in elems:
                elems.append(elem)
    return elems


def parse_obs_slurm(path):
    """Harvest SLURM resources and the environment block from the
    pipeline-generated run_cmfgen_OBS.slurm so ablation jobs run the
    same cmf_flux build with the same resources."""
    info = {'qos': None}
    env_lines = []
    in_env = False
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith('#SBATCH'):
            body = s[len('#SBATCH'):].strip()
            for key, opt in (('partition', '--partition'),
                             ('time', '--time'),
                             ('cpus', '--cpus-per-task'),
                             ('mem', '--mem'),
                             ('qos', '--qos')):
                if body.startswith(opt + '='):
                    info[key] = body.split('=', 1)[1].strip()
            continue
        if s.startswith('module purge'):
            in_env = True
        if in_env:
            env_lines.append(line)
            if s.startswith('export atomic='):
                in_env = False
    if not env_lines:
        raise ValueError(f'no module/export block found in {path}')
    info['env_block'] = '\n'.join(env_lines)
    for req in ('partition', 'time', 'cpus', 'mem'):
        if req not in info:
            raise ValueError(f'{path} lacks #SBATCH {req}; cannot inherit')
    m = re.search(r'export cmfdist="?([^"\n]+)"?', info['env_block'])
    if not m:
        raise ValueError(f'{path} lacks "export cmfdist="')
    info['cmfdist'] = m.group(1)
    return info


def config_notify_email(config_path):
    """notification_email from config_resolved.yaml (cmfgen section),
    matching the pipeline OBS job's behaviour."""
    import yaml
    with open(config_path) as f:
        cfg = yaml.safe_load(f) or {}
    return str((cfg.get('cmfgen', {}) or {}).get(
        'notification_email', '') or '').strip()


def find_config(cmfgen_dir):
    """Locate config_resolved.yaml for prepare_CMF_PARAM.py
    --vadat-overrides (needed so DO_NONMONO_KNOP / DO_NONMONO_SOB are
    stamped when the main run used them)."""
    for cand in (cmfgen_dir / 'config_resolved.yaml',
                 cmfgen_dir.parent / 'config_resolved.yaml',
                 cmfgen_dir.parent.parent / 'config_resolved.yaml'):
        if cand.is_file():
            return cand
    return None


def env_from_config(config_path):
    """Fallback when run_cmfgen_OBS.slurm is absent (SKIP_OBS configs):
    synthesise the OBS environment from config_resolved.yaml.

    Partition/time follow the hardcoded submission rules (AGENTS.md
    section 2.3.1): submission: debug -> interactive/01:59:00, anything
    else -> small/23:59:00.  install_dir gives the same cmf_flux build
    the pipeline used."""
    import yaml
    with open(config_path) as f:
        cfg = yaml.safe_load(f) or {}
    cmfgen_cfg = cfg.get('cmfgen', {}) or {}
    cmfdist = cmfgen_cfg.get('install_dir')
    if not cmfdist:
        raise ValueError(f'{config_path} lacks cmfgen.install_dir')
    if str(cfg.get('submission', 'normal')).lower() == 'debug':
        partition, time_lim = 'interactive', '01:59:00'
    else:
        partition, time_lim = 'small', '23:59:00'
    cpus = str(cmfgen_cfg.get('cpus_per_task', 8))
    mem = str(cmfgen_cfg.get('mem', '16G'))
    atomic = cmfgen_cfg.get('atomic_dir', '/u/sbrennan/atomic')
    project_dir = Path(__file__).resolve().parent.parent
    env_block = f"""module purge
module load intel/2024.0
module load mkl/2026.0
module load python-waterboa/2024.06
export cmfdist="{cmfdist}"
export ATOMIC="{atomic}"
export PROJECT_DIR="{project_dir}"
export CMFGEN_MODEL_TYPE="static_snapshot"
export SCRATCH_DIR="SCRATCH/"
mkdir -p SCRATCH/
export atomic="{atomic}\""""
    return {'partition': partition, 'time': time_lim, 'cpus': cpus,
            'mem': mem, 'qos': None, 'cmfdist': cmfdist,
            'env_block': env_block}


SLURM_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={model}_SDEC_{label}
#SBATCH --output=cmf_SDEC_{label}_%j.out
#SBATCH --error=cmf_SDEC_{label}_%j.err
#SBATCH --time={time}
#SBATCH --partition={partition}
{qos_line}#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --chdir={workdir}

set -euo pipefail

echo "SDEC {label} OBS job started on $(date)"

{env_block}

cp -f "$cmfdist/obs_template/batobs.sh" batobs.sh
cp -f "$cmfdist/obs_template/CMF_FLUX_PARAM_INIT" CMF_FLUX_PARAM_INIT
chmod +x batobs.sh

python3 "{scripts_dir}/prepare_CMF_PARAM.py" \\
  --template CMF_FLUX_PARAM_INIT --output CMF_FLUX_PARAM_INIT \\
  --vadat-overrides "{config}"

{ablation_block}

# Requested task size. SLURM_CPUS_ON_NODE is the whole node on general.
CPUS=${{SLURM_CPUS_PER_TASK:-{cpus}}}
export OMP_NUM_THREADS=$CPUS
export MKL_NUM_THREADS=$CPUS
export OMP_STACKSIZE=512M
export KMP_STACKSIZE=512M

for f in ../RVTJ ../VADAT ../MODEL_SPEC ../batch.sh; do
  [[ -e "$f" ]] || {{ echo "ERROR: missing $f (via SDEC/ links)" >&2; exit 2; }}
done
compgen -G "../POP*" >/dev/null || {{ echo "ERROR: no POP* files via ../" >&2; exit 2; }}

mkdir -p SCRATCH
for _f in POINT1 POINT2 SCRTEMP STEQ_VALS J_COMP JEW EWDATA CFDAT_OUT CONT_FREQ; do
  [[ -e "$_f" && ! -L "$_f" ]] && mv "$_f" SCRATCH/ 2>/dev/null || true
  ln -sf "SCRATCH/$_f" "$_f" 2>/dev/null || true
done

./batobs.sh

if [[ ! -s obs_fin_10 ]]; then
  echo "ERROR: batobs.sh did not produce obs_fin_10" >&2
  tail -n 80 OUT_FLUX batobs.log 2>/dev/null >&2 || true
  exit 3
fi

{verify_block}

# Descriptive copy at SDEC/ level: baseline has no suffix, ablations
# carry the removed element, the all-species run carries _CONT.
cp -f obs_fin_10 "{sdec_dir}/{spec_file}"
echo "spectrum copy: {spec_file}"

echo "SDEC {label} OBS job finished on $(date)"
"""

ABLATION_BLOCK = """# Ablation: scale {species} populations so its opacity vanishes.
python3 - <<'PYEOF'
key = 'SCL_{species}_ABUND'
lines = open('CMF_FLUX_PARAM_INIT').read().splitlines()
out = [l for l in lines if '[' + key + ']' not in l]
out.append('{scale}   [' + key + ']   !SDEC ablation of {elem}')
open('CMF_FLUX_PARAM_INIT', 'w').write('\\n'.join(out) + '\\n')
print('SDEC ablation injected:', key, '=', '{scale}')
PYEOF
"""

CONTINUUM_BLOCK = """# Continuum run: scale every species' abundance so all bound-bound
# and bound-free opacity vanishes, leaving electron scattering only.
python3 - <<'PYEOF'
species = '{species_list}'.split()
lines = open('CMF_FLUX_PARAM_INIT').read().splitlines()
out = [l for l in lines if '[SCL_' not in l]
for sp in species:
    out.append('{scale}   [SCL_' + sp + '_ABUND]   !SDEC continuum')
open('CMF_FLUX_PARAM_INIT', 'w').write('\\n'.join(out) + '\\n')
print('SDEC continuum: scaled', len(species), 'species')
PYEOF
"""

VERIFY_BLOCK = """# Verify the SCL key was honoured: cmf_flux writes a warning to LUER
# when it scales a species abundance.
if ! grep -qsi "Abundance of species.*scaled" OUT_FLUX batobs.log error.log OUT_PARAMS 2>/dev/null; then
  echo "ERROR: SCL_{species}_ABUND not honoured (no 'species scaled' warning)" >&2
  exit 4
fi
"""

VERIFY_CONT_BLOCK = """# Verify every SCL_*_ABUND was honoured (one warning per species).
# '|| true': cat exits 1 when some log files are absent and grep -c
# exits 1 on zero matches -- either would abort the script under
# 'set -e' before the check below runs.
n_scaled=$(cat OUT_FLUX batobs.log error.log OUT_PARAMS 2>/dev/null | grep -ci "Abundance of species.*scaled" || true)
if [[ "${{n_scaled:-0}}" -lt {nspec} ]]; then
  echo "ERROR: expected {nspec} scaled-species warnings, got ${{n_scaled:-0}}" >&2
  exit 4
fi
"""

PLOT_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={model}_SDEC_plot
#SBATCH --output=cmf_SDEC_plot_%j.out
#SBATCH --error=cmf_SDEC_plot_%j.err
#SBATCH --time=00:30:00
#SBATCH --partition={partition}
{qos_line}#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --chdir={workdir}

set -euo pipefail

{env_block}

python3 "{scripts_dir}/plot_sdec_decomposition.py" "{cmfgen_dir}" \\
  --out "{sdec_dir}/sdec_decomposition.png" --write-csv \\
  --species {species_args}

{notify_block}"""

NOTIFY_BLOCK = """
# Notify by email (same mailx pattern as the pipeline OBS job)
SDEC_PLOT="{sdec_dir}/sdec_decomposition.png"
SDEC_PLOT_LOG="{sdec_dir}/sdec_decomposition_log.png"
NOTIFY_EMAIL="{notify_email}"
_attach_args=()
[[ -f "$SDEC_PLOT" ]] && _attach_args+=(-a "$SDEC_PLOT")
[[ -f "$SDEC_PLOT_LOG" ]] && _attach_args+=(-a "$SDEC_PLOT_LOG")
_mail_subject="SDEC decomposition plot: {model}"
_mail_body=$(printf "SDEC species-ablation decomposition complete.\\nModel: %s\\nPlot: %s\\nSpecies: %s\\n" "{cmfgen_dir}" "$SDEC_PLOT" "{species_list}")
# Compute nodes cannot write to the local postfix maildrop (postdrop
# "Permission denied"), so sendmail-mode mailx hangs until the timeout.
# Try direct SMTP to the MPCDF relay first; fall back to sendmail mode.
_smtp_relay="${{CMFGEN_SMTP_RELAY:-smtp://hpcmailrelay.mpcdf.mpg.de:25}}"
printf '%s\\n' "$_mail_body" | timeout 30 mailx -S "mta=${{_smtp_relay}}" "${{_attach_args[@]}}" -s "$_mail_subject" "$NOTIFY_EMAIL" 2>/dev/null \\
    || printf '%s\\n' "$_mail_body" | timeout 30 mailx "${{_attach_args[@]}}" -s "$_mail_subject" "$NOTIFY_EMAIL" 2>/dev/null \\
    && echo "Notification email sent to $NOTIFY_EMAIL" \\
    || echo "WARNING: notification email failed or timed out (non-fatal)"
"""


def make_sdec_links(sdec_dir, cmfgen_dir, strict=True):
    """Populate SDEC/ with the ../-relative symlinks cmf_flux and
    batobs.sh resolve from inside each species workdir.

    strict=False links unconditionally -- used when SDEC jobs are
    submitted with a dependency before CMFGEN output exists."""
    linked = []
    for name in REQUIRED_LINKS:
        src = cmfgen_dir / name
        if strict and not src.exists():
            raise FileNotFoundError(f'required file missing: {src}')
        dst = sdec_dir / name
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(f'../{name}')
        linked.append(name)
    for pat in OPTIONAL_LINK_GLOBS:
        for src in sorted(cmfgen_dir.glob(pat)):
            if not src.is_file():
                continue
            dst = sdec_dir / src.name
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(f'../{src.name}')
            linked.append(src.name)
    return linked


def write_job_script(path, **kw):
    path.write_text(SLURM_TEMPLATE.format(**kw))


def main():
    ap = argparse.ArgumentParser(
        description='Submit SDEC species-ablation cmf_flux OBS jobs')
    ap.add_argument('cmfgen_dir',
                    help='converged CMFGEN run dir (contains RVTJ, VADAT, '
                         'batch.sh, MODEL_SPEC, POP*)')
    ap.add_argument('--species', nargs='+', default=None,
                    help='element symbols to ablate (default: all in '
                         'MODEL_SPEC)')
    ap.add_argument('--partition', default='interactive',
                    help='SLURM partition for the ablation jobs '
                         '(default: interactive -- a cmf_flux formal '
                         'solution takes ~30-45 min, well inside the '
                         '2h limit; pass "inherit" to use the OBS '
                         'job\'s partition)')
    ap.add_argument('--max-parallel', type=int, default=2,
                    help='max concurrent ablation jobs (default 2)')
    ap.add_argument('--scale', default='1.0D-10',
                    help='SCL_<species>_ABUND value (default 1.0D-10)')
    ap.add_argument('--continuum', action='store_true',
                    help='also submit a CONT job with every species in '
                         'MODEL_SPEC scaled -- the OBS-computed '
                         'continuum (electron scattering only)')
    ap.add_argument('--depends-on', default=None,
                    help='comma-separated job IDs the first jobs must '
                         'wait for (afterok)')
    ap.add_argument('--plot-depends-on', default=None,
                    help='comma-separated extra job IDs the plot job must '
                         'wait for (e.g. the pipeline OBS job that writes '
                         'the baseline obs_fin_10). Implies the baseline '
                         'arrives via those jobs, so no ALL job is added')
    ap.add_argument('--no-precheck', action='store_true',
                    help='skip input-file existence checks and link '
                         'unconditionally -- for chained submission '
                         'before CMFGEN output exists')
    ap.add_argument('--no-submit-all', action='store_true',
                    help='do not submit an unablated ALL baseline job; '
                         'the plot then falls back to obs/obs_fin_10, '
                         'which may carry a different I_NEXT generation '
                         'than the ablation runs')
    ap.add_argument('--skip-existing', action='store_true',
                    help='skip jobs whose SDEC/<label>/obs_fin_10 '
                         'already exists (idempotent resubmission)')
    ap.add_argument('--no-plot', action='store_true',
                    help='skip the final decomposition plot job')
    ap.add_argument('--notify-email', default=None, metavar='ADDR',
                    help='email the decomposition plot when the plot job '
                         'finishes (default: cmfgen.notification_email '
                         'from config_resolved.yaml; pass empty string '
                         'to disable)')
    ap.add_argument('--dry-run', action='store_true',
                    help='prepare directories and scripts but do not '
                         'call sbatch')
    args = ap.parse_args()

    cmfgen_dir = Path(args.cmfgen_dir).resolve()
    if not cmfgen_dir.is_dir():
        sys.exit(f'ERROR: not a directory: {cmfgen_dir}')

    if not args.no_precheck:
        for f in REQUIRED_LINKS:
            if not (cmfgen_dir / f).exists():
                sys.exit(f'ERROR: {f} not found in {cmfgen_dir} -- not a '
                         'converged CMFGEN run dir')
        if not (cmfgen_dir / 'RVTJ').stat().st_size:
            sys.exit(f'ERROR: {cmfgen_dir}/RVTJ is empty')
        pop_files = [p for p in sorted(cmfgen_dir.glob('POP*'))
                     if p.is_file()]
        if not pop_files:
            sys.exit('ERROR: no POP<species> files -- cmf_flux needs '
                     'them (zero species -> N_LINE_FREQ=0 segfault)')

    config = find_config(cmfgen_dir)
    if config is None:
        sys.exit('ERROR: config_resolved.yaml not found near '
                 f'{cmfgen_dir} -- required so DO_NONMONO_KNOP/'
                 'DO_NONMONO_SOB flags match the main run')

    obs_slurm = cmfgen_dir / 'run_cmfgen_OBS.slurm'
    if obs_slurm.is_file():
        try:
            obs_info = parse_obs_slurm(obs_slurm)
        except ValueError as e:
            sys.exit(f'ERROR: {e}')
    else:
        try:
            obs_info = env_from_config(config)
            print(f'  no run_cmfgen_OBS.slurm; env synthesised from '
                  f'{config} (cmfdist={obs_info["cmfdist"]}, '
                  f'{obs_info["partition"]})')
        except ValueError as e:
            sys.exit(f'ERROR: {obs_slurm} not found and config '
                     f'fallback failed: {e}')

    if args.partition != 'inherit':
        obs_info['partition'] = args.partition
        # A partition-bound QOS (e.g. n0004 on small) is invalid or
        # counterproductive on another partition -- drop it.
        obs_info['qos'] = None
        # The interactive partition has a hard 2h wall (AGENTS.md
        # section 2.3); a longer inherited limit would be rejected.
        if args.partition == 'interactive' and \
                obs_info['time'] > '01:59:00':
            obs_info['time'] = '01:59:00'

    model_name = cmfgen_dir.parent.parent.name \
        if cmfgen_dir.parent.name == 'cmfgen' else cmfgen_dir.name

    # ---- species selection ----
    if args.species:
        elems = []
        for s in args.species:
            canon = s[0].upper() + s[1:].lower()
            if canon not in ELEM_TO_SPECIES:
                sys.exit(f'ERROR: unknown element {s!r}')
            if canon not in elems:
                elems.append(canon)
    else:
        ms = cmfgen_dir / 'MODEL_SPEC'
        elems = parse_model_spec(ms) if ms.is_file() else []
        if not elems:
            sys.exit('ERROR: no *_ISF ions found in MODEL_SPEC '
                     + ('' if ms.is_file() else '(file missing -- pass '
                        '--species explicitly for chained submission)'))
    ms_file = cmfgen_dir / 'MODEL_SPEC'
    if ms_file.is_file():
        missing = [e for e in elems if e not in parse_model_spec(ms_file)]
        if missing:
            print(f'WARNING: {missing} not present in MODEL_SPEC; '
                  'ablation will be a no-op for them')

    scripts_dir = Path(__file__).resolve().parent
    sdec_dir = cmfgen_dir / 'SDEC'
    sdec_dir.mkdir(exist_ok=True)
    linked = make_sdec_links(sdec_dir, cmfgen_dir,
                             strict=not args.no_precheck)
    print(f'SDEC dir: {sdec_dir}')
    print(f'  ../ links: {", ".join(linked)}')

    # ---- job list: ALL baseline + one per element ----
    baseline = cmfgen_dir / 'obs' / 'obs_fin_10'
    jobs = []  # (label, species_name_or_None, elem_or_None)
    if args.no_submit_all:
        if not baseline.is_file() and not args.plot_depends_on:
            sys.exit('ERROR: obs/obs_fin_10 missing and --no-submit-all '
                     'given -- no baseline spectrum for the plot')
        print('  --no-submit-all: baseline falls back to obs/obs_fin_10')
    else:
        jobs.append(('ALL', None, None))
    if args.continuum:
        # CONT ablates every species in MODEL_SPEC regardless of any
        # --species subset -- the continuum must remove all opacity.
        ms = cmfgen_dir / 'MODEL_SPEC'
        cont_elems = (parse_model_spec(ms) if ms.is_file()
                      else list(elems))
        jobs.append(('CONT', 'CONT', 'CONT'))
        jobs_cont_species = [ELEM_TO_SPECIES[e] for e in cont_elems]
    for e in elems:
        jobs.append((e, ELEM_TO_SPECIES[e], e))

    # spec_tag identifies model + epoch so the copied spectrum file
    # carries its provenance: output_spectrum_<model>_<epoch>[_<label>].
    epoch = cmfgen_dir.name
    spec_tag = f'output_spectrum_{model_name}_{epoch}'

    def _spec_file(label):
        return spec_tag + ('' if label == 'ALL' else f'_{label}') + '.dat'

    if args.skip_existing:
        kept = []
        for label, species, elem in jobs:
            existing = sdec_dir / label / 'obs_fin_10'
            if existing.is_file() and existing.stat().st_size:
                # Backfill the descriptive copy for runs that predate
                # the naming convention.
                dst = sdec_dir / _spec_file(label)
                if not dst.exists():
                    dst.write_bytes(existing.read_bytes())
                    print(f'  {label}: wrote descriptive copy {dst.name}')
                print(f'  {label}: obs_fin_10 exists -- skipped '
                      '(--skip-existing)')
            else:
                kept.append((label, species, elem))
        jobs = kept

    qos_line = (f"#SBATCH --qos={obs_info['qos']}\n"
                if obs_info.get('qos') else '')

    # ---- write scripts ----
    for label, species, elem in jobs:
        wdir = sdec_dir / label
        wdir.mkdir(exist_ok=True)
        spec_file = _spec_file(label)
        if label == 'CONT':
            abl = CONTINUUM_BLOCK.format(
                species_list=' '.join(jobs_cont_species),
                scale=args.scale)
            ver = VERIFY_CONT_BLOCK.format(nspec=len(jobs_cont_species))
        elif species:
            abl = ABLATION_BLOCK.format(species=species, scale=args.scale,
                                        elem=elem)
            ver = VERIFY_BLOCK.format(species=species)
        else:
            abl = ver = ''
        write_job_script(
            wdir / f'run_sdec_{label}.slurm',
            model=model_name, label=label,
            time=obs_info['time'], partition=obs_info['partition'],
            qos_line=qos_line, cpus=obs_info['cpus'], mem=obs_info['mem'],
            workdir=wdir, env_block=obs_info['env_block'],
            scripts_dir=scripts_dir, config=config,
            ablation_block=abl, verify_block=ver,
            sdec_dir=sdec_dir, spec_file=spec_file)
    notify_email = (args.notify_email if args.notify_email is not None
                    else config_notify_email(config))
    notify_block = ''
    if notify_email:
        notify_block = NOTIFY_BLOCK.format(
            sdec_dir=sdec_dir, notify_email=notify_email,
            cmfgen_dir=cmfgen_dir, model=model_name,
            species_list=','.join(elems))

    plot_path = sdec_dir / 'run_sdec_plot.slurm'
    if not args.no_plot:
        plot_path.write_text(PLOT_TEMPLATE.format(
            model=model_name, partition=obs_info['partition'],
            qos_line=qos_line, workdir=sdec_dir,
            env_block=obs_info['env_block'], scripts_dir=scripts_dir,
            cmfgen_dir=cmfgen_dir, sdec_dir=sdec_dir,
            notify_block=notify_block, species_args=' '.join(elems)))

    print(f'  jobs prepared: {[j[0] for j in jobs]}'
          + ('' if args.no_plot else ' + plot'))
    if args.dry_run:
        print('DRY RUN: nothing submitted')
        return 0

    # ---- submit with i-N dependency so at most N run concurrently ----
    extra_dep = args.depends_on or ''
    job_ids = []
    for i, (label, _, _) in enumerate(jobs):
        deps = []
        if i >= args.max_parallel:
            deps.append(job_ids[i - args.max_parallel])
        elif extra_dep:
            deps.append(extra_dep)
        cmd = ['sbatch', '--parsable']
        if deps:
            cmd.append('--dependency=afterok:' + ':'.join(deps))
        cmd.append(str(sdec_dir / label / f'run_sdec_{label}.slurm'))
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f'ERROR: sbatch failed for {label}: {r.stderr.strip()}')
        jid = r.stdout.strip().split(';')[0]
        job_ids.append(jid)
        dep_str = f' (afterok:{":".join(deps)})' if deps else ''
        print(f'  {label}: job {jid}{dep_str}')
        time.sleep(0.2)

    plot_jid = None
    if not args.no_plot:
        plot_deps = list(job_ids)
        if args.plot_depends_on:
            plot_deps.extend(j for j in args.plot_depends_on.split(',')
                             if j.strip())
        cmd = ['sbatch', '--parsable']
        if plot_deps:
            cmd.append('--dependency=afterok:' + ':'.join(plot_deps))
        cmd.append(str(plot_path))
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f'WARNING: plot job sbatch failed: {r.stderr.strip()}',
                  file=sys.stderr)
        else:
            plot_jid = r.stdout.strip().split(';')[0]
            print(f'  plot: job {plot_jid} (afterok all)')

    log = sdec_dir / 'submit_log.txt'
    with open(log, 'w') as f:
        f.write(f'SDEC ablation submission {time.strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(f'cmfgen_dir: {cmfgen_dir}\n')
        f.write(f'cmfdist:    {obs_info["cmfdist"]}\n')
        f.write(f'config:     {config}\n')
        f.write(f'scale:      {args.scale}\n')
        f.write(f'max_parallel: {args.max_parallel}\n')
        f.write(f'notify_email: {notify_email or "(none)"}\n')
        for (label, species, _), jid in zip(jobs, job_ids):
            if label == 'CONT':
                key = 'all SCL_*_ABUND'
            elif species:
                key = f'SCL_{species}_ABUND'
            else:
                key = 'unablated'
            f.write(f'  {label}: job {jid}  ({key})\n')
        if plot_jid:
            f.write(f'  plot: job {plot_jid}\n')
    print(f'  log: {log}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
