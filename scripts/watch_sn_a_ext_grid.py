#!/usr/bin/env python3
"""Watch the SN_A_X_E07Ni05_R* STELLA grid, trigger SN_B resubmission,
and update the two-event A&A plot when everything completes.

Stage 1 - poll SLURM until the four SN_A extended-envelope jobs leave
the queue, then score each model against the SN 2024pba photometry
using the MJD < 60544 fit window (FIT_MAX_MJD in
plot_sn_a_pba_grid_vs_obs.py). The two-event plot's SN_A_MODELS list is
rewritten so the best model is the featured solid curve, and the plot
is regenerated.

Stage 2 - if the best new model beats the incumbent SN_A_S_E07Ni05,
the three SN_B Mdot-sequence configs (X_M10p0_V1e4_W3e-3 / W1e-2 /
W3e-2) are cloned onto the winning SN_A ejecta
(x_model_from_sn.source_dir) and submitted through
scripts/run_model_grid.py. The watcher then keeps polling until those
SN_B jobs finish, swaps them into the plot's SN_B_MODELS list, and
regenerates the plot.

All submission goes through the orchestrator - never raw sbatch.
State is kept in logs/sn_a_ext_grid_watch.state so a restart resumes
at the right stage without double-submitting.

Run:  nohup python scripts/watch_sn_a_ext_grid.py &
Log:  logs/sn_a_ext_grid_watch.log
"""
import os
import re
import subprocess
import sys
import time
from datetime import datetime

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
SCRIPTS = os.path.join(BASE, "scripts")
GRID_DIR = os.path.join(BASE, "SN_A_Y_grid")
SNB_GRID_DIR = os.path.join(BASE, "SN_B_v2_T6k_grid")
LOG_FILE = os.path.join(BASE, "logs", "sn_a_ext_grid_watch.log")
STATE_FILE = os.path.join(BASE, "logs", "sn_a_ext_grid_watch.state")
PLOT_SCRIPT = os.path.join(
    BASE, "A&A_paper", "plotting_scripts", "plot_sn2024pba_two_event.py")

NEW_MODELS = [
    "SN_A_X_E07Ni05_R400",
    "SN_A_X_E07Ni05_R400M20",
    "SN_A_X_E07Ni05_R600",
    "SN_A_X_E07Ni05_R600M20",
]
INCUMBENT = "SN_A_S_E07Ni05"
SNB_EJECTA_OLD = "SN_A_g_T20k_Ni06"
SNB_EJECTA_OLD_DIR = "archive/old_sn_a_models/SN_A/SN_A_g_T20k_Ni06"

# The SN_B Mdot sequence featured in the two-event figure.
# (template name, legend label, best flag)
SNB_TEMPLATES = [
    ("X_M10p0_V1e4_W3e-3", r"$\dot{M}_{\rm wind}=3\times10^{-3}$", False),
    ("X_M10p0_V1e4_W1e-2", r"$\dot{M}_{\rm wind}=1\times10^{-2}$", False),
    ("X_M10p0_V1e4_W3e-2", r"$\dot{M}_{\rm wind}=3\times10^{-2}$", True),
]
# SN_B entries kept verbatim from the current plot list.
SNB_TAIL = [
    ("T6k_M7.0_V10000", "SN_B_v2_T6k_grid/T6k_M7.0_V10000", False),
    ("T6k_M10.0_V10000", "SN_B_v2_T6k_grid/T6k_M10.0_V10000", False),
    ("T6k_M7.0_V8000", "SN_B_v2_T6k_grid/T6k_M7.0_V8000", False),
    ("T6k_M10.0_V8000", "SN_B_v2_T6k_grid/T6k_M10.0_V8000", False),
]

POLL_SEC = 600
MAX_WAIT_H = 40
# Model must reach at least this many days post-explosion for the
# MJD<60544 window (peak ~20 d + 36 d post-peak) to be fully covered.
MIN_COVERAGE_DAYS = 55.0

sys.path.insert(0, SCRIPTS)
from plot_sn_a_pba_grid_vs_obs import (  # noqa: E402
    OBS_FILE, OBS_PEAK_MJD, PEAK_BAND, FIT_MAX_MJD,
    load_model, load_observations, score_model,
)


def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as fh:
        fh.write(line + "\n")


def read_state():
    state = {}
    if os.path.exists(STATE_FILE):
        with open(STATE_FILE) as fh:
            for line in fh:
                if "=" in line:
                    k, v = line.strip().split("=", 1)
                    state[k] = v
    return state


def write_state(**kw):
    state = read_state()
    state.update(kw)
    with open(STATE_FILE, "w") as fh:
        for k, v in state.items():
            fh.write(f"{k}={v}\n")


def queued_jobs(prefixes):
    out = subprocess.run(
        ["squeue", "-u", os.environ.get("USER", "sbrennan"),
         "-h", "-o", "%j"],
        capture_output=True, text=True).stdout
    names = set(out.split())
    return sorted(n for n in names
                  if any(n == f"{p}_stella" for p in prefixes))


def wait_for_jobs(prefixes, stage):
    waited = 0.0
    while True:
        jobs = queued_jobs(prefixes)
        if not jobs:
            return True
        if waited >= MAX_WAIT_H * 3600:
            log(f"TIMEOUT after {MAX_WAIT_H} h in {stage}; "
                f"still queued: {jobs}")
            return False
        log(f"[{stage}] still queued/running: {jobs}; "
            f"sleeping {POLL_SEC}s")
        time.sleep(POLL_SEC)
        waited += POLL_SEC


def load_scored_obs():
    obs = load_observations(OBS_FILE)
    for band in ("g", "r", "i"):
        for k in ("det", "nd"):
            obs[band][k] = obs[band][k].copy()
            obs[band][k]["phase"] = obs[band][k]["MJD"] - OBS_PEAK_MJD
    return obs


def score_one(name, obs):
    """Return combined rms under the fit window, or None if unusable."""
    df, _ = load_model(os.path.join(GRID_DIR, name))
    if df is None or not len(df):
        log(f"  {name}: no .tt data - skipped")
        return None
    if df["time"].max() < MIN_COVERAGE_DAYS:
        log(f"  {name}: .tt only reaches {df['time'].max():.1f} d "
            f"(< {MIN_COVERAGE_DAYS} d) - incomplete, skipped")
        return None
    t_ipeak = df.loc[df[PEAK_BAND].idxmin(), "time"]
    sc = score_model(df, obs, -t_ipeak)
    n = sum(sc.get(f"n_{b}", 0) for b in ("g", "r", "i"))
    log(f"  {name}: rms={sc['combined']:.3f} "
        f"(g={sc.get('g', float('nan')):.3f} "
        f"r={sc.get('r', float('nan')):.3f} "
        f"i={sc.get('i', float('nan')):.3f}, n={n})")
    return sc["combined"]


def set_model_list(var, block):
    """Replace `var = [ ... ]` in the two-event plot script."""
    with open(PLOT_SCRIPT) as fh:
        txt = fh.read()
    new, nsub = re.subn(rf"(?ms)^{var} = \[.*?^\]",
                        lambda m: block.rstrip(), txt)
    if nsub != 1:
        log(f"  WARNING: could not rewrite {var} in {PLOT_SCRIPT}")
        return
    with open(PLOT_SCRIPT, "w") as fh:
        fh.write(new)
    log(f"  rewrote {var} in {os.path.basename(PLOT_SCRIPT)}")


def regen_plot():
    res = subprocess.run([sys.executable, PLOT_SCRIPT],
                         capture_output=True, text=True)
    if res.returncode != 0:
        log(f"  PLOT REGEN FAILED: {res.stderr[-1000:]}")
    else:
        log("  plot regenerated: "
            "A&A_paper/plots/sn2024pba_two_event.png/.svg")


def _entry(name, mdir, label, best):
    lab = f'r"{label}"' if label.startswith("$") else f'"{label}"'
    return (f'    dict(name="{name}", dir="{mdir}",\n'
            f'         label={lab}, best={best}),')


def update_plot_sn_a(best):
    """SN_A_MODELS: best model solid, others faded."""
    entries = [_entry(best, f"SN_A_Y_grid/{best}",
                      f"{best} (best fit)", True),
               _entry(SNB_EJECTA_OLD, SNB_EJECTA_OLD_DIR,
                      f"{SNB_EJECTA_OLD} (ejecta for SN_B)", True)]
    for name in NEW_MODELS:
        if name != best:
            entries.append(_entry(name, f"SN_A_Y_grid/{name}", name, False))
    for name in (INCUMBENT, "SN_A_X_E06", "SN_A_S_11dh", "SN_A_X_E08"):
        if name != best:
            entries.append(_entry(name, f"SN_A_Y_grid/{name}", name, False))
    set_model_list("SN_A_MODELS",
                   "SN_A_MODELS = [\n" + "\n".join(entries) + "\n]")


def update_plot_sn_b(tag):
    """SN_B_MODELS: Mdot sequence on the new ejecta + unchanged tail."""
    entries = []
    for tpl, label, best in SNB_TEMPLATES:
        new_name = f"{tpl}_on{tag}"
        entries.append(_entry(
            new_name, f"SN_B_v2_T6k_grid/{new_name}", label, best))
    for name, mdir, best in SNB_TAIL:
        entries.append(_entry(name, mdir, name, best))
    set_model_list("SN_B_MODELS",
                   "SN_B_MODELS = [\n" + "\n".join(entries) + "\n]")


def clone_snb_configs(tag, sn_a_name, sn_a_dir):
    """Clone the 3 SN_B Mdot configs onto the winning SN_A ejecta."""
    cfgs = []
    for tpl, _, _ in SNB_TEMPLATES:
        src_cfg = os.path.join(SNB_GRID_DIR, tpl, "config.yaml")
        new_name = f"{tpl}_on{tag}"
        new_dir = os.path.join(SNB_GRID_DIR, new_name)
        os.makedirs(new_dir, exist_ok=True)
        with open(src_cfg) as fh:
            txt = fh.read()
        txt = txt.replace(SNB_EJECTA_OLD, sn_a_name)
        txt = re.sub(r"^model_name:.*$",
                     f"model_name: {new_name}", txt, flags=re.M)
        txt = re.sub(r"^(\s*)source_dir:.*$",
                     rf"\1source_dir: {sn_a_dir}", txt, flags=re.M)
        txt = txt.replace(
            f"model_dir: ./SN_B_v2_T6k_grid/{tpl}",
            f"model_dir: ./SN_B_v2_T6k_grid/{new_name}")
        dst = os.path.join(new_dir, "config.yaml")
        with open(dst, "w") as fh:
            fh.write(txt)
        cfgs.append(dst)
        log(f"  wrote {dst}")
    return cfgs


def submit_snb(cfgs):
    rel = [os.path.relpath(c, BASE) for c in cfgs]
    cmd = [sys.executable, "scripts/run_model_grid.py"] + rel
    log(f"  submitting: {' '.join(cmd)}")
    res = subprocess.run(cmd, cwd=BASE, capture_output=True, text=True)
    log(res.stdout[-3000:])
    if res.returncode != 0:
        log(f"  SUBMISSION FAILED (rc={res.returncode}): "
            f"{res.stderr[-1500:]}")
        return False
    log("  SN_B resubmission complete")
    return True


def watch_snb_stage(tag):
    """Stage 2: wait for the cloned SN_B jobs, then update the plot."""
    prefixes = [f"{tpl}_on{tag}" for tpl, _, _ in SNB_TEMPLATES]
    if not wait_for_jobs(prefixes, "sn_b"):
        write_state(DONE="TIMEOUT_SNB")
        return
    log("all SN_B jobs left the queue; waiting 120 s for file flush")
    time.sleep(120)
    update_plot_sn_b(tag)
    regen_plot()
    write_state(DONE="COMPLETE")
    log("done: SN_B models on new ejecta are in the two-event plot")


def main():
    state = read_state()
    log(f"watcher start; state={state}; fit window MJD < {FIT_MAX_MJD}")

    if "SNB_TAG" in state:
        watch_snb_stage(state["SNB_TAG"])
        return

    if not wait_for_jobs(NEW_MODELS, "sn_a"):
        write_state(DONE="TIMEOUT_SNA")
        return

    log("all SN_A jobs left the queue; waiting 120 s for file flush")
    time.sleep(120)

    obs = load_scored_obs()
    scores = {}
    for name in NEW_MODELS:
        rms = score_one(name, obs)
        if rms is not None:
            scores[name] = rms
    incumbent_rms = score_one(INCUMBENT, obs)

    new_best = min(scores, key=scores.get) if scores else None
    ideal_found = (new_best is not None
                   and (incumbent_rms is None
                        or scores[new_best] < incumbent_rms))
    featured = new_best if ideal_found else (
        INCUMBENT if incumbent_rms is not None else new_best)
    log(f"best new model: {new_best} rms="
        f"{scores.get(new_best)}; incumbent rms={incumbent_rms}; "
        f"ideal_found={ideal_found}")
    if featured is None:
        write_state(DONE="FAIL_NO_SCORES")
        log("FAIL: no model produced a scorable lightcurve; "
            "plot left unchanged")
        return

    update_plot_sn_a(featured)
    regen_plot()

    if not ideal_found:
        write_state(DONE="NO_IDEAL")
        log("no new model beats the incumbent - "
            "no ideal model found, no SN_B resubmission")
        return

    best = new_best

    tag = best.replace("SN_A_X_E07Ni05_", "")
    sn_a_dir = os.path.join(GRID_DIR, best)
    log(f"ideal model found: {best} - cloning SN_B Mdot sequence "
        f"onto {sn_a_dir} (tag {tag})")
    cfgs = clone_snb_configs(tag, best, sn_a_dir)
    if submit_snb(cfgs):
        write_state(SNB_TAG=tag)
        watch_snb_stage(tag)
    else:
        write_state(DONE="SUBMIT_FAIL")


if __name__ == "__main__":
    main()
