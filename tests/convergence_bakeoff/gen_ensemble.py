#!/usr/bin/env python3
"""Generate the deterministic bake-off instance ensemble.

Writes instances/iNN_<regime>.in plus split.txt (tuning vs held-out).
The split is fixed here, before any tuning run; held-out instances must
not be touched until DECISION_RULE.md is committed.

Regimes (knobs on the toy model):
  easy      moderate stiffness/coupling, direct J solve
  stiff     gross-flux stiffness 1e11-1e12
  thin      tau << 1: strong non-local coupling, worst for AL*
  thick     tau >> 1: local coupling, AL* nearly exact
  noisy     tolerance-terminated Jacobi J solve (deterministic nonsmooth)
  singular  g0 near the 0.5 contraction limit on a thin grid
  lag       large lagged-store fraction rho
  mixed     stiff + thin + lag + noise combined

Usage: python3 gen_ensemble.py   (idempotent; overwrites instances/)
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))
IDIR = os.path.join(HERE, "instances")

# (name, seed, nd, nl, ntr, solveT, stiff, tau, g0, alstar, rho,
#  trscale, jsolve, jtol, jit, start, theta)
BASE = dict(nd=12, nl=4, ntr=3, solveT=1, stiff=1.0e9, tau=1.0,
            g0=0.45, alstar=0.4, rho=0.0, trscale=1.0e-15,
            jsolve=0, jtol=1.0e-8, jit=0.0, start=0.4, theta=0.1)

INSTANCES = []


def add(name, **kw):
    p = dict(BASE)
    p.update(kw)
    INSTANCES.append((name, p))


# --- easy (seeds 1-4): moderate everything -------------------------
add("i01_easy",  seed=1, stiff=1.0e8, tau=1.0,  g0=0.30, start=0.3)
add("i02_easy",  seed=2, stiff=3.0e8, tau=0.5,  g0=0.40, start=0.4)
add("i03_easy",  seed=3, stiff=1.0e8, tau=2.0,  g0=0.45, start=0.5)
add("i04_easy",  seed=4, stiff=3.0e8, tau=1.0,  g0=0.35, start=0.3)

# --- stiff (seeds 5-8): gross flux 1e11-1e12 ------------------------
add("i05_stiff", seed=5, stiff=1.0e11, tau=1.0, g0=0.40)
add("i06_stiff", seed=6, stiff=1.0e12, tau=0.5, g0=0.45)
add("i07_stiff", seed=7, stiff=3.0e11, tau=1.0, g0=0.35)
add("i08_stiff", seed=8, stiff=1.0e12, tau=2.0, g0=0.40, start=0.5)

# --- thin (seeds 1,2,3,5): tau << 1 ---------------------------------
add("i09_thin",  seed=1, stiff=1.0e9,  tau=1.0e-3, g0=0.48)
add("i10_thin",  seed=2, stiff=1.0e10, tau=1.0e-3, g0=0.45, rho=0.2)
add("i11_thin",  seed=3, stiff=1.0e9,  tau=5.0e-3, g0=0.49)
add("i12_thin",  seed=5, stiff=1.0e10, tau=1.0e-2, g0=0.48, jsolve=1,
    jtol=1.0e-7)

# --- thick (seeds 4,6,7,8): tau >> 1 --------------------------------
add("i13_thick", seed=4, stiff=1.0e9,  tau=10.0,  g0=0.45)
add("i14_thick", seed=6, stiff=1.0e10, tau=30.0,  g0=0.48)
add("i15_thick", seed=7, stiff=1.0e9,  tau=100.0, g0=0.49, rho=0.3)
add("i16_thick", seed=8, stiff=3.0e10, tau=30.0,  g0=0.47)

# --- noisy (seeds 2,3,5,7): Jacobi J solve --------------------------
add("i17_noisy", seed=2, stiff=1.0e9,  tau=1.0, g0=0.45, jsolve=1,
    jtol=1.0e-6)
add("i18_noisy", seed=3, stiff=1.0e10, tau=0.5, g0=0.40, jsolve=1,
    jtol=1.0e-7)
add("i19_noisy", seed=5, stiff=1.0e9,  tau=2.0, g0=0.45, jsolve=1,
    jtol=1.0e-9)
add("i20_noisy", seed=7, stiff=1.0e10, tau=1.0, g0=0.48, jsolve=1,
    jtol=1.0e-8, jit=1.0e-10)

# --- singular (seeds 1,4,6,8): g0 -> 0.5 on thin grids --------------
add("i21_sing",  seed=1, stiff=1.0e10, tau=1.0e-2, g0=0.499)
add("i22_sing",  seed=4, stiff=1.0e10, tau=5.0e-3, g0=0.495)
add("i23_sing",  seed=6, stiff=3.0e10, tau=1.0e-2, g0=0.497, rho=0.2)
add("i24_sing",  seed=8, stiff=1.0e10, tau=1.0e-3, g0=0.499)

# --- lag (seeds 2,3,5,6): strong lagged store -----------------------
add("i25_lag",   seed=2, stiff=1.0e9,  tau=0.1, g0=0.45, rho=0.5)
add("i26_lag",   seed=3, stiff=1.0e10, tau=1.0, g0=0.40, rho=0.8)
add("i27_lag",   seed=5, stiff=1.0e9,  tau=0.5, g0=0.45, rho=0.6,
    jsolve=1, jtol=1.0e-7)
add("i28_lag",   seed=6, stiff=1.0e10, tau=0.2, g0=0.48, rho=0.9)

# --- mixed (seeds 3,7,1,8): everything at once ----------------------
add("i29_mix",   seed=3, stiff=1.0e11, tau=1.0e-2, g0=0.48, rho=0.4,
    jsolve=1, jtol=1.0e-7)
add("i30_mix",   seed=7, stiff=1.0e11, tau=1.0e-3, g0=0.49, rho=0.5,
    jsolve=1, jtol=1.0e-8)
add("i31_mix",   seed=1, stiff=3.0e11, tau=1.0e-2, g0=0.47, rho=0.6,
    jsolve=1, jtol=1.0e-6)
add("i32_mix",   seed=8, stiff=1.0e12, tau=1.0e-3, g0=0.49, rho=0.7,
    jsolve=1, jtol=1.0e-7, jit=1.0e-11)


def fmt(p):
    return (
        "&toy\n"
        f" nml_seed = {p['seed']},\n"
        f" nml_nd = {p['nd']},\n"
        f" nml_nl = {p['nl']},\n"
        f" nml_ntr = {p['ntr']},\n"
        f" nml_solveT = {p['solveT']},\n"
        f" nml_stiff = {p['stiff']:.2e},\n"
        f" nml_tau = {p['tau']:.2e},\n"
        f" nml_g0 = {p['g0']:.3f},\n"
        f" nml_alstar = {p['alstar']:.2f},\n"
        f" nml_rho = {p['rho']:.2f},\n"
        f" nml_trscale = {p['trscale']:.1e},\n"
        f" nml_jsolve = {p['jsolve']},\n"
        f" nml_jtol = {p['jtol']:.1e},\n"
        f" nml_jit = {p['jit']:.1e},\n"
        f" nml_start = {p['start']:.2f},\n"
        f" nml_theta = {p['theta']:.2f},\n"
        " nml_budget = 200.0d0,\n"
        " nml_tfrac = 1.0d-4,\n"
        " nml_tol = 1.0d-6,\n"
        "/\n"
    )


# split BEFORE tuning: ~1/3 tuning (11), ~2/3 held-out (21).
# Fixed round-robin by index so every regime hits both sets and the
# choice cannot be re-picked after seeing results.
TUNING = {"i02", "i05", "i09", "i13", "i17", "i21",
          "i25", "i29", "i07", "i15", "i23"}


def main():
    os.makedirs(IDIR, exist_ok=True)
    names = []
    for name, p in INSTANCES:
        fname = name + ".in"
        with open(os.path.join(IDIR, fname), "w") as fh:
            fh.write(fmt(p))
        names.append(name.split("_")[0])
    seeds = sorted({p["seed"] for _, p in INSTANCES})
    with open(os.path.join(IDIR, "split.txt"), "w") as fh:
        fh.write("# instance split: TUNE = hyperparameter tuning set;\n")
        fh.write("# HELDOUT = pre-registered evaluation set.\n")
        fh.write("# Chosen in gen_ensemble.py before any tuning run.\n")
        for name, _ in INSTANCES:
            tag = "TUNE" if name.split("_")[0] in TUNING else "HELDOUT"
            fh.write(f"{name} {tag}\n")
    print(f"wrote {len(INSTANCES)} instances, seeds {seeds}, "
          f"tune={len(TUNING)}, heldout={len(INSTANCES)-len(TUNING)}")


if __name__ == "__main__":
    main()
