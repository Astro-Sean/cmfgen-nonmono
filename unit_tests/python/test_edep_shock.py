#!/usr/bin/env python3
"""Unit tests for new_main/subs/get_edep_shock_power.f.

The change: under DO_NONMONO_SOB or DO_NONMONO_KNOP the routine now
requires V > 0 and strictly decreasing V -- VLOC_SHOCK_POWER is a
velocity-space location, which is ambiguous when a velocity occurs on
more than one radial branch.  With both flags F the legacy path is
unchanged (checked bit-for-bit against the reference build).

Also covered:
  * deposition profile is the prescribed Gaussian in velocity space;
  * normalization: 4*pi*sum(shock_power*R^2*fcl)*1e30 == prescribed
    shock power (parsed back from the routine's own diagnostic);
  * upstream guards still fire: zero clump_fac under SCL_PWR_BY_FCL,
    and sum(wrk)=0 when VLOC is far outside the grid.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (
    BINDIR,
    F90DIR,
    PLOTDIR,
    REF,
    TREE,
    WORKDIR,
    Report,
    compile_objs,
    link,
    run,
)

HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "new_main/mod_subs/mod_cmfgen.f",
    "new_main/mod_subs/control_variable_mod.f",
    "new_main/mod_subs/shock_power_mod.f",
    "subs/mod_trap_derivatives.f",
    "subs/nordtau.f",
    "subs/lum_from_eta_v2.f",
]

ND = 30
VLOC, DVLOC, PSP = 4000.0, 500.0, 1.0e40


def build():
    bins = {}
    for tag, root, drv in (
        ("new", TREE, "drv_edep.f90"),
        ("ref", REF, "drv_edep_ref.f90"),
    ):
        objs = compile_objs(
            f"edep_{tag}",
            [root / h for h in HELPERS]
            + [root / "new_main" / "subs" / "get_edep_shock_power.f", F90DIR / drv],
        )
        bins[tag] = link(f"edep_{tag}", objs, moddir=BINDIR / f"obj_edep_{tag}")
    return bins


def write_input(d, v_mut=None, fcl_mut=None):
    d.mkdir(parents=True, exist_ok=True)
    r = np.geomspace(1.0e5, 1.0e3, ND)
    v = np.linspace(9000.0, 500.0, ND)
    fcl = np.full(ND, 0.2)
    dens = 1.0e-12 * (r / r[0]) ** -3.0
    if v_mut is not None:
        v = v_mut(v.copy())
    if fcl_mut is not None:
        fcl = fcl_mut(fcl.copy())
    with open(d / "edep_input.dat", "w") as f:
        f.write(f"{ND}\n")
        for i in range(ND):
            f.write(f"{r[i]:.8e} {v[i]:.8e} {dens[i]:.8e} " f"{fcl[i]:.8e}\n")
    return d, r, v, fcl


def parse(out):
    rows = {}
    cmg_pow = None
    for ln in out.splitlines():
        t = ln.split()
        if t and t[0].isdigit() and len(t) >= 3:
            rows[int(t[0]) - 1] = [float(x) for x in t[1:]]
        if "prescribed_shock_power from CMFGEN" in ln:
            cmg_pow = float(ln.split(":")[-1])
    arr = np.array([rows[i] for i in sorted(rows)]) if rows else None
    return arr, cmg_pow


def run_case(exe, d, mode, scl="F"):
    return run(exe, mode, VLOC, DVLOC, PSP, scl, cwd=d, check=False, timeout=30)


def main():
    rep = Report("GET_EDEP_SHOCK_POWER")
    bins = build()

    d, r, v, fcl = write_input(WORKDIR / "edep_mono")

    # -- legacy path: flags F -------------------------------------
    outs = {}
    for tag in ("new", "ref"):
        p = run_case(bins[tag], d, "F")
        assert "RESULT_OK" in p.stdout, p.stdout[-1500:]
        outs[tag] = parse(p.stdout)
    rep.check(
        "flags F: new == ref shock_power",
        np.array_equal(outs["new"][0], outs["ref"][0]),
        f"maxdiff={np.abs(outs['new'][0] - outs['ref'][0]).max():.2e}",
    )

    sp = outs["new"][0][:, 1]
    vv = outs["new"][0][:, 0]
    # Gaussian deposition profile: peak at node nearest VLOC
    wgt = np.exp(-(((vv - VLOC) / DVLOC) ** 2))
    rep.check(
        "flags F: argmax at |V-VLOC| min", np.argmax(sp) == np.argmin(np.abs(vv - VLOC))
    )
    # shock_power ∝ wgt everywhere wgt > 0
    m = wgt > 1e-20
    ratio = sp[m] / wgt[m]
    rep.check(
        "flags F: profile proportional to Gaussian",
        np.abs(ratio - ratio.mean()).max() / ratio.mean() < 1e-10,
    )
    rep.check("flags F: all shock_power non-negative", (sp >= 0).all())
    rep.check(
        "flags F: CMFGEN power == prescribed",
        outs["new"][1] is not None and abs(outs["new"][1] - PSP) / PSP < 1e-8,
        f"got {outs['new'][1]}",
    )

    # -- legacy path with non-monotonic V: still completes ---------
    def nonmono(x):
        x[10:14] = x[10] + np.array([400.0, 900.0, 300.0, -200.0])
        return x

    d, r, vn, _ = write_input(WORKDIR / "edep_nonmono", v_mut=nonmono)
    for tag in ("new", "ref"):
        p = run_case(bins[tag], d, "F")
        rep.check(
            f"flags F + nonmono V: {tag} still completes", "RESULT_OK" in p.stdout
        )

    # -- nonmono flags: guards ------------------------------------
    for mode, name in (
        ("S", "DO_NONMONO_SOB"),
        ("K", "DO_NONMONO_KNOP"),
        ("B", "both flags"),
    ):
        p = run_case(bins["new"], d, mode)
        rep.check(
            f"{name}: non-decreasing V aborts",
            "RESULT_OK" not in p.stdout and "decreasing V" in p.stdout,
        )
        # valid monotone V still runs under the flags
        p = run_case(bins["new"], WORKDIR / "edep_mono", mode)
        rep.check(f"{name}: monotone V completes", "RESULT_OK" in p.stdout)

    # V <= 0 guard
    def vzero(x):
        x[5] = 0.0
        return x

    d, *_ = write_input(WORKDIR / "edep_vzero", v_mut=vzero)
    p = run_case(bins["new"], d, "S")
    rep.check(
        "DO_NONMONO_SOB: V=0 aborts",
        "RESULT_OK" not in p.stdout and "positive V" in p.stdout,
    )

    # repeated (equal) V guard
    def vequal(x):
        x[7] = x[6]
        return x

    d, *_ = write_input(WORKDIR / "edep_veq", v_mut=vequal)
    p = run_case(bins["new"], d, "K")
    rep.check(
        "DO_NONMONO_KNOP: equal-V pair aborts",
        "RESULT_OK" not in p.stdout and "decreasing V" in p.stdout,
    )

    # flags-F path unaffected by same inputs (legacy tolerance)
    p = run_case(bins["new"], d, "F")
    rep.check("flags F: equal-V pair still completes (legacy)", "RESULT_OK" in p.stdout)

    # -- upstream guards still live --------------------------------
    d, *_ = write_input(WORKDIR / "edep_fcl0", fcl_mut=lambda x: x * 0.0)
    for tag in ("new", "ref"):
        p = run_case(bins[tag], d, "F", scl="T")
        rep.check(
            f"SCL_PWR_BY_FCL + zero fcl: {tag} aborts",
            "RESULT_OK" not in p.stdout and "clump_fac is zero" in p.stdout,
        )

    # VLOC far outside grid -> wgt all ~0 -> sum(wrk)=0 abort
    d, *_ = write_input(WORKDIR / "edep_far")
    for tag in ("new", "ref"):
        p = run(bins[tag], "F", 3.0e5, DVLOC, PSP, "F", cwd=d, check=False, timeout=30)
        rep.check(
            f"VLOC off grid: {tag} aborts on sum(wrk)=0",
            "RESULT_OK" not in p.stdout and "sum(wrk) is zero" in p.stdout,
        )

    # -- plot ------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(vv, wgt / wgt.max(), "k--", lw=0.9, label="Gaussian")
    ax.plot(vv, sp / sp.max(), "ro", ms=4, label="shock_power")
    ax.set_xlabel("V (km/s)")
    ax.set_ylabel("normalized deposition")
    ax.legend()
    ax.set_title("Gaussian shock-power deposition profile")
    fig.tight_layout()
    out = PLOTDIR / "edep_shock.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
