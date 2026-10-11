#!/usr/bin/env python3
"""Unit tests for the NG_ABSTAIN switch in
new_main/subs/do_ng_accel_band_v4.f (NG_MIT_BAND_OPT_V1).

The banded secant fit (NGACCEL) produces a single (A,B) pair shared by
all NT*ND components in a band.  With NG_ABSTAIN=F (upstream) one bad
component vetoes the whole depth and >3 bad depths cancel the whole
acceleration.  With NG_ABSTAIN=T a component is held at its current
value when:
  * NEWPOP <= 0               (secant crossed zero -- never physical)
  * POPS   <= 0               (acceptance ratio undefined)
  * RDPOPS(K,L,1) == RDPOPS(K,L,2)  (held/static row, e.g. R,V,SIGMA+1
                                   frozen by the caller's non-mono edit)
A vacuous-acceleration guard reports NG_DONE=F when every component
was held, so no identical SCRTEMP record is appended.

The driver links do_ng_accel_band_v4.f + subs/ngaccel.f.  The ref
build uses adapt_ng_ref.f, which drops the extra NG_ABSTAIN dummy
from the call.
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
    TREE / "tools" / "set_kind_module.f",
    TREE / "unix" / "term_io.f",
    TREE / "subs" / "ngaccel.f",
    TREE / "subs" / "equal.f",
    TREE / "subs" / "x02ajf.f",
    TREE / "new_main" / "mod_subs" / "control_variable_mod.f",
    F90DIR / "stub_scr_read.f",
]

NT, ND, NBAND = 4, 8, 0
C = 0.6  # shared contraction factor


def build():
    help_objs = compile_objs("ng_help", HELPERS)
    drv_obj = compile_objs(
        "ng_drv", [F90DIR / "drv_ng.f90"], moddir=BINDIR / "obj_ng_help"
    )
    bins = {}
    for tag, root in (("new", TREE), ("ref", REF)):
        objs = compile_objs(
            f"ng_{tag}",
            [
                root / "new_main" / "subs" / "do_ng_accel_band_v4.f",
                F90DIR / f"adapt_ng_{tag}.f",
            ],
            moddir=BINDIR / "obj_ng_help",
        )
        bins[tag] = link(
            f"ng_{tag}", help_objs + drv_obj + objs, moddir=BINDIR / "obj_ng_help"
        )
    return bins


def geometric(xstar, c=C, e=None):
    """History converging to xstar: J=1 newest, J=4 oldest."""
    e = 0.2 * xstar if e is None else e
    return np.array([xstar + c ** (4 - j) * e for j in range(1, 5)])


def write_case(name, rdpops, abstain, nband=NBAND):
    path = WORKDIR / name
    with open(path, "w") as f:
        f.write(f"{NT} {ND} {nband}\n")
        f.write("T\n" if abstain else "F\n")
        for k in range(NT):
            for l in range(ND):
                f.write(
                    f"{k + 1} {l + 1} "
                    + " ".join(f"{x:.17e}" for x in rdpops[k, l])
                    + "\n"
                )
    return path


def parse(out):
    ng_done = "NG_DONE= T" in out or "NG_DONE=T" in out
    pops = {}
    vetoes = set()
    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "POP":
            pops[int(t[1]), int(t[2])] = float(t[3])
        elif "No NG acceleration performed at depth" in ln:
            vetoes.add(int(ln.split()[-1].rstrip(",")))
        elif "Temperature correction has wrong sign" in ln:
            pass
    return ng_done, pops, vetoes, out


def all_good_history():
    """Every component a clean geometric contraction to x*."""
    rd = np.zeros((NT, ND, 4))
    xs = np.array([1.0, 10.0, 0.5, 5000.0])  # component "values"
    for k in range(NT):
        for l in range(ND):
            rd[k, l] = geometric(xs[k] * (1.0 + 0.05 * l))
    return rd


def main():
    rep = Report("NG_MIT_BAND_OPT_V1 NG_ABSTAIN")
    bins = build()

    # -- T1: clean contraction, F vs T and vs reference ----------
    rd = all_good_history()
    results = {}
    for tag, ab in (("new", True), ("new", False), ("ref", False)):
        path = write_case(f"ng_good_{tag}_{int(ab)}.txt", rd, ab)
        results[tag, ab] = parse(run(bins[tag], path).stdout)
    ng_t = results["new", True]
    ng_f = results["new", False]
    ng_r = results["ref", False]
    rep.check("clean history: NG_DONE=T under F and T", ng_t[0] and ng_f[0])
    same = all(abs(ng_t[1][k] - ng_f[1][k]) < 1e-8 for k in ng_t[1])
    rep.check("clean history: T and F give same POPS", same)
    same_r = all(abs(ng_f[1][k] - ng_r[1][k]) < 1e-8 for k in ng_f[1])
    rep.check("clean history: new(F) == reference", same_r)
    accel = any(
        abs(ng_t[1][(k + 1, l + 1)] - rd[k, l, 0]) > 1e-10
        for k in range(NT)
        for l in range(ND)
    )
    rep.check("clean history: acceleration actually moved POPS", accel)

    # -- T2: zero-diff component at one depth --------------------
    # K=1,L=3: newest two records identical (a held/static row).
    rd = all_good_history()
    x0 = 1.15
    rd[0, 2] = [x0, x0, x0 + 0.5, x0 + 1.0]
    path = write_case("ng_zdiff_T.txt", rd, True)
    ng_done_t, pops_t, vet_t, _ = parse(run(bins["new"], path).stdout)
    path = write_case("ng_zdiff_F.txt", rd, False)
    ng_done_f, pops_f, vet_f, _ = parse(run(bins["new"], path).stdout)
    rep.check(
        "zero-diff row held under T (unchanged)",
        ng_done_t and abs(pops_t[1, 3] - x0) < 1e-14,
    )
    moved_f = abs(pops_f[1, 3] - x0) > 1e-10
    rep.check(
        "zero-diff row moves under F (secant carries it)",
        moved_f or pops_f[1, 3] != x0,
        f"F value={pops_f[1, 3]:.6g}",
    )
    # its sibling components at that depth still accelerate under T
    sib = abs(pops_t[2, 3] - rd[1, 2, 0]) > 1e-10
    rep.check("zero-diff row: siblings still accelerate under T", sib)

    # -- T3: negative extrapolation holds one component ----------
    rd = all_good_history()
    # K=2,L=5 history crashes under the shared (A,B): secant crosses
    # zero, so NEWPOP <= 0 -- never physical for a population.
    rd[1, 4] = [1.0, 0.5, 5.0, 4.0]
    path = write_case("ng_neg_T.txt", rd, True)
    ng_done_t, pops_t, vet_t, _ = parse(run(bins["new"], path).stdout)
    path = write_case("ng_neg_F.txt", rd, False)
    ng_done_f, pops_f, vet_f, _ = parse(run(bins["new"], path).stdout)
    rep.check(
        "neg-extrapolation comp held under T",
        ng_done_t and abs(pops_t[2, 5] - 1.0) < 1e-10,
        f"POPS={pops_t.get((2, 5))}",
    )
    rep.check("neg-extrapolation vetoes depth under F", 5 in vet_f)
    rep.check("depth 5 not vetoed under T (held, not rejected)", 5 not in vet_t)
    # other depths still accelerate under T
    sib = abs(pops_t[1, 1] - rd[0, 0, 0]) > 1e-10
    rep.check("neg comp: other depths still accelerate under T", sib)

    # -- T4: >3 bad depths cancels under F, survives under T -----
    # NBAND=4 confines the corrupted (A,B) to depths 1-4; depths 5-8
    # form a clean band that still accelerates under T.
    rd = all_good_history()
    for l in (0, 1, 2, 3):
        rd[1, l] = [1.0, 0.5, 5.0, 4.0]
    path = write_case("ng_many_T.txt", rd, True, nband=4)
    ng_done_t, pops_t, vet_t, _ = parse(run(bins["new"], path).stdout)
    path = write_case("ng_many_F.txt", rd, False, nband=4)
    ng_done_f, pops_f, vet_f, _ = parse(run(bins["new"], path).stdout)
    rep.check("4 bad comps (band 1): NG_DONE=T under T", ng_done_t and len(vet_t) == 0)
    rep.check(
        "4 bad depths: NG_DONE=F under F (cancelled)", not ng_done_f and len(vet_f) >= 4
    )
    # under F all POPS restored to newest estimate
    restored = all(
        abs(pops_f[k + 1, l + 1] - rd[k, l, 0]) < 1e-14
        for k in range(NT)
        for l in range(ND)
    )
    rep.check("F cancellation restores RDPOPS(:,:,1)", restored)
    # band 2 (depths 5-8) accelerated under T
    band2 = any(
        abs(pops_t[k + 1, l + 1] - rd[k, l, 0]) > 1e-10
        for k in range(NT)
        for l in range(4, ND)
    )
    rep.check("clean band still accelerates under T", band2)

    # -- T5: vacuous acceleration guard --------------------------
    rd = np.zeros((NT, ND, 4))
    for k in range(NT):
        for l in range(ND):
            v = (k + 1) * 10.0 + l
            rd[k, l] = [v, v, v, v]
    path = write_case("ng_vac_T.txt", rd, True)
    ng_done_t, _, _, out_t = parse(run(bins["new"], path).stdout)
    path = write_case("ng_vac_F.txt", rd, False)
    ng_done_f, _, _, out_f = parse(run(bins["new"], path).stdout)
    rep.check(
        "all-static history: NG_DONE=F under T",
        not ng_done_t and "all components held" in out_t,
    )
    rep.check("all-static history: NG_DONE=T under F (upstream)", ng_done_f)

    # -- Plot: T2/T3 scenario ------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 4))
    rd = all_good_history()
    ax.plot(np.arange(1, 5), rd[1, 0, :] / rd[1, 0, 0], "ko-", label="contracting comp")
    ax.plot(
        np.arange(1, 5),
        [x0, x0, x0 + 0.5, x0 + 1.0],
        "rs-",
        label="zero-diff row (held)",
    )
    ax.plot(
        np.arange(1, 5),
        [0.5, 0.05, 0.4, 0.35],
        "b^-",
        label="counter-trend (neg extrap)",
    )
    ax.set_xlabel("history index J (1=newest)")
    ax.set_ylabel("value")
    ax.set_title("NG_ABSTAIN component classes")
    ax.legend()
    fig.tight_layout()
    out = PLOTDIR / "ng_abstain.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
