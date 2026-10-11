#!/usr/bin/env python3
"""Unit tests for the NONMONO_V VINF consistency check in
newsubs/rd_rv_file_v2.f.

Old rule (kept when NONMONO_V=F): abort when |V(1)-VINF|/V(1) > 0.1.
New rule (NONMONO_V=T): VINF is compared against MAXVAL(|V|) so an
interior velocity peak or a reversing profile no longer trips the
check, while a mismatched VINF still aborts.

The jfnk-tree source and the cmfgen_reference source differ in
signature (extra NONMONO_V dummy); a thin adapter wraps each so one
driver links both.  RVSIG_COL files are written to a scratch dir per
case; the driver chdirs implicitly by running with cwd=<case dir>.
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
    TREE / "unix" / "gen_asci_open.f",
    TREE / "subs" / "set_case.f",
]

ND = 20
RMAX, RP = 1.0e4, 1.0e2  # file radii already in final units


def build():
    help_objs = compile_objs("rv_help", HELPERS)
    drv_obj = compile_objs(
        "rv_drv", [F90DIR / "drv_rd_rv.f90"], moddir=BINDIR / "obj_rv_help"
    )
    bins = {}
    for tag, root in (("new", TREE), ("ref", REF)):
        objs = compile_objs(
            f"rv_{tag}",
            [root / "newsubs" / "rd_rv_file_v2.f", F90DIR / f"adapt_rd_rv_{tag}.f"],
            moddir=BINDIR / "obj_rv_help",
        )
        bins[tag] = link(
            f"rdrv_{tag}", help_objs + drv_obj + objs, moddir=BINDIR / "obj_rv_help"
        )
    return bins


def write_rvsig(case_dir, r, v, sig):
    case_dir.mkdir(parents=True, exist_ok=True)
    with open(case_dir / "RVSIG_COL", "w") as f:
        f.write("Synthetic RVSIG_COL for unit tests\n")
        f.write(f"{ND}  !Number of depth points\n")
        f.write("! R[1e10 cm]   V[km/s]   SIGMA\n")
        for i in range(ND):
            f.write(f"{r[i]:16.8e} {v[i]:16.8e} {sig[i]:16.8e}\n")


def case_dir(name, r, v, sig=None):
    if sig is None:
        sig = np.full(len(r), -0.5)
    d = WORKDIR / f"rv_{name}"
    write_rvsig(d, np.asarray(r, float), np.asarray(v, float), np.asarray(sig, float))
    return d


def nonmono_v(peak=6000.0, vout=2000.0, vin=500.0, rev=False):
    """V profile: low at both boundaries, peak interior at i=7."""
    v = np.linspace(vout, vin, ND)
    v[5:10] = [3000.0, 4500.0, peak, 4500.0, 3000.0]
    if rev:
        v[10:15] = -3000.0
    return v


def attempt(bins, tag, d, nonmono, vinf):
    r = run(bins[tag], "T" if nonmono else "F", vinf, RMAX, RP, cwd=d, check=False)
    ok = "RESULT_OK" in r.stdout
    rows = {}
    for ln in r.stdout.splitlines():
        t = ln.split()
        if t and t[0].isdigit():
            rows[int(t[0]) - 1] = [float(x) for x in t[1:]]
    arr = np.array([rows[i] for i in sorted(rows)]) if rows else None
    return ok, arr, r.stdout


def main():
    rep = Report("RD_RV_FILE_V2 NONMONO_V / VINF check")
    bins = build()
    r = np.geomspace(RMAX, RP, ND)

    # -- T1: monotonic V, V(1)=VINF: passes in both modes ----------
    v = np.geomspace(8000.0, 500.0, ND)
    d = case_dir("mono", r, v)
    for nm in (False, True):
        ok, a, _ = attempt(bins, "new", d, nm, 8000.0)
        rep.check(f"mono V, V(1)=VINF, nonmono={nm}: pass", ok and a is not None)

    # -- T2: interior peak, VINF=peak -----------------------------
    v = nonmono_v(peak=6000.0)
    d = case_dir("peak", r, v)
    ok_new, a_new, _ = attempt(bins, "new", d, True, 6000.0)
    rep.check(
        "interior peak, VINF=peak, NONMONO_V=T: pass", ok_new and a_new is not None
    )
    ok_new_f, _, _ = attempt(bins, "new", d, False, 6000.0)
    rep.check("interior peak, VINF=peak, NONMONO_V=F: abort", not ok_new_f)
    ok_ref, _, _ = attempt(bins, "ref", d, True, 6000.0)
    rep.check("interior peak, reference build: abort (old rule)", not ok_ref)

    # -- T3: interior peak, VINF=V(1) (mismatch vs peak) ----------
    ok_new, _, out = attempt(bins, "new", d, True, 2000.0)
    rep.check(
        "interior peak, VINF=V(1), NONMONO_V=T: abort",
        not ok_new and "inconsistent" in out,
    )
    ok_ref, _, _ = attempt(bins, "ref", d, True, 2000.0)
    rep.check("interior peak, VINF=V(1), reference: pass", ok_ref)

    # -- T4: VINF within 10% of peak -> pass ----------------------
    ok, _, _ = attempt(bins, "new", d, True, 5800.0)
    rep.check("peak, VINF 3.3% below peak: pass", ok)
    ok, _, _ = attempt(bins, "new", d, True, 5400.0)
    rep.check("peak, VINF 10% below peak: pass", ok)
    ok, _, _ = attempt(bins, "new", d, True, 5300.0)
    rep.check("peak, VINF 11.7% below peak: abort", not ok)

    # -- T5: VINF<=0.1 skips check --------------------------------
    ok, _, _ = attempt(bins, "new", d, True, 0.05)
    rep.check("VINF=0.05 (<0.1): check skipped, pass", ok)

    # -- T6: all-zero velocity, VINF>0.1 -> abort -----------------
    d0 = case_dir("zerov", r, np.zeros(ND))
    ok, _, out = attempt(bins, "new", d0, True, 6000.0)
    rep.check(
        "V=0 everywhere, VINF=6000, NONMONO_V=T: abort",
        not ok and "inconsistent" in out,
    )

    # -- T7: reversing velocity (negative interior) ---------------
    v = nonmono_v(peak=6000.0, rev=True)
    d = case_dir("rev", r, v)
    ok, _, _ = attempt(bins, "new", d, True, 6000.0)
    rep.check("reversing V, |V|peak=6000=VINF: pass", ok)
    ok, _, _ = attempt(bins, "new", d, True, 2000.0)
    rep.check("reversing V, VINF=2000 != peak: abort", not ok)

    # -- T8: grid errors still caught (pre-existing guards) -------
    rb = r.copy()
    rb[7] = rb[5]  # non-decreasing
    d = case_dir("rbad", rb, np.geomspace(8000.0, 500.0, ND))
    for tag in ("new", "ref"):
        ok, _, out = attempt(bins, tag, d, False, 8000.0)
        rep.check(f"{tag}: non-decreasing R aborts", not ok)
    rb = np.geomspace(RMAX, 90.0, ND)  # wrong RMAX/RP ratio
    d = case_dir("rscale", rb, np.geomspace(8000.0, 500.0, ND))
    ok, _, out = attempt(bins, "new", d, False, 8000.0)
    rep.check("new: RMAX/RP mismatch aborts", not ok)

    # -- Plot: the three velocity profiles used -------------------
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(np.geomspace(8000.0, 500.0, ND), "ko-", ms=4, label="monotonic")
    ax.plot(nonmono_v(peak=6000.0), "rs-", ms=4, label="interior peak")
    ax.plot(nonmono_v(peak=6000.0, rev=True), "b^-", ms=4, label="reversing")
    ax.axhline(6000.0, ls="--", c="gray", lw=0.8)
    ax.text(0, 6200.0, "VINF=peak", fontsize=8)
    ax.set_xlabel("depth index (1=outer)")
    ax.set_ylabel("V [km/s]")
    ax.set_title("RD_RV_FILE_V2 test profiles")
    ax.legend()
    fig.tight_layout()
    out = PLOTDIR / "rd_rv_vinf_profiles.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
