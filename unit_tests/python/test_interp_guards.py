#!/usr/bin/env python3
"""Unit tests for the strict-monotonicity guards added to the
interpolation routines.

Files under test:
  subs/lin_interp.f          LIN_INTERP(R,V,ND,U,W,NIN)
  subs/mon_interp.f          MON_INTERP
  subs/mon_interp_fast.f     MON_INTERP_FAST
  subs/mon_interp_fast_v2.f  MON_INTERP_FAST_V2
  subs/mon_int_funs_v2.f     MON_INT_FUNS_V2

Each routine is compiled twice: against the jfnk-tree source and
against ../cmfgen_reference/cur_cmf.  Valid inputs must give identical
results; invalid (non-strict / misdirected / out-of-range / too-small)
grids must abort cleanly in the new code with no RESULT_OK sentinel.

Note the deliberate asymmetry: MON_INTERP permits repeated output
abscissae (check is '<'), while the FAST variants reject them ('<=').
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
    Report,
    compile_objs,
    link,
    max_meaningful_diff,
    run,
)

HELPERS = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "unix" / "term_io.f",
]

ROUTINES = {
    "lin": "lin_interp.f",
    "mon": "mon_interp.f",
    "monf": "mon_interp_fast.f",
    "monf2": "mon_interp_fast_v2.f",
    "mif2": "mon_int_funs_v2.f",
}

# Expected behaviour of the NEW code per (routine, case):
#   True  -> must print RESULT_OK
#   False -> must abort (no RESULT_OK) and print an error
EXPECT = {
    "valid_dec": {"lin": True, "mon": True, "monf": True, "monf2": True, "mif2": True},
    "valid_inc": {"lin": True, "mon": True, "monf": True, "monf2": True, "mif2": True},
    "rnonmono": {
        "lin": False,
        "mon": False,
        "monf": False,
        "monf2": False,
        "mif2": False,
    },
    "rrepeat": {
        "lin": False,
        "mon": False,
        "monf": False,
        "monf2": False,
        "mif2": False,
    },
    "ureverse": {
        "lin": False,
        "mon": False,
        "monf": False,
        "monf2": False,
        "mif2": None,
    },
    "unonmono": {
        "lin": False,
        "mon": True,
        "monf": False,
        "monf2": False,
        "mif2": None,
    },
    "qzrrepeat": {"mon": True, "monf": False, "monf2": False},
    "urange": {"lin": False, "mon": False, "monf": False, "monf2": False},
    "ndsmall": {
        "lin": False,
        "mon": False,
        "monf": False,
        "monf2": False,
        "mif2": False,
    },
    "ninsmall": {"lin": False, "mon": False, "monf": False, "monf2": False},
}

ND, NIN, NX = 20, 30, 30


def build():
    help_objs = compile_objs("int_help", HELPERS)
    drv_obj = compile_objs(
        "int_drv", [F90DIR / "drv_interp.f90"], moddir=BINDIR / "obj_int_help"
    )
    bins = {}
    for tag, root in (("new", TREE), ("ref", REF)):
        objs = []
        for key, fname in ROUTINES.items():
            objs += compile_objs(
                f"int_{tag}_{key}",
                [root / "subs" / fname],
                moddir=BINDIR / "obj_int_help",
            )
        bins[tag] = link(
            f"interp_{tag}", help_objs + drv_obj + objs, moddir=BINDIR / "obj_int_help"
        )
    return bins


def parse_ok(out):
    """Return (ok, rows). ok = RESULT_OK sentinel present.

    A data row is any line in which every token parses as a float.
    Drivers either print an integer index plus values, or bare
    float pairs (mif2)."""
    rows = []
    ok = "RESULT_OK" in out
    for ln in out.splitlines():
        t = ln.split()
        if len(t) < 2:
            continue
        try:
            vals = [float(x) for x in t]
        except ValueError:
            continue
        rows.append(vals[1:] if t[0].isdigit() else vals)
    return ok, np.array(rows) if rows else None


def r_grid(dec=True):
    r = np.geomspace(1.0e5, 1.0e3, ND)
    return r if dec else r[::-1]


def main():
    rep = Report("Interpolation input guards")
    bins = build()

    # -- Guard matrix ---------------------------------------------------
    for case, per_rtn in EXPECT.items():
        for rtn, expect in per_rtn.items():
            if expect is None:
                continue
            r = run(bins["new"], rtn, case, check=False)
            ok_sent, _ = parse_ok(r.stdout)
            has_err = "Error" in r.stdout or "error" in r.stdout
            if expect:
                ok = ok_sent
                detail = "ran clean"
            else:
                ok = (not ok_sent) and has_err
                detail = f"abort={'yes' if has_err else 'NO'}"
            rep.check(
                f"{case:10s} {rtn:5s} -> " f"{'OK' if expect else 'abort'}", ok, detail
            )

    # -- New vs reference on valid inputs ------------------------------
    for case in ("valid_dec", "valid_inc"):
        for rtn in ROUTINES:
            ra = run(bins["new"], rtn, case)
            rb = run(bins["ref"], rtn, case)
            oka, a = parse_ok(ra.stdout)
            okb, b = parse_ok(rb.stdout)
            same = (
                oka and okb and a.shape == b.shape and max_meaningful_diff(a, b) == 0.0
            )
            rep.check(
                f"{case} {rtn}: new == reference",
                same,
                (
                    f"maxdiff={np.abs(a - b).max():.2e}"
                    if (oka and okb)
                    else "missing output"
                ),
            )

    # -- LIN_INTERP accuracy: linear W -> exact V ----------------------
    _, a = parse_ok(run(bins["new"], "lin", "valid_dec").stdout)
    # columns: R, V  (driver prints i, r, v)
    exact = 3.0 * a[:, 0] + 0.5
    err = np.abs(a[:, 1] - exact).max()
    rep.check("lin valid_dec: linear source exact", err < 1e-10, f"maxerr={err:.2e}")

    # -- MON_INTERP accuracy: exp(-R/3e4) smooth monotone --------------
    _, a = parse_ok(run(bins["new"], "mon", "valid_dec").stdout)
    qzr, q1, q2 = a[:, 0], a[:, 1], a[:, 2]
    err1 = np.abs(q1 - np.exp(-qzr / 3.0e4)).max()
    err2 = np.abs(q2 - qzr**2).max()
    scale2 = (qzr**2).max()
    rep.check(
        "mon valid_dec: exp interp <1% , R^2 interp <1%",
        err1 < 1e-2 and err2 / scale2 < 1e-2,
        f"err_exp={err1:.2e} err_r2={err2 / scale2:.2e}",
    )

    # -- MIF2 cubic evaluation vs analytic ------------------------------
    _, a = parse_ok(run(bins["new"], "mif2", "valid_dec").stdout)
    err = np.abs(a[:, 1] - np.exp(-a[:, 0] / 3.0e4)).max()
    rep.check(
        "mif2 valid_dec: cubic matches exp <5e-3", err < 5e-3, f"maxerr={err:.2e}"
    )

    # -- Plot: monotone cubic through a non-monotone-ish profile -------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    _, a = parse_ok(run(bins["new"], "mif2", "valid_dec").stdout)
    r = r_grid()
    ax[0].plot(r, np.exp(-r / 3.0e4), "ko", ms=4, label="nodes")
    ax[0].plot(a[:, 0], a[:, 1], "r-", lw=1, label="MON_INT_FUNS_V2")
    ax[0].plot(a[:, 0], np.exp(-a[:, 0] / 3.0e4), "k--", lw=0.8, label="exact")
    ax[0].set_xscale("log")
    ax[0].set_xlabel("R")
    ax[0].set_title("mif2 cubic vs exact")
    ax[0].legend()
    _, a = parse_ok(run(bins["new"], "mon", "valid_dec").stdout)
    ax[1].semilogy(a[:, 0], a[:, 1], "r-", label="MON_INTERP")
    _, b = parse_ok(run(bins["ref"], "mon", "valid_dec").stdout)
    ax[1].semilogy(b[:, 0], b[:, 1], "b--", lw=0.8, label="reference")
    ax[1].set_xlabel("QZR")
    ax[1].set_title("monotonic interp new vs ref")
    ax[1].legend()
    fig.tight_layout()
    out = PLOTDIR / "interp_guards.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
