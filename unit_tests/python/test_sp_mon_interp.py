#!/usr/bin/env python3
"""Guard checks for the single-precision monotone interpolators
(SP_MON_INTERP in pgplt/do_vec_op.f and MON_INTERP_SP in
pgplt/mon_interp_sp.f, both altered 05-Oct-2026).

Both routines gained fail-closed validation of the input radius grid
(R strictly monotonic), the output grid ordering, grid lengths, and
the QZR-in-range check.  The valid path is upstream monotone-cubic
interpolation and must be unchanged: column J=1 is VARRAY=R, which
any reasonable interpolant reproduces exactly.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, Report,
                    compile_objs, link_prod, prod_objs, run)

HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "subs/phys_con.f",
]

ND, NX, LIN = 10, 7, 3


def build():
    md = BINDIR / "obj_spmon"
    objs = prod_objs(
        [TREE / h for h in HELPERS]
        + [TREE / "pgplt" / "do_vec_op.f",
           TREE / "pgplt" / "mon_interp_sp.f"])
    objs += compile_objs(
        "spmon_drv", [F90DIR / "drv_sp_mon.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    # do_vec_op.o drags in the pgplt curve/pen module state; the
    # production archives resolve those intra-library references the
    # same way the shipped pgplt consumers do.
    libs = [TREE / "lib" / a
            for a in ("libpgplt.a", "libsubs.a", "libmod.a")]
    return link_prod("spmon", objs + libs, moddir=md)


ROUTINES = [(1, "SP_MON_INTERP"), (2, "MON_INTERP_SP")]


def main():
    rep = Report("SP monotone-interp guards")
    binary = build()

    for wh, tag in ROUTINES:
        # Case 1: valid increasing grid completes; J=1 (V=R)
        # reproduces QZR.
        r = run(binary, "1", wh, check=False)
        blob = r.stdout + r.stderr
        rep.check(f"{tag}: valid grid completes", r.returncode == 0,
                  detail=blob.strip().splitlines()[-1]
                  if blob.strip() else "")
        vals = {}
        for ln in r.stdout.splitlines():
            t = ln.split()
            if len(t) == 4 and t[0] == "SPM":
                vals[(int(t[1]), int(t[2]))] = float(t[3])
        qzr = [1.2 + 0.6 * i for i in range(NX)]
        if vals:
            err = max(abs(vals[(i + 1, 1)] - qzr[i])
                      for i in range(NX))
            rep.check(f"{tag}: linear column reproduced",
                      err < 2e-5, detail=f"max abs err {err:.3e}")
            finite = all(abs(v) < 1e6 for v in vals.values())
            rep.check(f"{tag}: all outputs finite", finite)
        else:
            rep.check(f"{tag}: linear column reproduced", False,
                      detail="no SPM output")
            rep.check(f"{tag}: all outputs finite", False,
                      detail="no SPM output")

        # Case 6: valid decreasing grid completes too.
        r6 = run(binary, "6", wh, check=False)
        rep.check(f"{tag}: decreasing grid completes",
                  r6.returncode == 0,
                  detail=(r6.stdout + r6.stderr).strip().splitlines()[-1]
                  if (r6.stdout + r6.stderr).strip() else "")

        # Cases 2-5: each malformed grid must STOP inside the routine.
        for case, what in [(2, "repeated R endpoint"),
                           (3, "reversed R grid"),
                           (4, "unordered QZR"),
                           (5, "QZR outside R range")]:
            rr = run(binary, str(case), wh, check=False)
            blob = rr.stdout + rr.stderr
            stopped = rr.returncode != 0 or "not monotonic" in blob \
                or "repeated" in blob or "not ordered" in blob \
                or "outside range" in blob
            rep.check(f"{tag}: {what} rejected", stopped,
                      detail=f"rc={rr.returncode} "
                             f"{blob.strip().splitlines()[-1] if blob.strip() else ''}")

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
