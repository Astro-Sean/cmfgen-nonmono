#!/usr/bin/env python3
"""Unit tests for NG_SHADOW_DIAG in new_main/subs/do_ng_accel_band_v4.f.

With NG_SHADOW_DIAG=T, DO_NG_BAND_ACCEL_V4 reads up to 16 SCRTEMP
records at each NG event and reports, per band and window size N,
the affine-fit coefficient L1 norm, the implied prediction move,
the predicted-residual ratio, and the Magritte-style settle ratio
(Ceulemans et al. 2024, Sec 4.2).  Diagnostics only: the committed
acceleration must be identical with the flag off or on.

The serving stub (stub_scr_read_ng.f) provides SCR_READ_V2 backed
by a case file, so the multi-record read path is exercised
end-to-end.  The reported |c|1/move/settle numbers are checked
against an independent numpy recompute of the same fit.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

from common import (
    BINDIR,
    F90DIR,
    PROD_MOD_INC,
    TREE,
    WORKDIR,
    Report,
    compile_objs,
    is_prod_source,
    link_prod,
    prod_obj,
    prod_objs,
    run,
)

HELPERS = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "unix" / "term_io.f",
    TREE / "subs" / "ngaccel.f",
    TREE / "subs" / "equal.f",
    TREE / "subs" / "x02ajf.f",
    TREE / "new_main" / "mod_subs" / "control_variable_mod.f",
]

NT, ND = 4, 8


def build():
    md = BINDIR / "obj_ng_shd"
    objs = prod_objs([h for h in HELPERS if is_prod_source(h)])
    objs += prod_objs(
        [TREE / "new_main" / "subs" / "do_ng_accel_band_v4.f"]
    )
    objs += compile_objs(
        "ng_shd",
        [F90DIR / "stub_scr_read_ng.f", F90DIR / "drv_ng_shadow.f90"],
        moddir=md,
        extra_flags=PROD_MOD_INC,
    )
    return link_prod("ng_shadow", objs, moddir=md)


def make_history(nrec, seed=7):
    """Contracting history: record i (0-based) is iterate i+1,
    deviation ~ q**(i+1) so the newest record is closest."""
    rng = np.random.default_rng(seed)
    xstar = np.array(
        [[(k + 1) * 10**k * (1 + 0.03 * l) for l in range(ND)] for k in range(NT)]
    )
    recs = np.zeros((NT, ND, nrec))
    for i in range(nrec):
        it = i + 1
        base = xstar * (1 + 0.02 * 0.75**it + 0.005 * 0.45**it)
        recs[:, :, i] = base * (1 + rng.normal(0, 2e-4, (NT, ND)))
    return recs


def write_case(name, recs, lstng=-1000):
    nt, nd, nrec = recs.shape
    path = WORKDIR / name
    with open(path, "w") as f:
        f.write(f"{nt} {nd} {nrec} {lstng}\n")
        for l in range(nd):
            f.write(f"{100.0 + l} {2000.0 + 10 * l} {0.1 + 0.001 * l}\n")
        for i in range(nrec):
            for l in range(nd):
                f.write(" ".join(f"{recs[k, l, i]:.17e}" for k in range(nt)) + "\n")
    return path


_HDR_RE = re.compile(
    r"ng shadow: irec=\s*(-?\d+)\s+lst_ng=\s*(-?\d+)\s+nwin=\s*(\d+)"
)
_LINE_RE = re.compile(
    r"ng shadow: band\s+(\d+)-\s*(\d+)\s+N=\s*(\d+)"
    r"(?:\s+\|c\|1=\s*(\S+)\s+move=\s*(\S+)\s+predres=\s*(\S+)\s+settle=\s*(\S+))?"
)


def _f(s):
    return float(s.replace("D", "E").replace("d", "e"))


def parse_shadow(out):
    """Return {'header': (irec,lst_ng,nwin), 'lines': [(lst,lend,n,cl1,move,predres,settle)], 'pops': {}, 'ng_done': bool}"""
    pops, lines, header, ng_done = {}, [], None, False
    for ln in out.splitlines():
        if ln.startswith("NG_DONE"):
            ng_done = "T" in ln.split("=")[1]
        elif ln.startswith("POP"):
            t = ln.split()
            pops[int(t[1]), int(t[2])] = float(t[3])
        elif ln.startswith("ng shadow: irec="):
            m = _HDR_RE.search(ln)
            if m:
                header = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
        elif ln.startswith("ng shadow: band"):
            m = _LINE_RE.search(ln)
            if not m:
                continue
            lst, lend, n = int(m.group(1)), int(m.group(2)), int(m.group(3))
            if m.group(4) is None:
                lines.append((lst, lend, n, None, None, None, None))
            else:
                lines.append(
                    (
                        lst,
                        lend,
                        n,
                        _f(m.group(4)),
                        _f(m.group(5)),
                        _f(m.group(6)),
                        _f(m.group(7)),
                    )
                )
    return {"header": header, "lines": lines, "pops": pops, "ng_done": ng_done}


def ref_fit(recs, idxs, band_l0, band_l1, nrec):
    """Recompute the constrained least-squares fit in numpy.

    idxs are 0-based residual indices: residual i is
    r_(i+1) = x_(i+1) - x_(i+2) with Fortran x_i = file record
    (nrec-i+1).  Gram weighted by 1/x_1; coefficients sum to 1
    and act on the post-iterates x_(i+2) (Fortran RDX(K,L,I+1)).
    The shifted "previous" prediction uses idxs=[1..N-2]."""
    x = [recs[:, :, nrec - i] for i in range(1, max(idxs) + 3)]
    w = 1.0 / x[0]
    nr = len(idxs)
    g = np.zeros((nr, nr))
    sl = slice(band_l0 - 1, band_l1)
    for a, i in enumerate(idxs):
        for b, j in enumerate(idxs):
            g[a, b] = np.sum(w[:, sl] * (x[i] - x[i + 1])[:, sl] * (x[j] - x[j + 1])[:, sl])
    c = np.linalg.solve(g, np.ones(nr))
    c /= c.sum()
    p = sum(c[t] * x[idxs[t] + 1] for t in range(nr))
    return c, p


def main():
    rep = Report("NG_SHADOW_DIAG")
    binary = build()
    recs = make_history(12)
    case = write_case("ng_shd_case.txt", recs)

    # -- T1: flag off = no diagnostics, identical outcome ----------
    r_off = run(binary, case, "F")
    r_on = run(binary, case, "T")
    off = parse_shadow(r_off.stdout)
    on = parse_shadow(r_on.stdout)
    rep.check("flag F: no 'ng shadow' output", off["header"] is None and not off["lines"])
    rep.check("flag T: header present", on["header"] is not None)

    # -- T2: diagnostics do not change the acceleration ------------
    rep.check("NG_DONE=T both modes", off["ng_done"] and on["ng_done"])
    same = off["pops"] == on["pops"]
    rep.check("POP output bit-identical F vs T", same)

    # -- T3: line coverage -----------------------------------------
    # nwin = min(12,16) with lst_ng<0 -> N=4..12 = 9 lines per band.
    ns = sorted(n for _, _, n, *_ in on["lines"])
    rep.check("N values cover 4..12", ns == sorted([4, 5, 6, 7, 8, 9, 10, 11, 12] * 2))
    rep.check("18 band lines (2 bands x 9 N)", len(on["lines"]) == 18)
    numeric = [ln for ln in on["lines"] if ln[3] is not None]
    singular = [ln for ln in on["lines"] if ln[3] is None]
    rep.check(
        "no singular/den~0 windows in the noisy case",
        not singular,
        f"got {len(singular)}",
    )
    rep.check(
        "all numeric fields finite",
        all(np.isfinite(ln[4:7]).all() for ln in numeric),
    )

    # -- T4: N=4 values match numpy recompute (both bands) ----------
    nrec = recs.shape[2]
    for lst, lend in ((5, 8), (1, 4)):
        c4, p4 = ref_fit(recs, [0, 1, 2], lst, lend, nrec)
        # shifted-window prev fit: residuals r_2..r_3 on x_3..x_4
        c3, p3 = ref_fit(recs, [1, 2], lst, lend, nrec)
        x = [recs[:, :, nrec - i] for i in range(1, 6)]
        sl = slice(lst - 1, lend)
        move = np.max(np.abs(p4 - x[0])[:, sl] / np.abs(x[0])[:, sl])
        setn = np.max(np.abs(p4 - p3)[:, sl] / np.abs(p3)[:, sl])
        setd = np.max(np.abs(x[0] - x[1])[:, sl] / np.abs(x[1])[:, sl])
        expect = {
            "cl1": float(np.abs(c4).sum()),
            "move": float(move),
            "settle": float(setn / setd),
        }
        row = next(
            ln for ln in on["lines"] if ln[0] == lst and ln[1] == lend and ln[2] == 4
        )
        for name, got, want in (
            ("|c|1", row[3], expect["cl1"]),
            ("move", row[4], expect["move"]),
            ("settle", row[6], expect["settle"]),
        ):
            rep.check(
                f"band {lst}-{lend} N=4 {name} matches numpy",
                # ES10.2 prints 3 significant figures (~0.5%)
                abs(got - want) < 6e-3 * abs(want) + 1e-6,
                f"got={got:.4e} want={want:.4e}",
            )

    # -- T5: LST_NG window cap ---------------------------------------
    case_ng = write_case("ng_shd_lstng.txt", recs, lstng=8)
    out_ng = parse_shadow(run(binary, case_ng, "T").stdout)
    ns_ng = sorted(n for _, _, n, *_ in out_ng["lines"])
    # irec=12, lst_ng=8 -> nwin=min(16,12,12-8)=4 -> only N=4.
    rep.check(
        "lst_ng=8 caps window at N=4 (2 band lines)",
        ns_ng == [4, 4],
        f"got Ns={ns_ng}",
    )

    # -- T6: short history (nrec=5 -> N=4 only) ----------------------
    recs5 = make_history(5)
    case5 = write_case("ng_shd_5.txt", recs5)
    out5 = parse_shadow(run(binary, case5, "T").stdout)
    ns5 = sorted(n for _, _, n, *_ in out5["lines"])
    rep.check("nrec=5 gives N=4,5 per band", ns5 == [4, 4, 5, 5], f"got Ns={ns5}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
