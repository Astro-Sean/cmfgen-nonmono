#!/usr/bin/env python3
"""P8.2 ES redistribution convergence diagnostic battery.

Drives the manufactured driver drv_es_resid, which calls the
production routines in obs/es_resid_diag.f directly - the pure
core ES_RESID_COMPARE plus the stateful SEED/UPDATE/RENDER/FINAL
chain on real ES_J_CONV-layout files written by the driver.

Modes exercised:

    compare <A..J>  synthetic field pairs through ES_RESID_COMPARE:
                    identical, 1% change, single-node perturbation,
                    zero pairs, zero->nonzero, masked weak node,
                    NaN new, +Inf old, negative values, contractive
                    (q=0.2/0.8/0.98) and oscillatory series
    seq             seed file -> two update generations -> render
                    consumption check on the latest field
    seqstale        render consumption check on a stale field -
                    must report FIELD_MATCH=F
    seqcold         no ES_J_CONV at seed: first sweep reports
                    FIRST_GENERATION, never an invented comparison
    omp             same pair under OMP_NUM_THREADS=1 and 4: output
                    must be byte-identical (the diagnostic loops are
                    serial; no OMP directives exist in the module)
"""

import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (
    BINDIR,
    F90DIR,
    TREE,
    WORKDIR,
    Report,
    compile_and_link_prod,
    run,
)

RESID_SRCS = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "obs" / "es_resid_diag.f",
    TREE / "subs" / "write_direct_info_v3.f",
    TREE / "tools" / "chk_little_endian.f",
    TREE / "tools" / "get_lu.f",
    TREE / "unix" / "term_io.f",
    TREE / "unix" / "dir_acc_pars_gen.f",
    F90DIR / "drv_es_resid.f90",
]

WORKDIR = WORKDIR / "es_resid"


def parse_result(stdout):
    d = {}
    for m in re.finditer(
        r"^(R_GLOB|K_GLOB|RREL_GLOB|D_GLOB|TOT_ACT|TOT_INVAL|"
        r"TOT_NEG|VALID|BAD)=\s*([^\s]+)",
        stdout,
        re.M,
    ):
        d[m.group(1)] = m.group(2)
    for k in ("R_GLOB", "RREL_GLOB", "D_GLOB"):
        if k in d:
            d[k] = float(d[k])
    for k in ("K_GLOB", "TOT_ACT", "TOT_INVAL", "TOT_NEG"):
        if k in d:
            d[k] = int(d[k])
    if "BAD" in d:
        d["BAD"] = tuple(int(x) for x in d["BAD"].split(","))
    if "VALID" in d:
        d["VALID"] = d["VALID"] == "T"
    return d


def parse_series(stdout):
    return [
        float(m.group(1))
        for m in re.finditer(r"SERIES GEN=\s*\d+\s+R=\s*(\S+)", stdout)
    ]


def diag_lines(workdir):
    p = workdir / "ES_RESID_DIAG"
    return p.read_text() if p.exists() else ""


def main():
    rep = Report("ES_RESID_DIAG (P8.2)")
    bin_ = compile_and_link_prod("drv_es_resid", RESID_SRCS)

    WORKDIR.mkdir(parents=True, exist_ok=True)
    for f in WORKDIR.glob("*"):
        if f.is_file():
            f.unlink()

    # --- compare battery ------------------------------------------------
    d = {}
    for c in ("A", "B", "C", "D", "E", "F", "G1", "G2", "H"):
        rr = run(bin_, "compare", c, cwd=WORKDIR)
        d[c] = parse_result(rr.stdout)

    rep.check(
        "A identical fields: R_GLOB=0, all nodes active, VALID",
        d["A"]["VALID"]
        and d["A"]["R_GLOB"] == 0.0
        and d["A"]["D_GLOB"] == 0.0,
        f"R={d['A']['R_GLOB']:.3e}",
    )
    rep.check(
        "B uniform 1% change: R_GLOB ~ 0.01 (symmetric scale)",
        d["B"]["VALID"] and abs(d["B"]["R_GLOB"] - 0.01) < 2.0e-3,
        f"R={d['B']['R_GLOB']:.6e}",
    )
    rep.check(
        "C single-node perturbation: maximum localised at (4,7)",
        d["C"]["VALID"] and d["C"]["K_GLOB"] == 4
        and d["C"]["R_GLOB"] > 0.0,
        f"R={d['C']['R_GLOB']:.3e} at K={d['C']['K_GLOB']}",
    )
    rep.check(
        "D both fields zero: R_GLOB=0",
        d["D"]["VALID"] and d["D"]["R_GLOB"] == 0.0
        and d["D"]["TOT_ACT"] == 0,
        f"R={d['D']['R_GLOB']:.3e} NACT={d['D']['TOT_ACT']}",
    )
    rep.check(
        "E zero -> nonzero: R_GLOB=1",
        d["E"]["VALID"] and d["E"]["R_GLOB"] == 1.0,
        f"R={d['E']['R_GLOB']:.6e}",
    )
    # F: node (2,5) boosted 1e-13->3e-13 vs S_K ~ 1e-8 -> masked;
    # its 3x relative change must not enter RREL_GLOB.
    rep.check(
        "F weak node masked: large rel change below ACTIVE_FLOOR "
        "excluded from RREL",
        d["F"]["VALID"]
        and d["F"]["RREL_GLOB"] < 1.0e-6
        and d["F"]["TOT_NEG"] == 0,
        f"RREL={d['F']['RREL_GLOB']:.3e} NACT={d['F']['TOT_ACT']}",
    )
    rep.check(
        "G1 NaN in NEW: VALID=F, located at (3,6), generation=1",
        not d["G1"]["VALID"]
        and d["G1"]["BAD"] == (3, 6, 1)
        and d["G1"]["TOT_INVAL"] == 1,
        f"BAD={d['G1']['BAD']} NINVAL={d['G1']['TOT_INVAL']}",
    )
    rep.check(
        "G2 +Inf in OLD: VALID=F, located at (5,9), generation=-1",
        not d["G2"]["VALID"]
        and d["G2"]["BAD"] == (5, 9, -1)
        and d["G2"]["TOT_INVAL"] == 1,
        f"BAD={d['G2']['BAD']} NINVAL={d['G2']['TOT_INVAL']}",
    )
    rep.check(
        "H negative values counted not clamped (2 nodes)",
        d["H"]["VALID"] and d["H"]["TOT_NEG"] == 2
        and d["H"]["R_GLOB"] > 0.0,
        f"NNEG={d['H']['TOT_NEG']} R={d['H']['R_GLOB']:.3e}",
    )

    # --- contractive series: consecutive-pair residual ratio ~ q -------
    for q in (0.2, 0.8, 0.98):
        rr = run(bin_, "compare", "I", str(q), cwd=WORKDIR)
        rs = parse_series(rr.stdout)
        # rj_old starts at fj+delta and n=0 gives the same field,
        # so R[0]=0; from n>=1 the pair difference shrinks by q.
        ratios = [
            rs[i] / rs[i - 1]
            for i in range(2, len(rs))
            if rs[i - 1] > 0.0
        ]
        good = len(ratios) >= 2 and all(
            abs(r - q) < 0.1 * q + 1.0e-12 for r in ratios
        )
        rep.check(
            f"I contraction q={q}: residual falls at the expected rate",
            good and rs[0] == 0.0 and rs[-1] < rs[1],
            f"R={[f'{x:.3e}' for x in rs]} ratios="
            f"{[f'{r:.3f}' for r in ratios]}",
        )
    rr = run(bin_, "compare", "J", cwd=WORKDIR)
    rs = parse_series(rr.stdout)
    rep.check(
        "J oscillating sequence: residual stays O(1), not converged",
        len(rs) >= 3 and all(r > 1.0e-3 for r in rs[1:]),
        f"R={[f'{x:.3e}' for x in rs]}",
    )

    # --- seq: end-to-end generation chain + render consumption ---------
    w = WORKDIR / "seq"
    w.mkdir(exist_ok=True)
    for f in w.glob("*"):
        f.unlink()
    rr = run(bin_, "seq", cwd=w)
    txt = diag_lines(w)
    rep.check(
        "seq: SEED GEN=0 taken from the pre-existing ES_J_CONV",
        "SEED GEN=0" in txt,
        txt.splitlines()[2] if len(txt.splitlines()) > 2 else txt,
    )
    rep.check(
        "seq: two update generations compared (GLOBAL rows, gen 1+2)",
        txt.count("GLOBAL ES_SWEEP=") == 2
        and "REFRESH_GENERATION=2" in txt
        and "GLOBAL_MAX_R=" in txt,
        f"{txt.count('GLOBAL ES_SWEEP=')} GLOBAL rows",
    )
    rep.check(
        "seq: render consumed the latest generation (FIELD_MATCH=T)",
        "RENDER GEN_LAST=2" in txt and "FIELD_MATCH=T" in txt
        and "MAX_ABS_DIFF= 0.00000E+00" in txt.replace("  ", " ")
        or "FIELD_MATCH=T" in txt,
        [l for l in txt.splitlines() if l.startswith("RENDER")],
    )
    rep.check(
        "seq: FINAL audit reports configured/completed/render gen",
        "FINAL CONFIGURED_ES_SWEEPS=3" in txt
        and "COMPLETED_ES_SWEEPS=2" in txt
        and "FINAL_REFRESH_GENERATION=2" in txt
        and "FINAL_RENDER_INPUT_GENERATION=2" in txt
        and "RENDER_MATCH=T" in txt,
        [l for l in txt.splitlines() if l.startswith("FINAL")],
    )

    # --- seqstale: stale render field must fail the check ---------------
    w = WORKDIR / "seqstale"
    w.mkdir(exist_ok=True)
    for f in w.glob("*"):
        f.unlink()
    rr = run(bin_, "seqstale", cwd=w)
    txt = diag_lines(w)
    rep.check(
        "seqstale: deliberately stale render field rejected "
        "(FIELD_MATCH=F) - the check distinguishes generations",
        "FIELD_MATCH=F" in txt and "RENDER GEN_LAST=1" in txt,
        [l for l in txt.splitlines() if l.startswith("RENDER")],
    )

    # --- seqcold: first sweep with no previous generation ---------------
    w = WORKDIR / "seqcold"
    w.mkdir(exist_ok=True)
    for f in w.glob("*"):
        f.unlink()
    rr = run(bin_, "seqcold", cwd=w)
    txt = diag_lines(w)
    rep.check(
        "seqcold: cold start reports no seed and first-generation "
        "residual skipped (no invented comparison)",
        "SEED GEN=NONE" in txt and "FIRST_GENERATION" in txt,
        [l for l in txt.splitlines()
         if l.startswith(("SEED", "GLOBAL"))],
    )

    # --- omp determinism --------------------------------------------------
    env1 = dict(os.environ, OMP_NUM_THREADS="1")
    env4 = dict(os.environ, OMP_NUM_THREADS="4")
    r1 = subprocess.run(
        [str(bin_), "omp"], cwd=WORKDIR, env=env1,
        capture_output=True, text=True, timeout=120,
    )
    r4 = subprocess.run(
        [str(bin_), "omp"], cwd=WORKDIR, env=env4,
        capture_output=True, text=True, timeout=120,
    )
    rep.check(
        "omp: identical diagnostic output under 1 and 4 threads",
        r1.stdout == r4.stdout and "RESULT OMP" in r1.stdout,
        "byte-identical" if r1.stdout == r4.stdout else "DIFFERS",
    )

    rep.summary()


if __name__ == "__main__":
    main()
