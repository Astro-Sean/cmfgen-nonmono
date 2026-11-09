#!/usr/bin/env python3
"""INS_LINE_V6 / INS_LINE_OBS_V5 blanketed-line scans
(subs/ins_line_v6.f, subs/ins_line_obs_v5.f).

Shared alteration: the DO WHILE scans that locate the next BLANK
transition.  Upstream combined the bounds test and the element
test with .AND. -- Fortran does not short-circuit, so an all-SOB
deck read TRANS_TYPE(N_LINES+1) out of bounds.  The split form
ends LN_INDX at N_LINES+1 when no BLANK exists.

Cases:
  v6_sob  all-SOB deck -> early return FREQ=NU_CONT, NFREQ=NCF
  v6_bla  BLANK+SOB -> merged grid; line indices set only for the
          BLANK line; LINES_THIS_FREQ populated inside the zone
  v5_sob  INCL_ALL=F, all-SOB -> pure continuum spacing
  v5_bla  INCL_ALL=F, one BLANK -> more frequencies than v5_sob
  v5_all  INCL_ALL=T, all-SOB -> lines included anyway
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)

NU_CONT = [40, 30, 25, 20, 18, 16, 14, 12, 10, 8, 6, 4, 3, 2, 1]


def build():
    md = BINDIR / "obj_ins"
    objs = [
        prod_obj(TREE / "subs" / "ins_line_v6.f"),
        prod_obj(TREE / "subs" / "ins_line_obs_v5.f"),
    ]
    objs += compile_objs(
        "ins_drv", [F90DIR / "drv_ins_line.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("ins", objs + libs, moddir=md)


def freqs(out):
    return [float(m.group(1))
            for m in re.finditer(r"FR\s+\d+\s+(\S+)", out)]


def casedir(tag):
    d = WORKDIR / f"ins_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def run_case(binary, tag, d):
    return run(binary, tag, cwd=d, check=False)


def main():
    rep = Report("INS_LINE_V6 / INS_LINE_OBS_V5 OOB-scan fixes")
    binary = build()

    d = casedir("v6_sob")
    r = run_case(binary, "v6_sob", d)
    fr = freqs(r.stdout)
    nf = int(re.search(r"NF\s+(\d+)", r.stdout).group(1))
    rep.check("v6_sob: NFREQ=NCF (early return)",
              nf == len(NU_CONT), detail=f"NFREQ={nf}")
    rep.check("v6_sob: FREQ == NU_CONT verbatim",
              fr == [float(x) for x in NU_CONT])
    rep.check("v6_sob: returned without OOB crash",
              r.returncode == 0)

    d = casedir("v6_bla")
    r = run_case(binary, "v6_bla", d)
    fr = freqs(r.stdout)
    nf = int(re.search(r"NF\s+(\d+)", r.stdout).group(1))
    m1 = re.search(r"LIX 1\s+(\d+)\s+(\d+)", r.stdout)
    m2 = re.search(r"LIX 2\s+(\d+)\s+(\d+)", r.stdout)
    nz = int(re.search(r"NZ\s+(\d+)", r.stdout).group(1))
    rep.check("v6_bla: NFREQ > NCF", nf > len(NU_CONT),
              detail=f"NFREQ={nf}")
    rep.check("v6_bla: strictly decreasing",
              all(fr[i] > fr[i + 1] for i in range(len(fr) - 1)))
    rep.check("v6_bla: BLANK line gets indices",
              m1 and 0 < int(m1.group(1)) <= int(m1.group(2)) <= nf,
              detail=str(m1.groups()) if m1 else "?")
    rep.check("v6_bla: SOB line skipped",
              m2 and m2.groups() == ("0", "0"))
    rep.check("v6_bla: LINES_THIS_FREQ populated", nz > 0,
              detail=f"NZ={nz}")
    rep.check("v6_bla: grid lands near line center",
              min(abs(f - 17.0) / 17.0 for f in fr) < 1e-5,
              detail=str(min(abs(f - 17.0) for f in fr)))

    d = casedir("v5_sob")
    r = run_case(binary, "v5_sob", d)
    fr_sob = freqs(r.stdout)
    nf_sob = len(fr_sob)
    rep.check("v5_sob: FREQ(1)=NU_MAX", fr_sob and fr_sob[0] == 40.0)
    rep.check("v5_sob: strictly decreasing",
              all(fr_sob[i] > fr_sob[i + 1]
                  for i in range(len(fr_sob) - 1)))
    rep.check("v5_sob: reaches ~1.1*NU_MIN",
              fr_sob and fr_sob[-1] <= 1.1 * 1.0 + 1e-9,
              detail=f"last={fr_sob[-1] if fr_sob else '?'}")

    d = casedir("v5_bla")
    r = run_case(binary, "v5_bla", d)
    fr_bla = freqs(r.stdout)
    rep.check("v5_bla: strictly decreasing",
              all(fr_bla[i] > fr_bla[i + 1]
                  for i in range(len(fr_bla) - 1)))
    rep.check("v5_bla: more freqs than all-SOB",
              len(fr_bla) > nf_sob,
              detail=f"{len(fr_bla)} vs {nf_sob}")

    d = casedir("v5_all")
    r = run_case(binary, "v5_all", d)
    fr_all = freqs(r.stdout)
    nf0 = int(re.search(r"NF0\s*(\d+)", r.stdout).group(1))
    rep.check("v5_all: strictly decreasing",
              all(fr_all[i] > fr_all[i + 1]
                  for i in range(len(fr_all) - 1)))
    rep.check("v5_all: all lines included",
              len(fr_all) > nf0,
              detail=f"INCL_ALL=T {len(fr_all)} vs F {nf0}")

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
