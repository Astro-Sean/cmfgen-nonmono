#!/usr/bin/env python3
"""AUTO_ADD_ION stale-index fix (new_main/subs/auto_add_ion.f).

Alteration under test: the high-stage and no-stage SUB_GUESS_DC_V2 calls
now pass ATM(NID)%NXzV_F.  Upstream passed ATM(ID)%NXzV_F where ID held
the stale index left over from the _IN scan loops - the last existing
ion in the high-stage case, or one below the species range when no
stage existed at all.  The generated XzV_IN header carries the passed
NLEV verbatim, so a stale ID wrote the wrong level count into the file
that seeds the new ion's departure coefficients.

The driver declares one species (carbon) with stages CIII/CIV/CV and
NXzV_F = 7/11/13, then calls production AUTO_ADD_ION.  Real model files
(RVTJ, EDDFACTOR, EDDFACTOR_INFO, *_F_OSCDAT) are symlinked from the
m15_2p82d_VAL_stock run directory into per-case work dirs; the _IN files
that must pre-exist are fabricated in the format SUB_GUESS_DC_V2 reads.
A second binary compiles the frozen-reference auto_add_ion.f so prod and
upstream behaviour are compared on identical fixtures.

  low     : only CIV_IN present -> CIII_IN written through the unaltered
            low-stage call (prod == ref, NLEV=7); CV_IN also written
            through the altered high-stage call (prod 13 vs ref 11)
  high    : CIII_IN+CIV_IN present -> only CV_IN generated; prod writes
            ATM(3)%NXzV_F=13, ref writes ATM(ID=2)%NXzV_F=11
  nostage : no _IN present -> all three generated through the altered
            no-stage call with self reference; prod writes 7/11/13,
            ref indexes ATM(0) (loop ran to completion below the range)
  off     : AUTO_ADD_ION_STAGES=F -> immediate RETURN, no files, prod
            and ref identical
"""
import os
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, REF, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)

MODEL = (Path("/ptmp/sbrennan/SN2024pba/cmfgen_validation")
         / "m15_2p82d_VAL_stock" / "cmfgen" / "2.82days")
ND = 70  # RVTJ / EDDFACTOR depth count in the fixture model
IONS = ("CIII", "CIV", "CV")
NLEV = {"CIII": 7, "CIV": 11, "CV": 13}


def build(prod=True):
    tag = "auto_ion" if prod else "auto_ion_ref"
    md = BINDIR / f"obj_{tag}_mods"
    md.mkdir(parents=True, exist_ok=True)
    if prod:
        objs = [prod_obj(TREE / "new_main" / "subs" / "auto_add_ion.f")]
    else:
        objs = compile_objs(
            f"{tag}_src", [REF / "new_main" / "subs" / "auto_add_ion.f"],
            moddir=md, extra_flags=PROD_MOD_INC, prod=True)
    objs += compile_objs(
        "auto_ion_drv", [F90DIR / "drv_auto_ion.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a", "libsubs.a")]
    return link_prod(tag, objs + libs, moddir=md)


def write_in_file(path, nlev=3, gs_pop=1.0e-2):
    """Minimal <ION>_IN fixture in the layout SUB_GUESS_DC_V2 reads:
    a '!Format date' line, then 'T1 T2 NLEV_RD ND_RD', then per depth a
    'T1 GS_ION_POP' row plus NLEV_RD departure coefficients."""
    lines = [" 07-Jul-1997                             !Format date\n",
             f" 1.0E+00  1.0E+00  {nlev}  {ND}\n"]
    dc_row = " " + "  ".join(["1.0"] * nlev) + "\n"
    for _ in range(ND):
        lines.append(f" 1.0E-02  {gs_pop:14.6E}\n")
        lines.append(dc_row)
    path.write_text("".join(lines))


def casedir(tag, present, oscdat):
    d = WORKDIR / f"aai_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    for name in ("RVTJ", "EDDFACTOR", "EDDFACTOR_INFO"):
        os.symlink(MODEL / name, d / name)
    for ion in oscdat:
        os.symlink(MODEL / f"{ion}_F_OSCDAT", d / f"{ion}_F_OSCDAT")
    for ion in present:
        (d / f"{ion}_IN").write_text("")  # INQUIRE-only placeholder
    if "CIV" in present:
        write_in_file(d / "CIV_IN")
    return d


def in_header_nlev(path):
    """NLEV (3rd field) from the generated XzV_IN header row."""
    lines = path.read_text().splitlines()
    for i, ln in enumerate(lines):
        if "!Format date" in ln:
            for nxt in lines[i + 1:]:
                flds = nxt.split()
                if len(flds) >= 4:
                    return int(flds[2]), int(flds[3])
    return None


def main():
    rep = Report("AUTO_ADD_ION stale-index fix")
    if not MODEL.exists():
        rep.summary()
        raise SystemExit(f"fixture model dir missing: {MODEL}")
    prod = build(prod=True)
    ref = build(prod=False)

    # ---- low: unaltered low-stage call + altered high-stage call ----
    d = casedir("low", present=("CIV",), oscdat=IONS)
    out = run(prod, "on", cwd=d).stdout
    rep.check("low: prod completes", "DONE_AUTO_ION" in out)
    got = in_header_nlev(d / "CIII_IN")
    rep.check("low: prod CIII_IN NLEV=7 (unaltered call)",
              got == (NLEV["CIII"], ND), f"hdr={got}")
    got = in_header_nlev(d / "CV_IN")
    rep.check("low: prod CV_IN NLEV=13 = ATM(NID)", got == (NLEV["CV"], ND),
              f"hdr={got}")

    d = casedir("low_ref", present=("CIV",), oscdat=IONS)
    out = run(ref, "on", cwd=d).stdout
    rep.check("low: ref completes", "DONE_AUTO_ION" in out)
    got = in_header_nlev(d / "CIII_IN")
    rep.check("low: ref CIII_IN NLEV=7 (call path unaltered upstream)",
              got == (NLEV["CIII"], ND), f"hdr={got}")
    got = in_header_nlev(d / "CV_IN")
    rep.check("low: ref CV_IN NLEV=11 (stale ID = last existing ion)",
              got == (NLEV["CIV"], ND), f"hdr={got}")

    # ---- high: isolated altered call ----
    d = casedir("high", present=("CIII", "CIV"), oscdat=IONS)
    run(prod, "on", cwd=d)
    rep.check("high: prod CIII_IN untouched (scan stop)",
              (d / "CIII_IN").read_text() == "")
    got = in_header_nlev(d / "CV_IN")
    rep.check("high: prod CV_IN NLEV=13", got == (NLEV["CV"], ND),
              f"hdr={got}")

    d = casedir("high_ref", present=("CIII", "CIV"), oscdat=IONS)
    run(ref, "on", cwd=d)
    got = in_header_nlev(d / "CV_IN")
    rep.check("high: ref CV_IN NLEV=11 (stale ID)", got == (NLEV["CIV"], ND),
              f"hdr={got}")

    # ---- nostage: altered no-stage call, all three generated ----
    d = casedir("nostage", present=(), oscdat=IONS)
    out = run(prod, "on", cwd=d).stdout
    rep.check("nostage: prod completes", "DONE_AUTO_ION" in out)
    ok = True
    for ion in IONS:
        got = in_header_nlev(d / f"{ion}_IN")
        if got != (NLEV[ion], ND):
            ok = False
            break
    rep.check("nostage: prod writes ATM(NID) NLEV for every stage",
              ok, f"{ion}_IN hdr={got}")

    d = casedir("nostage_ref", present=(), oscdat=IONS)
    r = run(ref, "on", cwd=d, check=False)
    ref_nlev = []
    for ion in IONS:
        p = d / f"{ion}_IN"
        got = in_header_nlev(p) if p.exists() and p.stat().st_size else None
        ref_nlev.append(got[0] if got else None)
    rep.check("nostage: ref diverges (ATM(0) stale index, or dies)",
              ref_nlev != [NLEV[i] for i in IONS],
              f"ref NLEVs={ref_nlev} rc={r.returncode}")

    # ---- off: flag-gated return identical ----
    d = casedir("off", present=(), oscdat=())
    out = run(prod, "off", cwd=d).stdout
    rep.check("off: returns immediately", "DONE_AUTO_ION" in out)
    rep.check("off: no _IN generated",
              not list(d.glob("*_IN")) and "Generating" not in out)

    rep.summary()


if __name__ == "__main__":
    main()
