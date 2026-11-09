#!/usr/bin/env python3
"""RD_CONT_J torn-record and interpolation coverage
(new_main/mod_subs/rd_cont_j.f).

Alteration under test: every direct-access EDDFACTOR record read
now carries IOSTAT with a zero-J fallback (upstream, a torn record
hard-crashed the LAMBDA warm-up). The RJ_HIGH read failure zeros
both bounds; the RJ_LOW read failure zeros only RJ_LOW; a mid-walk
failure clamps LOW_FREQ=HIGH_FREQ and keeps the last good field.

The driver builds a real unformatted direct-access file with the
production record layout ((RJ(1:ND),FREQ) per record), then calls
the routine with FIRST_FREQ=T/F sequences so the SAVE'd bound
arrays and the record walk both execute.

  torn_hi   : rec1 missing -> both bounds zeroed, RJ=0, ZETA-only
              source
  torn_lo   : rec2 missing -> RJ_LOW/LOW_FREQ zeroed, RJ=0
  first_ok  : FL between file freqs -> linear interp, ACCESS_F
              advanced past both records, P_OBS snapshot taken
  hi_clamp  : FL >= HIGH_FREQ -> adopt RJ_HIGH verbatim
  walk      : second call FL<LOW_FREQ -> record walk, interp in
              the new interval, ACCESS_F advanced
  walk_torn : torn record mid-walk -> LOW_FREQ=HIGH_FREQ fallback,
              RJ = last-good low field
  accurate  : ACCURATE=T -> STOP 'not yet implemented'
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_rd_cont_j"
    objs = [prod_obj(TREE / "new_main" / "mod_subs" / "rd_cont_j.f")]
    objs += compile_objs(
        "rcj_drv", [F90DIR / "drv_rd_cont_j.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in
            ("libsubs.a", "libmod.a", "libsubs.a")]
    return link_prod("rd_cont_j", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"rcj_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def val(out, key):
    m = re.search(rf"{key}=([TF])", out)
    return m.group(1) == "T" if m else None


def rj_vals(out):
    m = re.search(r"RJ=\s*([-+0-9.Ee\s]+)", out)
    if not m:
        return None
    return [float(x) for x in m.group(1).split()]


def run_case(binary, tag):
    d = casedir(tag)
    return run(binary, tag, cwd=d)


def side_checks(rep, tag, out):
    rep.check(f"{tag}: K_MOM = RJ/3", val(out, "KMOM_OK") is True)
    rep.check(f"{tag}: SOURCE = ZETA + THETA*RJ",
              val(out, "SRC_OK") is True)
    rep.check(f"{tag}: RSQHNU sentinel written",
              val(out, "RSQ_OK") is True)
    rep.check(f"{tag}: HBC_CMF set", val(out, "HBC_OK") is True)
    rep.check(f"{tag}: IPLUS seeded", val(out, "IPLUS_OK") is True)


def main():
    rep = Report("RD_CONT_J torn-record guards")
    binary = build()

    r = run_case(binary, "torn_hi")
    rep.check("torn_hi: read error reported",
              "Error reading EDDFACTOR in RD_CONT_J" in r.stdout)
    rep.check("torn_hi: run completes on fallback",
              "DONE_TORN_HI" in r.stdout)
    rep.check("torn_hi: RJ zeroed",
              rj_vals(r.stdout) == [0.0] * 4)

    r = run_case(binary, "torn_lo")
    rep.check("torn_lo: read error reported",
              "Error reading EDDFACTOR in RD_CONT_J" in r.stdout)
    rep.check("torn_lo: interp blends zeroed low bound (4.0)",
              rj_vals(r.stdout) == [4.0] * 4)

    r = run_case(binary, "first_ok")
    rep.check("first_ok: completes", "DONE_FIRST_OK" in r.stdout)
    rep.check("first_ok: linear interp FL=7 in [5,10]",
              rj_vals(r.stdout) == [5.6] * 4)
    rep.check("first_ok: ACCESS_F advanced past both records",
              "ACC=  3" in r.stdout or "ACC=3" in r.stdout)
    rep.check("first_ok: P_OBS snapshot on FIRST_FREQ",
              "P_OBS=" in r.stdout)
    side_checks(rep, "first_ok", r.stdout)

    r = run_case(binary, "hi_clamp")
    rep.check("hi_clamp: FL>=HIGH_FREQ adopts RJ_HIGH",
              rj_vals(r.stdout) == [8.0] * 4)

    r = run_case(binary, "walk")
    rep.check("walk: second call interpolates new interval",
              rj_vals(r.stdout) is not None
              and all(abs(x - 10.0/3.0) < 1e-5 for x in rj_vals(r.stdout)))
    rep.check("walk: ACCESS_F advanced through walk",
              "ACC=  4" in r.stdout or "ACC=4" in r.stdout)

    r = run_case(binary, "walk_torn")
    rep.check("walk_torn: mid-walk error reported",
              "Error reading EDDFACTOR in RD_CONT_J" in r.stdout)
    rep.check("walk_torn: LOW_FREQ=HIGH_FREQ fallback holds "
              "last good field",
              rj_vals(r.stdout) == [4.0] * 4)

    r = run_case(binary, "accurate")
    rep.check("accurate: unimplemented option STOPs",
              "Accurate option not yet implemented" in r.stdout
              and "ACC_NO_STOP" not in r.stdout)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
