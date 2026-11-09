#!/usr/bin/env python3
"""MOM_JNEXT Gauss-Seidel checkpoint (plane/mom_j_cmf_nonmono_v1.f).

Alterations under test (DO_NONMONO_KNOP core infrastructure):
  - DUMP_MOM_JNEXT_GS writes MOMJNXT_GS.tmp then atomically renames
    to MOMJNXT_GS - a kill mid-write cannot destroy the old file
  - MAGIC/VER envelope plus full validation on restore: iteration
    tag, dims, R/V endpoint fingerprint, complete expanded-grid
    fingerprint (R_NM + R_PNT_NM), then value checks (finite,
    non-negative J; finite H; non-negative drift/event metrics)
  - Rejection paths leave J/H/QUAR zeroed and flag the lambda
    rebuild; MOM_GS_TAG_EXP<=0 or a missing file short-circuit
    before the file is even opened
  - RESTORE_MOM_JNEXT_LATEST (in-run rollback): no tag check, same
    fingerprints, and an early rejection leaves the in-memory
    store untouched (it may still be better than nothing)

Driver allocates the non-monotonic module state directly
(J_NEXT_FRQ/H_NEXT_FRQ on a 6x4 grid, R_NM/R_PNT_NM expanded-grid
fingerprint, MOM_QUAR_AGE census) so the real production write/
read paths run end to end - not the no-op stub path used when the
arrays are unallocated.

  rt        : dump + restore round trip, bitwise J/H/QUAR, metrics
  no_file   : missing checkpoint -> cold-start fallback
  no_tag    : MOM_GS_TAG_EXP<=0 -> no file access at all
  stale     : file tag 7 vs expected 99 -> stale rejection
  bad_rv    : perturbed R(1) endpoint -> dims rejection
  grid_mm   : perturbed R_NM fingerprint -> dims rejection
  neg_j     : negative J in the file -> data rejection, J zeroed
  trunc     : truncated file -> read rejection
  latest    : LATEST rollback restores J/H, tag-independent
  latest_mm : LATEST fingerprint mismatch -> OK=F, store untouched
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_mom_gs"
    objs = [prod_obj(TREE / "plane" / "mom_j_cmf_nonmono_v1.f")]
    objs += compile_objs(
        "mom_gs_drv", [F90DIR / "drv_mom_gs.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("mom_gs", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"mom_gs_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def val(out, key):
    m = re.search(rf"{key}=([TF])", out)
    return m.group(1) == "T" if m else None


def fval(out, key):
    m = re.search(rf"{key}=\s*([-+0-9.Ee]+)", out)
    return float(m.group(1)) if m else None


def main():
    rep = Report("MOM_JNEXT GS checkpoint alterations")
    binary = build()

    # rt: full write/restore round trip.
    d = casedir("rt")
    r = run(binary, "rt", cwd=d)
    rep.check("rt: checkpoint file exists after dump",
              val(r.stdout, "FILE_OK") is True)
    rep.check("rt: .tmp renamed away (atomic write)",
              val(r.stdout, "TMP_GONE") is True)
    rep.check("rt: restore accepted (MOM_J_GS_OK)",
              val(r.stdout, "REST_OK") is True)
    rep.check("rt: J_NEXT_VALID set",
              val(r.stdout, "JVALID") is True)
    rep.check("rt: J_NEXT_FRQ bitwise round trip",
              val(r.stdout, "JBIT") is True)
    rep.check("rt: H_NEXT_FRQ bitwise round trip",
              val(r.stdout, "HBIT") is True)
    rep.check("rt: MOM_QUAR_AGE bitwise round trip",
              val(r.stdout, "QBIT") is True)
    rep.check("rt: drift metric returned",
              abs(fval(r.stdout, "CHG_FILE") - 1.0e-3) < 1e-12)
    rep.check("rt: event count returned",
              fval(r.stdout, "CAP_FILE") == 2.0)
    rep.check("rt: OUTGEN sweep report written",
              "MOM_J_CMF_NONMONO_V1: GS sweep" in r.stdout)
    rep.check("rt: acceptance sentence written",
              "present, accepted" in r.stdout)

    # no_file: nothing on disk -> fallback, no crash.
    d = casedir("no_file")
    r = run(binary, "no_file", cwd=d)
    rep.check("no_file: cold-start sentence",
              "absent - cold start" in r.stdout)
    rep.check("no_file: lambda-rebuild fallback",
              "FALLING BACK TO LAMBDA REBUILD" in r.stdout)
    rep.check("no_file: REST_OK=F", val(r.stdout, "REST_OK") is False)

    # no_tag: expected tag <=0 short-circuits before any open.
    d = casedir("no_tag")
    r = run(binary, "no_tag", cwd=d)
    rep.check("no_tag: cold-start sentence",
              "absent - cold start" in r.stdout)
    rep.check("no_tag: REST_OK=F", val(r.stdout, "REST_OK") is False)
    rep.check("no_tag: file never opened (no 'file present' line)",
              "file present" not in r.stdout)

    # stale: tag mismatch -> stale rejection.
    d = casedir("stale")
    r = run(binary, "stale", cwd=d)
    rep.check("stale: stale-tag rejection",
              "rejected because it is stale" in r.stdout)
    rep.check("stale: tag values reported",
              "tag 7 versus expected 99" in r.stdout)
    rep.check("stale: REST_OK=F", val(r.stdout, "REST_OK") is False)

    # bad_rv: perturbed endpoint -> dims rejection.
    d = casedir("bad_rv")
    r = run(binary, "bad_rv", cwd=d)
    rep.check("bad_rv: dims/grid rejection",
              "dimensions or the expanded grid do not match"
              in r.stdout)
    rep.check("bad_rv: REST_OK=F", val(r.stdout, "REST_OK") is False)

    # grid_mm: perturbed R_NM fingerprint -> dims rejection even
    # though endpoints/dims/tag all match.
    d = casedir("grid_mm")
    r = run(binary, "grid_mm", cwd=d)
    rep.check("grid_mm: expanded-grid fingerprint enforced",
              "dimensions or the expanded grid do not match"
              in r.stdout)
    rep.check("grid_mm: REST_OK=F", val(r.stdout, "REST_OK") is False)

    # neg_j: corrupted field -> data rejection + zeroed store.
    d = casedir("neg_j")
    r = run(binary, "neg_j", cwd=d)
    rep.check("neg_j: value validation rejects negative J",
              "data are non-finite or negative" in r.stdout)
    rep.check("neg_j: REST_OK=F", val(r.stdout, "REST_OK") is False)
    rep.check("neg_j: corrupt store zeroed",
              val(r.stdout, "J_ZEROED") is True)

    # trunc: file cut mid-record -> read rejection.
    d = casedir("trunc")
    run(binary, "dump_only", cwd=d)
    f = d / "MOMJNXT_GS"
    size = f.stat().st_size
    assert size > 200, f"checkpoint too small: {size}"
    with open(f, "r+b") as fh:
        fh.truncate(size // 2)
    r = run(binary, "restore_only", cwd=d)
    rep.check("trunc: read-failure rejection",
              "file could not be read or parsed" in r.stdout)
    rep.check("trunc: REST_OK=F", val(r.stdout, "REST_OK") is False)

    # latest: rollback restore is tag-independent.
    d = casedir("latest")
    r = run(binary, "latest", cwd=d)
    rep.check("latest: rollback restore OK",
              val(r.stdout, "LATEST_OK") is True)
    rep.check("latest: J bitwise restored",
              val(r.stdout, "JBIT") is True)
    rep.check("latest: H bitwise restored",
              val(r.stdout, "HBIT") is True)

    # latest_mm: early fingerprint failure must NOT touch the
    # in-memory store (rollback candidate is still the better seed).
    d = casedir("latest_mm")
    r = run(binary, "latest_mm", cwd=d)
    rep.check("latest_mm: rollback rejected",
              val(r.stdout, "LATEST_OK") is False)
    rep.check("latest_mm: pre-call store left untouched",
              val(r.stdout, "J_UNTOUCHED") is True)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
