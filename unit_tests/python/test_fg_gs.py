#!/usr/bin/env python3
"""FG_INEXT Gauss-Seidel checkpoint (plane/fg_j_cmf_nonmono_v1.f).

Alterations under test (DO_NONMONO_KNOP core infrastructure):
  - DUMP_FG_INEXT_GS writes FGINEXT_GS.tmp then atomically renames
    to FGINEXT_GS - a kill mid-write cannot destroy the old file
  - MAGIC/VER envelope plus full validation on restore: iteration
    tag, dims, FULL R/V/SIGMA fingerprint, COMP_DEPTH_IDX map,
    KN_FREQ grid sanity (positive, strictly decreasing), then
    value checks (finite, non-negative intensities and metrics)
  - GS_FOREIGN_SEED escape hatch: a tag-mismatched file may seed
    the field but is kept immature (I_NEXT_GS_OK stays F) so the
    warm-up re-converges before any BA step
  - Rejection paths leave the arrays zeroed and flag the lambda
    rebuild; FG_GS_TAG_EXP<=0 or a missing file short-circuit
    before the file is even opened
  - RESTORE_FG_INEXT_LATEST (in-run rollback): no tag check, same
    fingerprints, early rejection leaves the in-memory store
    untouched (stale field is no worse than pops-only rollback)

Driver allocates the non-monotonic module state directly
(I_P/I_M_DEPTH_NEXT on a 3x2x5 grid, KN_FREQ, COMP_DEPTH_IDX,
FG_RAY_SEEDED) so the real production write/read paths run end
to end - not the early-return path used when unallocated.

  rt        : dump + restore round trip, bitwise, metrics, seeded
  no_file   : missing checkpoint -> cold-start fallback
  no_tag    : FG_GS_TAG_EXP<=0 -> no file access at all
  stale     : file tag 7 vs expected 99 -> stale rejection
  foreign   : tag mismatch + GS_FOREIGN_SEED -> data read but
              I_NEXT_GS_OK stays F (immature seed)
  bad_rv    : perturbed R element -> grid fingerprint rejection
  bad_cmp   : perturbed COMP_DEPTH_IDX -> map mismatch rejection
  neg_i     : negative intensity in file -> data rejection, zeroed
  trunc     : truncated file -> read/format rejection
  latest    : LATEST rollback restores arrays, tag-independent
  latest_mm : LATEST SIGMA mismatch -> OK=F, store untouched
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_fg_gs"
    objs = [prod_obj(TREE / "plane" / "fg_j_cmf_nonmono_v1.f")]
    objs += compile_objs(
        "fg_gs_drv", [F90DIR / "drv_fg_gs.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("fg_gs", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"fg_gs_{tag}"
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
    rep = Report("FG_INEXT GS checkpoint alterations")
    binary = build()

    # rt: full write/restore round trip.
    d = casedir("rt")
    r = run(binary, "rt", cwd=d)
    rep.check("rt: checkpoint file exists after dump",
              val(r.stdout, "FILE_OK") is True)
    rep.check("rt: .tmp renamed away (atomic write)",
              val(r.stdout, "TMP_GONE") is True)
    rep.check("rt: restore accepted (I_NEXT_GS_OK)",
              val(r.stdout, "REST_OK") is True)
    rep.check("rt: I_P_DEPTH_NEXT bitwise round trip",
              val(r.stdout, "IPBIT") is True)
    rep.check("rt: I_M_DEPTH_NEXT bitwise round trip",
              val(r.stdout, "IMBIT") is True)
    rep.check("rt: FG_RAY_SEEDED committed on validate",
              val(r.stdout, "SEEDED") is True)
    rep.check("rt: drift metric returned",
              abs(fval(r.stdout, "CHG_FILE") - 2.0e-3) < 1e-12)
    rep.check("rt: cap count returned",
              fval(r.stdout, "CAP_FILE") == 3.0)
    rep.check("rt: OUTGEN Knop-store report written",
              "KNOP STORE" in r.stdout)
    rep.check("rt: acceptance line written",
              "RESTORED COMPLETE I_NEXT GENERATION" in r.stdout)

    # no_file / no_tag short-circuits.
    d = casedir("no_file")
    r = run(binary, "no_file", cwd=d)
    rep.check("no_file: missing-checkpoint report",
              "I_NEXT checkpoint present: NO" in r.stdout)
    rep.check("no_file: lambda-rebuild fallback",
              "FALLING BACK TO LAMBDA REBUILD" in r.stdout)
    rep.check("no_file: REST_OK=F", val(r.stdout, "REST_OK") is False)

    d = casedir("no_tag")
    r = run(binary, "no_tag", cwd=d)
    rep.check("no_tag: cold-start line",
              "cold start (NITSF<=0)" in r.stdout)
    rep.check("no_tag: file never opened",
              "checkpoint present" not in r.stdout)
    rep.check("no_tag: REST_OK=F", val(r.stdout, "REST_OK") is False)

    # stale: tag mismatch rejected when GS_FOREIGN_SEED=F.
    d = casedir("stale")
    r = run(binary, "stale", cwd=d)
    rep.check("stale: stale/foreign rejection",
              "valid: NO (stale/foreign checkpoint)" in r.stdout)
    rep.check("stale: REST_OK=F", val(r.stdout, "REST_OK") is False)

    # foreign: same mismatch accepted as immature seed.
    d = casedir("foreign")
    r = run(binary, "foreign", cwd=d)
    rep.check("foreign: foreign-seed notice",
              "accepting file as foreign seed" in r.stdout)
    rep.check("foreign: field restored bitwise",
              val(r.stdout, "IPBIT") is True)
    rep.check("foreign: immaturity marker written",
              "NOT treated as a validated generation" in r.stdout)
    rep.check("foreign: I_NEXT_GS_OK stays F",
              val(r.stdout, "REST_OK") is False)

    # bad_rv / bad_cmp: fingerprint enforcement.
    d = casedir("bad_rv")
    r = run(binary, "bad_rv", cwd=d)
    rep.check("bad_rv: grid fingerprint rejection",
              "valid: NO (grid fingerprint mismatch)" in r.stdout)
    rep.check("bad_rv: REST_OK=F", val(r.stdout, "REST_OK") is False)

    d = casedir("bad_cmp")
    r = run(binary, "bad_cmp", cwd=d)
    rep.check("bad_cmp: COMP_DEPTH map rejection",
              "valid: NO (COMP_DEPTH map mismatch)" in r.stdout)
    rep.check("bad_cmp: REST_OK=F", val(r.stdout, "REST_OK") is False)

    # neg_i: value validation + zeroed store.
    d = casedir("neg_i")
    r = run(binary, "neg_i", cwd=d)
    rep.check("neg_i: value validation rejects negative I",
              "valid: NO (non-finite/negative data)" in r.stdout)
    rep.check("neg_i: REST_OK=F", val(r.stdout, "REST_OK") is False)
    rep.check("neg_i: corrupt store zeroed",
              val(r.stdout, "IP_ZEROED") is True)

    # trunc: torn file -> read/format rejection.
    d = casedir("trunc")
    run(binary, "dump_only", cwd=d)
    f = d / "FGINEXT_GS"
    size = f.stat().st_size
    assert size > 200, f"checkpoint too small: {size}"
    with open(f, "r+b") as fh:
        fh.truncate(size // 2)
    r = run(binary, "restore_only", cwd=d)
    rep.check("trunc: read-failure rejection",
              "valid: NO (read/format failure)" in r.stdout)
    rep.check("trunc: REST_OK=F", val(r.stdout, "REST_OK") is False)

    # latest / latest_mm: in-run rollback semantics.
    d = casedir("latest")
    r = run(binary, "latest", cwd=d)
    rep.check("latest: rollback restore OK",
              val(r.stdout, "LATEST_OK") is True)
    rep.check("latest: I_P bitwise restored",
              val(r.stdout, "IPBIT") is True)
    rep.check("latest: I_M bitwise restored",
              val(r.stdout, "IMBIT") is True)

    d = casedir("latest_mm")
    r = run(binary, "latest_mm", cwd=d)
    rep.check("latest_mm: rollback rejected on fingerprint",
              val(r.stdout, "LATEST_OK") is False)
    rep.check("latest_mm: pre-call store left untouched",
              val(r.stdout, "IP_UNTOUCHED") is True)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
