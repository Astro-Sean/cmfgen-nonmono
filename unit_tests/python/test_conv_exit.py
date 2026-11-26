#!/usr/bin/env python3
"""test_conv_exit -- ten-gate ITER_SUM converged-exit replica.

Round-2 phase-2 patch A: the ITER_SUM replica was missing the ninth
conjunct (.NOT. LAM_FIELD_NOT_READY) that the real CONV_EXIT_MET in
CMFGEN_SUB carries, so an immature lagged field could print [MET]
while the true exit still held. The replica now calls
CONV_EXIT_EVAL (LAM_HANDOFF_MOD), the same evaluator both sides
share.

Tenth gate (25-Nov-2026): FCONV_CAP_IT_CNT -- the per-iteration
census of f-iteration cap exits, saved before the reset in
CMFGEN_SUB -- joins the gate. A frequency that left the
f/electron-scattering iteration at the pass cap carried an
unconverged ES source into the rates the SE solve used, so a
capped iteration is not a converged one. The suppressed
iteration-1 warning subset stays inside the census: the gate reads
the physical event, not the print flag.

Compiles the production lam_handoff_mod.f90 with drv_conv_exit and
asserts on the printed verdicts, then pins the production
integration points by source inspection:

  * exhaustive truth table: all 1024 input combinations must return
    ALL=AND and NPASS=popcount;
  * case A: nine old gates pass, field not ready -> HOLD;
  * case B: all ten pass -> MET, NPASS=10 (100%);
  * case G: a nonzero f-iteration cap census holds the exit and
    names 'f-cap-exits';
  * case C: each single failure prints exactly its own label;
  * case D: multiple failures print all labels;
  * source pins: the solve_for_pops call passes the live
    LAM_FIELD_NOT_READY, the percent denominator is 9.0, the real
    CONV_EXIT_MET blocks are unchanged, and CONV_DEM_OK keeps its
    BA_DEM_EXIT_OK mirror semantics.

Run: python3 unit_tests/python/test_conv_exit.py
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common

REP = common.Report("test_conv_exit (ten-gate exit replica)")

WORK = common.WORKDIR / "conv_exit"
WORK.mkdir(parents=True, exist_ok=True)

binpath = common.compile_and_link_prod(
    "drv_conv_exit",
    [
        common.SRC / "tools" / "set_kind_module.f",
        common.SRC / "new_main" / "mod_subs" / "lam_handoff_mod.f90",
        common.F90DIR / "drv_conv_exit.f90",
    ],
)
r = common.run(binpath, cwd=WORK)
REP.check("driver ran", r.returncode == 0)

vals = {}
for ln in r.stdout.splitlines():
    m = re.match(r"^\s*([A-Za-z0-9_]+)=(.+)$", ln.strip())
    if m:
        vals[m.group(1)] = m.group(2).strip()


def tf(name):
    return vals.get(name, "").upper().startswith("T")


def field(name, tag):
    m = re.search(tag + r"(\d+)", vals.get(name, ""))
    return m.group(1) if m else ""


def labels(name):
    return vals.get(name, "").split()


def fails_part(name):
    v = vals.get(name, "")
    i = v.find("fails=")
    return v[i + 6 :] if i >= 0 else ""


# ---- exhaustive truth table --------------------------------------
REP.check("exhaustive 1024-combination table has no violations",
          vals.get("EXH_nbad") == "0")

# ---- case A: nine old gates pass, field not ready ----------------
REP.check("A: field-not-ready alone holds the exit",
          not tf("A_holdfield"))
REP.check("A: npass is 9 of 10",
          field("A_holdfield", "npass=") == "9")
REP.check("A: labels name lam-field-not-ready",
          "lam-field-not-ready" in fails_part("A_holdfield"))

# ---- case B: all ten pass ----------------------------------------
REP.check("B: all ten pass -> MET",
          tf("B_allpass"))
REP.check("B: npass is 10 -> 100%",
          field("B_allpass", "npass=") == "10")
REP.check("B: no failures reported",
          "all pass" in vals.get("B_allpass", ""))

# ---- case C: each single failure names its own gate --------------
single = {
    "C_lamstate": "lambda-state",
    "C_jfnk": "jfnk-unconverged",
    "C_exm": "ex_maxch>=eps",
    "C_dem": "exempted-demand",
    "C_hydro": "hydro-revision",
    "C_held": "pops-held",
    "C_pend": "resacc-pending",
    "C_onng": "on-ng-iteration",
}
all_labels = set(single.values()) | {"lam-field-not-ready"}
for case, want in single.items():
    got = set(fails_part(case).split())
    REP.check(f"{case}: fails -> {want}",
              not tf(case) and want in got and len(got - {want}) == 0)

# ---- case D: multiple failures all listed ------------------------
multi = set(fails_part("D_multi").split())
REP.check(
    "D: every failing gate is listed",
    not tf("D_multi")
    and multi
    == {
        "lambda-state",
        "ex_maxch>=eps",
        "jfnk-unconverged",
        "exempted-demand",
        "pops-held",
        "resacc-pending",
        "on-ng-iteration",
        "lam-field-not-ready",
    },
)

# ---- case E: forced-lambda leg of the lambda-state gate ----------
REP.check("E: rd_lam without auto passes the lambda gate",
          tf("E_forcedlam"))

# ---- case G: f-iteration cap census holds the exit ---------------
REP.check("G: cap count > 0 blocks the exit",
          not tf("G_cap"))
REP.check("G: npass is 9 of 10",
          field("G_cap", "npass=") == "9")
REP.check("G: labels name f-cap-exits",
          "f-cap-exits" in fails_part("G_cap"))


# ---- CONV_MARGIN_EVAL --------------------------------------------
def mcase(name):
    """Parse a M_* margin line into (act_m, ch, rel_m, rch, why, trb)."""
    v = vals.get(name, "")
    m = re.match(
        r"(\S+)\s+ch=(\S+)\s+rel=\s*(\S+)\s+rch=(\S+)\s+why=(\S+)\s+trb=\s*(\S+)",
        v,
    )
    if not m:
        return None
    return (float(m.group(1)), m.group(2), float(m.group(3)),
            m.group(4), m.group(5), float(m.group(6)))


# SAT-released exm: the reported defect -- 59.8-style dex demand must
# not win the binding slot; ex_maxch is the worst ACTIVE margin.
c = mcase("M_sat")
REP.check("margin: SAT-released exm stays out of the active slot",
          c is not None and c[1] == "ex_maxch" and c[3] == "exm"
          and c[4] == "SAT" and c[2] > 5.0)

# Unsaturated material exm is binding and wins the active slot.
c = mcase("M_bind")
REP.check("margin: unsaturated exm binds the active slot",
          c is not None and c[1] == "exm"
          and abs(c[0] - 2.6990) < 1e-3)

# Small exempted demand is released ('small'), ex_maxch binds.
c = mcase("M_small")
REP.check("margin: small exm released with reason 'small'",
          c is not None and c[1] == "ex_maxch" and c[3] == "exm"
          and c[4] == "small" and c[2] < 0.0)

# Stage-B need/release are release legs but are NOT SAT.
c = mcase("M_trcneed")
REP.check("margin: trc-need releases exm without calling it SAT",
          c is not None and c[1] == "ex_maxch" and c[3] == "exm"
          and c[4] == "trc-need")
c = mcase("M_trcrel")
REP.check("margin: trc-rel releases exm without calling it SAT",
          c is not None and c[1] == "ex_maxch" and c[3] == "exm"
          and c[4] == "trc-rel")

# SAT outranks a simultaneous trc-need in the reason string.
c = mcase("M_sat_trc")
REP.check("margin: SAT dominates the release reason",
          c is not None and c[4] == "SAT")

# exr binding: unsaturated relative-floor demand wins the active slot.
c = mcase("M_exr_bind")
REP.check("margin: unsaturated exr binds the active slot",
          c is not None and c[1] == "exr"
          and abs(c[0] - 2.9031) < 1e-3)

# exr SAT-released reports on the released row.
c = mcase("M_exr_sat")
REP.check("margin: SAT-released exr reports on the released row",
          c is not None and c[3] == "exr" and c[4] == "SAT")

# BA_EX_RELFRAC disabled -> exr=0 is absent from both slots.
c = mcase("M_exr_off")
REP.check("margin: BA_EX_RELFRAC off -> exr contributes nothing",
          c is not None and c[1] == "ex_maxch" and c[3] == "exm"
          and c[4] == "small")

# BA_DEM_PFLOOR disabled -> exm=0 is absent from both slots.
c = mcase("M_pfl_off")
REP.check("margin: BA_DEM_PFLOOR off -> exm contributes nothing",
          c is not None and c[1] == "ex_maxch" and c[3] == "none")

# Stage-B armed reports its own margin; off reports -99.
c = mcase("M_trabs")
REP.check("margin: armed stage-B reports its own margin",
          c is not None and abs(c[5] - 1.0) < 1e-3
          and c[1] == "ex_maxch")
c = mcase("M_trabs_hi")
REP.check("margin: stage-B demand never takes the active slot",
          c is not None and c[1] == "ex_maxch" and c[5] == 5.0)

# ---- source pins --------------------------------------------------
POPS = (common.SRC / "new_main" / "mod_subs" / "solve_for_pops.f").read_text()
CMF = (common.SRC / "new_main" / "cmfgen_sub.f").read_text()
MOD = (common.SRC / "new_main" / "mod_subs" / "lam_handoff_mod.f90").read_text()

call = re.search(
    r"CALL CONV_EXIT_EVAL\(RD_LAMBDA,DO_LAMBDA_AUTO,LAST_LAMBDA,"
    r"[\s\S]{0,300}?LAM_FIELD_NOT_READY,LAM_FIELD_NOT_CONVERGED,"
    r"[\s\S]{0,80}?CONV_LAM_OK,CONV_ALL,[\s\S]{0,80}?"
    r"CONV_NPASS,CONV_FAILS\)",
    POPS,
)
REP.check(
    "replica calls CONV_EXIT_EVAL with the live strict/weak flags",
    call is not None,
)
REP.check(
    "percent-converged denominator is ten gates",
    "CONV_NPASS/10.0_LDP" in POPS,
)
REP.check(
    "evaluator carries the ninth conjunct on the strict flag",
    ".NOT. FIELD_NCONV" in MOD
    and "lam-field-not-ready" in MOD
    and "lam-field-not-converged" in MOD,
)
REP.check(
    "evaluator carries the tenth conjunct on the cap census",
    "CAP_CNT .EQ. 0" in MOD and "f-cap-exits" in MOD,
)
# Iteration-1 suppressed warnings do not gate the exit on their own:
# they are a subset of the same census, so the gate reads
# FCONV_CAP_IT_CNT and never the suppressed count.
REP.check(
    "gate reads the cap census, not the suppressed subset",
    "FCONV_SUPPRESSED_CNT" not in CMF.split("CONV_EXIT_MET=")[1]
    or "FCONV_SUPPRESSED_CNT .EQ. 0" in CMF.split("CONV_EXIT_MET=")[1],
)
# The census is saved before the per-iteration reset so the count
# survives to the exit gate.
save = re.search(
    r"FCONV_CAP_IT_CNT=FCONV_CAP_CNT[\s\S]{0,200}?"
    r"FCONV_CAP_CNT=0",
    CMF,
)
REP.check("cap census saved before the per-iteration reset",
          save is not None)
# Both real exit blocks keep all ten conjuncts (the patch must not
# touch the production gate).
blocks = re.findall(
    r"CONV_EXIT_MET=\([\s\S]{0,400}?\.NOT\. LAM_FIELD_NOT_CONVERGED"
    r"[\s\S]{0,80}?FCONV_CAP_IT_CNT \.EQ\. 0",
    CMF,
)
REP.check(
    "both CONV_EXIT_MET blocks keep the strict field verdict",
    len(blocks) >= 2,
)
for tok in (
    "JFNK_UNCONVERGED_COMMIT",
    "FCONV_CAP_IT_CNT .EQ. 0",
    "BA_EX_MAXCH .LT. EPS",
    "BA_DEM_EXIT_OK",
    "DONE_HYDRO_REVISION",
    "GS_POPS_HELD",
    "BAACC_PENDING",
    "MAIN_COUNTER .NE. LAST_NG",
):
    REP.check(
        f"exit blocks keep {tok}",
        all(tok in b for b in blocks),
    )
# CONV_DEM_OK mirror semantics preserved verbatim.
dem = re.search(
    r"CONV_DEM_OK=\(BA_DEM_EXM \.LT\. MAX\(10\.0_LDP\*EPS,1\.0_LDP\)"
    r"[\s\S]{0,300}?BA_DEM_EXR \.GE\. BA_SAT_RAW\)",
    POPS,
)
REP.check("CONV_DEM_OK mirrors BA_DEM_EXIT_OK", dem is not None)

# Margin reporter: the replica derives its classification from the
# shared evaluator, not from raw threshold comparisons.
REP.check(
    "margin reporter calls CONV_MARGIN_EVAL",
    "CALL CONV_MARGIN_EVAL(EPS,BA_EX_MAXCH,BA_DEM_EXM,BA_DEM_EXR" in POPS,
)
REP.check(
    "margin reporter prints the released demand on its own row",
    "'conv marg'" in POPS and "worst active margin" in POPS,
)
REP.check(
    "old unclassified margin loop is gone",
    "LOG10(BA_DEM_EXM/MAX(10.0_LDP*EPS" not in POPS,
)
REP.check(
    "classifier mirrors the release predicates",
    "EXM .GE. SAT_RAW" in MOD and "TRC_NEED" in MOD
    and "TRC_REL" in MOD and "EXM .LT. SMALL" in MOD,
)


def main():
    return REP.summary()


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)
