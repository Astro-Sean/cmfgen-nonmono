#!/usr/bin/env python3
"""test_conv_exit -- nine-gate ITER_SUM converged-exit replica.

Round-2 phase-2 patch A: the ITER_SUM replica was missing the ninth
conjunct (.NOT. LAM_FIELD_NOT_READY) that the real CONV_EXIT_MET in
CMFGEN_SUB carries, so an immature lagged field could print [MET]
while the true exit still held. The replica now calls
CONV_EXIT_EVAL (LAM_HANDOFF_MOD), the same evaluator both sides
share.

Compiles the production lam_handoff_mod.f90 with drv_conv_exit and
asserts on the printed verdicts, then pins the production
integration points by source inspection:

  * exhaustive truth table: all 512 input combinations must return
    ALL=AND and NPASS=popcount;
  * case A: eight old gates pass, field not ready -> HOLD;
  * case B: all nine pass -> MET, NPASS=9 (100%);
  * case C: each single failure prints exactly its own label;
  * case D: multiple failures print all labels;
  * cases F/G: the ES_NRDY veto holds the exit and names itself
    without joining the nine-gate NPASS;
  * source pins: the solve_for_pops call passes the live
    LAM_FIELD_NOT_READY and the resolved ES conjunct, the percent
    denominator is 9.0, the real CONV_EXIT_MET blocks are
    unchanged, and CONV_DEM_OK keeps its BA_DEM_EXIT_OK mirror
    semantics.

Run: python3 unit_tests/python/test_conv_exit.py
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common

REP = common.Report("test_conv_exit (nine-gate exit replica)")

WORK = common.WORKDIR / "conv_exit"
WORK.mkdir(parents=True, exist_ok=True)

binpath = common.compile_and_link_prod(
    "drv_conv_exit",
    [
        common.TREE / "tools" / "set_kind_module.f",
        common.TREE / "new_main" / "mod_subs" / "lam_handoff_mod.f90",
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
REP.check("exhaustive 512-combination table has no violations",
          vals.get("EXH_nbad") == "0")

# ---- case A: eight old gates pass, field not ready ---------------
REP.check("A: field-not-ready alone holds the exit",
          not tf("A_holdfield"))
REP.check("A: npass is 8 of 9",
          field("A_holdfield", "npass=") == "8")
REP.check("A: labels name lam-field-not-ready",
          "lam-field-not-ready" in fails_part("A_holdfield"))

# ---- case B: all nine pass ---------------------------------------
REP.check("B: all nine pass -> MET",
          tf("B_allpass"))
REP.check("B: npass is 9 -> 100%",
          field("B_allpass", "npass=") == "9")
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

# ---- cases F/G: the ES_NRDY veto (W3, 24-Nov-2026) ---------------
# All nine gates pass but the resolved ES veto is armed: the exit
# holds, the label names the ES census, and NPASS stays the nine
# physical gates -- the conditional gate is not a tenth point of
# the percent figure.
REP.check("F: es_nrdy alone holds the exit", not tf("F_esnrdy"))
REP.check("F: npass stays the nine physical gates",
          field("F_esnrdy", "npass=") == "9")
REP.check("F: label names es-ali-unconverged",
          "es-ali-unconverged" in fails_part("F_esnrdy"))
got = set(fails_part("G_esheld").split())
REP.check(
    "G: es_nrdy + pops-held lists both labels",
    not tf("G_esheld")
    and "es-ali-unconverged" in got
    and "pops-held" in got,
)


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
POPS = (common.TREE / "new_main" / "mod_subs" / "solve_for_pops.f").read_text()
CMF = (common.TREE / "new_main" / "cmfgen_sub.f").read_text()
MOD = (common.TREE / "new_main" / "mod_subs" / "lam_handoff_mod.f90").read_text()

call = re.search(
    r"CALL CONV_EXIT_EVAL\(RD_LAMBDA,DO_LAMBDA_AUTO,LAST_LAMBDA,"
    r"[\s\S]{0,300}?LAM_FIELD_NOT_READY,"
    r"[\s\S]{0,120}?CONV_LAM_OK,CONV_ALL,"
    r"[\s\S]{0,80}?CONV_NPASS,CONV_FAILS\)",
    POPS,
)
REP.check(
    "replica calls CONV_EXIT_EVAL with the live LAM_FIELD_NOT_READY",
    call is not None,
)
REP.check(
    "replica passes the resolved ES_NRDY (ES_ALI_GATE=HARD arm)",
    "ES_FIELD_NOT_READY .AND. ES_ALI_GATE .EQ. 'HARD'" in POPS,
)
REP.check(
    "percent-converged denominator is nine gates",
    "CONV_NPASS/9.0_LDP" in POPS,
)
REP.check(
    "evaluator carries the ninth conjunct",
    ".NOT. FIELD_NRDY" in MOD and "lam-field-not-ready" in MOD,
)
REP.check(
    "evaluator carries the ES_NRDY veto and label",
    ".NOT. ES_NRDY" in MOD and "es-ali-unconverged" in MOD,
)
# Both real exit blocks keep all nine conjuncts plus the ES_NRDY
# clause (W3): capture through the ES_ALI_GATE conjunct so a block
# missing it fails the count.
blocks = re.findall(
    r"CONV_EXIT_MET=\([\s\S]{0,500}?ES_ALI_GATE \.EQ\. 'HARD'\)", CMF
)
REP.check(
    "both CONV_EXIT_MET blocks keep LAM_FIELD_NOT_READY + ES gate",
    len(blocks) >= 2
    and all(".NOT. LAM_FIELD_NOT_READY" in b for b in blocks),
)
for tok in (
    "JFNK_UNCONVERGED_COMMIT",
    "BA_EX_MAXCH .LT. EPS",
    "BA_DEM_EXIT_OK",
    "DONE_HYDRO_REVISION",
    "GS_POPS_HELD",
    "BAACC_PENDING",
    "MAIN_COUNTER .NE. LAST_NG",
    "ES_FIELD_NOT_READY .AND. ES_ALI_GATE .EQ. 'HARD'",
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
