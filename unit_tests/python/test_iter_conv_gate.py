#!/usr/bin/env python3
"""Static checks on the iteration-exit gate and field-maturity coverage.

NUMERICAL_VERIFICATION_TEST_PLAN TEST 11/14 (source level): the
convergence exit must not fire on small population demands while the
lagged radiation field is still evolving.

Checks pin the current gate structure. Both audit defects this file
tracked are now landed and asserted positively:

  T-01 F-01: FIXED. Both CONV_EXIT_MET blocks conjunct .NOT.
             LAM_FIELD_NOT_READY -- the per-call field-readiness
             verdict (held/frozen/raw-immature/stability window)
             published by SOLVE_FOR_POPS; the deferred-exit
             diagnostic mirrors it. Forced F on stock and plain
             SOB+LAM_ES=T, so those gates are unchanged.
  T-02 F-16: FIXED. The census block opens under HAS_LAGGED_FIELD,
             which covers DO_NONMONO_SOB + USE_LAM_ES=F alongside
             KNOP and DO_SOB_FG_KNOP, so the moment field's maturity
             now reaches the demote gate and the handoff field gate.

Static only: CONV_EXIT_MET and the census gate cannot be exercised
standalone.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import TREE, Report

CMF = (TREE / "new_main" / "cmfgen_sub.f").read_text()
POPS = (TREE / "new_main" / "mod_subs" / "solve_for_pops.f").read_text()

REP = Report("test_iter_conv_gate (exit-gate maturity coverage)")


def conv_blocks(text):
    """Each CONV_EXIT_MET assignment with its continuation lines.

    The source is tab-indented fixed form: continuation markers are a
    single non-zero digit as the first non-blank character, followed
    by whitespace.  A label is a digit run immediately followed by
    statement text, so `100 IF(` does not match.
    """
    cont = re.compile(r"^\s*[1-9]\s")
    blocks = []
    for m in re.finditer(r"^\s+CONV_EXIT_MET\s*=", text, re.M):
        lines = text[m.start() :].split("\n")
        body = [lines[0]]
        for ln in lines[1:]:
            if cont.match(ln):
                body.append(ln)
            else:
                break
        blocks.append("\n".join(body))
    return blocks


BLOCKS = conv_blocks(CMF)
JOINED = "\n".join(BLOCKS)


def test_conv_exit_structure():
    REP.check(
        "CONV_EXIT_MET assigned at both 9999 exit paths",
        len(BLOCKS) >= 2,
        f"{len(BLOCKS)} blocks",
    )
    REP.check(
        "every exit block prices the raw demanded correction",
        all("BA_EX_MAXCH" in b for b in BLOCKS),
    )
    REP.check(
        "existing hold conjuncts intact (held field already blocks exit)",
        "GS_POPS_HELD" in JOINED,
    )


def test_maturity_predicates_exist():
    for tok in ("GS_FG_IMMATURE", "MOM_UNSAFE_B", "MOM_FREEZE_B"):
        REP.check(f"maturity predicate {tok} computed", tok in POPS)


def test_T01_exit_maturity_conjunct():
    # F-01 landed: both exit blocks carry .NOT. LAM_FIELD_NOT_READY,
    # the per-call verdict from LAM_FIELD_GATE (held .OR. frozen
    # .OR. HAS_LAGGED_FIELD .AND. (raw-immature .OR. under the
    # N_MATURE/GS_REL_STABLE window)). Forced F off the lagged-field
    # routes, so stock and SOB+LAM_ES=T are unchanged.
    REP.check(
        "T-01: every CONV_EXIT_MET block gates on LAM_FIELD_NOT_READY",
        len(BLOCKS) >= 2
        and all(".NOT. LAM_FIELD_NOT_READY" in b for b in BLOCKS),
    )
    # The deferred-exit diagnostic mirrors the term so its "only
    # resacc pending" claim stays true under an immature field.
    REP.check(
        "T-01: deferred-exit diagnostic mirrors the maturity term",
        re.search(
            r"CONV_EXIT_MET \.AND\. BAACC_PENDING[\s\S]{0,300}?"
            r"\.NOT\. LAM_FIELD_NOT_READY",
            CMF,
        )
        is not None,
    )
    # The verdict is published per call, unconditionally (before the
    # RD_LAMBDA branches), so the gate never reads a stale value.
    REP.check(
        "T-01: verdict assigned unconditionally each SOLVE_FOR_POPS call",
        re.search(
            r"LAM_FIELD_NOT_READY=LAM_FIELD_GATE\(HAS_LAGGED_FIELD,"
            r"GS_POPS_HELD,\s*\n\s*1\s+GS_FREEZE,GS_RAW_IMMATURE,"
            r"N_MATURE,GS_REL_STABLE\)",
            POPS,
        )
        is not None,
    )
    # On stock the gate reduces to held/frozen flags -- both forced F
    # off KNOP -- so the new conjunct cannot alter the stock exit.
    REP.check(
        "T-01: LAM_FIELD_GATE reduces to held/frozen without lagged field",
        re.search(
            r"LAM_FIELD_GATE=HELD \.OR\. FREEZE",
            (TREE / "new_main" / "mod_subs" / "lam_handoff_mod.f90").read_text(),
        )
        is not None,
    )


def test_T02_sob_mom_census_fixed():
    # F-16 landed: the census block opens under HAS_LAGGED_FIELD,
    # which RD_CONTROL_VARIABLES derives as KNOP .OR. SOB_FG_KNOP
    # .OR. (SOB .AND. .NOT. USE_LAM_ES) -- the SOB moment route is
    # covered. The getters are pure read-backs so the widened gate
    # only admits the metrics; GS_FREEZE stays forced-F off KNOP.
    gate = re.search(
        r"IF\(HAS_LAGGED_FIELD\)THEN\s*\n\s*CALL GET_FG_GS_MATURITY",
        POPS,
    )
    REP.check(
        "T-02: maturity census gated on HAS_LAGGED_FIELD "
        "(SOB+LAM_ES=F covered)",
        gate is not None,
    )
    # The demote gate judges the moment leg on the widened block:
    # MOM_RAN .AND. MOM_UNSAFE_B demotes a scheduled BA step.
    sob_demote = re.search(
        r"MOM_RAN \.AND\. MOM_UNSAFE_B", POPS
    )
    REP.check(
        "T-02: moment-field demote gate live on the census block",
        sob_demote is not None,
    )
    # The lagged-field route flag must include the SOB moment route.
    rcv = (TREE / "new_main" / "mod_subs" / "rd_control_variables.f").read_text()
    REP.check(
        "T-02: HAS_LAGGED_FIELD includes SOB .AND. .NOT. USE_LAM_ES",
        re.search(
            r"HAS_LAGGED_FIELD=DO_NONMONO_KNOP[\s\S]{0,200}?"
            r"DO_NONMONO_SOB \.AND\. \.NOT\. USE_LAM_ES",
            rcv,
        )
        is not None,
    )


def main():
    test_conv_exit_structure()
    test_maturity_predicates_exist()
    test_T01_exit_maturity_conjunct()
    test_T02_sob_mom_census_fixed()
    return REP.summary()


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)
