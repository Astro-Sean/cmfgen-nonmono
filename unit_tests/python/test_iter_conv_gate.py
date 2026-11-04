#!/usr/bin/env python3
"""Static checks on the iteration-exit gate and field-maturity coverage.

NUMERICAL_VERIFICATION_TEST_PLAN TEST 11/14 (source level): the
convergence exit must not fire on small population demands while the
lagged radiation field is still evolving.

Checks pin the CURRENT gate structure (PASS) and the documented
defects (XFAIL -- each asserts the defect is still present, so it
fails loudly the moment T-01/T-02 lands and must then be flipped to
the positive assertion):

  T-01 F-01: CONV_EXIT_MET on nonmono paths carries no radiation-
             maturity conjunct.
  T-02 F-16: the GS/MOM maturity census in SOLVE_FOR_POPS opens only
             under DO_NONMONO_KNOP .OR. DO_SOB_FG_KNOP, so the
             SOB+USE_LAM_ES=F moment path runs with no immature-
             field protection.

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

MATURITY_RE = re.compile(
    r"FIELD_MATURE|FIELD_IMMATURE|GS_FG_IMMATURE|MOM_FREEZE_B|"
    r"MOM_UNSAFE_B|RF_MATURE|GS_FIELD"
)


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


def test_T01_exit_lacks_maturity_xfail():
    # Defect F-01 present => PASS; flips FAIL when T-01 lands, at
    # which point this check is rewritten to assert the conjunct.
    REP.check(
        "XFAIL T-01: no CONV_EXIT_MET block references a field-"
        "maturity predicate",
        not any(MATURITY_RE.search(b) for b in BLOCKS),
    )
    m = re.search(
        r"CONV_EXIT_MET[^!]*?DO_NONMONO[^!]*?(MATURE|IMMATURE)", CMF, re.S
    )
    REP.check(
        "XFAIL T-01: no nonmono-scoped maturity term in the exit gate",
        m is None,
    )


def test_T02_sob_mom_outside_census_xfail():
    # F-16: the census block opens under KNOP/DO_SOB_FG_KNOP only.
    gate = re.search(
        r"IF\(DO_NONMONO_KNOP \.OR\. DO_SOB_FG_KNOP\)THEN\s*\n\s*"
        r"CALL GET_FG_GS_MATURITY",
        POPS,
    )
    REP.check(
        "XFAIL T-02: maturity census still KNOP-scoped "
        "(SOB+LAM_ES=F uncovered)",
        gate is not None,
    )
    sob_demote = "mom_immature_demote" in POPS or re.search(
        r"DO_NONMONO_SOB[^\n]*\n[^!]*?MOM_(UNSAFE_B|FREEZE_B)", POPS
    )
    REP.check(
        "XFAIL T-02: no SOB-path demote/census on the moment field",
        not sob_demote,
    )


def main():
    test_conv_exit_structure()
    test_maturity_predicates_exist()
    test_T01_exit_lacks_maturity_xfail()
    test_T02_sob_mom_outside_census_xfail()
    return REP.summary()


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)
