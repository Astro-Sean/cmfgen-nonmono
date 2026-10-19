#!/usr/bin/env python3
"""Static regression checks for the KNOP lagged-field drift fixes.

Reads the production Fortran sources as text and asserts the guards
landed: fail-closed refresh commit defaults, the NEW_FREQ sweep-write
gate, the fully-poisoned-row heal, the rollback-restore value audit,
the FG_TRIAL hold exclusion, and the raw-release counter wiring.
"""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FG = (ROOT / "plane/fg_j_cmf_nonmono_v1.f").read_text()
POPS = (ROOT / "new_main/mod_subs/solve_for_pops.f").read_text()
RD_CONTROL = (ROOT / "new_main/mod_subs/rd_control_variables.f").read_text()


def section(source, start, end):
    first = source.index(start)
    last = source.index(end, first)
    return source[first:last]


def check(name, passed):
    if not passed:
        raise AssertionError(name)
    print(f"[PASS] {name}")


def test_heal_store_zeroes_fully_poisoned_rows():
    body = section(FG, "SUBROUTINE HEAL_STORE", "RETURN\n\tEND")
    check(
        "HEAL_STORE no longer skips fully poisoned rows",
        "IF(NV .EQ. 0)CYCLE" not in body,
    )
    check(
        "HEAL_STORE zeroes the row when no sane slot exists",
        re.search(r"NV \.EQ\. 0\)THEN[\s\S]*?STORE\(KC,LS,LS2\)=0\.0_LDP", body)
        is not None,
    )


def test_sweep_store_write_gated_on_new_freq():
    block = section(FG, "Store intensities at depth grid points", "DO K=1,ND")
    check(
        "in-place I_NEXT write block is gated on NEW_FREQ",
        re.search(
            r"IF\(HAS_COMPRESSING \.AND\. NEW_FREQ \.AND\.", block
        )
        is not None,
    )


def test_rollback_restore_validates_field_values():
    body = section(
        FG,
        "SUBROUTINE RESTORE_FG_INEXT_LATEST",
        "END SUBROUTINE RESTORE_FG_INEXT_LATEST"
        if "END SUBROUTINE RESTORE_FG_INEXT_LATEST" in FG
        else "Observation-frame restore",
    )
    accept = body.index("OK=.TRUE.")
    audit = body.index("FG_GS_CAP_FILE .LT. 0")
    check(
        "rollback restore audits values before accepting",
        audit < accept
        and "I_P_DEPTH_NEXT .LT. 0.0_LDP" in body[:accept]
        and "I_M_DEPTH_NEXT .LT. 0.0_LDP" in body[:accept],
    )


def test_fail_closed_commit_defaults():
    check(
        "FG_COMMIT_CONV defaults to TRUE",
        "FG_COMMIT_CONV_V=.TRUE." in RD_CONTROL
        and "FG_COMMIT_CONV=.TRUE." in FG,
    )
    check(
        "FG_COMMIT_NC defaults to 0.05",
        "FG_COMMIT_NC_V=0.05_LDP" in RD_CONTROL
        and "FG_COMMIT_NC=0.05_LDP" in FG,
    )
    check(
        "FG_BULK_CONTR defaults to TRUE",
        "FG_BULK_CONTR_V=.TRUE." in RD_CONTROL
        and "FG_BULK_CONTR=.TRUE." in FG,
    )


def test_fg_trial_excluded_from_population_hold():
    check(
        "GS_HOLD_POPS condition excludes FG_TRIAL",
        re.search(
            r"GS_HOLD_POPS \.AND\. GS_FREEZE \.AND\. \.NOT\. MOM_TRIAL"
            r" \.AND\.\s*\n\s*\d?\s*\.NOT\. FG_TRIAL",
            POPS,
        )
        is not None,
    )


def test_rel_damp_window_default_and_raw_release_counter():
    check(
        "GS_REL_N defaults to 2",
        re.search(r"GS_REL_N=2\s+CALL RD_STORE_INT", RD_CONTROL) is not None,
    )
    check(
        "raw-cure releases counted into GS_STATE",
        "GS_RAW_RELEASES_TOT=GS_RAW_RELEASES_TOT+1" in POPS
        and "'n_raw_releases='" in POPS,
    )


if __name__ == "__main__":
    test_heal_store_zeroes_fully_poisoned_rows()
    test_sweep_store_write_gated_on_new_freq()
    test_rollback_restore_validates_field_values()
    test_fail_closed_commit_defaults()
    test_fg_trial_excluded_from_population_hold()
    test_rel_damp_window_default_and_raw_release_counter()
    print("all KNOP lagged-field checks passed")
