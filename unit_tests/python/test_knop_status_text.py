#!/usr/bin/env python3
"""Static regression checks for the SOB/moment status-text rewrite.

Reads the production Fortran sources as text and asserts the
sentence-based narration landed: the shared KNOP_STATUS_MOD helpers
exist with the binding sentences, the old column codes and bare
booleans are gone from the human OUTGEN blocks, and the machine
records (GS_STATE keys, MOMJNXT_GS layout, ITER_SUM/ITER_HISTORY
fields) are untouched.
"""

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MOD = (ROOT / "new_main/mod_subs/knop_status_mod.f90").read_text()
MOM = (ROOT / "plane/mom_j_cmf_nonmono_v1.f").read_text()
POPS = (ROOT / "new_main/mod_subs/solve_for_pops.f").read_text()
RD_CONTROL = (ROOT / "new_main/mod_subs/rd_control_variables.f").read_text()
CJB = (ROOT / "new_main/mod_subs/comp_j_blank.f").read_text()
CMF = (ROOT / "new_main/cmfgen_sub.f").read_text()


FAILED = []


def check(msg, cond):
    if cond:
        print(f"  ok  {msg}")
    else:
        print(f"  FAIL {msg}")
        FAILED.append(msg)


def section(source, start, end):
    first = source.index(start)
    last = source.index(end, first)
    return source[first:last]


def test_status_module_and_binding_sentences():
    for fn in (
        "mom_state_sentence",
        "mom_sweep_sentence",
        "mom_refresh_sentence",
        "mom_pass_sentence",
        "mom_restart_sentence",
        "mom_direction_sentence",
    ):
        check(f"module exports {fn}", f"public :: {fn}" in MOD)
    for frag in (
        "holding populations while the lagged J_NEXT field converges",
        "populations may move this step (hold disabled)",
        "catastrophic drift; holding populations",
        "drifting over most cells - holding populations",
        "trying a full step with the rollback net armed",
        "taking a live LAMBDA step",
        "Moment refresh converged at pass",
        "bulk-converged at pass",
        "Moment refresh stalled at pass",
        "Moment refresh diverged at pass",
        "retained the best pass",
        "retained the pre-refresh field",
        "Moment refresh finished:",
        "Moment pass ",
        "J_NEXT checkpoint: absent - cold start",
        "present, accepted",
        "stale or from another model",
        "dimensions or the expanded grid do not match",
        "non-finite or negative; rebuilding from zero",
        "compressing direction only",
        "expanding direction only",
        "alternating compressing and expanding directions",
    ):
        check(f"binding sentence '{frag[:50]}...'", frag in MOD)


def test_mom_sweep_table_is_prose():
    body = section(MOM, "MOM_J_CMF_NONMONO_V1: GS sweep", "MOMJNXT_GS.tmp")
    for token in (
        "n50=",
        "n5=",
        "negJ=",
        "tiny=",
        "pquar=",
        "nonMM=",
        "negEPS=",
        "eps_cap=",
        "dom_cmp=",
        "dom_exp=",
        "conv=",
        "dir=",
        "cmp_frac=",
        "nsupp=",
        "cols:",
    ):
        check(f"sweep table drops '{token}'", token not in body)
    check("sweep table keeps SW_ST", "SW_ST" in body)
    check("sweep table keeps SV_ST", "SV_ST" in body)
    check("direction sentence used", "MOM_DIRECTION_SENTENCE" in body)
    check("refresh delta sentence used", "'delta'" in body)


def test_mom_pass_and_exit_prose():
    check(
        "pass header names directions",
        "DIR_W='compressing direction'" in MOM
        and "DIR_W='expanding direction'" in MOM,
    )
    for token in ("'n5='", "'n50='", "'match='", "'off='", "'cap='"):
        check(f"pass table drops {token}", token not in MOM)
    check("compressing stall sentence", "Compressing direction" in MOM)
    check("expanding stall sentence", "Expanding direction" in MOM)
    check("both-directions stall", "Both directions stalled" in MOM)
    check("converged sentence", "'converged'" in MOM)
    check("bulk sentence", "'bulk'" in MOM)
    check("diverge sentence", "'diverge'" in MOM)
    check("keep-best sentence", "'keep_best'" in MOM)
    check("keep-prior sentence", "'keep_prior'" in MOM)
    check("Aitken summary prose", "inner Aitken acceleration" in MOM)


def test_solve_for_pops_prose():
    check("module imported", "USE KNOP_STATUS_MOD" in POPS)
    check(
        "held sentence wired",
        "MOM_STATE_SENTENCE(.TRUE.,GS_HOLD_POPS" in POPS,
    )
    check(
        "pops census sentence",
        "MOM_SWEEP_SENTENCE('pops'" in POPS,
    )
    check(
        "sweep census sentence",
        "MOM_SWEEP_SENTENCE('sweep'" in POPS,
    )
    check(
        "refresh delta sentence",
        "MOM_REFRESH_SENTENCE(MOM_RF_CONV" in POPS,
    )
    body = section(POPS, "KNOP GS SWEEP", "marginal-release BA trial")
    for token in ("nbad5", "negJ / quar", "converged", "n5 / n50", "MOM bulk"):
        check(f"GS sweep drops '{token}'", token not in body)
    check(
        "trial terms use yes/no clauses",
        "moment-trial gate terms" in POPS and "MERGE('yes','no '" in POPS,
    )
    check("no 'MOM_TRIAL terms' blob", "MOM_TRIAL terms" not in POPS)


def test_mode_and_format_lines():
    check(
        "LAM_ES narration fixed",
        "Field solver: ray solutions with lambda" in CJB
        and "equatons" not in CJB,
    )
    check(
        "moment-solver mode line",
        "Field solver: moment equations (MOM_J_CMF_NONMONO_V1)" in CJB,
    )
    check(
        "SOB auto-set prose",
        "using the SOB line transfer" in RD_CONTROL
        and "lambda iteration with ray solutions enabled" in RD_CONTROL,
    )
    check(
        "iteration count widened to I8",
        "'(A,I8)')' Current great iteration count is" in CMF,
    )
    check(
        "no narrow I4 on the count line",
        not re.search(
            r"\(A,I4\)\)' Current great iteration count is", CMF
        ),
    )


def test_machine_contracts_unchanged():
    check(
        "GS_STATE mom_refresh_conv kept",
        "'mom_refresh_conv='" in POPS,
    )
    check(
        "GS_STATE pops_held kept",
        "'pops_held='" in POPS,
    )
    check(
        "ITER_SUM header keys kept",
        "' ITER_SUM: it=',MAIN_COUNTER,'  mode='" in POPS
        and "'  held=',GS_POPS_HELD" in POPS,
    )
    check(
        "audit demand row kept",
        "MAXCH demand / applied" in POPS,
    )
    check(
        "ITER_HISTORY mom columns kept",
        "mom_n5f mom_n50f mom_negjf mom_quarf" in POPS,
    )
    check("checkpoint MAGIC kept", "MAGIC=20260916" in MOM)
    check("checkpoint VER records kept", "VER .NE. 3" in MOM)
    check(
        "VADAT echo keys kept",
        "RD_STORE_LOG" in RD_CONTROL,
    )


def test_synthetic_block_still_parses():
    """Feed a synthetic SOB-flavoured OUTGEN through the audit parser."""
    audit = ROOT / "tests/python/jfnk_outgen_audit.py"
    if not audit.exists():
        check("audit script present", False)
        return
    sys.path.insert(0, str(audit.parent))
    import jfnk_outgen_audit

    block = "\n".join(
        [
            " T [DO_JFNK]",
            " JFNK_MAX_EVALS  40",
            " Current great iteration count is    1057",
            " ITER_TYPE it= 1057  Running LAMBDA",
            " ------------- Moment refresh pass 4 - compressing direction ---",
            "  [OK]   Moment refresh converged at pass 7 (residual 1.8E-03).",
            " ITER_SUM: it= 1057  mode=BA  frozen=T  held=F",
            "    MAXCH demand / applied        :  2.0E+03   1.5E+02",
            " Total F evaluations: 40",
        ]
    )
    import tempfile

    with tempfile.NamedTemporaryFile(
        "w", suffix="_OUTGEN", delete=False
    ) as fh:
        fh.write(block)
        tmp = Path(fh.name)
    try:
        out = jfnk_outgen_audit.parse_outgen(tmp)
    finally:
        tmp.unlink()
    check("iteration count parsed", 1057 in [it for it, _ in out["iter_types"]])
    check("modctl row parsed", len(out["modctl"]) == 1)
    check("DO_JFNK parsed", out["do_jfnk"])
    check("JFNK_MAX_EVALS parsed", out["max_evals"] == 40)


if __name__ == "__main__":
    test_status_module_and_binding_sentences()
    test_mom_sweep_table_is_prose()
    test_mom_pass_and_exit_prose()
    test_solve_for_pops_prose()
    test_mode_and_format_lines()
    test_machine_contracts_unchanged()
    test_synthetic_block_still_parses()
    if FAILED:
        print(f"{len(FAILED)} checks failed")
        sys.exit(1)
    print("all KNOP status-text checks passed")
