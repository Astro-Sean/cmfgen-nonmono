#!/usr/bin/env python3
"""Static checks for the POP_AUDIT population diagnostic wiring.

Verifies that the VADAT keys are declared in control_variable_mod.f
and parsed with correct defaults in rd_control_variables.f, that
SOLVE_FOR_POPS calls POP_AUDIT_V1 under the flag at the end of the
call, and that the audit routine itself is general: it iterates
SPECIES_BEG_ID..SPECIES_END_ID (never a hard-coded species), treats
the terminal XzV_PRES=F ion slot as a single density row, sums every
stage into the conservation check regardless of the ion filter, and
compares against POP_SPECIES -- the quantity EQ_SPECIES enforces.
"""

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TREE = HERE.parent.parent
MODSUBS = TREE / "new_main" / "mod_subs"

SFP = (MODSUBS / "solve_for_pops.f").read_text()
CVM = (MODSUBS / "control_variable_mod.f").read_text()
RCV = (MODSUBS / "rd_control_variables.f").read_text()

FAILED = []


def check(name, ok):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        FAILED.append(name)


def test_vadat_keys():
    for key in ("POP_AUDIT", "POP_AUDIT_LEV"):
        check(f"{key} declared logical in control_variable_mod",
              re.search(rf"^\s*LOGICAL {key}\b", CVM, re.M) is not None)
        check(f"{key} parsed with FALSE default",
              re.search(rf"{key}=\.FALSE\.\s*\n\s*CALL RD_STORE_LOG\({key}",
                        CVM) is None
              and f"{key}=.FALSE." in RCV and f"RD_STORE_LOG({key}" in RCV)
    for key in ("POP_AUDIT_SPEC", "POP_AUDIT_ION", "POP_AUDIT_D"):
        check(f"{key} declared integer",
              re.search(rf"^\s*INTEGER {key}\b", CVM, re.M) is not None)
        check(f"{key} parsed with 0 default and clamped",
              f"{key}=0" in RCV and f"RD_STORE_INT({key}" in RCV
              and f"IF({key} .LT. 0){key}=0" in RCV)


def test_call_site():
    check("POP_AUDIT call guarded at end of SOLVE_FOR_POPS",
          "IF(POP_AUDIT)CALL POP_AUDIT_V1(POPS,NT,ND,MAIN_COUNTER,LUER)"
          in SFP)
    pos_call = SFP.find("IF(POP_AUDIT)CALL POP_AUDIT_V1")
    pos_end = SFP.find("\n\tRETURN\n\tEND", pos_call)
    check("call placed before final RETURN", 0 < pos_call < pos_end)


def test_audit_routine_generality():
    body = SFP[SFP.find("SUBROUTINE POP_AUDIT_V1"):]
    check("iterates SPECIES_BEG_ID..SPECIES_END_ID",
          "SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)" in body)
    check("uses SPECIES_PRES presence flag", "SPECIES_PRES(ISPEC)" in body)
    check("terminal ion handled via XzV_PRES branch",
          ".NOT." not in body or "XzV_PRES" in body)
    check("sums SL rows from EQXzV", "POPS(R0+L-1,K)" in body)
    check("terminal slot reads single density row",
          "SUM_ION=POPS(R0,K)" in body)
    check("conservation vs POP_SPECIES", "POP_SPECIES(K,ISPEC)" in body)
    check("ion filter does not affect conservation sum",
          body.find("TOT=TOT+SUM_ION") < body.find("IF(ION_SHOW)"))
    check("departure coefficient from XzVLTE", "XzVLTE" in body)
    check("DXzV aliased-row staleness check",
          "POPS(R0+NLVL,K)" in body)
    check("negative level census", "JNEG=JNEG+1" in body)
    check("non-finite census (x-x != 0)", "PLEV-PLEV .NE. 0.0_LDP" in body)
    check("no hard-coded species/ion strings",
          "'OV'" not in body and "'OXY'" not in body
          and "ZXzV .EQ." not in body)


def test_audit_comment_documents_layout():
    body = SFP[SFP.find("! POP_AUDIT_V1"):SFP.find("SUBROUTINE POP_AUDIT_V1")]
    check("comment documents DXzV = next-ion aliased row",
          "aliased" in body and "next ion" in body)


def main():
    test_vadat_keys()
    test_call_site()
    test_audit_routine_generality()
    test_audit_comment_documents_layout()
    if FAILED:
        print(f"test_pop_audit: {len(FAILED)} checks failed")
        sys.exit(1)
    print("test_pop_audit: all checks passed")


if __name__ == "__main__":
    main()
