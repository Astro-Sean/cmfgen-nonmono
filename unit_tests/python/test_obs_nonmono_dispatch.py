#!/usr/bin/env python3
"""Static regression checks for the OBS (CMF_FLUX) non-mono fixes.

Pins the OBS-side changes made 09-Nov-2026:

1. INC_GRID/ACCURATE hard rejection on a genuinely non-monotonic
   velocity field (N_NONMONO_OBS>0), inside cmf_flux_sub_v5.f.
2. COMP_JCONT_V4.INC dispatch: under DO_NONMONO_SOB the moment solve
   is MOM_J_CMF_NONMONO_V1 (the population solver), not the clamped
   MOM_J_CMF_V8; stock V8 retained when the flag is off.
3. The FIRST_FREQ one-sided dLOG_NU fix that MOM_J_CMF_NONMONO_V1
   needs at the bluest frequency (chi_hat consumes dLOG_NU even
   under INIT).
4. DO_NONMONO_FREQ_MATCH obs-side mirror of the VADAT auto-set
   (absent -> T, explicit F -> overridden with a banner) under
   DO_NONMONO_KNOP or DO_NONMONO_SOB+USE_LAM_ES=F, forwarded to
   SET_FG_FREQ_MATCH / SET_MOM_FREQ_MATCH.
5. OBS_FRAME_SUB_V9 edge-clamp census (N_EDGE_CLAMP etc.) at the
   three clamp sites and the post-loop report.
6. CHK_OBS_FLUX guards before the OBSFLUX/OBSFRAME writes and the
   OBS_FREQ strictly-decreasing check after INS_LINE_OBS_V4.
7. COMP_JCONT_V4.INC FG dispatch (10-Nov-2026): FG_J_CMF_NONMONO_V1
   is gated on DO_NONMONO_KNOP alone -- a DO_NONMONO_SOB population
   was converged on the clamped-SIGMA V13 formal J, so its OBS render
   must use the clamped V12, not the Knop solver (cold I_NEXT on the
   fine OBS grid distorted the spectrum).
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUB = (ROOT / "obs/cmf_flux_sub_v5.f").read_text()
INC = (ROOT / "obs/COMP_JCONT_V4.INC").read_text()
OBSF = (ROOT / "obs/obs_frame_sub_v9.f").read_text()
RDC = (ROOT / "obs/rd_cmf_flux_controls.f").read_text()
MODF = (ROOT / "obs/cmf_flux_cntrl_var_mod.f").read_text()


def section(source, start, end):
    first = source.index(start)
    last = source.index(end, first)
    return source[first:last]


def check(name, passed):
    if not passed:
        raise AssertionError(name)
    print(f"[PASS] {name}")


def test_inc_grid_rejection():
    """ACCURATE=T must abort on a genuinely non-monotonic field."""
    body = section(
        SUB,
        "IF(ACCURATE .AND. N_NONMONO_OBS .GT. 0)THEN",
        "END IF",
    )
    check(
        "INC_GRID rejection gated on the actual field (N_NONMONO_OBS)",
        "N_NONMONO_OBS .GT. 0" in body,
    )
    check(
        "INC_GRID rejection aborts",
        "STOP" in body,
    )
    check(
        "rejection diagnostic names EXTEND_VTSIGMA log(V) hazard",
        "EXTEND_VTSIGMA" in body and "log(V)" in body,
    )
    # The check must sit AFTER the detector that fills N_NONMONO_OBS,
    # i.e. after the detector's closing block, before the frequency
    # grid is built.
    det_end = SUB.index("INC_GRID=T (ACCURATE) is not supported")
    det_pos = SUB.index("N_NONMONO_OBS=N_NONMONO")
    check(
        "INC_GRID check runs after the velocity detector",
        det_end > det_pos,
    )


def test_moment_dispatch():
    """SOB+LAM_ES=F routes to MOM_J_CMF_NONMONO_V1; stock keeps V8."""
    body = section(
        INC,
        "ELSE IF(USE_LAM_ES)THEN",
        "CALL TUNE(ITWO,'MOM_J_CMF')",
    )
    check(
        "nonmono dispatch gated on DO_NONMONO_SOB",
        "IF(DO_NONMONO_SOB)THEN" in body,
    )
    check(
        "SET_MOM_NU_NEXT called before the solver",
        body.index("CALL SET_MOM_NU_NEXT")
        < body.index("CALL MOM_J_CMF_NONMONO_V1"),
    )
    call = section(
        body, "CALL MOM_J_CMF_NONMONO_V1", "\n\t       ELSE"
    )
    check(
        "SIGMA_TRUE passed to the nonmono moment solver",
        "SIGMA,SIGMA_TRUE" in call or "SIGMA_TRUE" in call,
    )
    check(
        "boundary controls match the population call",
        "INNER_BND_METH" in call
        and "IB_STAB_FACTOR" in call
        and "H_CHK_OPTION" in call,
    )
    check(
        "FREQ_INDX and NCF passed (lagged store indexing)",
        "FREQ_INDX" in call and "NCF" in call,
    )
    v8 = section(body, "ELSE\n", "END IF")
    check(
        "stock V8 path preserved when DO_NONMONO_SOB=F",
        "CALL MOM_J_CMF_V8" in v8,
    )


def test_first_freq_dlog_nu():
    """The nonmono moment solver needs a valid dLOG_NU at ML=1."""
    # FG_COUNT=0.0_LDP appears in both the extended and non-extended
    # FIRST_FREQ blocks; the fix lives in the non-extended (last) one.
    start = INC.rindex("FG_COUNT=0.0_LDP")
    body = INC[start:INC.index("dLOG_NU=LOG(FL_OLD/FL)", start)]
    check(
        "one-sided dLOG_NU gated on SOB+LAM_ES=F",
        "DO_NONMONO_SOB .AND. .NOT. USE_LAM_ES" in body,
    )
    check(
        "forward spacing ln(nu_1/nu_2) used",
        "LOG(FL/NU(MIN(ML+1,NCF)))" in body,
    )
    check(
        "guard requires a positive, smaller neighbour frequency",
        "NU(MIN(ML+1,NCF)) .GT. 0.0_LDP" in body
        and "NU(MIN(ML+1,NCF)) .LT. FL" in body,
    )


def test_freq_match_mirror():
    """Obs reader mirrors the population FREQ_MATCH enforcement."""
    check(
        "DO_NONMONO_FREQ_MATCH declared in obs control module",
        "LOGICAL :: DO_NONMONO_FREQ_MATCH" in MODF,
    )
    body = section(
        RDC,
        "DO_NONMONO_FREQ_MATCH=.FALSE.",
        "CALL SET_MOM_FREQ_MATCH",
    )
    check(
        "auto-set applies under KNOP or SOB+LAM_ES=F",
        "DO_NONMONO_KNOP .OR." in body
        and "DO_NONMONO_SOB .AND. .NOT. USE_LAM_ES" in body,
    )
    check(
        "explicit F overridden with a banner (not silently)",
        "OBS override: explicit DO_NONMONO_FREQ_MATCH=F -> T" in body,
    )
    check(
        "absent key auto-set message",
        "OBS auto-set: DO_NONMONO_FREQ_MATCH=T" in body,
    )
    check(
        "both solver modules receive the flag",
        "CALL SET_FG_FREQ_MATCH(DO_NONMONO_FREQ_MATCH)" in RDC
        and "CALL SET_MOM_FREQ_MATCH(DO_NONMONO_FREQ_MATCH)" in RDC,
    )


def test_edge_clamp_census():
    """All three clamp sites update the census; report is printed."""
    check(
        "census counters declared",
        "N_EDGE_CLAMP" in OBSF and "CLAMP_GRID_MAX" in OBSF,
    )
    check(
        "counters zeroed before the LS loop",
        "N_EDGE_CLAMP=0" in OBSF,
    )
    check(
        "all three clamp sites update the census",
        OBSF.count("N_EDGE_CLAMP=N_EDGE_CLAMP+1") == 3,
    )
    check(
        "updates are OpenMP-safe",
        OBSF.count("$OMP CRITICAL (CLAMP_CENSUS)") == 3,
    )
    check(
        "window and full-grid overshoot tracked separately",
        "CLAMP_FRAC_MAX=MAX(CLAMP_FRAC_MAX,T3)" in OBSF
        and "T4 .GT. CLAMP_GRID_MAX" in OBSF,
    )
    report = section(OBSF, "IF(N_EDGE_CLAMP .GT. 0)THEN", "END IF")
    check(
        "report prints count and worst location",
        "N_EDGE_CLAMP" in report and "CLAMP_I_WORST" in report,
    )
    check(
        "out-of-grid clamp raises a warning to widen MIN_CF/MAX_CF",
        "CLAMP_GRID_MAX .GT. 0.0_LDP" in report
        and "MIN_CF/MAX_CF" in report,
    )


def test_output_checks():
    """OBS_FREQ monotonicity and OBS_FLUX finiteness are enforced."""
    body = section(
        SUB,
        "Finished calling INS_LINE_OBS_V4",
        "CALL TUNE(1,'SET_PROF_STORE')",
    )
    check(
        "OBS_FREQ strictly-decreasing check aborts",
        "OBS_FREQ(ML) .GE. OBS_FREQ(ML-1)" in body and "STOP" in body,
    )
    check(
        "OBS_FREQ finiteness check aborts",
        "OBS_FREQ not finite" in body,
    )
    check(
        "CHK_OBS_FLUX called before OBSFLUX write",
        "CALL CHK_OBS_FLUX(OBS_FLUX,N_OBS,'OBSFLUX',LUER)" in SUB,
    )
    check(
        "CHK_OBS_FLUX called before both OBSFRAME writes",
        "CALL CHK_OBS_FLUX(OBS_FLUX,N_OBS,'OBSFRAME phase 1',LUER)"
        in SUB
        and "CALL CHK_OBS_FLUX(OBS_FLUX,N_OBS,'OBSFRAME phase 2',LUER)"
        in SUB,
    )
    chk = section(SUB, "SUBROUTINE CHK_OBS_FLUX", "RETURN\n\tEND")
    check(
        "NaN/Inf flux aborts",
        "OBS_FLUX(I) .NE. OBS_FLUX(I)" in chk and "STOP" in chk,
    )
    check(
        "negative flux censused not repaired",
        "N_NEG=N_NEG+1" in chk and "OBS_FLUX(I)=ABS" not in chk,
    )


def test_sob_phase_guards():
    """Two-phase structure invariants remain in place."""
    check(
        "phase 2 forces coherent ES and fresh EDDFAC",
        "IF(SOB_PHASE .EQ. 2)COHERENT_ES=.TRUE." in SUB
        and "IF(SOB_PHASE .EQ. 2)COMPUTE_EDDFAC=.TRUE." in SUB,
    )
    check(
        "phase-2 reload branch preserves deposits",
        "ES_COUNTER .GT. 1 .OR. SOB_PHASE .GT. 1" in SUB,
    )
    check(
        "OBSFLUX written in phase 1 only",
        "IF(SOB_PHASE .EQ. 1)THEN\n\tCALL CHK_OBS_FLUX(OBS_FLUX,N_OBS,"
        "'OBSFLUX',LUER)" in SUB,
    )
    check(
        "ACCESS_F restarted per phase",
        "ACCESS_F=5" in SUB,
    )
    check(
        "maser-drop diagnostic retained",
        "SOB_MASER_DROP" in SUB
        and "SOB_DROP_EMIS" in SUB
        and "SOB_INJ_EMIS" in SUB,
    )


def test_fg_dispatch_knop_only():
    """FG_J_CMF_NONMONO_V1 runs only under DO_NONMONO_KNOP."""
    # Anchor on the dispatch call itself: an earlier guard block also
    # opens with IF(DO_NONMONO_KNOP)THEN (the INC_GRID rejection).
    call_pos = INC.index("'Calling FG_J_CMF_NONMONO_V1")
    gate_pos = INC.rindex("IF(DO_NONMONO_KNOP)THEN", 0, call_pos)
    v12_pos = INC.index("CALL FG_J_CMF_V12", call_pos)
    body = INC[gate_pos:v12_pos + len("CALL FG_J_CMF_V12")]
    check(
        "Knop formal solver called inside the KNOP gate",
        "CALL FG_J_CMF_NONMONO_V1" in body,
    )
    gate_line = body.splitlines()[0]
    check(
        "gate does not include DO_NONMONO_SOB",
        "DO_NONMONO_SOB" not in gate_line,
    )
    check(
        "SOB-alone falls through to clamped V12",
        "CALL FG_J_CMF_V12" in body,
    )
    # The SET_FG_GS_OBS_GRID registration and the OBS refresh gates
    # must match the same routing or the lagged store is never used.
    check(
        "OBS grid registration gated on KNOP only",
        "IF(DO_NONMONO_KNOP)THEN\n\t  CALL SET_FG_GS_OBS_GRID"
        in SUB,
    )
    check(
        "FG refresh gated on KNOP only",
        "IF(DO_NONMONO_KNOP .AND.\n\t1    ES_COUNTER"
        in SUB,
    )
    check(
        "MOM refresh gated on SOB + not USE_LAM_ES",
        "IF(DO_NONMONO_SOB .AND. .NOT. USE_LAM_ES .AND." in SUB,
    )


def main():
    test_inc_grid_rejection()
    test_moment_dispatch()
    test_first_freq_dlog_nu()
    test_freq_match_mirror()
    test_edge_clamp_census()
    test_output_checks()
    test_sob_phase_guards()
    test_fg_dispatch_knop_only()
    print("all obs nonmono dispatch checks passed")


if __name__ == "__main__":
    main()
