! drv_knop_spec_gate.f - production-object harness for
! KNOP_SPEC_GATE_EVAL (obs/knop_spec_gate.f, P8.3).  The python
! wrapper stages GS_STATE / .CMFGEN_KNOP_STALLED / manifest.json
! fixtures in a scratch model dir and runs this binary from its
! obs/ child, exactly where batobs.sh puts cmf_flux.exe.
! Argument 'sob' exercises the non-KNOP gate-off path.
!
      PROGRAM DRV_KNOP_SPEC_GATE
      USE CMF_FLUX_CNTRL_VAR_MOD
      IMPLICIT NONE
      INTEGER LU,NARG
      CHARACTER(LEN=8) ARG
      LU=6
      NARG=COMMAND_ARGUMENT_COUNT()
      DO_NONMONO_KNOP=.TRUE.
      IF(NARG .GE. 1)THEN
        CALL GET_COMMAND_ARGUMENT(1,ARG)
        IF(ARG .EQ. 'sob')DO_NONMONO_KNOP=.FALSE.
      END IF
      CALL KNOP_SPEC_GATE_EVAL(LU)
      WRITE(LU,'(A,I0,1X,A,1X,L1)')'RESULT ',KNOP_SPEC_STATE,
     1    TRIM(KNOP_SPEC_CLASS),KNOP_SPEC_ALLOW
      END
