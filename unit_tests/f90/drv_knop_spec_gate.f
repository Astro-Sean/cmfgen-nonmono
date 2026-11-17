! drv_knop_spec_gate.f - production-object harness for
! KNOP_SPEC_GATE_EVAL (obs/knop_spec_gate.f, P8.3).  The python
! wrapper stages GS_STATE / .CMFGEN_KNOP_STALLED / manifest.json
! fixtures in a scratch model dir and runs this binary from its
! obs/ child, exactly where batobs.sh puts cmf_flux.exe.
!
      PROGRAM DRV_KNOP_SPEC_GATE
      USE CMF_FLUX_CNTRL_VAR_MOD
      IMPLICIT NONE
      INTEGER LU
      LU=6
      DO_NONMONO_KNOP=.TRUE.
      CALL KNOP_SPEC_GATE_EVAL(LU)
      WRITE(LU,'(A,I0,1X,A,1X,L1)')'RESULT ',KNOP_SPEC_STATE,
     1    TRIM(KNOP_SPEC_CLASS),KNOP_SPEC_ALLOW
      END
