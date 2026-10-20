! Driver for the DUMP_MOM_JNEXT_GS sweep report: fills the module
! counters with distinctive values and calls the production routine
! so the test can verify the 'matrix'/'eq resid' rows print the
! real arguments on numeric descriptors (the pre-fix formats left
! one descriptor too few, so MOM_MIN_DRA/MOM_MAX_COH landed on A
! descriptors and printed garbage).
      PROGRAM DRV_MOM_SWEEP
      USE SET_KIND_MODULE
      USE MOD_MOM_J_NONMONO_V1
      USE CONTROL_VARIABLE_MOD, ONLY: CONV_DIAG_LVL
      IMPLICIT NONE
      REAL(KIND=LDP) R(4),V(4),SIGMA(4)
!
      ALLOCATE(J_NEXT_FRQ(4,2))
      J_NEXT_FRQ=0.0_LDP
      R=1.0_LDP
      V=1.0_LDP
      SIGMA=0.5_LDP
!
! CONV_DIAG_LVL=2 turns on the extended report rows (per-pass,
! bulk-share, argmax-history, boundary tag) so this run exercises
! their formats.
!
      CONV_DIAG_LVL=2
      MOM_GS_NTOT=100_8
      MOM_GS_NCF=2
      MOM_SW_NCROW=60_8
      MOM_SW_NROW=100_8
      MOM_GS_CHG_P1=8.0E-2_LDP
      MOM_GS_CHG_P2=3.0E-2_LDP
      MOM_MIN_DR=1.1E-3_LDP
      MOM_MIN_DR_I=2
      MOM_MIN_DR_L=3
      MOM_MIN_DRA=4.5E-6_LDP
      MOM_MIN_DRA_I=7
      MOM_MIN_DRA_L=8
      MOM_EQ_RES=9.9E-2_LDP
      MOM_EQ_RES_I=11
      MOM_EQ_RES_L=12
      MOM_MAX_COH=1.3E+1_LDP
      MOM_RF_DIR='AUTO'
      CALL DUMP_MOM_JNEXT_GS(0,R,V,SIGMA,4,97)
      END
