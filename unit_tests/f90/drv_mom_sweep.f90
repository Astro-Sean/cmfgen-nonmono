! Driver for the DUMP_MOM_JNEXT_GS sweep report: fills the module
! counters with distinctive values and calls the production routine
! so the test can verify the 'matrix'/'eq resid' rows print the
! real arguments on numeric descriptors (the pre-fix formats left
! one descriptor too few, so MOM_MIN_DRA/MOM_MAX_COH landed on A
! descriptors and printed garbage).
      PROGRAM DRV_MOM_SWEEP
      USE SET_KIND_MODULE
      USE MOD_MOM_J_NONMONO_V1
      USE FG_J_CMF_NONMONO_MOD_V1
      USE CONTROL_VARIABLE_MOD, ONLY: CONV_DIAG_LVL
      IMPLICIT NONE
      REAL(KIND=LDP) R(4),V(4),SIGMA(4)
!
      ALLOCATE(J_NEXT_FRQ(8,2))
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
      MOM_GS_NLOW=2_8
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
!
! Call 1: interior argmax with alternating last-4 values -> 'cycle'.
!
      MOM_ARG_HI=(/9,7,5,5/)
      MOM_ARG_HL=(/3,1,2,2/)
      MOM_ARG_HV=(/0.10_LDP,0.05_LDP,0.10_LDP,0.05_LDP/)
      MOM_GS_ARGI=5
      MOM_GS_ARGL=2
      MOM_GS_MAXCHG=8.0E-2_LDP
      CALL DUMP_MOM_JNEXT_GS(0,R,V,SIGMA,4,97)
!
! Calls 2-3: argmax pinned on a boundary row; the stable position
! classifies 'bndry' and tags the drift line bnd_row.
!
      MOM_GS_ARGI=1
      MOM_GS_ARGL=4
      MOM_GS_MAXCHG=0.15_LDP
      CALL DUMP_MOM_JNEXT_GS(0,R,V,SIGMA,4,97)
      MOM_GS_MAXCHG=0.20_LDP
      CALL DUMP_MOM_JNEXT_GS(0,R,V,SIGMA,4,97)
!
! FG variant: the same alternating seed on the (depth,ray,freq)
! history classifies 'cycle'.  N_COMP=0 keeps the checkpoint write
! from streaming unallocated intensity arrays.
!
      FG_ARG_HK=(/3,3,2,2/)
      FG_ARG_HR=(/4,5,1,1/)
      FG_ARG_HL=(/6,1,2,2/)
      FG_ARG_HV=(/0.10_LDP,0.05_LDP,0.10_LDP,0.05_LDP/)
      FG_GS_ARGK=2
      FG_GS_ARGR=1
      FG_GS_ARGL=2
      FG_GS_MAXCHG=8.0E-2_LDP
      N_COMP=0
      ALLOCATE(I_P_DEPTH_NEXT(0,1,1))
      ALLOCATE(I_M_DEPTH_NEXT(0,1,1))
      CALL DUMP_FG_INEXT_GS(0,R,V,SIGMA,4,98)
      END
