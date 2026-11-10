! Driver for FG_J_CMF_V12 (plane/fg_j_cmf_v12.f).
!
! Alterations under test (bug fixes, not flag-gated):
!   1. I_CORE floor at the inner boundary: under DIF the boundary
!      intensity (ETA_RAY + DBB*mu_z)/CHI_RAY can go negative where the
!      diffusion term |DBB| exceeds ETA_RAY on an optically thin base.
!      Prod floors it at 1e-30; upstream seeds the outward sweep with the
!      negative value.
!   2. Negative I_M/I_P sweep values are floored to 1e-30 and counted
!      (diagnostic 'FG_J_CMF_V12: N negative I_X repaired on ray LS').
!   3. JNU_STORE negatives below the underflow noise floor are repaired
!      in place before FEDD=K/J divides.
!
! The driver builds a small expanding atmosphere (ND=12, NC=6, NP=18),
! the IMPAR impact parameters and the NORDANG/GENANG quadrature weights,
! then calls the solver once with INIT=T, NEW_FREQ=T, DIF=T on the
! INTEGRAL path (where the guards live).  Case 'neg' passes DBB=-1e6 so
! I_CORE is deeply negative; 'ctl' passes DBB=0 for an identical-path
! control.
!
      PROGRAM DRV_FG12
      USE SET_KIND_MODULE
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: ND=12, NC=6, NP=ND+NC
      REAL(KIND=LDP) :: ETA(ND),CHI(ND),ESEC(ND),V(ND),SIGMA(ND),R(ND)
      REAL(KIND=LDP) :: P(NP)
      REAL(KIND=LDP) :: JNU(ND),FEDD(ND),IPLUS_P(NP)
      REAL(KIND=LDP) :: JQW(ND,NP),HQW(ND,NP),KQW(ND,NP),NQW(ND,NP)
      REAL(KIND=LDP) :: HMIDQW(ND,NP),NMIDQW(ND,NP)
      REAL(KIND=LDP) :: VDOP_VEC(ND)
      REAL(KIND=LDP) :: RP,DELMU,DBB,IC,FREQ,dLOG_NU
      REAL(KIND=LDP) :: REXT_FAC,VDOP_FRAC
      REAL(KIND=LDP) :: RETURNED_IN_HBC,RETURNED_OUT_HBC
      CHARACTER(LEN=32) :: CASE
      INTEGER :: I
      LOGICAL :: MID
      EXTERNAL JWEIGHT_V2,HWEIGHT_V2,KWEIGHT_V2,NWEIGHT_V2
!
      CALL GETARG(1,CASE)
!
      DO I=1,ND
        R(I)=20.0_LDP-18.0_LDP*(I-1)/REAL(ND-1,KIND=LDP)
        V(I)=2000.0_LDP-1800.0_LDP*(I-1)/REAL(ND-1,KIND=LDP)
        SIGMA(I)=-0.5_LDP
        ETA(I)=1.0E-6_LDP
        CHI(I)=1.0E-8_LDP
        ESEC(I)=1.0E-10_LDP
        VDOP_VEC(I)=30.0_LDP
        JNU(I)=1.0E-8_LDP
        FEDD(I)=1.0_LDP/3.0_LDP
      END DO
!
! jump puts a decade-scale step in the source function S=ETA/CHI in
! the innermost interval.  The -dS correction terms in the sweep
! recursion then outweigh S itself where DTAU is small, so I_M and
! I_P can undershoot below zero: the branch the fork's repair guards
! cover.  A smooth ETA never gets |dS| above ~S per interval.
!
      IF(CASE .EQ. 'jump')THEN
        CHI(:)=1.0E-4_LDP
        ETA(ND)=1.0E-2_LDP
        ETA(ND-1)=1.0E-7_LDP
      ELSE IF(CASE .EQ. 'wall')THEN
        CHI(ND-1:ND)=1.0E6_LDP
        ETA(ND-1:ND)=1.0E-12_LDP
      ELSE IF(CASE .EQ. 'thin')THEN
        ETA(:)=1.0E-30_LDP
        CHI(:)=1.0E-12_LDP
        ESEC(:)=1.0E-30_LDP
      END IF
      IPLUS_P(:)=0.0_LDP
!
! IMPAR layout: core rays equally spaced in mu inside R(ND), tangent
! rays on the remaining grid radii.
!
      RP=R(ND)
      P(1)=0.0_LDP
      DELMU=1.0_LDP/NC
      DO I=2,NC
        P(I)=RP*SQRT(1.0_LDP-(DELMU*(NC-I+1))**2)
      END DO
      DO I=NC+1,NP
        P(I)=R(ND-I+NC+1)
      END DO
!
      CALL NORDANGQW_V2(JQW,R,P,NC,ND,NP,JWEIGHT_V2)
      CALL NORDANGQW_V2(HQW,R,P,NC,ND,NP,HWEIGHT_V2)
      CALL NORDANGQW_V2(KQW,R,P,NC,ND,NP,KWEIGHT_V2)
      CALL NORDANGQW_V2(NQW,R,P,NC,ND,NP,NWEIGHT_V2)
      MID=.TRUE.
      CALL GENANGQW_V2(HMIDQW,R,P,NC,ND,NP,HWEIGHT_V2,MID)
      CALL GENANGQW_V2(NMIDQW,R,P,NC,ND,NP,NWEIGHT_V2,MID)
!
      IC=1.0_LDP
      FREQ=1.0_LDP
      dLOG_NU=1.0E-6_LDP
      REXT_FAC=1.0_LDP
      VDOP_FRAC=0.5_LDP
      DBB=0.0_LDP
      IF(CASE .EQ. 'neg')DBB=-1.0E6_LDP
!
      CALL FG_J_CMF_V12(ETA,CHI,ESEC,V,SIGMA,R,P, &
         JNU,FEDD,JQW,HQW,KQW,NQW,HMIDQW,NMIDQW, &
         RETURNED_IN_HBC,RETURNED_OUT_HBC,IPLUS_P, &
         FREQ,dLOG_NU,.TRUE.,DBB,IC, &
         VDOP_VEC,VDOP_FRAC,REXT_FAC, &
         'ZERO','INTEGRAL',.FALSE.,.TRUE.,.TRUE.,NC,NP,ND)
!
      WRITE(6,'(A,*(ES14.5))')'JNU=',JNU
      WRITE(6,'(A,*(ES14.5))')'FEDD=',FEDD
      WRITE(6,'(A,2ES14.5)')'HBC=',RETURNED_IN_HBC,RETURNED_OUT_HBC
      WRITE(6,*)'DONE_FG12'
!
      END
