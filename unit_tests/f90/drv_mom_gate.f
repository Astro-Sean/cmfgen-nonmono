! Standalone driver for the phase-2 FIX 1 / FIX 2 units in
! plane/mom_j_cmf_nonmono_v1.f:
!
!   MOM_SAME_DIR_RHO -- same-direction residual-ratio estimator for
!     the adaptive JNEXT_RELAX update (audit finding F). The driver
!     fills RES_PASS(2,NP) with manufactured per-direction residual
!     sequences and prints RHO / NR / LOK.
!
!   MOM_RF_COMMIT_GATE -- the single post-AA acceptance gate (audit
!     finding G). The driver builds small manufactured
!     seed / pre-AA / candidate fields and prints KEEP plus the gate
!     warn stream, which lists each failing check.
!
! Usage: drv_mom_gate <case>
!
      PROGRAM DRV_MOM_GATE
      USE SET_KIND_MODULE
      IMPLICIT NONE
!
      INTEGER ND,NCF,NP
      PARAMETER (ND=8,NCF=4,NP=12)
!
      CHARACTER(LEN=32) CASE
      INTEGER I,L,K
      REAL(KIND=LDP) RES_PASS(2,NP)
      REAL(KIND=LDP) JF(ND,NCF),HF(ND,NCF)
      REAL(KIND=LDP) JS(ND,NCF),HS(ND,NCF)
      REAL(KIND=LDP) JP(ND,NCF),HP(ND,NCF)
      REAL(KIND=LDP) RHO,X,XN,R1,R2
      INTEGER NR
      LOGICAL LOK,KEEP
!
      CALL GETARG(1,CASE)
      IF(LEN_TRIM(CASE) .EQ. 0)STOP 'DRV_MOM_GATE: need <case>'
!
      RES_PASS=0.0_LDP
      RHO=0.0_LDP
      NR=0
      LOK=.FALSE.
!
      IF(CASE .EQ. 'rhoalt' .OR. CASE .EQ. 'rhosingle' .OR.
     1   CASE .EQ. 'rhodiv')THEN
!
! Scalar contraction maps with a shared fixed point: pass k applies
! map rho_d to the iterate, the recorded residual is the iterate
! step |x_k - x_{k-1}|. rhoalt alternates expanding rho1 and
! compressing rho2, so every same-direction pair spans a full
! two-pass cycle and the estimator must return the composite
! factor rho1*rho2. rhosingle runs one direction only, rhodiv
! diverges to check the routine reports rho>1 unclamped.
!
        R1=0.5_LDP;  R2=0.9_LDP
        IF(CASE .EQ. 'rhosingle')R2=-1.0_LDP
        IF(CASE .EQ. 'rhodiv')THEN
          R1=1.3_LDP;  R2=1.3_LDP
        END IF
        X=1.0_LDP
        DO K=1,NP
          IF(R2 .LT. 0.0_LDP)THEN
            XN=R1*X
            RES_PASS(1,K)=ABS(XN-X)
          ELSE IF(MOD(K,2) .EQ. 1)THEN
            XN=R1*X
            RES_PASS(1,K)=ABS(XN-X)
          ELSE
            XN=R2*X
            RES_PASS(2,K)=ABS(XN-X)
          END IF
          X=XN
        END DO
        CALL MOM_SAME_DIR_RHO(RES_PASS,NP,RHO,NR,LOK)
        WRITE(*,'(A,ES16.8E3)')'RHO=',RHO
        WRITE(*,'(A,I4)')'NRAT=',NR
        WRITE(*,'(A,L2)')'LOK=',LOK
!
      ELSE IF(CASE .EQ. 'rhobeat')THEN
!
! Beat case: each direction's own residual sequence contracts at
! 0.5 per pass, but the two sequences sit orders apart, so every
! stride-1 cross-map ratio alternates between ~1E3 and ~5E-4 --
! the old pass-ordered estimator read this as divergence/garbage.
!
        DO K=1,NP
          IF(MOD(K,2) .EQ. 1)THEN
            RES_PASS(1,K)=1.0_LDP*0.5_LDP**((K+1)/2-1)
          ELSE
            RES_PASS(2,K)=1.0E+3_LDP*0.5_LDP**(K/2-1)
          END IF
        END DO
        CALL MOM_SAME_DIR_RHO(RES_PASS,NP,RHO,NR,LOK)
        WRITE(*,'(A,ES16.8E3)')'RHO=',RHO
        WRITE(*,'(A,I4)')'NRAT=',NR
        WRITE(*,'(A,L2)')'LOK=',LOK
!
      ELSE IF(CASE .EQ. 'rhoskip')THEN
!
! Dropped direction: expanding passes run twice then stall out
! (SKIP_D), compressing passes keep going. Slot 1 yields one
! ratio, slot 2 two; the pooled median must read 0.8.
!
        RES_PASS(1,1)=1.0_LDP
        RES_PASS(1,3)=0.5_LDP
        RES_PASS(2,2)=0.8_LDP
        RES_PASS(2,4)=0.64_LDP
        RES_PASS(2,5)=0.512_LDP
        CALL MOM_SAME_DIR_RHO(RES_PASS,NP,RHO,NR,LOK)
        WRITE(*,'(A,ES16.8E3)')'RHO=',RHO
        WRITE(*,'(A,I4)')'NRAT=',NR
        WRITE(*,'(A,L2)')'LOK=',LOK
!
      ELSE IF(CASE .EQ. 'rhofew')THEN
!
! Under-sampled history: fewer than 3 ratios -> LOK must come back
! false so the caller falls back to the cross-call estimate.
!
        RES_PASS(1,1)=1.0_LDP
        RES_PASS(1,2)=0.7_LDP
        CALL MOM_SAME_DIR_RHO(RES_PASS,NP,RHO,NR,LOK)
        WRITE(*,'(A,ES16.8E3)')'RHO=',RHO
        WRITE(*,'(A,I4)')'NRAT=',NR
        WRITE(*,'(A,L2)')'LOK=',LOK
!
      ELSE
!
! Gate cases. Manufactured fields on the small grid: a smooth
! positive seed JS and signed HS, a pre-AA field JP/HP that passed
! the guards, and a candidate JF/HF modified to violate one check.
! JP is also the restore target, so each case prints KEEP plus the
! gate's own warn line naming the failing check.
!
        DO L=1,NCF
          DO I=1,ND
            JS(I,L)=1.0_LDP+0.1_LDP*DBLE(I)
            HS(I,L)=0.02_LDP*(DBLE(I)-4.5_LDP)
            JP(I,L)=JS(I,L)
            HP(I,L)=HS(I,L)
            JF(I,L)=JP(I,L)
            HF(I,L)=HP(I,L)
          END DO
        END DO
!
        IF(CASE .EQ. 'clean')THEN
!
! Nothing past the guards moved: candidate identical to pre-AA,
! which equals the seed. All checks must pass.
!
        ELSE IF(CASE .EQ. 'aafix')THEN
!
! Legitimate AA step: the pre-AA field regressed a little vs the
! seed (a few cells shifted >5%), the candidate holds that field
! unchanged. A field identical to its restore target must always
! pass -- the gate only bites writes past the snapshot.
!
          JP(3,1)=JS(3,1)*1.2_LDP
          JP(5,2)=JS(5,2)*0.8_LDP
          HP(2,3)=HS(2,3)*1.3_LDP
          JF(1:ND,1:NCF)=JP(1:ND,1:NCF)
          HF(1:ND,1:NCF)=HP(1:ND,1:NCF)
        ELSE IF(CASE .EQ. 'aaimp')THEN
!
! AA improves on the pre-AA field: pre-AA regressed vs the seed in
! three >50% cells; the candidate repairs two of them, so the
! worst-cell count strictly improves and the gate must commit.
!
          JP(2,1)=JS(2,1)*1.7_LDP
          JP(4,2)=JS(4,2)*0.4_LDP
          JP(6,3)=JS(6,3)*1.8_LDP
          JF(1:ND,1:NCF)=JP(1:ND,1:NCF)
          HF(1:ND,1:NCF)=HP(1:ND,1:NCF)
          JF(2,1)=JS(2,1)*1.05_LDP
          JF(4,2)=JS(4,2)*0.95_LDP
        ELSE IF(CASE .EQ. 'neg')THEN
!
! AA write drives a positive cell non-positive.
!
          JF(4,2)=-0.5_LDP
        ELSE IF(CASE .EQ. 'nonfin')THEN
!
! Poison entry past the 1E40 bound.
!
          JF(3,1)=1.0E+41_LDP
        ELSE IF(CASE .EQ. 'jump')THEN
!
! One cell committed more than 50x over the retained store.
!
          JF(5,1)=JS(5,1)+51.0_LDP*
     1      MAX(ABS(JS(5,1)),MAXVAL(ABS(JS))/50.0_LDP)
        ELSE IF(CASE .EQ. 'bulk')THEN
!
! Bulk churn past the 1.5x bound: pre-AA moved two cells >5%, the
! candidate stirs eight without moving any past 50%.
!
          JP(2,1)=JS(2,1)*1.2_LDP
          JP(6,3)=JS(6,3)*1.2_LDP
          JF(1:ND,1:NCF)=JP(1:ND,1:NCF)
          HF(1:ND,1:NCF)=HP(1:ND,1:NCF)
          DO I=1,ND
            JF(I,4)=JS(I,4)*1.1_LDP
          END DO
          JF(3,1)=JS(3,1)*1.1_LDP
          JF(7,2)=JS(7,2)*1.1_LDP
        ELSE IF(CASE .EQ. 'seedgate')THEN
!
! Reverted-store semantics: the pre-AA field IS the seed (a
! keep_prior refresh), so its n5/n50 counts are zero and any >5%
! move the candidate makes is a regression.
!
          JF(4,2)=JS(4,2)*1.1_LDP
        ELSE
          STOP 'DRV_MOM_GATE: unknown case'
        END IF
!
        CALL MOM_RF_COMMIT_GATE(JF,HF,JS,HS,JP,HP,ND,NCF,6,KEEP)
        WRITE(*,'(A,L2)')'KEEP=',KEEP
!
      END IF
!
      STOP
      END
