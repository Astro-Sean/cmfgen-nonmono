! Driver for MOM_NU_CHAIN_SOLVE (plane/mom_j_cmf_nonmono_v1.f,
! [MOM_NUSOLV]). Seeds the NU_* capture arrays and J_NEXT_FRQ with
! manufactured row equations whose exact solutions are closed-form,
! calls the chain solve directly, and dumps the solved field. The
! python side recomputes the analytic answers.
!
! Modes:
!   exp    : all-expanding chain (kup=k-1), geometric solution
!   cmp    : all-compressing chain (kup=k+1)
!   ring2  : one expanding/compressing 2-cycle between boundaries
!   ring4  : stride-3 compressing edge closing a 4-member ring
!            (DO_NONMONO_FREQ_MATCH strides >1 make these)
!   mix    : chains feeding two 2-cycles
!   depth  : nonzero TA/TC depth terms (lagged reads)
!   fresh  : NU_FRESH=F columns frozen, fresh columns solved
!   sing   : ring with cyclic product A=1 -> frozen (incompatible)
!   neg    : chain solving to negative values -> rejected
!   nonfin : TB=0 and nonfinite RHS rows -> frozen
!
! Usage: drv_mom_nu <mode> <tag>
!
      PROGRAM DRV_MOM_NU
      USE SET_KIND_MODULE
      USE MOD_MOM_J_NONMONO_V1
      IMPLICIT NONE
!
      INTEGER ND,NCF
      PARAMETER (ND=4,NCF=10)
!
      INTEGER I,K,IOS
      REAL(KIND=LDP) CHG,T3,AV0,RV0
      INTEGER*8 N5,N50
      CHARACTER(LEN=32) MODE,OUTTAG
      REAL(KIND=LDP), ALLOCATABLE :: J0(:,:)
!
      CALL GETARG(1,MODE)
      CALL GETARG(2,OUTTAG)
      IF(LEN_TRIM(OUTTAG) .EQ. 0)OUTTAG='nu'
!
! No damping: the checks want the exact chain solution, not a relaxed
! step toward it.
!
      JNEXT_RELAX=1.0_LDP
      MOM_KN_REFRESHING=.FALSE.
!
      ND_NM=ND
      ALLOCATE(J_NEXT_FRQ(ND,NCF),NU_TA(ND,NCF),NU_TB(ND,NCF),
     1    NU_TC(ND,NCF),NU_X0(ND,NCF),NU_CUP(ND,NCF),
     2    NU_KUP(ND,NCF),NU_FRESH(NCF),STAT=IOS)
      IF(IOS .NE. 0)STOP 'DRV_MOM_NU: alloc failed'
      J_NEXT_FRQ=1.0_LDP
      NU_TA=0.0_LDP; NU_TC=0.0_LDP
      NU_TB=1.0_LDP; NU_X0=0.5_LDP
      NU_CUP=0.0_LDP; NU_KUP=0
      NU_FRESH(1:NCF)=.TRUE.
      AV0=0.5_LDP
      RV0=0.5_LDP
!
      SELECT CASE(MODE)
      CASE('exp')
! All expanding: x(k) = r + a*x(k-1), x(1)=r -> geometric.
        DO K=2,NCF
          DO I=1,ND-1
            NU_CUP(I,K)=AV0
            NU_KUP(I,K)=K-1
          END DO
        END DO
        NU_X0(1:ND-1,1:NCF)=RV0
      CASE('cmp')
! All compressing: x(k) = r + a*x(k+1), x(NCF)=r.
        DO K=1,NCF-1
          DO I=1,ND-1
            NU_CUP(I,K)=AV0
            NU_KUP(I,K)=K+1
          END DO
        END DO
        NU_X0(1:ND-1,1:NCF)=RV0
      CASE('ring2')
! Nodes 2,3 form a 2-cycle (2 compressing reads 3; 3 expanding
! reads 2); nodes 1 and 4+ are uncoupled boundaries.
        NU_X0(1:ND-1,1:NCF)=RV0
        DO I=1,ND-1
          NU_CUP(I,2)=AV0; NU_KUP(I,2)=3
          NU_CUP(I,3)=0.3_LDP; NU_KUP(I,3)=2
        END DO
      CASE('ring4')
! Stride-3 compressing edge at node 1 closes a 4-member ring:
! 1->4 (cmp), 4->3, 3->2, 2->1 (all expanding).
        NU_X0(1:ND-1,1:NCF)=RV0
        DO I=1,ND-1
          NU_CUP(I,1)=AV0; NU_KUP(I,1)=4
          NU_CUP(I,4)=0.4_LDP; NU_KUP(I,4)=3
          NU_CUP(I,3)=0.4_LDP; NU_KUP(I,3)=2
          NU_CUP(I,2)=0.4_LDP; NU_KUP(I,2)=1
        END DO
      CASE('mix')
! Chains feeding rings: 1,2 expand; {3,4} 2-cycle; 5,6 compress into
! 2-cycle {6,7}; 8 is a red-edge outflow boundary.
        NU_X0(1:ND-1,1:NCF)=RV0
        DO I=1,ND-1
          NU_CUP(I,2)=0.4_LDP; NU_KUP(I,2)=1
          NU_CUP(I,3)=0.5_LDP; NU_KUP(I,3)=4
          NU_CUP(I,4)=0.3_LDP; NU_KUP(I,4)=3
          NU_CUP(I,5)=0.4_LDP; NU_KUP(I,5)=6
          NU_CUP(I,6)=0.5_LDP; NU_KUP(I,6)=7
          NU_CUP(I,7)=0.3_LDP; NU_KUP(I,7)=6
        END DO
      CASE('depth')
! Expanding chain plus TA/TC reads on the other depths: depth 2's
! equation sees J(1),J(3) -- the lagged depth terms.
        DO K=2,NCF
          DO I=1,ND-1
            NU_CUP(I,K)=AV0
            NU_KUP(I,K)=K-1
          END DO
        END DO
        NU_X0(1:ND-1,1:NCF)=RV0
        NU_TA(2,1:NCF)=0.1_LDP
        NU_TC(2,1:NCF)=0.1_LDP
        J_NEXT_FRQ(1,1:NCF)=0.7_LDP
        J_NEXT_FRQ(3,1:NCF)=0.9_LDP
      CASE('fresh')
! Only the first half of the band was captured this pass; the rest
! must stay at the seeded field.
        DO K=2,NCF
          DO I=1,ND-1
            NU_CUP(I,K)=AV0
            NU_KUP(I,K)=K-1
          END DO
        END DO
        NU_X0(1:ND-1,1:NCF)=RV0
        NU_FRESH(6:NCF)=.FALSE.
        J_NEXT_FRQ(1:ND,6:NCF)=3.0_LDP
      CASE('sing')
! All-a=1 ring: cyclic product exactly 1 -- the incompatible cycle
! freezes its members.
        NU_X0(1:ND-1,1:NCF)=RV0
        DO I=1,ND-1
          NU_CUP(I,2)=1.0_LDP; NU_KUP(I,2)=3
          NU_CUP(I,3)=1.0_LDP; NU_KUP(I,3)=2
        END DO
        J_NEXT_FRQ(1:ND,2)=4.0_LDP
        J_NEXT_FRQ(1:ND,3)=5.0_LDP
      CASE('neg')
! Negative sources produce a negative exact solution -- every cell
! must be rejected to the previous iterate.
        DO K=2,NCF
          DO I=1,ND-1
            NU_CUP(I,K)=AV0
            NU_KUP(I,K)=K-1
          END DO
        END DO
        NU_X0(1:ND-1,1:NCF)=-2.0_LDP
        J_NEXT_FRQ(1:ND,1:NCF)=0.7_LDP
      CASE('nonfin')
! Row (I=2,K=4): dead diagonal; (I=1,K=6): nonfinite RHS. Freeze
! them; downwind readers consume the old value.
        DO K=2,NCF
          DO I=1,ND-1
            NU_CUP(I,K)=AV0
            NU_KUP(I,K)=K-1
          END DO
        END DO
        NU_X0(1:ND-1,1:NCF)=RV0
        NU_TB(2,4)=0.0_LDP
        T3=0.0_LDP
        NU_X0(1,6)=1.0_LDP/T3
        J_NEXT_FRQ(2,4)=2.5_LDP
        J_NEXT_FRQ(1,6)=2.5_LDP
      CASE DEFAULT
        WRITE(*,'(A)')'RESULT unknown_mode=1'
        STOP 1
      END SELECT
!
      ALLOCATE(J0(ND,NCF))
      J0=J_NEXT_FRQ
      CALL MOM_NU_CHAIN_SOLVE(NCF,CHG,N5,N50)
!
      WRITE(*,'(A,ES16.8)')'NU_CHG=',CHG
      WRITE(*,'(A,I10)')'NU_N5=',N5
      WRITE(*,'(A,I10)')'NU_N50=',N50
      WRITE(*,'(A,I10)')'NU_NRING=',MOM_NU_NRING
      WRITE(*,'(A,I10)')'NU_NSING=',MOM_NU_NSING
      WRITE(*,'(A,I10)')'NU_NFRZ=',MOM_NU_NFRZ
      WRITE(*,'(A,I10)')'NU_NNEG=',MOM_NU_NNEG
!
      OPEN(UNIT=301,FILE=TRIM(OUTTAG)//'_xv.txt',STATUS='REPLACE')
      WRITE(301,'(A)')'# I K JOLD XNEW'
      DO I=1,ND
        DO K=1,NCF
          WRITE(301,'(2I6,2ES25.16)')I,K,J0(I,K),J_NEXT_FRQ(I,K)
        END DO
      END DO
      CLOSE(301)
!
      WRITE(*,'(A)')'RESULT_OK'
      END
