! Driver for GFF / GFF_VEC (subs/gff.f).
!
! The alteration under test is the GFF_VEC bad-T diagnostic: a
! non-positive (or NaN) temperature is flagged to stdout before
! LOG10 turns it into a NaN, so a corrupt temperature can be
! traced to its depth index.
!
! Usage: drv CASE
!   ok   : all-positive T, valid Gaunt factors, no diagnostic
!   bad  : T(3) <= 0 -> exactly one 'bad T' line, I=3
!   nan  : T(2) = NaN -> .NOT.(T.GT.0) still catches it
!   gff  : scalar GFF vs GFF_VEC element agreement on valid T
!
      PROGRAM DRV_GFF
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=5
      REAL(KIND=LDP) T(ND),GFF_VAL(ND),RNU,ZION,GFF
      CHARACTER(LEN=32) CASE
      INTEGER I
      CALL GETARG(1,CASE)
      RNU=3.29_LDP
      ZION=1.0_LDP
      T=(/1.0_LDP,2.0_LDP,4.0_LDP,8.0_LDP,16.0_LDP/)
      IF(CASE .EQ. 'bad')T(3)=-1.0_LDP
      IF(CASE .EQ. 'nan')T(2)=0.0_LDP/0.0_LDP
      CALL GFF_VEC(GFF_VAL,RNU,T,ZION,ND)
      DO I=1,ND
        WRITE(*,'(A,I2,1X,ES16.9)')'GV ',I,GFF_VAL(I)
      END DO
      IF(CASE .EQ. 'gff')THEN
        WRITE(*,'(A,ES16.9)')'SC ',GFF(RNU,T(2),ZION)
      END IF
      END
