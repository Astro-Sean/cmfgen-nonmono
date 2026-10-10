C
C
	  SUBROUTINE COMPION(HYD,NHYD,DHYD,NION,YV,ND)
	  USE SET_KIND_MODULE
	  IMPLICIT NONE
C
C Altered 20-Mar-1997: NB: YV is now REAL(KIND=LDP)
C
	  INTEGER NHYD,NION,ND,I,J
	  REAL(KIND=LDP) HYD(NHYD,ND),DHYD(NION,ND),T1
	  REAL(KIND=LDP) YV(ND)
C
	  IF(HYD(1,ND) .NE. 0
	1              .AND. DHYD(1,ND) .NE. 0)THEN
	    DO J=1,ND
	      YV(J)=0.0_LDP
	      DO I=1,NHYD
	        YV(J)=YV(J)+HYD(I,J)
	      END DO
	      T1=0.0_LDP
	      DO I=1,NION
	        T1=T1+DHYD(I,J)
	      END DO
	      YV(J)=LOG10(YV(J)/T1)
	    END DO
	  ELSE
	    WRITE(6,*)'This species unavailable - populations zero'
	    WRITE(6,*)HYD(1,ND),DHYD(1,ND)
	  END IF
C
	  RETURN
	  END
