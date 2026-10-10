C
C Subroutine to change the color settings of the various color indexes
C
      SUBROUTINE CHANGE_COLOR(RED,BLUE,GREEN)
	USE SET_KIND_MODULE
      USE GEN_IN_INTERFACE
C
C Altered 15-Nov-2021 : Cleaned up real variables in calls.
C Altered 01-Jul-1997 : Call tah utilizes GEN_IN_MULTI_? changed.
C Created September 1996 by Gregson Vaux
C
      IMPLICIT NONE
C
      INTEGER I,NC,JUST,AXIS,IRET
      REAL*4 RED(0:15),BLUE(0:15),GREEN(0:15)
      REAL*4 XBEG,XEND,YBEG,YEND,HEIGHT
      REAL*4 COL(3),T1,T2
      INTEGER ROUT,BOUT,GOUT
      CHARACTER*8 ICH(1:15),RCH(0:15),BCH(0:15),GCH(0:15)
C
      INTEGER, PARAMETER :: ITHREE=3
      REAL, PARAMETER :: RZERO=0.0
      REAL, PARAMETER :: RTHREE=3.0
C
C Clear screen and ouput sample of pens to plot device
C
      XBEG=0.0
      XEND=100.0
      YBEG=0.0
      YEND=34.0
      JUST=0
      AXIS=-2
      HEIGHT=1.0
      CALL PGENV(XBEG,XEND,YBEG,YEND,JUST,AXIS)
      DO I=0,15
        CALL PGQCR(I,RED(I),GREEN(I),BLUE(I))
      END DO
      CALL PGSCH(HEIGHT)
C
C Loop to change pen colors
C
100   CALL PGSCI(1)
      CALL PGERAS
      CALL PGTEXT(RZERO,RZERO,'PEN        ,  RED,       GREEN,         BLUE')
      ROUT=(RED(0)*1000)
      GOUT=(GREEN(0)*1000)
      BOUT=(BLUE(0)*1000)
      CALL PGNUMB(ROUT,-3,1,RCH,NC)
      CALL PGNUMB(GOUT,-3,1,GCH,NC)
      CALL PGNUMB(BOUT,-3,1,BCH,NC)
      CALL PGTEXT(RZERO,RTHREE,'0')
      T1=20.0; CALL PGTEXT(T1,RTHREE,RCH)
      T1=60.0; CALL PGTEXT(T1,RTHREE,GCH)
      T1=40.0; CALL PGTEXT(T1,RTHREE,BCH)
      DO I=1,15
        CALL PGSCI(I)
        ROUT=(RED(I)*1000)
        BOUT=(BLUE(I)*1000)
        GOUT=(GREEN(I)*1000)
        CALL PGNUMB(I,0,1,ICH,NC)
        CALL PGNUMB(ROUT,-3,1,RCH,NC)
        CALL PGNUMB(GOUT,-3,1,GCH,NC)
        CALL PGNUMB(BOUT,-3,1,BCH,NC)
!
	T2=(I*2.0)+3
        CALL PGTEXT(RZERO,T2,ICH)
        T1=20.0; CALL PGTEXT(T1,T2,RCH)
        T1=40.0; CALL PGTEXT(T1,T2,GCH)
        T1=60.0; CALL PGTEXT(T1,T2,BCH)
      END DO
      CALL GEN_IN(I,'Pen 0-15 (>15 exits)')
      IF (I .GT. 15) RETURN
C
C Set defaults and get new values.
C
      COL(1)=RED(I)
      COL(2)=GREEN(I)
      COL(3)=BLUE(I)
      CALL GEN_IN(COL,IRET,ITHREE,'RGB color values (0-1)')
      RED(I)=COL(1)
      GREEN(I)=COL(2)
      BLUE(I)=COL(3)
      CALL PGSCR(I,RED(I),GREEN(I),BLUE(I))
      GOTO 100
      END
