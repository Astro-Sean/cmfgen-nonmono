C-----------------------------------------------------------------------
C   Purpose             - EVALUATE THE COMPLEMENTED ERROR FUNCTION OF
C                           A DOUBLE PRECISION ARGUMENT
C
C   Usage               - RESULT = S15ADF(Y)
C
C   Arguments    Y      - Input double precision argument of the
C                           complemented error function.
C                ERFC  - Output double precision value of the
C                           complemented error function. erfc must be
C                           typed double precision in the calling
C                           program.
C
C-----------------------------------------------------------------------
C
      FUNCTION S15ADF(Y,IFAIL)
	USE SET_KIND_MODULE
C
C Specifications for arguments
C
      REAL(KIND=LDP) S15ADF,Y
C
C Local variables
C
      REAL(KIND=LDP) P(5),Q(4),P1(9),Q1(8),P2(6),Q2(5)
      REAL(KIND=LDP) XMIN,XLARGE,SQRPI,X,
     *                   RES,XSQ,XNUM,XDEN,XI,XBIG
      INTEGER            ISW,I,IFAIL
C
C Coefficients for 0.0 .LE. Y .LT.
C                                  .477
      DATA               P(1)/113.8641541510502_LDP/,
     *                   P(2)/377.4852376853020_LDP/,
     *                   P(3)/3209.377589138469_LDP/,
     *                   P(4)/.1857777061846032_LDP/,
     *                   P(5)/3.161123743870566_LDP/
      DATA               Q(1)/244.0246379344442_LDP/,
     *                   Q(2)/1282.616526077372_LDP/,
     *                   Q(3)/2844.236833439171_LDP/,
     *                   Q(4)/23.60129095234412_LDP/
C                                  COEFFICIENTS FOR .477 .LE. Y
C                                  .LE. 4.0
      DATA               P1(1)/8.883149794388376_LDP/,
     *                   P1(2)/66.11919063714163_LDP/,
     *                   P1(3)/298.6351381974001_LDP/,
     *                   P1(4)/881.9522212417691_LDP/,
     *                   P1(5)/1712.047612634071_LDP/,
     *                   P1(6)/2051.078377826071_LDP/,
     *                   P1(7)/1230.339354797997_LDP/,
     *                   P1(8)/2.153115354744038E-8_LDP/,
     *                   P1(9)/.5641884969886701_LDP/
      DATA               Q1(1)/117.6939508913125_LDP/,
     *                   Q1(2)/537.1811018620099_LDP/,
     *                   Q1(3)/1621.389574566690_LDP/,
     *                   Q1(4)/3290.799235733460_LDP/,
     *                   Q1(5)/4362.619090143247_LDP/,
     *                   Q1(6)/3439.367674143722_LDP/,
     *                   Q1(7)/1230.339354803749_LDP/,
     *                   Q1(8)/15.74492611070983_LDP/
C                                  COEFFICIENTS FOR 4.0 .LT. Y
      DATA               P2(1)/-3.603448999498044E-01_LDP/,
     *                   P2(2)/-1.257817261112292E-01_LDP/,
     *                   P2(3)/-1.608378514874228E-02_LDP/,
     *                   P2(4)/-6.587491615298378E-04_LDP/,
     *                   P2(5)/-1.631538713730210E-02_LDP/,
     *                   P2(6)/-3.053266349612323E-01_LDP/
      DATA               Q2(1)/1.872952849923460_LDP/,
     *                   Q2(2)/5.279051029514284E-01_LDP/,
     *                   Q2(3)/6.051834131244132E-02_LDP/,
     *                   Q2(4)/2.335204976268692E-03_LDP/,
     *                   Q2(5)/2.568520192289822_LDP/
C                                  CONSTANTS
      DATA               XMIN/1.0E-10_LDP/,XLARGE/6.375_LDP/
C                                  ERFC(XBIG) .APPROX. DETAP
      DATA               XBIG/13.3_LDP/
      DATA               SQRPI/.5641895835477563_LDP/
C                                  FIRST EXECUTABLE STATEMENT
      IFAIL=0
      X = Y
      ISW = 1
      IF (X.GE.0.0_LDP) GO TO 5
      ISW = -1
      X = -X
    5 IF (X.LT..477_LDP) GO TO 10
      IF (X.LE.4.0_LDP) GO TO 30
      IF (ISW .GT. 0) GO TO 40
      IF (X.LT.XLARGE) GO TO 45
      RES = 2.0_LDP
      GO TO 70
C                                  ABS(Y) .LT. .477, EVALUATE
C                                  APPROXIMATION FOR ERFC
   10 IF (X.LT.XMIN) GO TO 20
      XSQ = X*X
      XNUM = P(4)*XSQ+P(5)
      XDEN = XSQ+Q(4)
      DO 15 I = 1,3
         XNUM = XNUM*XSQ+P(I)
         XDEN = XDEN*XSQ+Q(I)
   15 CONTINUE
      RES = X*XNUM/XDEN
      GO TO 25
   20 RES = X*P(3)/Q(3)
   25 IF (ISW.EQ.-1) RES = -RES
      RES = 1.0_LDP-RES
      GO TO 70
C                                  .477 .LE. ABS(Y) .LE. 4.0
C                                  EVALUATE APPROXIMATION FOR ERFC
   30 XSQ = X*X
      XNUM = P1(8)*X+P1(9)
      XDEN = X+Q1(8)
      DO 35 I=1,7
         XNUM = XNUM*X+P1(I)
         XDEN = XDEN*X+Q1(I)
   35 CONTINUE
      RES = XNUM/XDEN
      GO TO 60
C                                  4.0 .LT. ABS(Y), EVALUATE
C                                  MINIMAX APPROXIMATION FOR ERFC
   40 IF (X.GT.XBIG) GO TO 65
   45 XSQ = X*X
      XI = 1.0_LDP/XSQ
      XNUM= P2(5)*XI+P2(6)
      XDEN = XI+Q2(5)
      DO 50 I = 1,4
         XNUM = XNUM*XI+P2(I)
         XDEN = XDEN*XI+Q2(I)
   50 CONTINUE
      RES = (SQRPI+XI*XNUM/XDEN)/X
   60 RES = RES*EXP(-XSQ)
      IF (ISW.EQ.-1) RES = 2.0_LDP-RES
      GO TO 70
   65 RES = 0.0_LDP
   70 S15ADF = RES
C
      RETURN
      END
