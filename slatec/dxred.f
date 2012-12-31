*DECK DXRED
      SUBROUTINE DXRED (X, IX, IERROR)
C***BEGIN PROLOGUE  DXRED
C***PURPOSE  To provide double-precision floating-point arithmetic
C            with an extended exponent range.
C***LIBRARY   SLATEC
C***CATEGORY  A3D
C***TYPE      DOUBLE PRECISION (XRED-S, DXRED-D)
C***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
C***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
C           Smith, John M., (NBS and George Mason University)
C***DESCRIPTION
C     DOUBLE PRECISION X
C     INTEGER IX
C
C                  IF
C                  RADIX**(-2L) .LE. (ABS(X),IX) .LE. RADIX**(2L)
C                  THEN DXRED TRANSFORMS (X,IX) SO THAT IX=0.
C                  IF (X,IX) IS OUTSIDE THE ABOVE RANGE,
C                  THEN DXRED TAKES NO ACTION.
C                  THIS SUBROUTINE IS USEFUL IF THE
C                  RESULTS OF EXTENDED-RANGE CALCULATIONS
C                  ARE TO BE USED IN SUBSEQUENT ORDINARY
C                  DOUBLE-PRECISION CALCULATIONS.
C
C***SEE ALSO  DXSET
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DXBLK2
C***REVISION HISTORY  (YYMMDD)
C   820712  DATE WRITTEN
C   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  DXRED
      DOUBLE PRECISION X
      INTEGER IX
      DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R, XA
      INTEGER L, L2, KMAX
      COMMON /DXBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE /DXBLK2/
C
C***FIRST EXECUTABLE STATEMENT  DXRED
      IERROR=0
      IF (X.EQ.0.0D0) GO TO 90
      XA = ABS(X)
      IF (IX.EQ.0) GO TO 70
      IXA = ABS(IX)
      IXA1 = IXA/L2
      IXA2 = MOD(IXA,L2)
      IF (IX.GT.0) GO TO 40
   10 CONTINUE
      IF (XA.GT.1.0D0) GO TO 20
      XA = XA*RAD2L
      IXA1 = IXA1 + 1
      GO TO 10
   20 XA = XA/RADIX**IXA2
      IF (IXA1.EQ.0) GO TO 70
      DO 30 I=1,IXA1
        IF (XA.LT.1.0D0) GO TO 100
        XA = XA/RAD2L
   30 CONTINUE
      GO TO 70
C
   40 CONTINUE
      IF (XA.LT.1.0D0) GO TO 50
      XA = XA/RAD2L
      IXA1 = IXA1 + 1
      GO TO 40
   50 XA = XA*RADIX**IXA2
      IF (IXA1.EQ.0) GO TO 70
      DO 60 I=1,IXA1
        IF (XA.GT.1.0D0) GO TO 100
        XA = XA*RAD2L
   60 CONTINUE
   70 IF (XA.GT.RAD2L) GO TO 100
      IF (XA.GT.1.0D0) GO TO 80
      IF (RAD2L*XA.LT.1.0D0) GO TO 100
   80 X = SIGN(XA,X)
   90 IX = 0
  100 RETURN
      END
