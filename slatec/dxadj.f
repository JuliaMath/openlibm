*DECK DXADJ
      SUBROUTINE DXADJ (X, IX, IERROR)
C***BEGIN PROLOGUE  DXADJ
C***PURPOSE  To provide double-precision floating-point arithmetic
C            with an extended exponent range.
C***LIBRARY   SLATEC
C***CATEGORY  A3D
C***TYPE      DOUBLE PRECISION (XADJ-S, DXADJ-D)
C***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
C***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
C           Smith, John M., (NBS and George Mason University)
C***DESCRIPTION
C     DOUBLE PRECISION X
C     INTEGER IX
C
C                  TRANSFORMS (X,IX) SO THAT
C                  RADIX**(-L) .LE. ABS(X) .LT. RADIX**L.
C                  ON MOST COMPUTERS THIS TRANSFORMATION DOES
C                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS
C                  THE NUMBER BASE OF DOUBLE-PRECISION ARITHMETIC.
C
C***SEE ALSO  DXSET
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***COMMON BLOCKS    DXBLK2
C***REVISION HISTORY  (YYMMDD)
C   820712  DATE WRITTEN
C   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  DXADJ
      DOUBLE PRECISION X
      INTEGER IX
      DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R
      INTEGER L, L2, KMAX
      COMMON /DXBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE /DXBLK2/
C
C   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE
C IS
C     2*L .LE. KMAX
C
C THIS CONDITION MUST BE MET BY APPROPRIATE CODING
C IN SUBROUTINE DXSET.
C
C***FIRST EXECUTABLE STATEMENT  DXADJ
      IERROR=0
      IF (X.EQ.0.0D0) GO TO 50
      IF (ABS(X).GE.1.0D0) GO TO 20
      IF (RADIXL*ABS(X).GE.1.0D0) GO TO 60
      X = X*RAD2L
      IF (IX.LT.0) GO TO 10
      IX = IX - L2
      GO TO 70
   10 IF (IX.LT.-KMAX+L2) GO TO 40
      IX = IX - L2
      GO TO 70
   20 IF (ABS(X).LT.RADIXL) GO TO 60
      X = X/RAD2L
      IF (IX.GT.0) GO TO 30
      IX = IX + L2
      GO TO 70
   30 IF (IX.GT.KMAX-L2) GO TO 40
      IX = IX + L2
      GO TO 70
   40 CALL XERMSG ('SLATEC', 'DXADJ', 'overflow in auxiliary index',
     +             207, 1)
      IERROR=207
      RETURN
   50 IX = 0
   60 IF (ABS(IX).GT.KMAX) GO TO 40
   70 RETURN
      END
