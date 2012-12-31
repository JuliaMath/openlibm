*DECK R9PAK
      FUNCTION R9PAK (Y, N)
C***BEGIN PROLOGUE  R9PAK
C***PURPOSE  Pack a base 2 exponent into a floating point number.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  A6B
C***TYPE      SINGLE PRECISION (R9PAK-S, D9PAK-D)
C***KEYWORDS  FNLIB, PACK
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Pack a base 2 exponent into floating point number Y.  This
C routine is almost the inverse of R9UPAK.  It is not exactly
C the inverse, because ABS(X) need not be between 0.5 and
C 1.0.  If both R9PAK and 2.0**N were known to be in range, we
C could compute
C       R9PAK = Y * 2.0**N .
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  I1MACH, R1MACH, R9UPAK, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901009  Routine used I1MACH(7) where it should use I1MACH(10),
C           Corrected (RWC)
C***END PROLOGUE  R9PAK
      LOGICAL FIRST
      SAVE NMIN, NMAX, A1N210, FIRST
      DATA A1N210 / 3.321928094 887362 E0/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  R9PAK
      IF (FIRST) THEN
         A1N2B = 1.0
         IF (I1MACH(10).NE.2) A1N2B = R1MACH(5)*A1N210
         NMIN = A1N2B*I1MACH(12)
         NMAX = A1N2B*I1MACH(13)
      ENDIF
      FIRST = .FALSE.
C
      CALL R9UPAK(Y,R9PAK,NY)
C
      NSUM = N + NY
      IF (NSUM.LT.NMIN) GO TO 40
      IF (NSUM .GT. NMAX) CALL XERMSG ('SLATEC', 'R9PAK',
     +   'PACKED NUMBER OVERFLOWS', 2, 2)
C
      IF (NSUM.EQ.0) RETURN
      IF (NSUM.GT.0) GO TO 30
C
 20   R9PAK = 0.5*R9PAK
      NSUM = NSUM + 1
      IF(NSUM.NE.0) GO TO 20
      RETURN
C
30    R9PAK = 2.0*R9PAK
      NSUM = NSUM - 1
      IF(NSUM.NE.0) GO TO 30
      RETURN
C
40    CALL XERMSG ('SLATEC', 'R9PAK', 'PACKED NUMBER UNDERFLOWS', 1, 1)
      R9PAK = 0.0
      RETURN
C
      END
