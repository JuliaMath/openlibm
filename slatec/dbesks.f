*DECK DBESKS
      SUBROUTINE DBESKS (XNU, X, NIN, BK)
C***BEGIN PROLOGUE  DBESKS
C***PURPOSE  Compute a sequence of modified Bessel functions of the
C            third kind of fractional order.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B3
C***TYPE      DOUBLE PRECISION (BESKS-S, DBESKS-D)
C***KEYWORDS  FNLIB, FRACTIONAL ORDER, MODIFIED BESSEL FUNCTION,
C             SEQUENCE OF BESSEL FUNCTIONS, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESKS computes a sequence of modified Bessel functions of the third
C kind of order XNU + I at X, where X .GT. 0, XNU lies in (-1,1),
C and I = 0, 1, ... , NIN - 1, if NIN is positive and I = 0, 1, ... ,
C NIN + 1, if NIN is negative.  On return, the vector BK(.) contains
C the results at X for order starting at XNU.  XNU, X, and BK are
C double precision.  NIN is an integer.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBSKES, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESKS
      DOUBLE PRECISION XNU, X, BK(*), EXPXI, XMAX, D1MACH
      SAVE XMAX
      DATA XMAX / 0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBESKS
      IF (XMAX.EQ.0.D0) XMAX = -LOG (D1MACH(1))
C
      IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DBESKS',
     +   'X SO BIG BESSEL K UNDERFLOWS', 1, 2)
C
      CALL DBSKES (XNU, X, NIN, BK)
C
      EXPXI = EXP (-X)
      N = ABS (NIN)
      DO 20 I=1,N
        BK(I) = EXPXI * BK(I)
 20   CONTINUE
C
      RETURN
      END
