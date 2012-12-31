*DECK BESKS
      SUBROUTINE BESKS (XNU, X, NIN, BK)
C***BEGIN PROLOGUE  BESKS
C***PURPOSE  Compute a sequence of modified Bessel functions of the
C            third kind of fractional order.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B3
C***TYPE      SINGLE PRECISION (BESKS-S, DBESKS-D)
C***KEYWORDS  FNLIB, FRACTIONAL ORDER, MODIFIED BESSEL FUNCTION,
C             SEQUENCE OF BESSEL FUNCTIONS, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESKS computes a sequence of modified Bessel functions of the third
C kind of order XNU + I at X, where X .GT. 0, XNU lies in (-1,1),
C and I = 0, 1, ... , NIN - 1, if NIN is positive and I = 0, 1, ... ,
C NIN + 1, if NIN is negative.  On return, the vector BK(.) Contains
C the results at X for order starting at XNU.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BESKES, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  BESKS
      DIMENSION BK(*)
      SAVE XMAX
      DATA XMAX / 0.0 /
C***FIRST EXECUTABLE STATEMENT  BESKS
      IF (XMAX.EQ.0.0) XMAX = -LOG (R1MACH(1))
C
      IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'BESKS',
     +   'X SO BIG BESSEL K UNDERFLOWS', 1, 2)
C
      CALL BESKES (XNU, X, NIN, BK)
C
      EXPXI = EXP (-X)
      N = ABS (NIN)
      DO 20 I=1,N
        BK(I) = EXPXI * BK(I)
 20   CONTINUE
C
      RETURN
      END
