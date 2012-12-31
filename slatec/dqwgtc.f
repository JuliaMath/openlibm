*DECK DQWGTC
      DOUBLE PRECISION FUNCTION DQWGTC (X, C, P2, P3, P4, KP)
C***BEGIN PROLOGUE  DQWGTC
C***SUBSIDIARY
C***PURPOSE  This function subprogram is used together with the
C            routine DQAWC and defines the WEIGHT function.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (QWGTC-S, DQWGTC-D)
C***KEYWORDS  CAUCHY PRINCIPAL VALUE, WEIGHT FUNCTION
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***SEE ALSO  DQK15W
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810101  DATE WRITTEN
C   830518  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DQWGTC
C
      DOUBLE PRECISION C,P2,P3,P4,X
      INTEGER KP
C***FIRST EXECUTABLE STATEMENT  DQWGTC
      DQWGTC = 0.1D+01/(X-C)
      RETURN
      END
