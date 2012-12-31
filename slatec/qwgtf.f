*DECK QWGTF
      REAL FUNCTION QWGTF (X, OMEGA, P2, P3, P4, INTEGR)
C***BEGIN PROLOGUE  QWGTF
C***SUBSIDIARY
C***PURPOSE  This function subprogram is used together with the
C            routine QAWF and defines the WEIGHT function.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (QWGTF-S, DQWGTF-D)
C***KEYWORDS  COS OR SIN IN WEIGHT FUNCTION
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***SEE ALSO  QK15W
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810101  DATE WRITTEN
C   830518  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  QWGTF
C
      REAL OMEGA,OMX,P2,P3,P4,X
      INTEGER INTEGR
C***FIRST EXECUTABLE STATEMENT  QWGTF
      OMX = OMEGA*X
      GO TO(10,20),INTEGR
   10 QWGTF = COS(OMX)
      GO TO 30
   20 QWGTF = SIN(OMX)
   30 RETURN
      END
