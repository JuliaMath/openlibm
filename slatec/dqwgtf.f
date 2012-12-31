*DECK DQWGTF
      DOUBLE PRECISION FUNCTION DQWGTF (X, OMEGA, P2, P3, P4, INTEGR)
C***BEGIN PROLOGUE  DQWGTF
C***SUBSIDIARY
C***PURPOSE  This function subprogram is used together with the
C            routine DQAWF and defines the WEIGHT function.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (QWGTF-S, DQWGTF-D)
C***KEYWORDS  COS OR SIN IN WEIGHT FUNCTION
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
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DQWGTF
C
      DOUBLE PRECISION OMEGA,OMX,P2,P3,P4,X
      INTEGER INTEGR
C***FIRST EXECUTABLE STATEMENT  DQWGTF
      OMX = OMEGA*X
      GO TO(10,20),INTEGR
   10 DQWGTF = COS(OMX)
      GO TO 30
   20 DQWGTF = SIN(OMX)
   30 RETURN
      END
