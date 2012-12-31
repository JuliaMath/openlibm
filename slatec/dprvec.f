*DECK DPRVEC
      DOUBLE PRECISION FUNCTION DPRVEC (M, U, V)
C***BEGIN PROLOGUE  DPRVEC
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (PRVEC-S, DPRVEC-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C  This subroutine computes the inner product of a vector U
C  with the imaginary product or mate vector corresponding to V.
C
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  DDOT
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DPRVEC
C
      DOUBLE PRECISION DDOT
      INTEGER M, N, NP
      DOUBLE PRECISION U(*), V(*), VP
C***FIRST EXECUTABLE STATEMENT  DPRVEC
      N = M/2
      NP = N + 1
      VP = DDOT(N,U(1),1,V(NP),1)
      DPRVEC = DDOT(N,U(NP),1,V(1),1) - VP
      RETURN
      END
