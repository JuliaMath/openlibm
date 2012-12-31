*DECK PRVEC
      FUNCTION PRVEC (M, U, V)
C***BEGIN PROLOGUE  PRVEC
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (PRVEC-S, DPRVEC-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C  This subroutine computes the inner product of a vector U
C  with the imaginary product or mate vector corresponding to V
C
C***SEE ALSO  BVSUP
C***ROUTINES CALLED  SDOT
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  PRVEC
C
      DIMENSION U(*),V(*)
C***FIRST EXECUTABLE STATEMENT  PRVEC
      N=M/2
      NP=N+1
      VP=SDOT(N,U(1),1,V(NP),1)
      PRVEC=SDOT(N,U(NP),1,V(1),1) - VP
      RETURN
      END
