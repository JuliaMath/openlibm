*DECK MPMLP
      SUBROUTINE MPMLP (U, V, W, J)
C***BEGIN PROLOGUE  MPMLP
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPMLP-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C Performs inner multiplication loop for MPMUL. Carries are not pro-
C pagated in inner loop, which saves time at the expense of space.
C
C***SEE ALSO  DQDOTA, DQDOTI
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  MPMLP
      INTEGER U(*), V(*), W
C***FIRST EXECUTABLE STATEMENT  MPMLP
      DO 10 I = 1, J
   10 U(I) = U(I) + W*V(I)
      RETURN
      END
