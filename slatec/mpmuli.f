*DECK MPMULI
      SUBROUTINE MPMULI (X, IY, Z)
C***BEGIN PROLOGUE  MPMULI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPMULI-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C Multiplies 'mp' X by single-precision integer IY giving 'mp' Z.
C This is faster than using MPMUL.  Result is ROUNDED.
C Multiplication by 1 may be used to normalize a number
C even if the last digit is B.
C
C***SEE ALSO  DQDOTA, DQDOTI
C***ROUTINES CALLED  MPMUL2
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  MPMULI
      INTEGER X(*), Z(*)
C***FIRST EXECUTABLE STATEMENT  MPMULI
      CALL MPMUL2 (X, IY, Z, 0)
      RETURN
      END
