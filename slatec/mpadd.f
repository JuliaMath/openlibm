*DECK MPADD
      SUBROUTINE MPADD (X, Y, Z)
C***BEGIN PROLOGUE  MPADD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPADD-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C Adds X and Y, forming result in Z, where X, Y and Z are 'mp'
C  (multiple precision) numbers.  Four guard digits are used,
C  and then R*-rounding.
C
C***SEE ALSO  DQDOTA, DQDOTI
C***ROUTINES CALLED  MPADD2
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  MPADD
      INTEGER X(*), Y(*), Z(*)
C***FIRST EXECUTABLE STATEMENT  MPADD
      CALL MPADD2 (X, Y, Z, Y, 0)
      RETURN
      END
