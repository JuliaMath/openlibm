*DECK NUMXER
      FUNCTION NUMXER (NERR)
C***BEGIN PROLOGUE  NUMXER
C***PURPOSE  Return the most recent error number.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      INTEGER (NUMXER-I)
C***KEYWORDS  ERROR NUMBER, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        NUMXER returns the most recent error number,
C        in both NUMXER and the parameter NERR.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   910411  Made user-callable and added KEYWORDS section.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  NUMXER
C***FIRST EXECUTABLE STATEMENT  NUMXER
      NERR = J4SAVE(1,0,.FALSE.)
      NUMXER = NERR
      RETURN
      END
