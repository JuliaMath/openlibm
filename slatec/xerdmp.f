*DECK XERDMP
      SUBROUTINE XERDMP
C***BEGIN PROLOGUE  XERDMP
C***PURPOSE  Print the error tables and then clear them.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERDMP-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XERDMP prints the error tables, then clears them.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  XERSVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900510  Changed call of XERSAV to XERSVE.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERDMP
C***FIRST EXECUTABLE STATEMENT  XERDMP
      CALL XERSVE (' ',' ',' ',0,0,0,KOUNT)
      RETURN
      END
