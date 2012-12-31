*DECK XERMAX
      SUBROUTINE XERMAX (MAX)
C***BEGIN PROLOGUE  XERMAX
C***PURPOSE  Set maximum number of times any error message is to be
C            printed.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERMAX-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XERMAX sets the maximum number of times any message
C        is to be printed.  That is, non-fatal messages are
C        not to be printed after they have occurred MAX times.
C        Such non-fatal messages may be printed less than
C        MAX times even if they occur MAX times, if error
C        suppression mode (KONTRL=0) is ever in effect.
C
C     Description of Parameter
C      --Input--
C        MAX - the maximum number of times any one message
C              is to be printed.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERMAX
C***FIRST EXECUTABLE STATEMENT  XERMAX
      JUNK = J4SAVE(4,MAX,.TRUE.)
      RETURN
      END
