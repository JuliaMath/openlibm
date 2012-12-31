*DECK XSETUN
      SUBROUTINE XSETUN (IUNIT)
C***BEGIN PROLOGUE  XSETUN
C***PURPOSE  Set output file to which error messages are to be sent.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3B
C***TYPE      ALL (XSETUN-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XSETUN sets the output file to which error messages are to
C        be sent.  Only one file will be used.  See XSETUA for
C        how to declare more than one file.
C
C     Description of Parameter
C      --Input--
C        IUNIT - an input parameter giving the logical unit number
C                to which error messages are to be sent.
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
C***END PROLOGUE  XSETUN
C***FIRST EXECUTABLE STATEMENT  XSETUN
      JUNK = J4SAVE(3,IUNIT,.TRUE.)
      JUNK = J4SAVE(5,1,.TRUE.)
      RETURN
      END
