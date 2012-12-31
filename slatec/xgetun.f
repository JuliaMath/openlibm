*DECK XGETUN
      SUBROUTINE XGETUN (IUNIT)
C***BEGIN PROLOGUE  XGETUN
C***PURPOSE  Return the (first) output file to which error messages
C            are being sent.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XGETUN-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XGETUN gets the (first) output file to which error messages
C        are being sent.  To find out if more than one file is being
C        used, one must use the XGETUA routine.
C
C     Description of Parameter
C      --Output--
C        IUNIT - the logical unit number of the  (first) unit to
C                which error messages are being sent.
C                A value of zero means that the default file, as
C                defined by the I1MACH routine, is being used.
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
C***END PROLOGUE  XGETUN
C***FIRST EXECUTABLE STATEMENT  XGETUN
      IUNIT = J4SAVE(3,0,.FALSE.)
      RETURN
      END
