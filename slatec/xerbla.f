*DECK XERBLA
      SUBROUTINE XERBLA (SRNAME, INFO)
C***BEGIN PROLOGUE  XERBLA
C***SUBSIDIARY
C***PURPOSE  Error handler for the Level 2 and Level 3 BLAS Routines.
C***LIBRARY   SLATEC
C***CATEGORY  R3
C***TYPE      ALL (XERBLA-A)
C***KEYWORDS  ERROR MESSAGE
C***AUTHOR  Dongarra, J. J., (ANL)
C***DESCRIPTION
C
C  Purpose
C  =======
C
C  It is called by Level 2 and 3 BLAS routines if an input parameter
C  is invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*6.
C           On entry, SRNAME specifies the name of the routine which
C           called XERBLA.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   860720  DATE WRITTEN
C   910610  Routine rewritten to serve as an interface between the
C           Level 2 and Level 3 BLAS routines and the SLATEC error
C           handler XERMSG.  (BKS)
C***END PROLOGUE  XERBLA
C
C     ..    Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*6        SRNAME
      CHARACTER*2        XERN1
C
C***FIRST EXECUTABLE STATEMENT  XERBLA
C
      WRITE (XERN1, '(I2)') INFO
      CALL XERMSG ('SLATEC', SRNAME, 'On entry to '//SRNAME//
     $             ' parameter number '//XERN1//' had an illegal value',
     $             INFO,1)
C
      RETURN
C
C     End of XERBLA.
C
      END
