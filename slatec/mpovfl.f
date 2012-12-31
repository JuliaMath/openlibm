*DECK MPOVFL
      SUBROUTINE MPOVFL (X)
C***BEGIN PROLOGUE  MPOVFL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPOVFL-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  Called on multiple-precision overflow, i.e. when the
C  exponent of 'mp' number X would exceed M.  At present execution is
C  terminated with an error message after calling MPMAXR(X), but it
C  would be possible to return, possibly updating a counter and
C  terminating execution after a preset number of overflows.  Action
C  could easily be determined by a flag in labelled common.
C
C  The argument X(*) is an INTEGER array of size 30.  See the comments
C  in the routine MPBLAS for the reason for this choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***ROUTINES CALLED  MPCHK, MPERR, MPMAXR
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   ??????  Modified for use with BLAS.  Blank COMMON changed to named
C           COMMON.  R given dimension 12.
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPOVFL
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, X(*)
C***FIRST EXECUTABLE STATEMENT  MPOVFL
      CALL MPCHK (1, 4)
C SET X TO LARGEST POSSIBLE POSITIVE NUMBER
      CALL MPMAXR (X)
      WRITE (LUN, 10)
   10 FORMAT (' *** CALL TO MPOVFL, MP OVERFLOW OCCURRED ***')
C TERMINATE EXECUTION BY CALLING MPERR
      CALL MPERR
      RETURN
      END
