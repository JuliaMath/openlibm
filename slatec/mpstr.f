*DECK MPSTR
      SUBROUTINE MPSTR (X, Y)
C***BEGIN PROLOGUE  MPSTR
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPSTR-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  Sets Y = X for 'mp' X and Y.
C
C  The arguments X(*) and Y(*) are INTEGER arrays of size 30.  See the
C  comments in the routine MPBLAS for the reason for this choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   ??????  Modified for use with BLAS.  Blank COMMON changed to named
C           COMMON.  R given dimension 12.
C   890206  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPSTR
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, X(*), Y(*)
C***FIRST EXECUTABLE STATEMENT  MPSTR
      DO 10 I = 1, T+2
         Y(I) = X(I)
   10 CONTINUE
      RETURN
      END
