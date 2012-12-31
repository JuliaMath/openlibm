*DECK MPMAXR
      SUBROUTINE MPMAXR (X)
C***BEGIN PROLOGUE  MPMAXR
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPMAXR-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  Sets X to the largest possible positive 'mp' number.
C
C  The argument X(*) is an INTEGER arrays of size 30.  See the comments
C  in the routine MPBLAS for the reason for this choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***ROUTINES CALLED  MPCHK
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   ??????  Modified for use with BLAS.  Blank COMMON changed to named
C           COMMON.  R given dimension 12.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPMAXR
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, X(*)
C***FIRST EXECUTABLE STATEMENT  MPMAXR
      CALL MPCHK (1, 4)
      IT = B - 1
C SET FRACTION DIGITS TO B-1
      DO 10 I = 1, T
   10 X(I+2) = IT
C SET SIGN AND EXPONENT
      X(1) = 1
      X(2) = M
      RETURN
      END
