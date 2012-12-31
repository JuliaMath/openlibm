*DECK MPCMD
      SUBROUTINE MPCMD (X, DZ)
C***BEGIN PROLOGUE  MPCMD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPCMD-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  Converts multiple-precision X to double-precision DZ. Assumes
C  X is in allowable range for double-precision numbers. There is
C  some loss of accuracy if the exponent is large.
C
C  The argument X(*) is INTEGER array of size 30.  See the comments in
C  the routine MPBLAS for the reason for this choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***ROUTINES CALLED  MPCHK, MPERR
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   ??????  Modified for use with BLAS.  Blank COMMON changed to named
C           COMMON.  R given dimension 12.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPCMD
      DOUBLE PRECISION DB, DZ, DZ2
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, X(*), TM
C***FIRST EXECUTABLE STATEMENT  MPCMD
      CALL MPCHK (1, 4)
      DZ = 0D0
      IF (X(1).EQ.0) RETURN
      DB = DBLE(B)
      DO 10 I = 1, T
      DZ = DB*DZ + DBLE(X(I+2))
      TM = I
C CHECK IF FULL DOUBLE-PRECISION ACCURACY ATTAINED
      DZ2 = DZ + 1D0
C TEST BELOW NOT ALWAYS EQUIVALENT TO - IF (DZ2.LE.DZ) GO TO 20,
C FOR EXAMPLE ON CYBER 76.
      IF ((DZ2-DZ).LE.0D0) GO TO 20
   10 CONTINUE
C NOW ALLOW FOR EXPONENT
   20 DZ = DZ*(DB**(X(2)-TM))
C CHECK REASONABLENESS OF RESULT.
      IF (DZ.LE.0D0) GO TO 30
C LHS SHOULD BE .LE. 0.5 BUT ALLOW FOR SOME ERROR IN LOG
      IF (ABS(DBLE(X(2))-(LOG(DZ)/
     1    LOG(DBLE(B))+0.5D0)).GT.0.6D0) GO TO 30
      IF (X(1).LT.0) DZ = -DZ
      RETURN
C FOLLOWING MESSAGE INDICATES THAT X IS TOO LARGE OR SMALL -
C TRY USING MPCMDE INSTEAD.
   30 WRITE (LUN, 40)
   40 FORMAT (' *** FLOATING-POINT OVER/UNDER-FLOW IN MPCMD ***')
      CALL MPERR
      RETURN
      END
