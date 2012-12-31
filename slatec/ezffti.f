*DECK EZFFTI
      SUBROUTINE EZFFTI (N, WSAVE)
C***BEGIN PROLOGUE  EZFFTI
C***PURPOSE  Initialize a work array for EZFFTF and EZFFTB.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (EZFFTI-S)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine EZFFTI initializes the work array WSAVE which is used in
C  both EZFFTF and EZFFTB.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15.
C          The same work array can be used for both EZFFTF and EZFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  EZFFT1
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*).
C   861211  REVISION DATE from Version 3.2
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  EZFFTI
      DIMENSION WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  EZFFTI
      IF (N .EQ. 1) RETURN
      CALL EZFFT1 (N,WSAVE(2*N+1),WSAVE(3*N+1))
      RETURN
      END
