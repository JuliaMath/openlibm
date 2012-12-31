*DECK SINQI
      SUBROUTINE SINQI (N, WSAVE)
C***BEGIN PROLOGUE  SINQI
C***PURPOSE  Initialize a work array for SINQF and SINQB.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A3
C***TYPE      SINGLE PRECISION (SINQI-S)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine SINQI initializes the array WSAVE which is used in
C  both SINQF and SINQB.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15.
C          The same work array can be used for both SINQF and SINQB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of SINQF or SINQB.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  COSQI
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*)
C   861211  REVISION DATE from Version 3.2
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SINQI
      DIMENSION WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINQI
      CALL COSQI (N,WSAVE)
      RETURN
      END
