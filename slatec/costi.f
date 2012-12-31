*DECK COSTI
      SUBROUTINE COSTI (N, WSAVE)
C***BEGIN PROLOGUE  COSTI
C***PURPOSE  Initialize a work array for COST.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A3
C***TYPE      SINGLE PRECISION (COSTI-S)
C***KEYWORDS  COSINE FOURIER TRANSFORM, FFTPACK
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine COSTI initializes the array WSAVE which is used in
C  subroutine COST.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N-1 is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15.
C          Different WSAVE arrays are required for different values
C          of N.  The contents of WSAVE must not be changed between
C          calls of COST.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTI
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing references to intrinsic function FLOAT
C               to REAL, and
C           (c) changing definition of variable PI by using
C               FORTRAN intrinsic function ATAN instead of a DATA
C               statement.
C   881128  Modified by Dick Valent to meet prologue standards.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  COSTI
      DIMENSION WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  COSTI
      IF (N .LE. 3) RETURN
      PI = 4.*ATAN(1.)
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      DT = PI/NM1
      FK = 0.
      DO 101 K=2,NS2
         KC = NP1-K
         FK = FK+1.
         WSAVE(K) = 2.*SIN(FK*DT)
         WSAVE(KC) = 2.*COS(FK*DT)
  101 CONTINUE
      CALL RFFTI (NM1,WSAVE(N+1))
      RETURN
      END
