*DECK SINTI
      SUBROUTINE SINTI (N, WSAVE)
C***BEGIN PROLOGUE  SINTI
C***PURPOSE  Initialize a work array for SINT.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A3
C***TYPE      SINGLE PRECISION (SINTI-S)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine SINTI initializes the array WSAVE which is used in
C  subroutine SINT.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N+1 is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array with at least INT(3.5*N+16) locations.
C          Different WSAVE arrays are required for different values
C          of N.  The contents of WSAVE must not be changed between
C          calls of SINT.
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
C***END PROLOGUE  SINTI
      DIMENSION WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINTI
      IF (N .LE. 1) RETURN
      PI = 4.*ATAN(1.)
      NP1 = N+1
      NS2 = N/2
      DT = PI/NP1
      KS = N+2
      KF = KS+NS2-1
      FK = 0.
      DO 101 K=KS,KF
         FK = FK+1.
         WSAVE(K) = 2.*SIN(FK*DT)
  101 CONTINUE
      CALL RFFTI (NP1,WSAVE(KF+1))
      RETURN
      END
