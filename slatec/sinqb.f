*DECK SINQB
      SUBROUTINE SINQB (N, X, WSAVE)
C***BEGIN PROLOGUE  SINQB
C***PURPOSE  Compute the unnormalized inverse of SINQF.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A3
C***TYPE      SINGLE PRECISION (SINQB-S)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine SINQB computes the fast Fourier transform of quarter
C  wave data.  That is, SINQB computes a sequence from its
C  representation in terms of a sine series with odd wave numbers.
C  the transform is defined below at output parameter X.
C
C  SINQF is the unnormalized inverse of SINQB since a call of SINQB
C  followed by a call of SINQF will multiply the input sequence X
C  by 4*N.
C
C  The array WSAVE which is used by subroutine SINQB must be
C  initialized by calling subroutine SINQI(N,WSAVE).
C
C  Input Parameters
C
C  N       the length of the array X to be transformed.  The method
C          is most efficient when N is a product of small primes.
C
C  X       an array which contains the sequence to be transformed
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15
C          in the program that calls SINQB.  The WSAVE array must be
C          initialized by calling subroutine SINQI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C
C  Output Parameters
C
C  X       For I=1,...,N
C
C               X(I)= the sum from K=1 to K=N of
C
C                 4*X(K)*SIN((2*K-1)*I*PI/(2*N))
C
C               a call of SINQB followed by a call of
C               SINQF will multiply the sequence X by 4*N.
C               Therefore SINQF is the unnormalized inverse
C               of SINQB.
C
C  WSAVE   contains initialization calculations which must not
C          be destroyed between calls of SINQB or SINQF.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  COSQB
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*).
C   861211  REVISION DATE from Version 3.2
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SINQB
      DIMENSION X(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINQB
      IF (N .GT. 1) GO TO 101
      X(1) = 4.*X(1)
      RETURN
  101 NS2 = N/2
      DO 102 K=2,N,2
         X(K) = -X(K)
  102 CONTINUE
      CALL COSQB (N,X,WSAVE)
      DO 103 K=1,NS2
         KC = N-K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
  103 CONTINUE
      RETURN
      END
