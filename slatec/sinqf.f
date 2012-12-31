*DECK SINQF
      SUBROUTINE SINQF (N, X, WSAVE)
C***BEGIN PROLOGUE  SINQF
C***PURPOSE  Compute the forward sine transform with odd wave numbers.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A3
C***TYPE      SINGLE PRECISION (SINQF-S)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine SINQF computes the fast Fourier transform of quarter
C  wave data.  That is, SINQF computes the coefficients in a sine
C  series representation with only odd wave numbers.  The transform
C  is defined below at output parameter X.
C
C  SINQB is the unnormalized inverse of SINQF since a call of SINQF
C  followed by a call of SINQB will multiply the input sequence X
C  by 4*N.
C
C  The array WSAVE which is used by subroutine SINQF must be
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
C          in the program that calls SINQF.  The WSAVE array must be
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
C               X(I) = (-1)**(I-1)*X(N)
C
C                  + the sum from K=1 to K=N-1 of
C
C                  2*X(K)*SIN((2*I-1)*K*PI/(2*N))
C
C               A call of SINQF followed by a call of
C               SINQB will multiply the sequence X by 4*N.
C               Therefore SINQB is the unnormalized inverse
C               of SINQF.
C
C  WSAVE   contains initialization calculations which must not
C          be destroyed between calls of SINQF or SINQB.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  COSQF
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*)
C   861211  REVISION DATE from Version 3.2
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SINQF
      DIMENSION X(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINQF
      IF (N .EQ. 1) RETURN
      NS2 = N/2
      DO 101 K=1,NS2
         KC = N-K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
  101 CONTINUE
      CALL COSQF (N,X,WSAVE)
      DO 102 K=2,N,2
         X(K) = -X(K)
  102 CONTINUE
      RETURN
      END
