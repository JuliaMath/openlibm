*DECK COSQF
      SUBROUTINE COSQF (N, X, WSAVE)
C***BEGIN PROLOGUE  COSQF
C***PURPOSE  Compute the forward cosine transform with odd wave numbers.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A3
C***TYPE      SINGLE PRECISION (COSQF-S)
C***KEYWORDS  COSINE FOURIER TRANSFORM, FFTPACK
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine COSQF computes the fast Fourier transform of quarter
C  wave data. That is, COSQF computes the coefficients in a cosine
C  series representation with only odd wave numbers.  The transform
C  is defined below at Output Parameter X
C
C  COSQF is the unnormalized inverse of COSQB since a call of COSQF
C  followed by a call of COSQB will multiply the input sequence X
C  by 4*N.
C
C  The array WSAVE which is used by subroutine COSQF must be
C  initialized by calling subroutine COSQI(N,WSAVE).
C
C
C  Input Parameters
C
C  N       the length of the array X to be transformed.  The method
C          is most efficient when N is a product of small primes.
C
C  X       an array which contains the sequence to be transformed
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15
C          in the program that calls COSQF.  The WSAVE array must be
C          initialized by calling subroutine COSQI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C
C  Output Parameters
C
C  X       For I=1,...,N
C
C               X(I) = X(1) plus the sum from K=2 to K=N of
C
C                  2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N))
C
C               A call of COSQF followed by a call of
C               COSQB will multiply the sequence X by 4*N.
C               Therefore COSQB is the unnormalized inverse
C               of COSQF.
C
C  WSAVE   contains initialization calculations which must not
C          be destroyed between calls of COSQF or COSQB.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  COSQF1
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing definition of variable SQRT2 by using
C               FORTRAN intrinsic function SQRT instead of a DATA
C               statement.
C   861211  REVISION DATE from Version 3.2
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  COSQF
      DIMENSION X(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  COSQF
      SQRT2 = SQRT(2.)
      IF (N-2) 102,101,103
  101 TSQX = SQRT2*X(2)
      X(2) = X(1)-TSQX
      X(1) = X(1)+TSQX
  102 RETURN
  103 CALL COSQF1 (N,X,WSAVE,WSAVE(N+1))
      RETURN
      END
