*DECK EZFFTF
      SUBROUTINE EZFFTF (N, R, AZERO, A, B, WSAVE)
C***BEGIN PROLOGUE  EZFFTF
C***PURPOSE  Compute a simplified real, periodic, fast Fourier forward
C            transform.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (EZFFTF-S)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine EZFFTF computes the Fourier coefficients of a real
C  periodic sequence (Fourier analysis).  The transform is defined
C  below at Output Parameters AZERO, A and B.  EZFFTF is a simplified
C  but slower version of RFFTF.
C
C  Input Parameters
C
C  N       the length of the array R to be transformed.  The method
C          is most efficient when N is the product of small primes.
C
C  R       a real array of length N which contains the sequence
C          to be transformed.  R is not destroyed.
C
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15
C          in the program that calls EZFFTF.  The WSAVE array must be
C          initialized by calling subroutine EZFFTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C          The same WSAVE array can be used by EZFFTF and EZFFTB.
C
C  Output Parameters
C
C  AZERO   the sum from I=1 to I=N of R(I)/N
C
C  A,B     for N even B(N/2)=0. and A(N/2) is the sum from I=1 to
C          I=N of (-1)**(I-1)*R(I)/N
C
C          for N even define KMAX=N/2-1
C          for N odd  define KMAX=(N-1)/2
C
C          then for  K=1,...,KMAX
C
C               A(K) equals the sum from I=1 to I=N of
C
C                    2./N*R(I)*COS(K*(I-1)*2*PI/N)
C
C               B(K) equals the sum from I=1 to I=N of
C
C                    2./N*R(I)*SIN(K*(I-1)*2*PI/N)
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTF
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing references to intrinsic function FLOAT
C               to REAL.
C   881128  Modified by Dick Valent to meet prologue standards.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  EZFFTF
      DIMENSION R(*), A(*), B(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  EZFFTF
      IF (N-2) 101,102,103
  101 AZERO = R(1)
      RETURN
  102 AZERO = .5*(R(1)+R(2))
      A(1) = .5*(R(1)-R(2))
      RETURN
  103 DO 104 I=1,N
         WSAVE(I) = R(I)
  104 CONTINUE
      CALL RFFTF (N,WSAVE,WSAVE(N+1))
      CF = 2./N
      CFM = -CF
      AZERO = .5*CF*WSAVE(1)
      NS2 = (N+1)/2
      NS2M = NS2-1
      DO 105 I=1,NS2M
         A(I) = CF*WSAVE(2*I)
         B(I) = CFM*WSAVE(2*I+1)
  105 CONTINUE
      IF (MOD(N,2) .EQ. 0) A(NS2) = .5*CF*WSAVE(N)
      RETURN
      END
