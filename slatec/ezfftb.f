*DECK EZFFTB
      SUBROUTINE EZFFTB (N, R, AZERO, A, B, WSAVE)
C***BEGIN PROLOGUE  EZFFTB
C***PURPOSE  A simplified real, periodic, backward fast Fourier
C            transform.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (EZFFTB-S)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine EZFFTB computes a real periodic sequence from its
C  Fourier coefficients (Fourier synthesis).  The transform is
C  defined below at Output Parameter R.  EZFFTB is a simplified
C  but slower version of RFFTB.
C
C  Input Parameters
C
C  N       the length of the output array R.  The method is most
C          efficient when N is the product of small primes.
C
C  AZERO   the constant Fourier coefficient
C
C  A,B     arrays which contain the remaining Fourier coefficients.
C          These arrays are not destroyed.
C
C          The length of these arrays depends on whether N is even or
C          odd.
C
C          If N is even, N/2    locations are required.
C          If N is odd, (N-1)/2 locations are required
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15
C          in the program that calls EZFFTB.  The WSAVE array must be
C          initialized by calling subroutine EZFFTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C          The same WSAVE array can be used by EZFFTF and EZFFTB.
C
C  Output Parameters
C
C  R       if N is even, define KMAX=N/2
C          if N is odd,  define KMAX=(N-1)/2
C
C          Then for I=1,...,N
C
C               R(I)=AZERO plus the sum from K=1 to K=KMAX of
C
C               A(K)*COS(K*(I-1)*2*PI/N)+B(K)*SIN(K*(I-1)*2*PI/N)
C
C  ********************* Complex Notation **************************
C
C          For J=1,...,N
C
C          R(J) equals the sum from K=-KMAX to K=KMAX of
C
C               C(K)*EXP(I*K*(J-1)*2*PI/N)
C
C          where
C
C               C(K) = .5*CMPLX(A(K),-B(K))   for K=1,...,KMAX
C
C               C(-K) = CONJG(C(K))
C
C               C(0) = AZERO
C
C                    and I=SQRT(-1)
C
C  *************** Amplitude - Phase Notation ***********************
C
C          For I=1,...,N
C
C          R(I) equals AZERO plus the sum from K=1 to K=KMAX of
C
C               ALPHA(K)*COS(K*(I-1)*2*PI/N+BETA(K))
C
C          where
C
C               ALPHA(K) = SQRT(A(K)*A(K)+B(K)*B(K))
C
C               COS(BETA(K))=A(K)/ALPHA(K)
C
C               SIN(BETA(K))=-B(K)/ALPHA(K)
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTB
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*)
C   861211  REVISION DATE from Version 3.2
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  EZFFTB
      DIMENSION R(*), A(*), B(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  EZFFTB
      IF (N-2) 101,102,103
  101 R(1) = AZERO
      RETURN
  102 R(1) = AZERO+A(1)
      R(2) = AZERO-A(1)
      RETURN
  103 NS2 = (N-1)/2
      DO 104 I=1,NS2
         R(2*I) = .5*A(I)
         R(2*I+1) = -.5*B(I)
  104 CONTINUE
      R(1) = AZERO
      IF (MOD(N,2) .EQ. 0) R(N) = A(NS2+1)
      CALL RFFTB (N,R,WSAVE(N+1))
      RETURN
      END
