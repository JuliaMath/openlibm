*DECK COST
      SUBROUTINE COST (N, X, WSAVE)
C***BEGIN PROLOGUE  COST
C***PURPOSE  Compute the cosine transform of a real, even sequence.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A3
C***TYPE      SINGLE PRECISION (COST-S)
C***KEYWORDS  COSINE FOURIER TRANSFORM, FFTPACK
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine COST computes the discrete Fourier cosine transform
C  of an even sequence X(I).  The transform is defined below at output
C  parameter X.
C
C  COST is the unnormalized inverse of itself since a call of COST
C  followed by another call of COST will multiply the input sequence
C  X by 2*(N-1).  The transform is defined below at output parameter X.
C
C  The array WSAVE which is used by subroutine COST must be
C  initialized by calling subroutine COSTI(N,WSAVE).
C
C  Input Parameters
C
C  N       the length of the sequence X.  N must be greater than 1.
C          The method is most efficient when N-1 is a product of
C          small primes.
C
C  X       an array which contains the sequence to be transformed
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15
C          in the program that calls COST.  The WSAVE array must be
C          initialized by calling subroutine COSTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C
C  Output Parameters
C
C  X       For I=1,...,N
C
C             X(I) = X(1)+(-1)**(I-1)*X(N)
C
C               + the sum from K=2 to K=N-1
C
C                 2*X(K)*COS((K-1)*(I-1)*PI/(N-1))
C
C               A call of COST followed by another call of
C               COST will multiply the sequence X by 2*(N-1).
C               Hence COST is the unnormalized inverse
C               of itself.
C
C  WSAVE   contains initialization calculations which must not be
C          destroyed between calls of COST.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTF
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*)
C   861211  REVISION DATE from Version 3.2
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  COST
      DIMENSION X(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  COST
      NM1 = N-1
      NP1 = N+1
      NS2 = N/2
      IF (N-2) 106,101,102
  101 X1H = X(1)+X(2)
      X(2) = X(1)-X(2)
      X(1) = X1H
      RETURN
  102 IF (N .GT. 3) GO TO 103
      X1P3 = X(1)+X(3)
      TX2 = X(2)+X(2)
      X(2) = X(1)-X(3)
      X(1) = X1P3+TX2
      X(3) = X1P3-TX2
      RETURN
  103 C1 = X(1)-X(N)
      X(1) = X(1)+X(N)
      DO 104 K=2,NS2
         KC = NP1-K
         T1 = X(K)+X(KC)
         T2 = X(K)-X(KC)
         C1 = C1+WSAVE(KC)*T2
         T2 = WSAVE(K)*T2
         X(K) = T1-T2
         X(KC) = T1+T2
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) X(NS2+1) = X(NS2+1)+X(NS2+1)
      CALL RFFTF (NM1,X,WSAVE(N+1))
      XIM2 = X(2)
      X(2) = C1
      DO 105 I=4,N,2
         XI = X(I)
         X(I) = X(I-2)-X(I-1)
         X(I-1) = XIM2
         XIM2 = XI
  105 CONTINUE
      IF (MODN .NE. 0) X(N) = XIM2
  106 RETURN
      END
