*DECK SINT
      SUBROUTINE SINT (N, X, WSAVE)
C***BEGIN PROLOGUE  SINT
C***PURPOSE  Compute the sine transform of a real, odd sequence.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A3
C***TYPE      SINGLE PRECISION (SINT-S)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine SINT computes the discrete Fourier sine transform
C  of an odd sequence X(I).  The transform is defined below at
C  output parameter X.
C
C  SINT is the unnormalized inverse of itself since a call of SINT
C  followed by another call of SINT will multiply the input sequence
C  X by 2*(N+1).
C
C  The array WSAVE which is used by subroutine SINT must be
C  initialized by calling subroutine SINTI(N,WSAVE).
C
C  Input Parameters
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N+1 is the product of small primes.
C
C  X       an array which contains the sequence to be transformed
C
C
C  WSAVE   a work array with dimension at least INT(3.5*N+16)
C          in the program that calls SINT.  The WSAVE array must be
C          initialized by calling subroutine SINTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C
C  Output Parameters
C
C  X       For I=1,...,N
C
C               X(I)= the sum from K=1 to K=N
C
C                    2*X(K)*SIN(K*I*PI/(N+1))
C
C               A call of SINT followed by another call of
C               SINT will multiply the sequence X by 2*(N+1).
C               Hence SINT is the unnormalized inverse
C               of itself.
C
C  WSAVE   contains initialization calculations which must not be
C          destroyed between calls of SINT.
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
C           (b) changing definition of variable SQRT3 by using
C               FORTRAN intrinsic function SQRT instead of a DATA
C               statement.
C   881128  Modified by Dick Valent to meet prologue standards.
C   891009  Removed unreferenced statement label.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SINT
      DIMENSION X(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINT
      SQRT3 = SQRT(3.)
      IF (N-2) 101,102,103
  101 X(1) = X(1)+X(1)
      RETURN
  102 XH = SQRT3*(X(1)+X(2))
      X(2) = SQRT3*(X(1)-X(2))
      X(1) = XH
      RETURN
  103 NP1 = N+1
      NS2 = N/2
      WSAVE(1) = 0.
      KW = NP1
      DO 104 K=1,NS2
         KW = KW+1
         KC = NP1-K
         T1 = X(K)-X(KC)
         T2 = WSAVE(KW)*(X(K)+X(KC))
         WSAVE(K+1) = T1+T2
         WSAVE(KC+1) = T2-T1
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) WSAVE(NS2+2) = 4.*X(NS2+1)
      NF = NP1+NS2+1
      CALL RFFTF (NP1,WSAVE,WSAVE(NF))
      X(1) = .5*WSAVE(1)
      DO 105 I=3,N,2
         X(I-1) = -WSAVE(I)
         X(I) = X(I-2)+WSAVE(I-1)
  105 CONTINUE
      IF (MODN .NE. 0) RETURN
      X(N) = -WSAVE(N+1)
      RETURN
      END
