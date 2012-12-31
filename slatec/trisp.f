*DECK TRISP
      SUBROUTINE TRISP (N, A, B, C, D, U, Z)
C***BEGIN PROLOGUE  TRISP
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPELI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (TRISP-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine solves for a non-zero eigenvector corresponding
C     to the zero eigenvalue of the transpose of the rank
C     deficient ONE matrix with subdiagonal A, diagonal B, and
C     superdiagonal C , with A(1) in the (1,N) position, with
C     C(N) in the (N,1) position, and all other elements zero.
C
C***SEE ALSO  SEPELI
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  TRISP
C
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,D(*)       ,
     1                U(*)       ,Z(*)
C***FIRST EXECUTABLE STATEMENT  TRISP
      BN = B(N)
      D(1) = A(2)/B(1)
      V = A(1)
      U(1) = C(N)/B(1)
      NM2 = N-2
      DO  10 J=2,NM2
         DEN = B(J)-C(J-1)*D(J-1)
         D(J) = A(J+1)/DEN
         U(J) = -C(J-1)*U(J-1)/DEN
         BN = BN-V*U(J-1)
         V = -V*D(J-1)
   10 CONTINUE
      DEN = B(N-1)-C(N-2)*D(N-2)
      D(N-1) = (A(N)-C(N-2)*U(N-2))/DEN
      AN = C(N-1)-V*D(N-2)
      BN = BN-V*U(N-2)
      DEN = BN-AN*D(N-1)
C
C     SET LAST COMPONENT EQUAL TO ONE
C
      Z(N) = 1.0
      Z(N-1) = -D(N-1)
      NM1 = N-1
      DO  20 J=2,NM1
         K = N-J
         Z(K) = -D(K)*Z(K+1)-U(K)*Z(N)
   20 CONTINUE
      RETURN
      END
