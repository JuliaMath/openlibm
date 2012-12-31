*DECK D1UPDT
      SUBROUTINE D1UPDT (M, N, S, LS, U, V, W, SING)
C***BEGIN PROLOGUE  D1UPDT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DNSQ and DNSQE
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (R1UPDT-S, D1UPDT-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Given an M by N lower trapezoidal matrix S, an M-vector U,
C     and an N-vector V, the problem is to determine an
C     orthogonal matrix Q such that
C
C                   t
C           (S + U*V )*Q
C
C     is again lower trapezoidal.
C
C     This subroutine determines Q as the product of 2*(N - 1)
C     transformations
C
C           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     where GV(I), GW(I) are Givens rotations in the (I,N) plane
C     which eliminate elements in the I-th and N-th planes,
C     respectively. Q itself is not accumulated, rather the
C     information to recover the GV, GW rotations is returned.
C
C     The SUBROUTINE statement is
C
C       SUBROUTINE D1UPDT(M,N,S,LS,U,V,W,SING)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of S.
C
C       N is a positive integer input variable set to the number
C         of columns of S. N must not exceed M.
C
C       S is an array of length LS. On input S must contain the lower
C         trapezoidal matrix S stored by columns. On output S contains
C         the lower trapezoidal matrix produced as described above.
C
C       LS is a positive integer input variable not less than
C         (N*(2*M-N+1))/2.
C
C       U is an input array of length M which must contain the
C         vector U.
C
C       V is an array of length N. On input V must contain the vector
C         V. On output V(I) contains the information necessary to
C         recover the Givens rotation GV(I) described above.
C
C       W is an output array of length M. W(I) contains information
C         necessary to recover the Givens rotation GW(I) described
C         above.
C
C       SING is a LOGICAL output variable. SING is set TRUE if any
C         of the diagonal elements of the output S are zero. Otherwise
C         SING is set FALSE.
C
C***SEE ALSO  DNSQ, DNSQE
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  D1UPDT
      DOUBLE PRECISION D1MACH
      INTEGER I, J, JJ, L, LS, M, N, NM1, NMJ
      DOUBLE PRECISION COS, COTAN, GIANT, ONE, P25, P5, S(*),
     1     SIN, TAN, TAU, TEMP, U(*), V(*), W(*), ZERO
      LOGICAL SING
      SAVE ONE, P5, P25, ZERO
      DATA ONE,P5,P25,ZERO /1.0D0,5.0D-1,2.5D-1,0.0D0/
C
C     GIANT IS THE LARGEST MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  D1UPDT
      GIANT = D1MACH(2)
C
C     INITIALIZE THE DIAGONAL ELEMENT POINTER.
C
      JJ = (N*(2*M - N + 1))/2 - (M - N)
C
C     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.
C
      L = JJ
      DO 10 I = N, M
         W(I) = S(L)
         L = L + 1
   10    CONTINUE
C
C     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
C     IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W.
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 NMJ = 1, NM1
         J = N - NMJ
         JJ = JJ - (M - J + 1)
         W(J) = ZERO
         IF (V(J) .EQ. ZERO) GO TO 50
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C        J-TH ELEMENT OF V.
C
         IF (ABS(V(N)) .GE. ABS(V(J))) GO TO 20
            COTAN = V(N)/V(J)
            SIN = P5/SQRT(P25+P25*COTAN**2)
            COS = SIN*COTAN
            TAU = ONE
            IF (ABS(COS)*GIANT .GT. ONE) TAU = ONE/COS
            GO TO 30
   20    CONTINUE
            TAN = V(J)/V(N)
            COS = P5/SQRT(P25+P25*TAN**2)
            SIN = COS*TAN
            TAU = SIN
   30    CONTINUE
C
C        APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION
C        NECESSARY TO RECOVER THE GIVENS ROTATION.
C
         V(N) = SIN*V(J) + COS*V(N)
         V(J) = TAU
C
C        APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W.
C
         L = JJ
         DO 40 I = J, M
            TEMP = COS*S(L) - SIN*W(I)
            W(I) = SIN*S(L) + COS*W(I)
            S(L) = TEMP
            L = L + 1
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.
C
      DO 80 I = 1, M
         W(I) = W(I) + V(N)*U(I)
   80    CONTINUE
C
C     ELIMINATE THE SPIKE.
C
      SING = .FALSE.
      IF (NM1 .LT. 1) GO TO 140
      DO 130 J = 1, NM1
         IF (W(J) .EQ. ZERO) GO TO 120
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C        J-TH ELEMENT OF THE SPIKE.
C
         IF (ABS(S(JJ)) .GE. ABS(W(J))) GO TO 90
            COTAN = S(JJ)/W(J)
            SIN = P5/SQRT(P25+P25*COTAN**2)
            COS = SIN*COTAN
            TAU = ONE
            IF (ABS(COS)*GIANT .GT. ONE) TAU = ONE/COS
            GO TO 100
   90    CONTINUE
            TAN = W(J)/S(JJ)
            COS = P5/SQRT(P25+P25*TAN**2)
            SIN = COS*TAN
            TAU = SIN
  100    CONTINUE
C
C        APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W.
C
         L = JJ
         DO 110 I = J, M
            TEMP = COS*S(L) + SIN*W(I)
            W(I) = -SIN*S(L) + COS*W(I)
            S(L) = TEMP
            L = L + 1
  110       CONTINUE
C
C        STORE THE INFORMATION NECESSARY TO RECOVER THE
C        GIVENS ROTATION.
C
         W(J) = TAU
  120    CONTINUE
C
C        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.
C
         IF (S(JJ) .EQ. ZERO) SING = .TRUE.
         JJ = JJ + (M - J + 1)
  130    CONTINUE
  140 CONTINUE
C
C     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.
C
      L = JJ
      DO 150 I = N, M
         S(L) = W(I)
         L = L + 1
  150    CONTINUE
      IF (S(JJ) .EQ. ZERO) SING = .TRUE.
      RETURN
C
C     LAST CARD OF SUBROUTINE D1UPDT.
C
      END
