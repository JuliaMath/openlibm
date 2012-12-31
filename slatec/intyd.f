*DECK INTYD
      SUBROUTINE INTYD (T, K, YH, NYH, DKY, IFLAG)
C***BEGIN PROLOGUE  INTYD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DEBDF
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (INTYD-S, DINTYD-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   INTYD approximates the solution and derivatives at T by polynomial
C   interpolation. Must be used in conjunction with the integrator
C   package DEBDF.
C ----------------------------------------------------------------------
C INTYD computes interpolated values of the K-th derivative of the
C dependent variable vector Y, and stores it in DKY.
C This routine is called by DEBDF with K = 0,1 and T = TOUT, but may
C also be called by the user for any K up to the current order.
C (see detailed instructions in LSODE usage documentation.)
C ----------------------------------------------------------------------
C The computed values in DKY are gotten by interpolation using the
C Nordsieck history array YH.  This array corresponds uniquely to a
C vector-valued polynomial of degree NQCUR or less, and DKY is set
C to the K-th derivative of this polynomial at T.
C The formula for DKY is..
C              Q
C  DKY(I)  =  sum  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1)
C             J=K
C where  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR.
C The quantities  NQ = NQCUR, L = NQ+1, N = NEQ, TN, and H are
C communicated by common.  The above sum is done in reverse order.
C IFLAG is returned negative if either K or T is out of bounds.
C ----------------------------------------------------------------------
C
C***SEE ALSO  DEBDF
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DEBDF1
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  INTYD
C
CLLL. OPTIMIZE
      INTEGER K, NYH, IFLAG, I, IC, IER, IOWND, IOWNS, J, JB, JB2,
     1   JJ, JJ1, JP1, JSTART, KFLAG, L, MAXORD, METH, MITER, N, NFE,
     2   NJE, NQ, NQU, NST
      REAL T, YH, DKY,
     1   ROWND, ROWNS, EL0, H, HMIN, HMXI, HU, TN, UROUND,
     2   C, R, S, TP
      DIMENSION YH(NYH,*), DKY(*)
      COMMON /DEBDF1/ ROWND, ROWNS(210),
     1   EL0, H, HMIN, HMXI, HU, TN, UROUND, IOWND(14), IOWNS(6),
     2   IER, JSTART, KFLAG, L, METH, MITER, MAXORD, N, NQ, NST, NFE,
     3   NJE, NQU
C
C***FIRST EXECUTABLE STATEMENT  INTYD
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
      IF ((T-TP)*(T-TN) .GT. 0.0E0) GO TO 90
C
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1,NQ
 10     IC = IC*JJ
 15   C = IC
      DO 20 I = 1,N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1,JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1,J
 30       IC = IC*JJ
 35     C = IC
        DO 40 I = 1,N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      DO 60 I = 1,N
 60     DKY(I) = R*DKY(I)
      RETURN
C
 80   IFLAG = -1
      RETURN
 90   IFLAG = -2
      RETURN
C----------------------- END OF SUBROUTINE INTYD -----------------------
      END
