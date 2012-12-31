*DECK DINTYD
      SUBROUTINE DINTYD (T, K, YH, NYH, DKY, IFLAG)
C***BEGIN PROLOGUE  DINTYD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (INTYD-S, DINTYD-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   DINTYD approximates the solution and derivatives at T by polynomial
C   interpolation. Must be used in conjunction with the integrator
C   package DDEBDF.
C ----------------------------------------------------------------------
C DINTYD computes interpolated values of the K-th derivative of the
C dependent variable vector Y, and stores it in DKY.
C This routine is called by DDEBDF with K = 0,1 and T = TOUT, but may
C also be called by the user for any K up to the current order.
C (see detailed instructions in LSODE usage documentation.)
C ----------------------------------------------------------------------
C The computed values in DKY are gotten by interpolation using the
C Nordsieck history array YH.  This array corresponds uniquely to a
C vector-valued polynomial of degree NQCUR or less, and DKY is set
C to the K-th derivative of this polynomial at T.
C The formula for DKY is..
C              Q
C  DKY(I)  =  Sum  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1)
C             J=K
C where  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR.
C The quantities  NQ = NQCUR, L = NQ+1, N = NEQ, TN, and H are
C communicated by common.  The above sum is done in reverse order.
C IFLAG is returned negative if either K or T is out of bounds.
C ----------------------------------------------------------------------
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DINTYD
C
      INTEGER I, IC, IER, IFLAG, IOWND, IOWNS, J, JB, JB2, JJ, JJ1,
     1      JP1, JSTART, K, KFLAG, L, MAXORD, METH, MITER, N, NFE,
     2      NJE, NQ, NQU, NST, NYH
      DOUBLE PRECISION C, DKY, EL0, H, HMIN, HMXI, HU, R, ROWND,
     1      ROWNS, S, T, TN, TP, UROUND, YH
      DIMENSION YH(NYH,*),DKY(*)
      COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND,
     1                IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER,
     2                MAXORD,N,NQ,NST,NFE,NJE,NQU
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 130
C***FIRST EXECUTABLE STATEMENT  DINTYD
         IFLAG = 0
         IF (K .LT. 0 .OR. K .GT. NQ) GO TO 110
            TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
            IF ((T - TP)*(T - TN) .LE. 0.0D0) GO TO 10
               IFLAG = -2
C     .........EXIT
               GO TO 130
   10       CONTINUE
C
            S = (T - TN)/H
            IC = 1
            IF (K .EQ. 0) GO TO 30
               JJ1 = L - K
               DO 20 JJ = JJ1, NQ
                  IC = IC*JJ
   20          CONTINUE
   30       CONTINUE
            C = IC
            DO 40 I = 1, N
               DKY(I) = C*YH(I,L)
   40       CONTINUE
            IF (K .EQ. NQ) GO TO 90
               JB2 = NQ - K
               DO 80 JB = 1, JB2
                  J = NQ - JB
                  JP1 = J + 1
                  IC = 1
                  IF (K .EQ. 0) GO TO 60
                     JJ1 = JP1 - K
                     DO 50 JJ = JJ1, J
                        IC = IC*JJ
   50                CONTINUE
   60             CONTINUE
                  C = IC
                  DO 70 I = 1, N
                     DKY(I) = C*YH(I,JP1) + S*DKY(I)
   70             CONTINUE
   80          CONTINUE
C     .........EXIT
               IF (K .EQ. 0) GO TO 130
   90       CONTINUE
            R = H**(-K)
            DO 100 I = 1, N
               DKY(I) = R*DKY(I)
  100       CONTINUE
         GO TO 120
  110    CONTINUE
C
            IFLAG = -1
  120    CONTINUE
  130 CONTINUE
      RETURN
C     ----------------------- END OF SUBROUTINE DINTYD
C     -----------------------
      END
