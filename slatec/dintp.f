*DECK DINTP
      SUBROUTINE DINTP (X, Y, XOUT, YOUT, YPOUT, NEQN, KOLD, PHI, IVC,
     +   IV, KGI, GI, ALPHA, OG, OW, OX, OY)
C***BEGIN PROLOGUE  DINTP
C***PURPOSE  Approximate the solution at XOUT by evaluating the
C            polynomial computed in DSTEPS at XOUT.  Must be used in
C            conjunction with DSTEPS.
C***LIBRARY   SLATEC (DEPAC)
C***CATEGORY  I1A1B
C***TYPE      DOUBLE PRECISION (SINTRP-S, DINTP-D)
C***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
C             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR,
C             SMOOTH INTERPOLANT
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   The methods in subroutine  DSTEPS  approximate the solution near  X
C   by a polynomial.  Subroutine  DINTP  approximates the solution at
C   XOUT  by evaluating the polynomial there.  Information defining this
C   polynomial is passed from  DSTEPS  so  DINTP  cannot be used alone.
C
C   Subroutine DSTEPS is completely explained and documented in the text
C   "Computer Solution of Ordinary Differential Equations, the Initial
C   Value Problem"  by L. F. Shampine and M. K. Gordon.
C
C   Input to DINTP --
C
C   The user provides storage in the calling program for the arrays in
C   the call list
C      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),OY(NEQN)
C                AND ALPHA(12),OG(13),OW(12),GI(11),IV(10)
C   and defines
C      XOUT -- point at which solution is desired.
C   The remaining parameters are defined in  DSTEPS  and passed to
C   DINTP  from that subroutine
C
C   Output from  DINTP --
C
C      YOUT(*) -- solution at  XOUT
C      YPOUT(*) -- derivative of solution at  XOUT
C   The remaining parameters are returned unaltered from their input
C   values.  Integration with  DSTEPS  may be continued.
C
C***REFERENCES  H. A. Watts, A smoother interpolant for DE/STEP, INTRP
C                 II, Report SAND84-0293, Sandia Laboratories, 1984.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   840201  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DINTP
C
      INTEGER I, IQ, IV, IVC, IW, J, JQ, KGI, KOLD, KP1, KP2,
     1        L, M, NEQN
      DOUBLE PRECISION ALP, ALPHA, C, G, GDI, GDIF, GI, GAMMA, H, HI,
     1       HMU, OG, OW, OX, OY, PHI, RMU, SIGMA, TEMP1, TEMP2, TEMP3,
     2       W, X, XI, XIM1, XIQ, XOUT, Y, YOUT, YPOUT
C
      DIMENSION Y(*),YOUT(*),YPOUT(*),PHI(NEQN,16),OY(*)
      DIMENSION G(13),C(13),W(13),OG(13),OW(12),ALPHA(12),GI(11),IV(10)
C
C***FIRST EXECUTABLE STATEMENT  DINTP
      KP1 = KOLD + 1
      KP2 = KOLD + 2
C
      HI = XOUT - OX
      H = X - OX
      XI = HI/H
      XIM1 = XI - 1.D0
C
C   INITIALIZE W(*) FOR COMPUTING G(*)
C
      XIQ = XI
      DO 10 IQ = 1,KP1
        XIQ = XI*XIQ
        TEMP1 = IQ*(IQ+1)
 10     W(IQ) = XIQ/TEMP1
C
C   COMPUTE THE DOUBLE INTEGRAL TERM GDI
C
      IF (KOLD .LE. KGI) GO TO 50
      IF (IVC .GT. 0) GO TO 20
      GDI = 1.0D0/TEMP1
      M = 2
      GO TO 30
 20   IW = IV(IVC)
      GDI = OW(IW)
      M = KOLD - IW + 3
 30   IF (M .GT. KOLD) GO TO 60
      DO 40 I = M,KOLD
 40     GDI = OW(KP2-I) - ALPHA(I)*GDI
      GO TO 60
 50   GDI = GI(KOLD)
C
C   COMPUTE G(*) AND C(*)
C
 60   G(1) = XI
      G(2) = 0.5D0*XI*XI
      C(1) = 1.0D0
      C(2) = XI
      IF (KOLD .LT. 2) GO TO 90
      DO 80 I = 2,KOLD
        ALP = ALPHA(I)
        GAMMA = 1.0D0 + XIM1*ALP
        L = KP2 - I
        DO 70 JQ = 1,L
 70       W(JQ) = GAMMA*W(JQ) - ALP*W(JQ+1)
        G(I+1) = W(1)
 80     C(I+1) = GAMMA*C(I)
C
C   DEFINE INTERPOLATION PARAMETERS
C
 90   SIGMA = (W(2) - XIM1*W(1))/GDI
      RMU = XIM1*C(KP1)/GDI
      HMU = RMU/H
C
C   INTERPOLATE FOR THE SOLUTION -- YOUT
C   AND FOR THE DERIVATIVE OF THE SOLUTION -- YPOUT
C
      DO 100 L = 1,NEQN
        YOUT(L) = 0.0D0
 100    YPOUT(L) = 0.0D0
      DO 120 J = 1,KOLD
        I = KP2 - J
        GDIF = OG(I) - OG(I-1)
        TEMP2 = (G(I) - G(I-1)) - SIGMA*GDIF
        TEMP3 = (C(I) - C(I-1)) + RMU*GDIF
        DO 110 L = 1,NEQN
          YOUT(L) = YOUT(L) + TEMP2*PHI(L,I)
 110      YPOUT(L) = YPOUT(L) + TEMP3*PHI(L,I)
 120    CONTINUE
      DO 130 L = 1,NEQN
        YOUT(L) = ((1.0D0 - SIGMA)*OY(L) + SIGMA*Y(L)) +
     1             H*(YOUT(L) + (G(1) - SIGMA*OG(1))*PHI(L,1))
 130    YPOUT(L) = HMU*(OY(L) - Y(L)) +
     1                (YPOUT(L) + (C(1) + RMU*OG(1))*PHI(L,1))
C
      RETURN
      END
