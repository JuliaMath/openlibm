*DECK CAIRY
      SUBROUTINE CAIRY (Z, ID, KODE, AI, NZ, IERR)
C***BEGIN PROLOGUE  CAIRY
C***PURPOSE  Compute the Airy function Ai(z) or its derivative dAi/dz
C            for complex argument z.  A scaling option is available
C            to help avoid underflow and overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10D
C***TYPE      COMPLEX (CAIRY-C, ZAIRY-C)
C***KEYWORDS  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD,
C             BESSEL FUNCTION OF ORDER TWO THIRDS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C         On KODE=1, CAIRY computes the complex Airy function Ai(z)
C         or its derivative dAi/dz on ID=0 or ID=1 respectively. On
C         KODE=2, a scaling option exp(zeta)*Ai(z) or exp(zeta)*dAi/dz
C         is provided to remove the exponential decay in -pi/3<arg(z)
C         <pi/3 and the exponential growth in pi/3<abs(arg(z))<pi where
C         zeta=(2/3)*z**(3/2).
C
C         While the Airy functions Ai(z) and dAi/dz are analytic in
C         the whole z-plane, the corresponding scaled functions defined
C         for KODE=2 have a cut along the negative real axis.
C
C         Input
C           Z      - Argument of type COMPLEX
C           ID     - Order of derivative, ID=0 or ID=1
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            AI=Ai(z)  on ID=0
C                            AI=dAi/dz on ID=1
C                            at z=Z
C                        =2  returns
C                            AI=exp(zeta)*Ai(z)  on ID=0
C                            AI=exp(zeta)*dAi/dz on ID=1
C                            at z=Z where zeta=(2/3)*z**(3/2)
C
C         Output
C           AI     - Result of type COMPLEX
C           NZ     - Underflow indicator
C                    NZ=0    Normal return
C                    NZ=1    AI=0 due to underflow in
C                            -pi/3<arg(Z)<pi/3 on KODE=1
C           IERR   - Error flag
C                    IERR=0  Normal return     - COMPUTATION COMPLETED
C                    IERR=1  Input error       - NO COMPUTATION
C                    IERR=2  Overflow          - NO COMPUTATION
C                            (Re(Z) too large with KODE=1)
C                    IERR=3  Precision warning - COMPUTATION COMPLETED
C                            (Result has less than half precision)
C                    IERR=4  Precision error   - NO COMPUTATION
C                            (Result has no precision)
C                    IERR=5  Algorithmic error - NO COMPUTATION
C                            (Termination condition not met)
C
C *Long Description:
C
C         Ai(z) and dAi/dz are computed from K Bessel functions by
C
C                Ai(z) =  c*sqrt(z)*K(1/3,zeta)
C               dAi/dz = -c*   z   *K(2/3,zeta)
C                    c =  1/(pi*sqrt(3))
C                 zeta =  (2/3)*z**(3/2)
C
C         when abs(z)>1 and from power series when abs(z)<=1.
C
C         In most complex variable computation, one must evaluate ele-
C         mentary functions.  When the magnitude of Z is large, losses
C         of significance by argument reduction occur.  Consequently, if
C         the magnitude of ZETA=(2/3)*Z**(3/2) exceeds U1=SQRT(0.5/UR),
C         then losses exceeding half precision are likely and an error
C         flag IERR=3 is triggered where UR=R1MACH(4)=UNIT ROUNDOFF.
C         Also, if the magnitude of ZETA is larger than U2=0.5/UR, then
C         all significance is lost and IERR=4.  In order to use the INT
C         function, ZETA must be further restricted not to exceed
C         U3=I1MACH(9)=LARGEST INTEGER.  Thus, the magnitude of ZETA
C         must be restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2,
C         and U3 are approximately 2.0E+3, 4.2E+6, 2.1E+9 in single
C         precision and 4.7E+7, 2.3E+15, 2.1E+9 in double precision.
C         This makes U2 limiting is single precision and U3 limiting
C         in double precision.  This means that the magnitude of Z
C         cannot exceed approximately 3.4E+4 in single precision and
C         2.1E+6 in double precision.  This also means that one can
C         expect to retain, in the worst cases on 32-bit machines,
C         no digits in single precision and only 6 digits in double
C         precision.
C
C         The approximate relative error in the magnitude of a complex
C         Bessel function can be expressed as P*10**S where P=MAX(UNIT
C         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
C         sents the increase in error due to argument reduction in the
C         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
C         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
C         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
C         have only absolute accuracy.  This is most likely to occur
C         when one component (in magnitude) is larger than the other by
C         several orders of magnitude.  If one component is 10**K larger
C         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
C         0) significant digits; or, stated another way, when K exceeds
C         the exponent of P, no significant digits remain in the smaller
C         component.  However, the phase angle retains absolute accuracy
C         because, in complex arithmetic with precision P, the smaller
C         component will not (as a rule) decrease below P times the
C         magnitude of the larger component. In these extreme cases,
C         the principal phase angle is on the order of +P, -P, PI/2-P,
C         or -PI/2+P.
C
C***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
C                 matical Functions, National Bureau of Standards
C                 Applied Mathematics Series 55, U. S. Department
C                 of Commerce, Tenth Printing (1972) or later.
C               2. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument and Large Order, Report SAND83-0643,
C                 Sandia National Laboratories, Albuquerque, NM, May
C                 1983.
C               3. D. E. Amos, A Subroutine Package for Bessel Functions
C                 of a Complex Argument and Nonnegative Order, Report
C                 SAND85-1018, Sandia National Laboratory, Albuquerque,
C                 NM, May 1985.
C               4. D. E. Amos, A portable package for Bessel functions
C                 of a complex argument and nonnegative order, ACM
C                 Transactions on Mathematical Software, 12 (September
C                 1986), pp. 265-273.
C
C***ROUTINES CALLED  CACAI, CBKNU, I1MACH, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C***END PROLOGUE  CAIRY
      COMPLEX AI, CONE, CSQ, CY, S1, S2, TRM1, TRM2, Z, ZTA, Z3
      REAL AA, AD, AK, ALIM, ATRM, AZ, AZ3, BK, CK, COEF, C1, C2, DIG,
     * DK, D1, D2, ELIM, FID, FNU, RL, R1M5, SFAC, TOL, TTH, ZI, ZR,
     * Z3I, Z3R, R1MACH, BB, ALAZ
      INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
      DIMENSION CY(1)
      DATA TTH, C1, C2, COEF /6.66666666666666667E-01,
     * 3.55028053887817240E-01,2.58819403792806799E-01,
     * 1.83776298473930683E-01/
      DATA  CONE / (1.0E0,0.0E0) /
C***FIRST EXECUTABLE STATEMENT  CAIRY
      IERR = 0
      NZ=0
      IF (ID.LT.0 .OR. ID.GT.1) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (IERR.NE.0) RETURN
      AZ = ABS(Z)
      TOL = MAX(R1MACH(4),1.0E-18)
      FID = ID
      IF (AZ.GT.1.0E0) GO TO 60
C-----------------------------------------------------------------------
C     POWER SERIES FOR ABS(Z).LE.1.
C-----------------------------------------------------------------------
      S1 = CONE
      S2 = CONE
      IF (AZ.LT.TOL) GO TO 160
      AA = AZ*AZ
      IF (AA.LT.TOL/AZ) GO TO 40
      TRM1 = CONE
      TRM2 = CONE
      ATRM = 1.0E0
      Z3 = Z*Z*Z
      AZ3 = AZ*AA
      AK = 2.0E0 + FID
      BK = 3.0E0 - FID - FID
      CK = 4.0E0 - FID
      DK = 3.0E0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = MIN(D1,D2)
      AK = 24.0E0 + 9.0E0*FID
      BK = 30.0E0 - 9.0E0*FID
      Z3R = REAL(Z3)
      Z3I = AIMAG(Z3)
      DO 30 K=1,25
        TRM1 = TRM1*CMPLX(Z3R/D1,Z3I/D1)
        S1 = S1 + TRM1
        TRM2 = TRM2*CMPLX(Z3R/D2,Z3I/D2)
        S2 = S2 + TRM2
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = MIN(D1,D2)
        IF (ATRM.LT.TOL*AD) GO TO 40
        AK = AK + 18.0E0
        BK = BK + 18.0E0
   30 CONTINUE
   40 CONTINUE
      IF (ID.EQ.1) GO TO 50
      AI = S1*CMPLX(C1,0.0E0) - Z*S2*CMPLX(C2,0.0E0)
      IF (KODE.EQ.1) RETURN
      ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
      AI = AI*CEXP(ZTA)
      RETURN
   50 CONTINUE
      AI = -S2*CMPLX(C2,0.0E0)
      IF (AZ.GT.TOL) AI = AI + Z*Z*S1*CMPLX(C1/(1.0E0+FID),0.0E0)
      IF (KODE.EQ.1) RETURN
      ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
      AI = AI*CEXP(ZTA)
      RETURN
C-----------------------------------------------------------------------
C     CASE FOR ABS(Z).GT.1.0
C-----------------------------------------------------------------------
   60 CONTINUE
      FNU = (1.0E0+FID)/3.0E0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C-----------------------------------------------------------------------
      K1 = I1MACH(12)
      K2 = I1MACH(13)
      R1M5 = R1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303E0*(K*R1M5-3.0E0)
      K1 = I1MACH(11) - 1
      AA = R1M5*K1
      DIG = MIN(AA,18.0E0)
      AA = AA*2.303E0
      ALIM = ELIM + MAX(-AA,-41.45E0)
      RL = 1.2E0*DIG + 3.0E0
      ALAZ=ALOG(AZ)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA=0.5E0/TOL
      BB=I1MACH(9)*0.5E0
      AA=MIN(AA,BB)
      AA=AA**TTH
      IF (AZ.GT.AA) GO TO 260
      AA=SQRT(AA)
      IF (AZ.GT.AA) IERR=3
      CSQ=CSQRT(Z)
      ZTA=Z*CSQ*CMPLX(TTH,0.0E0)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
      IFLAG = 0
      SFAC = 1.0E0
      ZI = AIMAG(Z)
      ZR = REAL(Z)
      AK = AIMAG(ZTA)
      IF (ZR.GE.0.0E0) GO TO 70
      BK = REAL(ZTA)
      CK = -ABS(BK)
      ZTA = CMPLX(CK,AK)
   70 CONTINUE
      IF (ZI.NE.0.0E0) GO TO 80
      IF (ZR.GT.0.0E0) GO TO 80
      ZTA = CMPLX(0.0E0,AK)
   80 CONTINUE
      AA = REAL(ZTA)
      IF (AA.GE.0.0E0 .AND. ZR.GT.0.0E0) GO TO 100
      IF (KODE.EQ.2) GO TO 90
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.GT.(-ALIM)) GO TO 90
      AA = -AA + 0.25E0*ALAZ
      IFLAG = 1
      SFAC = TOL
      IF (AA.GT.ELIM) GO TO 240
   90 CONTINUE
C-----------------------------------------------------------------------
C     CBKNU AND CACAI RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
C-----------------------------------------------------------------------
      MR = 1
      IF (ZI.LT.0.0E0) MR = -1
      CALL CACAI(ZTA, FNU, KODE, MR, 1, CY, NN, RL, TOL, ELIM, ALIM)
      IF (NN.LT.0) GO TO 250
      NZ = NZ + NN
      GO TO 120
  100 CONTINUE
      IF (KODE.EQ.2) GO TO 110
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.LT.ALIM) GO TO 110
      AA = -AA - 0.25E0*ALAZ
      IFLAG = 2
      SFAC = 1.0E0/TOL
      IF (AA.LT.(-ELIM)) GO TO 180
  110 CONTINUE
      CALL CBKNU(ZTA, FNU, KODE, 1, CY, NZ, TOL, ELIM, ALIM)
  120 CONTINUE
      S1 = CY(1)*CMPLX(COEF,0.0E0)
      IF (IFLAG.NE.0) GO TO 140
      IF (ID.EQ.1) GO TO 130
      AI = CSQ*S1
      RETURN
  130 AI = -Z*S1
      RETURN
  140 CONTINUE
      S1 = S1*CMPLX(SFAC,0.0E0)
      IF (ID.EQ.1) GO TO 150
      S1 = S1*CSQ
      AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
      RETURN
  150 CONTINUE
      S1 = -S1*Z
      AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
      RETURN
  160 CONTINUE
      AA = 1.0E+3*R1MACH(1)
      S1 = CMPLX(0.0E0,0.0E0)
      IF (ID.EQ.1) GO TO 170
      IF (AZ.GT.AA) S1 = CMPLX(C2,0.0E0)*Z
      AI = CMPLX(C1,0.0E0) - S1
      RETURN
  170 CONTINUE
      AI = -CMPLX(C2,0.0E0)
      AA = SQRT(AA)
      IF (AZ.GT.AA) S1 = Z*Z*CMPLX(0.5E0,0.0E0)
      AI = AI + S1*CMPLX(C1,0.0E0)
      RETURN
  180 CONTINUE
      NZ = 1
      AI = CMPLX(0.0E0,0.0E0)
      RETURN
  240 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  250 CONTINUE
      IF(NN.EQ.(-1)) GO TO 240
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END
