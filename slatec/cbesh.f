*DECK CBESH
      SUBROUTINE CBESH (Z, FNU, KODE, M, N, CY, NZ, IERR)
C***BEGIN PROLOGUE  CBESH
C***PURPOSE  Compute a sequence of the Hankel functions H(m,a,z)
C            for superscript m=1 or 2, real nonnegative orders a=b,
C            b+1,... where b>0, and nonzero complex argument z.  A
C            scaling option is available to help avoid overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10A4
C***TYPE      COMPLEX (CBESH-C, ZBESH-C)
C***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF THE THIRD KIND, H BESSEL FUNCTIONS,
C             HANKEL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C         On KODE=1, CBESH computes an N member sequence of complex
C         Hankel (Bessel) functions CY(L)=H(M,FNU+L-1,Z) for super-
C         script M=1 or 2, real nonnegative orders FNU+L-1, L=1,...,
C         N, and complex nonzero Z in the cut plane -pi<arg(Z)<=pi.
C         On KODE=2, CBESH returns the scaled functions
C
C            CY(L) = H(M,FNU+L-1,Z)*exp(-(3-2*M)*Z*i),  i**2=-1
C
C         which removes the exponential behavior in both the upper
C         and lower half planes.  Definitions and notation are found
C         in the NBS Handbook of Mathematical Functions (Ref. 1).
C
C         Input
C           Z      - Nonzero argument of type COMPLEX
C           FNU    - Initial order of type REAL, FNU>=0
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            CY(L)=H(M,FNU+L-1,Z), L=1,...,N
C                        =2  returns
C                            CY(L)=H(M,FNU+L-1,Z)*exp(-(3-2M)*Z*i),
C                            L=1,...,N
C           M      - Superscript of Hankel function, M=1 or 2
C           N      - Number of terms in the sequence, N>=1
C
C         Output
C           CY     - Result vector of type COMPLEX
C           NZ     - Number of underflows set to zero
C                    NZ=0    Normal return
C                    NZ>0    CY(L)=0 for NZ values of L (if M=1 and
C                            Im(Z)>0 or if M=2 and Im(Z)<0, then
C                            CY(L)=0 for L=1,...,NZ; in the com-
C                            plementary half planes, the underflows
C                            may not be in an uninterrupted sequence)
C           IERR   - Error flag
C                    IERR=0  Normal return     - COMPUTATION COMPLETED
C                    IERR=1  Input error       - NO COMPUTATION
C                    IERR=2  Overflow          - NO COMPUTATION
C                            (abs(Z) too small and/or FNU+N-1
C                            too large)
C                    IERR=3  Precision warning - COMPUTATION COMPLETED
C                            (Result has half precision or less
C                            because abs(Z) or FNU+N-1 is large)
C                    IERR=4  Precision error   - NO COMPUTATION
C                            (Result has no precision because
C                            abs(Z) or FNU+N-1 is too large)
C                    IERR=5  Algorithmic error - NO COMPUTATION
C                            (Termination condition not met)
C
C *Long Description:
C
C         The computation is carried out by the formula
C
C            H(m,a,z) = (1/t)*exp(-a*t)*K(a,z*exp(-t))
C                   t = (3-2*m)*i*pi/2
C
C         where the K Bessel function is computed as described in the
C         prologue to CBESK.
C
C         Exponential decay of H(m,a,z) occurs in the upper half z
C         plane for m=1 and the lower half z plane for m=2.  Exponential
C         growth occurs in the complementary half planes.  Scaling
C         by exp(-(3-2*m)*z*i) removes the exponential behavior in the
C         whole z plane as z goes to infinity.
C
C         For negative orders, the formula
C
C            H(m,-a,z) = H(m,a,z)*exp((3-2*m)*a*pi*i)
C
C         can be used.
C
C         In most complex variable computation, one must evaluate ele-
C         mentary functions.  When the magnitude of Z or FNU+N-1 is
C         large, losses of significance by argument reduction occur.
C         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
C         losses exceeding half precision are likely and an error flag
C         IERR=3 is triggered where UR=R1MACH(4)=UNIT ROUNDOFF.  Also,
C         if either is larger than U2=0.5/UR, then all significance is
C         lost and IERR=4.  In order to use the INT function, arguments
C         must be further restricted not to exceed the largest machine
C         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1
C         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and
C         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision
C         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This
C         makes U2 limiting in single precision and U3 limiting in
C         double precision.  This means that one can expect to retain,
C         in the worst cases on IEEE machines, no digits in single pre-
C         cision and only 6 digits in double precision.  Similar con-
C         siderations hold for other machines.
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
C         magnitude of the larger component.  In these extreme cases,
C         the principal phase angle is on the order of +P, -P, PI/2-P,
C         or -PI/2+P.
C
C***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
C                 matical Functions, National Bureau of Standards
C                 Applied Mathematics Series 55, U. S. Department
C                 of Commerce, Tenth Printing (1972) or later.
C               2. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument, Report SAND83-0086, Sandia National
C                 Laboratories, Albuquerque, NM, May 1983.
C               3. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument and Large Order, Report SAND83-0643,
C                 Sandia National Laboratories, Albuquerque, NM, May
C                 1983.
C               4. D. E. Amos, A Subroutine Package for Bessel Functions
C                 of a Complex Argument and Nonnegative Order, Report
C                 SAND85-1018, Sandia National Laboratory, Albuquerque,
C                 NM, May 1985.
C               5. D. E. Amos, A portable package for Bessel functions
C                 of a complex argument and nonnegative order, ACM
C                 Transactions on Mathematical Software, 12 (September
C                 1986), pp. 265-273.
C
C***ROUTINES CALLED  CACON, CBKNU, CBUNK, CUOIK, I1MACH, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C***END PROLOGUE  CBESH
C
      COMPLEX CY, Z, ZN, ZT, CSGN
      REAL AA, ALIM, ALN, ARG, AZ, CPN, DIG, ELIM, FMM, FN, FNU, FNUL,
     * HPI, RHPI, RL, R1M5, SGN, SPN, TOL, UFL, XN, XX, YN, YY, R1MACH,
     * BB, ASCLE, RTOL, ATOL
      INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, M,
     * MM, MR, N, NN, NUF, NW, NZ, I1MACH
      DIMENSION CY(N)
C
      DATA HPI /1.57079632679489662E0/
C
C***FIRST EXECUTABLE STATEMENT  CBESH
      NZ=0
      XX = REAL(Z)
      YY = AIMAG(Z)
      IERR = 0
      IF (XX.EQ.0.0E0 .AND. YY.EQ.0.0E0) IERR=1
      IF (FNU.LT.0.0E0) IERR=1
      IF (M.LT.1 .OR. M.GT.2) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      TOL = MAX(R1MACH(4),1.0E-18)
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
      FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
      RL = 1.2E0*DIG + 3.0E0
      FN = FNU + (NN-1)
      MM = 3 - M - M
      FMM = MM
      ZN = Z*CMPLX(0.0E0,-FMM)
      XN = REAL(ZN)
      YN = AIMAG(ZN)
      AZ = ABS(Z)
C-----------------------------------------------------------------------
C     TEST FOR RANGE
C-----------------------------------------------------------------------
      AA = 0.5E0/TOL
      BB=I1MACH(9)*0.5E0
      AA=MIN(AA,BB)
      IF(AZ.GT.AA) GO TO 240
      IF(FN.GT.AA) GO TO 240
      AA=SQRT(AA)
      IF(AZ.GT.AA) IERR=3
      IF(FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
      UFL = R1MACH(1)*1.0E+3
      IF (AZ.LT.UFL) GO TO 220
      IF (FNU.GT.FNUL) GO TO 90
      IF (FN.LE.1.0E0) GO TO 70
      IF (FN.GT.2.0E0) GO TO 60
      IF (AZ.GT.TOL) GO TO 70
      ARG = 0.5E0*AZ
      ALN = -FN*ALOG(ARG)
      IF (ALN.GT.ELIM) GO TO 220
      GO TO 70
   60 CONTINUE
      CALL CUOIK(ZN, FNU, KODE, 2, NN, CY, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 220
      NZ = NZ + NUF
      NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
      IF (NN.EQ.0) GO TO 130
   70 CONTINUE
      IF ((XN.LT.0.0E0) .OR. (XN.EQ.0.0E0 .AND. YN.LT.0.0E0 .AND.
     * M.EQ.2)) GO TO 80
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR.
C     YN.GE.0. .OR. M=1)
C-----------------------------------------------------------------------
      CALL CBKNU(ZN, FNU, KODE, NN, CY, NZ, TOL, ELIM, ALIM)
      GO TO 110
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C-----------------------------------------------------------------------
   80 CONTINUE
      MR = -MM
      CALL CACON(ZN, FNU, KODE, MR, NN, CY, NW, RL, FNUL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 230
      NZ=NW
      GO TO 110
   90 CONTINUE
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
      MR = 0
      IF ((XN.GE.0.0E0) .AND. (XN.NE.0.0E0 .OR. YN.GE.0.0E0 .OR.
     * M.NE.2)) GO TO 100
      MR = -MM
      IF (XN.EQ.0.0E0 .AND. YN.LT.0.0E0) ZN = -ZN
  100 CONTINUE
      CALL CBUNK(ZN, FNU, KODE, MR, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 230
      NZ = NZ + NW
  110 CONTINUE
C-----------------------------------------------------------------------
C     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
C
C     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
C-----------------------------------------------------------------------
      SGN = SIGN(HPI,-FMM)
C-----------------------------------------------------------------------
C     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = FNU
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-(INU-IR))*SGN
      RHPI = 1.0E0/SGN
      CPN = RHPI*COS(ARG)
      SPN = RHPI*SIN(ARG)
C     ZN = CMPLX(-SPN,CPN)
      CSGN = CMPLX(-SPN,CPN)
C     IF (MOD(INUH,2).EQ.1) ZN = -ZN
      IF (MOD(INUH,2).EQ.1) CSGN = -CSGN
      ZT = CMPLX(0.0E0,-FMM)
      RTOL = 1.0E0/TOL
      ASCLE = UFL*RTOL
      DO 120 I=1,NN
C       CY(I) = CY(I)*ZN
C       ZN = ZN*ZT
        ZN=CY(I)
        AA=REAL(ZN)
        BB=AIMAG(ZN)
        ATOL=1.0E0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 125
          ZN = ZN*CMPLX(RTOL,0.0E0)
          ATOL = TOL
  125   CONTINUE
        ZN = ZN*CSGN
        CY(I) = ZN*CMPLX(ATOL,0.0E0)
        CSGN = CSGN*ZT
  120 CONTINUE
      RETURN
  130 CONTINUE
      IF (XN.LT.0.0E0) GO TO 220
      RETURN
  220 CONTINUE
      IERR=2
      NZ=0
      RETURN
  230 CONTINUE
      IF(NW.EQ.(-1)) GO TO 220
      NZ=0
      IERR=5
      RETURN
  240 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
