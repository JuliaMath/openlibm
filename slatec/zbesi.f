*DECK ZBESI
      SUBROUTINE ZBESI (ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
C***BEGIN PROLOGUE  ZBESI
C***PURPOSE  Compute a sequence of the Bessel functions I(a,z) for
C            complex argument z and real nonnegative orders a=b,b+1,
C            b+2,... where b>0.  A scaling option is available to
C            help avoid overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10B4
C***TYPE      COMPLEX (CBESI-C, ZBESI-C)
C***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT, I BESSEL FUNCTIONS,
C             MODIFIED BESSEL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C                    ***A DOUBLE PRECISION ROUTINE***
C         On KODE=1, ZBESI computes an N-member sequence of complex
C         Bessel functions CY(L)=I(FNU+L-1,Z) for real nonnegative
C         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
C         -pi<arg(Z)<=pi where Z=ZR+i*ZI.  On KODE=2, CBESI returns
C         the scaled functions
C
C            CY(L) = exp(-abs(X))*I(FNU+L-1,Z), L=1,...,N and X=Re(Z)
C
C         which removes the exponential growth in both the left and
C         right half-planes as Z goes to infinity.
C
C         Input
C           ZR     - DOUBLE PRECISION real part of argument Z
C           ZI     - DOUBLE PRECISION imag part of argument Z
C           FNU    - DOUBLE PRECISION initial order, FNU>=0
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            CY(L)=I(FNU+L-1,Z), L=1,...,N
C                        =2  returns
C                            CY(L)=exp(-abs(X))*I(FNU+L-1,Z), L=1,...,N
C                            where X=Re(Z)
C           N      - Number of terms in the sequence, N>=1
C
C         Output
C           CYR    - DOUBLE PRECISION real part of result vector
C           CYI    - DOUBLE PRECISION imag part of result vector
C           NZ     - Number of underflows set to zero
C                    NZ=0    Normal return
C                    NZ>0    CY(L)=0, L=N-NZ+1,...,N
C           IERR   - Error flag
C                    IERR=0  Normal return     - COMPUTATION COMPLETED
C                    IERR=1  Input error       - NO COMPUTATION
C                    IERR=2  Overflow          - NO COMPUTATION
C                            (Re(Z) too large on KODE=1)
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
C         The computation of I(a,z) is carried out by the power series
C         for small abs(z), the asymptotic expansion for large abs(z),
C         the Miller algorithm normalized by the Wronskian and a
C         Neumann series for intermediate magnitudes of z, and the
C         uniform asymptotic expansions for I(a,z) and J(a,z) for
C         large orders a.  Backward recurrence is used to generate
C         sequences or reduce orders when necessary.
C
C         The calculations above are done in the right half plane and
C         continued into the left half plane by the formula
C
C            I(a,z*exp(t)) = exp(t*a)*I(a,z), Re(z)>0
C                        t = i*pi or -i*pi
C
C         For negative orders, the formula
C
C            I(-a,z) = I(a,z) + (2/pi)*sin(pi*a)*K(a,z)
C
C         can be used.  However, for large orders close to integers the
C         the function changes radically.  When a is a large positive
C         integer, the magnitude of I(-a,z)=I(a,z) is a large
C         negative power of ten. But when a is not an integer,
C         K(a,z) dominates in magnitude with a large positive power of
C         ten and the most that the second term can be reduced is by
C         unit roundoff from the coefficient. Thus, wide changes can
C         occur within unit roundoff of a large integer for a. Here,
C         large means a>abs(z).
C
C         In most complex variable computation, one must evaluate ele-
C         mentary functions.  When the magnitude of Z or FNU+N-1 is
C         large, losses of significance by argument reduction occur.
C         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
C         losses exceeding half precision are likely and an error flag
C         IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is double
C         precision unit roundoff limited to 18 digits precision.  Also,
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
C***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZBINU
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C***END PROLOGUE  ZBESI
C     COMPLEX CONE,CSGN,CW,CY,CZERO,Z,ZN
      DOUBLE PRECISION AA, ALIM, ARG, CONEI, CONER, CSGNI, CSGNR, CYI,
     * CYR, DIG, ELIM, FNU, FNUL, PI, RL, R1M5, STR, TOL, ZI, ZNI, ZNR,
     * ZR, D1MACH, AZ, BB, FN, ZABS, ASCLE, RTOL, ATOL, STI
      INTEGER I, IERR, INU, K, KODE, K1,K2,N,NZ,NN, I1MACH
      DIMENSION CYR(N), CYI(N)
      EXTERNAL ZABS
      DATA PI /3.14159265358979324D0/
      DATA CONER, CONEI /1.0D0,0.0D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESI
      IERR = 0
      NZ=0
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
      TOL = MAX(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303D0*(K*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*K1
      DIG = MIN(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + MAX(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      FN = FNU+(N-1)
      AA = 0.5D0/TOL
      BB=I1MACH(9)*0.5D0
      AA = MIN(AA,BB)
      IF (AZ.GT.AA) GO TO 260
      IF (FN.GT.AA) GO TO 260
      AA = SQRT(AA)
      IF (AZ.GT.AA) IERR=3
      IF (FN.GT.AA) IERR=3
      ZNR = ZR
      ZNI = ZI
      CSGNR = CONER
      CSGNI = CONEI
      IF (ZR.GE.0.0D0) GO TO 40
      ZNR = -ZR
      ZNI = -ZI
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = FNU
      ARG = (FNU-INU)*PI
      IF (ZI.LT.0.0D0) ARG = -ARG
      CSGNR = COS(ARG)
      CSGNI = SIN(ARG)
      IF (MOD(INU,2).EQ.0) GO TO 40
      CSGNR = -CSGNR
      CSGNI = -CSGNI
   40 CONTINUE
      CALL ZBINU(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL,
     * ELIM, ALIM)
      IF (NZ.LT.0) GO TO 120
      IF (ZR.GE.0.0D0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
C-----------------------------------------------------------------------
      NN = N - NZ
      IF (NN.EQ.0) RETURN
      RTOL = 1.0D0/TOL
      ASCLE = D1MACH(1)*RTOL*1.0D+3
      DO 50 I=1,NN
C       STR = CYR(I)*CSGNR - CYI(I)*CSGNI
C       CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR
C       CYR(I) = STR
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 55
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   55   CONTINUE
        STR = AA*CSGNR - BB*CSGNI
        STI = AA*CSGNI + BB*CSGNR
        CYR(I) = STR*ATOL
        CYI(I) = STI*ATOL
        CSGNR = -CSGNR
        CSGNI = -CSGNI
   50 CONTINUE
      RETURN
  120 CONTINUE
      IF(NZ.EQ.(-2)) GO TO 130
      NZ = 0
      IERR=2
      RETURN
  130 CONTINUE
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
