*DECK ZBESY
      SUBROUTINE ZBESY (ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, CWRKR,
     +   CWRKI, IERR)
C***BEGIN PROLOGUE  ZBESY
C***PURPOSE  Compute a sequence of the Bessel functions Y(a,z) for
C            complex argument z and real nonnegative orders a=b,b+1,
C            b+2,... where b>0.  A scaling option is available to
C            help avoid overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10A4
C***TYPE      COMPLEX (CBESY-C, ZBESY-C)
C***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF SECOND KIND, WEBER'S FUNCTION,
C             Y BESSEL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         On KODE=1, ZBESY computes an N member sequence of complex
C         Bessel functions CY(L)=Y(FNU+L-1,Z) for real nonnegative
C         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
C         -pi<arg(Z)<=pi where Z=ZR+i*ZI.  On KODE=2, CBESY returns
C         the scaled functions
C
C            CY(L) = exp(-abs(Y))*Y(FNU+L-1,Z),  L=1,...,N, Y=Im(Z)
C
C         which remove the exponential growth in both the upper and
C         lower half planes as Z goes to infinity.  Definitions and
C         notation are found in the NBS Handbook of Mathematical
C         Functions (Ref. 1).
C
C         Input
C           ZR     - DOUBLE PRECISION real part of nonzero argument Z
C           ZI     - DOUBLE PRECISION imag part of nonzero argument Z
C           FNU    - DOUBLE PRECISION initial order, FNU>=0
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            CY(L)=Y(FNU+L-1,Z), L=1,...,N
C                        =2  returns
C                            CY(L)=Y(FNU+L-1,Z)*exp(-abs(Y)), L=1,...,N
C                            where Y=Im(Z)
C           N      - Number of terms in the sequence, N>=1
C           CWRKR  - DOUBLE PRECISION work vector of dimension N
C           CWRKI  - DOUBLE PRECISION work vector of dimension N
C
C         Output
C           CYR    - DOUBLE PRECISION real part of result vector
C           CYI    - DOUBLE PRECISION imag part of result vector
C           NZ     - Number of underflows set to zero
C                    NZ=0    Normal return
C                    NZ>0    CY(L)=0 for NZ values of L, usually on
C                            KODE=2 (the underflows may not be in an
C                            uninterrupted sequence)
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
C            Y(a,z) = (H(1,a,z) - H(2,a,z))/(2*i)
C
C         where the Hankel functions are computed as described in CBESH.
C
C         For negative orders, the formula
C
C            Y(-a,z) = Y(a,z)*cos(a*pi) + J(a,z)*sin(a*pi)
C
C         can be used.  However, for large orders close to half odd
C         integers the function changes radically.  When a is a large
C         positive half odd integer, the magnitude of Y(-a,z)=J(a,z)*
C         sin(a*pi) is a large negative power of ten.  But when a is
C         not a half odd integer, Y(a,z) dominates in magnitude with a
C         large positive power of ten and the most that the second term
C         can be reduced is by unit roundoff from the coefficient.
C         Thus,  wide changes can occur within unit roundoff of a large
C         half odd integer.  Here, large means a>abs(z).
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
C***ROUTINES CALLED  D1MACH, I1MACH, ZBESH
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C***END PROLOGUE  ZBESY
C
C     COMPLEX CWRK,CY,C1,C2,EX,HCI,Z,ZU,ZV
      DOUBLE PRECISION CWRKI, CWRKR, CYI, CYR, C1I, C1R, C2I, C2R,
     * ELIM, EXI, EXR, EY, FNU, HCII, STI, STR, TAY, ZI, ZR,
     * D1MACH, ASCLE, RTOL, ATOL, AA, BB, TOL, R1M5
      INTEGER I, IERR, K, KODE, K1, K2, N, NZ, NZ1, NZ2, I1MACH
      DIMENSION CYR(N), CYI(N), CWRKR(N), CWRKI(N)
C***FIRST EXECUTABLE STATEMENT  ZBESY
      IERR = 0
      NZ=0
      IF (ZR.EQ.0.0D0 .AND. ZI.EQ.0.0D0) IERR=1
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      HCII = 0.5D0
      CALL ZBESH(ZR, ZI, FNU, KODE, 1, N, CYR, CYI, NZ1, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 170
      CALL ZBESH(ZR, ZI, FNU, KODE, 2, N, CWRKR, CWRKI, NZ2, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 170
      NZ = MIN(NZ1,NZ2)
      IF (KODE.EQ.2) GO TO 60
      DO 50 I=1,N
        STR = CWRKR(I) - CYR(I)
        STI = CWRKI(I) - CYI(I)
        CYR(I) = -STI*HCII
        CYI(I) = STR*HCII
   50 CONTINUE
      RETURN
   60 CONTINUE
      TOL = MAX(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      K = MIN(ABS(K1),ABS(K2))
      R1M5 = D1MACH(5)
C-----------------------------------------------------------------------
C     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
C-----------------------------------------------------------------------
      ELIM = 2.303D0*(K*R1M5-3.0D0)
      EXR = COS(ZR)
      EXI = SIN(ZR)
      EY = 0.0D0
      TAY = ABS(ZI+ZI)
      IF (TAY.LT.ELIM) EY = EXP(-TAY)
      IF (ZI.LT.0.0D0) GO TO 90
      C1R = EXR*EY
      C1I = EXI*EY
      C2R = EXR
      C2I = -EXI
   70 CONTINUE
      NZ = 0
      RTOL = 1.0D0/TOL
      ASCLE = D1MACH(1)*RTOL*1.0D+3
      DO 80 I=1,N
C       STR = C1R*CYR(I) - C1I*CYI(I)
C       STI = C1R*CYI(I) + C1I*CYR(I)
C       STR = -STR + C2R*CWRKR(I) - C2I*CWRKI(I)
C       STI = -STI + C2R*CWRKI(I) + C2I*CWRKR(I)
C       CYR(I) = -STI*HCII
C       CYI(I) = STR*HCII
        AA = CWRKR(I)
        BB = CWRKI(I)
        ATOL = 1.0D0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 75
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   75   CONTINUE
        STR = (AA*C2R - BB*C2I)*ATOL
        STI = (AA*C2I + BB*C2R)*ATOL
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 85
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   85   CONTINUE
        STR = STR - (AA*C1R - BB*C1I)*ATOL
        STI = STI - (AA*C1I + BB*C1R)*ATOL
        CYR(I) = -STI*HCII
        CYI(I) =  STR*HCII
        IF (STR.EQ.0.0D0 .AND. STI.EQ.0.0D0 .AND. EY.EQ.0.0D0) NZ = NZ
     *   + 1
   80 CONTINUE
      RETURN
   90 CONTINUE
      C1R = EXR
      C1I = EXI
      C2R = EXR*EY
      C2I = -EXI*EY
      GO TO 70
  170 CONTINUE
      NZ = 0
      RETURN
      END
