*DECK DGAUS8
      SUBROUTINE DGAUS8 (FUN, A, B, ERR, ANS, IERR)
C***BEGIN PROLOGUE  DGAUS8
C***PURPOSE  Integrate a real function of one variable over a finite
C            interval using an adaptive 8-point Legendre-Gauss
C            algorithm.  Intended primarily for high accuracy
C            integration or integration of smooth functions.
C***LIBRARY   SLATEC
C***CATEGORY  H2A1A1
C***TYPE      DOUBLE PRECISION (GAUS8-S, DGAUS8-D)
C***KEYWORDS  ADAPTIVE QUADRATURE, AUTOMATIC INTEGRATOR,
C             GAUSS QUADRATURE, NUMERICAL INTEGRATION
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract  *** a DOUBLE PRECISION routine ***
C        DGAUS8 integrates real functions of one variable over finite
C        intervals using an adaptive 8-point Legendre-Gauss algorithm.
C        DGAUS8 is intended primarily for high accuracy integration
C        or integration of smooth functions.
C
C        The maximum number of significant digits obtainable in ANS
C        is the smaller of 18 and the number of digits carried in
C        double precision arithmetic.
C
C     Description of Arguments
C
C        Input--* FUN, A, B, ERR are DOUBLE PRECISION *
C        FUN - name of external function to be integrated.  This name
C              must be in an EXTERNAL statement in the calling program.
C              FUN must be a DOUBLE PRECISION function of one DOUBLE
C              PRECISION argument.  The value of the argument to FUN
C              is the variable of integration which ranges from A to B.
C        A   - lower limit of integration
C        B   - upper limit of integration (may be less than A)
C        ERR - is a requested pseudorelative error tolerance.  Normally
C              pick a value of ABS(ERR) so that DTOL .LT. ABS(ERR) .LE.
C              1.0D-3 where DTOL is the larger of 1.0D-18 and the
C              double precision unit roundoff D1MACH(4).  ANS will
C              normally have no more error than ABS(ERR) times the
C              integral of the absolute value of FUN(X).  Usually,
C              smaller values of ERR yield more accuracy and require
C              more function evaluations.
C
C              A negative value for ERR causes an estimate of the
C              absolute error in ANS to be returned in ERR.  Note that
C              ERR must be a variable (not a constant) in this case.
C              Note also that the user must reset the value of ERR
C              before making any more calls that use the variable ERR.
C
C        Output--* ERR,ANS are double precision *
C        ERR - will be an estimate of the absolute error in ANS if the
C              input value of ERR was negative.  (ERR is unchanged if
C              the input value of ERR was non-negative.)  The estimated
C              error is solely for information to the user and should
C              not be used as a correction to the computed integral.
C        ANS - computed value of integral
C        IERR- a status code
C            --Normal codes
C               1 ANS most likely meets requested error tolerance,
C                 or A=B.
C              -1 A and B are too nearly equal to allow normal
C                 integration.  ANS is set to zero.
C            --Abnormal code
C               2 ANS probably does not meet requested error tolerance.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, I1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   810223  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  DGAUS8
      INTEGER IERR, K, KML, KMX, L, LMN, LMX, LR, MXL, NBITS,
     1 NIB, NLMN, NLMX
      INTEGER I1MACH
      DOUBLE PRECISION A,AA,AE,ANIB,ANS,AREA,B,C,CE,EE,EF,
     1 EPS, ERR, EST, GL, GLR, GR, HH, SQ2, TOL, VL, VR, W1, W2, W3,
     2 W4, X1, X2, X3, X4, X, H
      DOUBLE PRECISION D1MACH, G8, FUN
      DIMENSION AA(60), HH(60), LR(60), VL(60), GR(60)
      SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2,
     1 NLMN, KMX, KML
      DATA X1, X2, X3, X4/
     1     1.83434642495649805D-01,     5.25532409916328986D-01,
     2     7.96666477413626740D-01,     9.60289856497536232D-01/
      DATA W1, W2, W3, W4/
     1     3.62683783378361983D-01,     3.13706645877887287D-01,
     2     2.22381034453374471D-01,     1.01228536290376259D-01/
      DATA SQ2/1.41421356D0/
      DATA NLMN/1/,KMX/5000/,KML/6/
      G8(X,H)=H*((W1*(FUN(X-X1*H) + FUN(X+X1*H))
     1           +W2*(FUN(X-X2*H) + FUN(X+X2*H)))
     2          +(W3*(FUN(X-X3*H) + FUN(X+X3*H))
     3           +W4*(FUN(X-X4*H) + FUN(X+X4*H))))
C***FIRST EXECUTABLE STATEMENT  DGAUS8
C
C     Initialize
C
      K = I1MACH(14)
      ANIB = D1MACH(5)*K/0.30102000D0
      NBITS = ANIB
      NLMX = MIN(60,(NBITS*5)/8)
      ANS = 0.0D0
      IERR = 1
      CE = 0.0D0
      IF (A .EQ. B) GO TO 140
      LMX = NLMX
      LMN = NLMN
      IF (B .EQ. 0.0D0) GO TO 10
      IF (SIGN(1.0D0,B)*A .LE. 0.0D0) GO TO 10
      C = ABS(1.0D0-A/B)
      IF (C .GT. 0.1D0) GO TO 10
      IF (C .LE. 0.0D0) GO TO 140
      ANIB = 0.5D0 - LOG(C)/0.69314718D0
      NIB = ANIB
      LMX = MIN(NLMX,NBITS-NIB-7)
      IF (LMX .LT. 1) GO TO 130
      LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0D0**(5-NBITS))/2.0D0
      IF (ERR .EQ. 0.0D0) TOL = SQRT(D1MACH(4))
      EPS = TOL
      HH(1) = (B-A)/4.0D0
      AA(1) = A
      LR(1) = 1
      L = 1
      EST = G8(AA(L)+2.0D0*HH(L),2.0D0*HH(L))
      K = 8
      AREA = ABS(EST)
      EF = 0.5D0
      MXL = 0
C
C     Compute refined estimates, estimate the error, etc.
C
   20 GL = G8(AA(L)+HH(L),HH(L))
      GR(L) = G8(AA(L)+3.0D0*HH(L),HH(L))
      K = K + 16
      AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
C     IF (L .LT .LMN) GO TO 11
      GLR = GL + GR(L)
      EE = ABS(EST-GLR)*EF
      AE = MAX(EPS*AREA,TOL*ABS(GLR))
      IF (EE-AE) 40, 40, 50
   30 MXL = 1
   40 CE = CE + (EST-GLR)
      IF (LR(L)) 60, 60, 80
C
C     Consider the left half of this level
C
   50 IF (K .GT. KMX) LMX = KML
      IF (L .GE. LMX) GO TO 30
      L = L + 1
      EPS = EPS*0.5D0
      EF = EF/SQ2
      HH(L) = HH(L-1)*0.5D0
      LR(L) = -1
      AA(L) = AA(L-1)
      EST = GL
      GO TO 20
C
C     Proceed to right half at this level
C
   60 VL(L) = GLR
   70 EST = GR(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 4.0D0*HH(L)
      GO TO 20
C
C     Return one level
C
   80 VR = GLR
   90 IF (L .LE. 1) GO TO 120
      L = L - 1
      EPS = EPS*2.0D0
      EF = EF*SQ2
      IF (LR(L)) 100, 100, 110
  100 VL(L) = VL(L+1) + VR
      GO TO 70
  110 VR = VL(L+1) + VR
      GO TO 90
C
C     Exit
C
  120 ANS = VR
      IF ((MXL.EQ.0) .OR. (ABS(CE).LE.2.0D0*TOL*AREA)) GO TO 140
      IERR = 2
      CALL XERMSG ('SLATEC', 'DGAUS8',
     +   'ANS is probably insufficiently accurate.', 3, 1)
      GO TO 140
  130 IERR = -1
      CALL XERMSG ('SLATEC', 'DGAUS8',
     +   'A and B are too nearly equal to allow normal integration. $$'
     +   // 'ANS is set to zero and IERR to -1.', 1, -1)
  140 IF (ERR .LT. 0.0D0) ERR = CE
      RETURN
      END
