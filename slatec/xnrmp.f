*DECK XNRMP
      SUBROUTINE XNRMP (NU, MU1, MU2, SARG, MODE, SPN, IPN, ISIG,
     1   IERROR)
C***BEGIN PROLOGUE  XNRMP
C***PURPOSE  Compute normalized Legendre polynomials.
C***LIBRARY   SLATEC
C***CATEGORY  C3A2, C9
C***TYPE      SINGLE PRECISION (XNRMP-S, DXNRMP-D)
C***KEYWORDS  LEGENDRE FUNCTIONS
C***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
C           Smith, John M., (NBS and George Mason University)
C***DESCRIPTION
C
C        SUBROUTINE TO CALCULATE NORMALIZED LEGENDRE POLYNOMIALS
C        (DXNRMP is double-precision version)
C        XNRMP calculates normalized Legendre polynomials of varying
C        order and fixed argument and degree. The order MU and degree
C        NU are non-negative integers and the argument is real. Because
C        the algorithm requires the use of numbers outside the normal
C        machine range, this subroutine employs a special arithmetic
C        called extended-range arithmetic. See J.M. Smith, F.W.J. Olver,
C        and D.W. Lozier, Extended-Range Arithmetic and Normalized
C        Legendre Polynomials, ACM Transactions on Mathematical Soft-
C        ware, 93-105, March 1981, for a complete description of the
C        algorithm and special arithmetic. Also see program comments
C        in XSET.
C
C        The normalized Legendre polynomials are multiples of the
C        associated Legendre polynomials of the first kind where the
C        normalizing coefficients are chosen so as to make the integral
C        from -1 to 1 of the square of each function equal to 1. See
C        E. Jahnke, F. Emde and F. Losch, Tables of Higher Functions,
C        McGraw-Hill, New York, 1960, p. 121.
C
C        The input values to XNRMP are NU, MU1, MU2, SARG, and MODE.
C        These must satisfy
C          1. NU .GE. 0 specifies the degree of the normalized Legendre
C             polynomial that is wanted.
C          2. MU1 .GE. 0 specifies the lowest-order normalized Legendre
C             polynomial that is wanted.
C          3. MU2 .GE. MU1 specifies the highest-order normalized Leg-
C             endre polynomial that is wanted.
C         4a. MODE = 1 and -1.0 .LE. SARG .LE. 1.0 specifies that
C             Normalized Legendre(NU, MU, SARG) is wanted for MU = MU1,
C             MU1 + 1, ..., MU2.
C         4b. MODE = 2 and -3.14159... .LT. SARG .LT. 3.14159... spec-
C             ifies that Normalized Legendre(NU, MU, COS(SARG)) is want-
C             ed for MU = MU1, MU1 + 1, ..., MU2.
C
C        The output of XNRMP consists of the two vectors SPN and IPN
C        and the error estimate ISIG. The computed values are stored as
C        extended-range numbers such that
C             (SPN(1),IPN(1))=NORMALIZED LEGENDRE(NU,MU1,X)
C             (SPN(2),IPN(2))=NORMALIZED LEGENDRE(NU,MU1+1,X)
C                .
C                .
C             (SPN(K),IPN(K))=NORMALIZED LEGENDRE(NU,MU2,X)
C        where K = MU2 - MU1 + 1 and X = SARG or COS(SARG) according
C        to whether MODE = 1 or 2. Finally, ISIG is an estimate of the
C        number of decimal digits lost through rounding errors in the
C        computation. For example if SARG is accurate to 12 significant
C        decimals, then the computed function values are accurate to
C        12 - ISIG significant decimals (except in neighborhoods of
C        zeros).
C
C        The interpretation of (SPN(I),IPN(I)) is SPN(I)*(IR**IPN(I))
C        where IR is the internal radix of the computer arithmetic. When
C        IPN(I) = 0 the value of the normalized Legendre polynomial is
C        contained entirely in SPN(I) and subsequent single-precision
C        computations can be performed without further consideration of
C        extended-range arithmetic. However, if IPN(I) .NE. 0 the corre-
C        sponding value of the normalized Legendre polynomial cannot be
C        represented in single-precision because of overflow or under-
C        flow. THE USER MUST TEST IPN(I) IN HIS/HER PROGRAM. In the case
C        that IPN(I) is nonzero, the user should try using double pre-
C        cision if it has a wider exponent range. If double precision
C        fails, the user could rewrite his/her program to use extended-
C        range arithmetic.
C
C        The interpretation of (SPN(I),IPN(I)) can be changed to
C        SPN(I)*(10**IPN(I)) by calling the extended-range subroutine
C        XCON. This should be done before printing the computed values.
C        As an example of usage, the Fortran coding
C              J = K
C              DO 20 I = 1, K
C              CALL XCON(SPN(I), IPN(I),IERROR)
C              IF (IERROR.NE.0) RETURN
C              PRINT 10, SPN(I), IPN(I)
C           10 FORMAT(1X, E30.8 , I15)
C              IF ((IPN(I) .EQ. 0) .OR. (J .LT. K)) GO TO 20
C              J = I - 1
C           20 CONTINUE
C        will print all computed values and determine the largest J
C        such that IPN(1) = IPN(2) = ... = IPN(J) = 0. Because of the
C        change of representation caused by calling XCON, (SPN(I),
C        IPN(I)) for I = J+1, J+2, ... cannot be used in subsequent
C        extended-range computations.
C
C        IERROR is an error indicator. If no errors are detected,
C        IERROR=0 when control returns to the calling routine. If
C        an error is detected, IERROR is returned as nonzero. The
C        calling routine must check the value of IERROR.
C
C        If IERROR=112 or 113, invalid input was provided to XNRMP.
C        If IERROR=101,102,103, or 104, invalid input was provided
C        to XSET.
C        If IERROR=105 or 106, an internal consistency error occurred
C        in XSET (probably due to a software malfunction in the
C        library routine I1MACH).
C        If IERROR=107, an overflow or underflow of an extended-range
C        number was detected in XADJ.
C        If IERROR=108, an overflow or underflow of an extended-range
C        number was detected in XC210.
C
C***SEE ALSO  XSET
C***REFERENCES  Smith, Olver and Lozier, Extended-Range Arithmetic and
C                 Normalized Legendre Polynomials, ACM Trans on Math
C                 Softw, v 7, n 1, March 1981, pp 93--105.
C***ROUTINES CALLED  XADD, XADJ, XERMSG, XRED, XSET
C***REVISION HISTORY  (YYMMDD)
C   820712  DATE WRITTEN
C   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  XNRMP
      INTEGER NU, MU1, MU2, MODE, IPN, ISIG
      REAL SARG, SPN
      DIMENSION SPN(*), IPN(*)
      REAL C1,C2,P,P1,P2,P3,S,SX,T,TX,X,RK
C CALL XSET TO INITIALIZE EXTENDED-RANGE ARITHMETIC (SEE XSET
C LISTING FOR DETAILS)
C***FIRST EXECUTABLE STATEMENT  XNRMP
      IERROR=0
      CALL XSET (0, 0, 0.0, 0,IERROR)
      IF (IERROR.NE.0) RETURN
C
C        TEST FOR PROPER INPUT VALUES.
C
      IF (NU.LT.0) GO TO 110
      IF (MU1.LT.0) GO TO 110
      IF (MU1.GT.MU2) GO TO 110
      IF (NU.EQ.0) GO TO 90
      IF (MODE.LT.1 .OR. MODE.GT.2) GO TO 110
      GO TO (10, 20), MODE
   10 IF (ABS(SARG).GT.1.0) GO TO 120
      IF (ABS(SARG).EQ.1.0) GO TO 90
      X = SARG
      SX = SQRT((1.0+ABS(X))*((0.5-ABS(X))+0.5))
      TX = X/SX
      ISIG = LOG10(2.0*NU*(5.0+TX**2))
      GO TO 30
   20 IF (ABS(SARG).GT.4.0*ATAN(1.0)) GO TO 120
      IF (SARG.EQ.0.0) GO TO 90
      X = COS(SARG)
      SX = ABS(SIN(SARG))
      TX = X/SX
      ISIG = LOG10(2.0*NU*(5.0+ABS(SARG*TX)))
C
C        BEGIN CALCULATION
C
   30 MU = MU2
      I = MU2 - MU1 + 1
C
C        IF MU.GT.NU, NORMALIZED LEGENDRE(NU,MU,X)=0.
C
   40 IF (MU.LE.NU) GO TO 50
      SPN(I) = 0.0
      IPN(I) = 0
      I = I - 1
      MU = MU - 1
      IF (I .GT. 0) GO TO 40
      ISIG = 0
      GO TO 160
   50 MU = NU
C
C        P1 = 0. = NORMALIZED LEGENDRE(NU,NU+1,X)
C
      P1 = 0.0
      IP1 = 0
C
C        CALCULATE P2 = NORMALIZED LEGENDRE(NU,NU,X)
C
      P2 = 1.0
      IP2 = 0
      P3 = 0.5
      RK = 2.0
      DO 60 J=1,NU
        P3 = ((RK+1.0)/RK)*P3
        P2 = P2*SX
        CALL XADJ(P2, IP2,IERROR)
        IF (IERROR.NE.0) RETURN
        RK = RK + 2.0
   60 CONTINUE
      P2 = P2*SQRT(P3)
      CALL XADJ(P2, IP2,IERROR)
      IF (IERROR.NE.0) RETURN
      S = 2.0*TX
      T = 1.0/NU
      IF (MU2.LT.NU) GO TO 70
      SPN(I) = P2
      IPN(I) = IP2
      I = I - 1
      IF (I .EQ. 0) GO TO 140
C
C        RECURRENCE PROCESS
C
   70 P = MU*T
      C1 = 1.0/SQRT((1.0-P+T)*(1.0+P))
      C2 = S*P*C1*P2
      C1 = -SQRT((1.0+P+T)*(1.0-P))*C1*P1
      CALL XADD(C2, IP2, C1, IP1, P, IP,IERROR)
      IF (IERROR.NE.0) RETURN
      MU = MU - 1
      IF (MU.GT.MU2) GO TO 80
C
C        STORE IN ARRAY SPN FOR RETURN TO CALLING ROUTINE.
C
      SPN(I) = P
      IPN(I) = IP
      I = I - 1
      IF (I .EQ. 0) GO TO 140
   80 P1 = P2
      IP1 = IP2
      P2 = P
      IP2 = IP
      IF (MU.LE.MU1) GO TO 140
      GO TO 70
C
C        SPECIAL CASE WHEN X=-1 OR +1, OR NU=0.
C
   90 K = MU2 - MU1 + 1
      DO 100 I=1,K
        SPN(I) = 0.0
        IPN(I) = 0
  100 CONTINUE
      ISIG = 0
      IF (MU1.GT.0) GO TO 160
      ISIG = 1
      SPN(1) = SQRT(NU+0.5)
      IPN(1) = 0
      IF (MOD(NU,2).EQ.0) GO TO 160
      IF (MODE.EQ.1 .AND. SARG.EQ.1.0) GO TO 160
      IF (MODE.EQ.2) GO TO 160
      SPN(1) = -SPN(1)
      GO TO 160
C
C          ERROR PRINTOUTS AND TERMINATION.
C
  110 CALL XERMSG ('SLATEC', 'XNRMP', 'NU, MU1, MU2 or MODE not valid',
     +             112, 1)
      IERROR=112
      RETURN
  120 CALL XERMSG ('SLATEC', 'XNRMP', 'SARG out of range', 113, 1)
      IERROR=113
      RETURN
C
C        RETURN TO CALLING PROGRAM
C
  140 K = MU2 - MU1 + 1
      DO 150 I=1,K
        CALL XRED(SPN(I),IPN(I),IERROR)
        IF (IERROR.NE.0) RETURN
  150 CONTINUE
  160 RETURN
      END
