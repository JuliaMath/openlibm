*DECK DXNRMP
      SUBROUTINE DXNRMP (NU, MU1, MU2, DARG, MODE, DPN, IPN, ISIG,
     1   IERROR)
C***BEGIN PROLOGUE  DXNRMP
C***PURPOSE  Compute normalized Legendre polynomials.
C***LIBRARY   SLATEC
C***CATEGORY  C3A2, C9
C***TYPE      DOUBLE PRECISION (XNRMP-S, DXNRMP-D)
C***KEYWORDS  LEGENDRE FUNCTIONS
C***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
C           Smith, John M., (NBS and George Mason University)
C***DESCRIPTION
C
C        SUBROUTINE TO CALCULATE NORMALIZED LEGENDRE POLYNOMIALS
C        (XNRMP is single-precision version)
C        DXNRMP calculates normalized Legendre polynomials of varying
C        order and fixed argument and degree. The order MU and degree
C        NU are non-negative integers and the argument is real. Because
C        the algorithm requires the use of numbers outside the normal
C        machine range, this subroutine employs a special arithmetic
C        called extended-range arithmetic. See J.M. Smith, F.W.J. Olver,
C        and D.W. Lozier, Extended-Range Arithmetic and Normalized
C        Legendre Polynomials, ACM Transactions on Mathematical Soft-
C        ware, 93-105, March 1981, for a complete description of the
C        algorithm and special arithmetic. Also see program comments
C        in DXSET.
C
C        The normalized Legendre polynomials are multiples of the
C        associated Legendre polynomials of the first kind where the
C        normalizing coefficients are chosen so as to make the integral
C        from -1 to 1 of the square of each function equal to 1. See
C        E. Jahnke, F. Emde and F. Losch, Tables of Higher Functions,
C        McGraw-Hill, New York, 1960, p. 121.
C
C        The input values to DXNRMP are NU, MU1, MU2, DARG, and MODE.
C        These must satisfy
C          1. NU .GE. 0 specifies the degree of the normalized Legendre
C             polynomial that is wanted.
C          2. MU1 .GE. 0 specifies the lowest-order normalized Legendre
C             polynomial that is wanted.
C          3. MU2 .GE. MU1 specifies the highest-order normalized Leg-
C             endre polynomial that is wanted.
C         4a. MODE = 1 and -1.0D0 .LE. DARG .LE. 1.0D0 specifies that
C             Normalized Legendre(NU, MU, DARG) is wanted for MU = MU1,
C             MU1 + 1, ..., MU2.
C         4b. MODE = 2 and -3.14159... .LT. DARG .LT. 3.14159... spec-
C             ifies that Normalized Legendre(NU, MU, COS(DARG)) is
C             wanted for MU = MU1, MU1 + 1, ..., MU2.
C
C        The output of DXNRMP consists of the two vectors DPN and IPN
C        and the error estimate ISIG. The computed values are stored as
C        extended-range numbers such that
C             (DPN(1),IPN(1))=NORMALIZED LEGENDRE(NU,MU1,DX)
C             (DPN(2),IPN(2))=NORMALIZED LEGENDRE(NU,MU1+1,DX)
C                .
C                .
C             (DPN(K),IPN(K))=NORMALIZED LEGENDRE(NU,MU2,DX)
C        where K = MU2 - MU1 + 1 and DX = DARG or COS(DARG) according
C        to whether MODE = 1 or 2. Finally, ISIG is an estimate of the
C        number of decimal digits lost through rounding errors in the
C        computation. For example if DARG is accurate to 12 significant
C        decimals, then the computed function values are accurate to
C        12 - ISIG significant decimals (except in neighborhoods of
C        zeros).
C
C        The interpretation of (DPN(I),IPN(I)) is DPN(I)*(IR**IPN(I))
C        where IR is the internal radix of the computer arithmetic. When
C        IPN(I) = 0 the value of the normalized Legendre polynomial is
C        contained entirely in DPN(I) and subsequent double-precision
C        computations can be performed without further consideration of
C        extended-range arithmetic. However, if IPN(I) .NE. 0 the corre-
C        sponding value of the normalized Legendre polynomial cannot be
C        represented in double-precision because of overflow or under-
C        flow. THE USER MUST TEST IPN(I) IN HIS/HER PROGRAM. In the case
C        that IPN(I) is nonzero, the user could rewrite his/her program
C        to use extended range arithmetic.
C
C
C
C        The interpretation of (DPN(I),IPN(I)) can be changed to
C        DPN(I)*(10**IPN(I)) by calling the extended-range subroutine
C        DXCON. This should be done before printing the computed values.
C        As an example of usage, the Fortran coding
C              J = K
C              DO 20 I = 1, K
C              CALL DXCON(DPN(I), IPN(I),IERROR)
C              IF (IERROR.NE.0) RETURN
C              PRINT 10, DPN(I), IPN(I)
C           10 FORMAT(1X, D30.18 , I15)
C              IF ((IPN(I) .EQ. 0) .OR. (J .LT. K)) GO TO 20
C              J = I - 1
C           20 CONTINUE
C        will print all computed values and determine the largest J
C        such that IPN(1) = IPN(2) = ... = IPN(J) = 0. Because of the
C        change of representation caused by calling DXCON, (DPN(I),
C        IPN(I)) for I = J+1, J+2, ... cannot be used in subsequent
C        extended-range computations.
C
C        IERROR is an error indicator. If no errors are detected,
C        IERROR=0 when control returns to the calling routine. If
C        an error is detected, IERROR is returned as nonzero. The
C        calling routine must check the value of IERROR.
C
C        If IERROR=212 or 213, invalid input was provided to DXNRMP.
C        If IERROR=201,202,203, or 204, invalid input was provided
C        to DXSET.
C        If IERROR=205 or 206, an internal consistency error occurred
C        in DXSET (probably due to a software malfunction in the
C        library routine I1MACH).
C        If IERROR=207, an overflow or underflow of an extended-range
C        number was detected in DXADJ.
C        If IERROR=208, an overflow or underflow of an extended-range
C        number was detected in DXC210.
C
C***SEE ALSO  DXSET
C***REFERENCES  Smith, Olver and Lozier, Extended-Range Arithmetic and
C                 Normalized Legendre Polynomials, ACM Trans on Math
C                 Softw, v 7, n 1, March 1981, pp 93--105.
C***ROUTINES CALLED  DXADD, DXADJ, DXRED, DXSET, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820712  DATE WRITTEN
C   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  DXNRMP
      INTEGER NU, MU1, MU2, MODE, IPN, ISIG
      DOUBLE PRECISION DARG, DPN
      DIMENSION DPN(*), IPN(*)
      DOUBLE PRECISION C1,C2,P,P1,P2,P3,S,SX,T,TX,X,DK
C CALL DXSET TO INITIALIZE EXTENDED-RANGE ARITHMETIC (SEE DXSET
C LISTING FOR DETAILS)
C***FIRST EXECUTABLE STATEMENT  DXNRMP
      IERROR=0
      CALL DXSET (0, 0, 0.0D0, 0,IERROR)
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
   10 IF (ABS(DARG).GT.1.0D0) GO TO 120
      IF (ABS(DARG).EQ.1.0D0) GO TO 90
      X = DARG
      SX = SQRT((1.0D0+ABS(X))*((0.5D0-ABS(X))+0.5D0))
      TX = X/SX
      ISIG = LOG10(2.0D0*NU*(5.0D0+TX**2))
      GO TO 30
   20 IF (ABS(DARG).GT.4.0D0*ATAN(1.0D0)) GO TO 120
      IF (DARG.EQ.0.0D0) GO TO 90
      X = COS(DARG)
      SX = ABS(SIN(DARG))
      TX = X/SX
      ISIG = LOG10(2.0D0*NU*(5.0D0+ABS(DARG*TX)))
C
C        BEGIN CALCULATION
C
   30 MU = MU2
      I = MU2 - MU1 + 1
C
C        IF MU.GT.NU, NORMALIZED LEGENDRE(NU,MU,X)=0.
C
   40 IF (MU.LE.NU) GO TO 50
      DPN(I) = 0.0D0
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
      P1 = 0.0D0
      IP1 = 0
C
C        CALCULATE P2 = NORMALIZED LEGENDRE(NU,NU,X)
C
      P2 = 1.0D0
      IP2 = 0
      P3 = 0.5D0
      DK = 2.0D0
      DO 60 J=1,NU
        P3 = ((DK+1.0D0)/DK)*P3
        P2 = P2*SX
        CALL DXADJ(P2, IP2,IERROR)
        IF (IERROR.NE.0) RETURN
        DK = DK + 2.0D0
   60 CONTINUE
      P2 = P2*SQRT(P3)
      CALL DXADJ(P2, IP2,IERROR)
      IF (IERROR.NE.0) RETURN
      S = 2.0D0*TX
      T = 1.0D0/NU
      IF (MU2.LT.NU) GO TO 70
      DPN(I) = P2
      IPN(I) = IP2
      I = I - 1
      IF (I .EQ. 0) GO TO 140
C
C        RECURRENCE PROCESS
C
   70 P = MU*T
      C1 = 1.0D0/SQRT((1.0D0-P+T)*(1.0D0+P))
      C2 = S*P*C1*P2
      C1 = -SQRT((1.0D0+P+T)*(1.0D0-P))*C1*P1
      CALL DXADD(C2, IP2, C1, IP1, P, IP,IERROR)
      IF (IERROR.NE.0) RETURN
      MU = MU - 1
      IF (MU.GT.MU2) GO TO 80
C
C        STORE IN ARRAY DPN FOR RETURN TO CALLING ROUTINE.
C
      DPN(I) = P
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
        DPN(I) = 0.0D0
        IPN(I) = 0
  100 CONTINUE
      ISIG = 0
      IF (MU1.GT.0) GO TO 160
      ISIG = 1
      DPN(1) = SQRT(NU+0.5D0)
      IPN(1) = 0
      IF (MOD(NU,2).EQ.0) GO TO 160
      IF (MODE.EQ.1 .AND. DARG.EQ.1.0D0) GO TO 160
      IF (MODE.EQ.2) GO TO 160
      DPN(1) = -DPN(1)
      GO TO 160
C
C          ERROR PRINTOUTS AND TERMINATION.
C
  110 CALL XERMSG ('SLATEC', 'DXNRMP', 'NU, MU1, MU2 or MODE not valid',
     +             212, 1)
      IERROR=212
      RETURN
  120 CALL XERMSG ('SLATEC', 'DXNRMP', 'DARG out of range', 213, 1)
      IERROR=213
      RETURN
C
C        RETURN TO CALLING PROGRAM
C
  140 K = MU2 - MU1 + 1
      DO 150 I=1,K
        CALL DXRED(DPN(I),IPN(I),IERROR)
        IF (IERROR.NE.0) RETURN
  150 CONTINUE
  160 RETURN
      END
