*DECK BINT4
      SUBROUTINE BINT4 (X, Y, NDATA, IBCL, IBCR, FBCL, FBCR, KNTOPT, T,
     +   BCOEF, N, K, W)
C***BEGIN PROLOGUE  BINT4
C***PURPOSE  Compute the B-representation of a cubic spline
C            which interpolates given data.
C***LIBRARY   SLATEC
C***CATEGORY  E1A
C***TYPE      SINGLE PRECISION (BINT4-S, DBINT4-D)
C***KEYWORDS  B-SPLINE, CUBIC SPLINES, DATA FITTING, INTERPOLATION
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C         BINT4 computes the B representation (T,BCOEF,N,K) of a
C         cubic spline (K=4) which interpolates data (X(I)),Y(I))),
C         I=1,NDATA.  Parameters IBCL, IBCR, FBCL, FBCR allow the
C         specification of the spline first or second derivative at
C         both X(1) and X(NDATA).  When this data is not specified
C         by the problem, it is common practice to use a natural
C         spline by setting second derivatives at X(1) and X(NDATA)
C         to zero (IBCL=IBCR=2,FBCL=FBCR=0.0).  The spline is defined on
C         T(4) .LE. X .LE. T(N+1) with (ordered) interior knots at X(I))
C         values where N=NDATA+2.  The knots T(1), T(2), T(3) lie to
C         the left of T(4)=X(1) and the knots T(N+2), T(N+3), T(N+4)
C         lie to the right of T(N+1)=X(NDATA) in increasing order.  If
C         no extrapolation outside (X(1),X(NDATA)) is anticipated, the
C         knots T(1)=T(2)=T(3)=T(4)=X(1) and T(N+2)=T(N+3)=T(N+4)=
C         T(N+1)=X(NDATA) can be specified by KNTOPT=1.  KNTOPT=2
C         selects a knot placement for T(1), T(2), T(3) to make the
C         first 7 knots symmetric about T(4)=X(1) and similarly for
C         T(N+2), T(N+3), T(N+4) about T(N+1)=X(NDATA).  KNTOPT=3
C         allows the user to make his own selection, in increasing
C         order, for T(1), T(2), T(3) to the left of X(1) and T(N+2),
C         T(N+3), T(N+4) to the right of X(NDATA) in the work array
C         W(1) through W(6).  In any case, the interpolation on
C         T(4) .LE. X .LE. T(N+1) by using function BVALU is unique
C         for given boundary conditions.
C
C     Description of Arguments
C         Input
C           X      - X vector of abscissae of length NDATA, distinct
C                    and in increasing order
C           Y      - Y vector of ordinates of length NDATA
C           NDATA  - number of data points, NDATA .GE. 2
C           IBCL   - selection parameter for left boundary condition
C                    IBCL = 1 constrain the first derivative at
C                             X(1) to FBCL
C                         = 2 constrain the second derivative at
C                             X(1) to FBCL
C           IBCR   - selection parameter for right boundary condition
C                    IBCR = 1 constrain first derivative at
C                             X(NDATA) to FBCR
C                    IBCR = 2 constrain second derivative at
C                             X(NDATA) to FBCR
C           FBCL   - left boundary values governed by IBCL
C           FBCR   - right boundary values governed by IBCR
C           KNTOPT - knot selection parameter
C                    KNTOPT = 1 sets knot multiplicity at T(4) and
C                               T(N+1) to 4
C                           = 2 sets a symmetric placement of knots
C                               about T(4) and T(N+1)
C                           = 3 sets TNP)=WNP) and T(N+1+I)=w(3+I),I=1,3
C                               where WNP),I=1,6 is supplied by the user
C           W      - work array of dimension at least 5*(NDATA+2)
C                    if KNTOPT=3, then W(1),W(2),W(3) are knot values to
C                    the left of X(1) and W(4),W(5),W(6) are knot
C                    values to the right of X(NDATA) in increasing
C                    order to be supplied by the user
C
C         Output
C           T      - knot array of length N+4
C           BCOEF  - B-spline coefficient array of length N
C           N      - number of coefficients, N=NDATA+2
C           K      - order of spline, K=4
C
C     Error Conditions
C         Improper  input is a fatal error
C         Singular system of equations is a fatal error
C
C***REFERENCES  D. E. Amos, Computation with splines and B-splines,
C                 Report SAND78-1968, Sandia Laboratories, March 1979.
C               Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C               Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  BNFAC, BNSLV, BSPVD, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BINT4
C
      INTEGER I, IBCL, IBCR, IFLAG, ILB, ILEFT, IT, IUB, IW, IWP, J,
     1 JW, K, KNTOPT, N, NDATA, NDM, NP, NWROW
      REAL BCOEF,FBCL,FBCR,T, TOL,TXN,TX1,VNIKX,W,WDTOL,WORK,X, XL,
     1 Y
      REAL R1MACH
      DIMENSION X(*), Y(*), T(*), BCOEF(*), W(5,*), VNIKX(4,4), WORK(15)
C***FIRST EXECUTABLE STATEMENT  BINT4
      WDTOL = R1MACH(4)
      TOL = SQRT(WDTOL)
      IF (NDATA.LT.2) GO TO 200
      NDM = NDATA - 1
      DO 10 I=1,NDM
        IF (X(I).GE.X(I+1)) GO TO 210
   10 CONTINUE
      IF (IBCL.LT.1 .OR. IBCL.GT.2) GO TO 220
      IF (IBCR.LT.1 .OR. IBCR.GT.2) GO TO 230
      IF (KNTOPT.LT.1 .OR. KNTOPT.GT.3) GO TO 240
      K = 4
      N = NDATA + 2
      NP = N + 1
      DO 20 I=1,NDATA
        T(I+3) = X(I)
   20 CONTINUE
      GO TO (30, 50, 90), KNTOPT
C     SET UP KNOT ARRAY WITH MULTIPLICITY 4 AT X(1) AND X(NDATA)
   30 CONTINUE
      DO 40 I=1,3
        T(4-I) = X(1)
        T(NP+I) = X(NDATA)
   40 CONTINUE
      GO TO 110
C     SET UP KNOT ARRAY WITH SYMMETRIC PLACEMENT ABOUT END POINTS
   50 CONTINUE
      IF (NDATA.GT.3) GO TO 70
      XL = (X(NDATA)-X(1))/3.0E0
      DO 60 I=1,3
        T(4-I) = T(5-I) - XL
        T(NP+I) = T(NP+I-1) + XL
   60 CONTINUE
      GO TO 110
   70 CONTINUE
      TX1 = X(1) + X(1)
      TXN = X(NDATA) + X(NDATA)
      DO 80 I=1,3
        T(4-I) = TX1 - X(I+1)
        T(NP+I) = TXN - X(NDATA-I)
   80 CONTINUE
      GO TO 110
C     SET UP KNOT ARRAY LESS THAN X(1) AND GREATER THAN X(NDATA) TO BE
C     SUPPLIED BY USER IN WORK LOCATIONS W(1) THROUGH W(6) WHEN KNTOPT=3
   90 CONTINUE
      DO 100 I=1,3
        T(4-I) = W(4-I,1)
        JW = MAX(1,I-1)
        IW = MOD(I+2,5)+1
        T(NP+I) = W(IW,JW)
        IF (T(4-I).GT.T(5-I)) GO TO 250
        IF (T(NP+I).LT.T(NP+I-1)) GO TO 250
  100 CONTINUE
  110 CONTINUE
C
      DO 130 I=1,5
        DO 120 J=1,N
          W(I,J) = 0.0E0
  120   CONTINUE
  130 CONTINUE
C     SET UP LEFT INTERPOLATION POINT AND LEFT BOUNDARY CONDITION FOR
C     RIGHT LIMITS
      IT = IBCL + 1
      CALL BSPVD(T, K, IT, X(1), K, 4, VNIKX, WORK)
      IW = 0
      IF (ABS(VNIKX(3,1)).LT.TOL) IW = 1
      DO 140 J=1,3
        W(J+1,4-J) = VNIKX(4-J,IT)
        W(J,4-J) = VNIKX(4-J,1)
  140 CONTINUE
      BCOEF(1) = Y(1)
      BCOEF(2) = FBCL
C     SET UP INTERPOLATION EQUATIONS FOR POINTS I=2 TO I=NDATA-1
      ILEFT = 4
      IF (NDM.LT.2) GO TO 170
      DO 160 I=2,NDM
        ILEFT = ILEFT + 1
        CALL BSPVD(T, K, 1, X(I), ILEFT, 4, VNIKX, WORK)
        DO 150 J=1,3
          W(J+1,3+I-J) = VNIKX(4-J,1)
  150   CONTINUE
        BCOEF(I+1) = Y(I)
  160 CONTINUE
C     SET UP RIGHT INTERPOLATION POINT AND RIGHT BOUNDARY CONDITION FOR
C     LEFT LIMITS(ILEFT IS ASSOCIATED WITH T(N)=X(NDATA-1))
  170 CONTINUE
      IT = IBCR + 1
      CALL BSPVD(T, K, IT, X(NDATA), ILEFT, 4, VNIKX, WORK)
      JW = 0
      IF (ABS(VNIKX(2,1)).LT.TOL) JW = 1
      DO 180 J=1,3
        W(J+1,3+NDATA-J) = VNIKX(5-J,IT)
        W(J+2,3+NDATA-J) = VNIKX(5-J,1)
  180 CONTINUE
      BCOEF(N-1) = FBCR
      BCOEF(N) = Y(NDATA)
C     SOLVE SYSTEM OF EQUATIONS
      ILB = 2 - JW
      IUB = 2 - IW
      NWROW = 5
      IWP = IW + 1
      CALL BNFAC(W(IWP,1), NWROW, N, ILB, IUB, IFLAG)
      IF (IFLAG.EQ.2) GO TO 190
      CALL BNSLV(W(IWP,1), NWROW, N, ILB, IUB, BCOEF)
      RETURN
C
C
  190 CONTINUE
      CALL XERMSG ('SLATEC', 'BINT4',
     +   'THE SYSTEM OF EQUATIONS IS SINGULAR', 2, 1)
      RETURN
  200 CONTINUE
      CALL XERMSG ('SLATEC', 'BINT4', 'NDATA IS LESS THAN 2', 2, 1)
      RETURN
  210 CONTINUE
      CALL XERMSG ('SLATEC', 'BINT4',
     +   'X VALUES ARE NOT DISTINCT OR NOT ORDERED', 2, 1)
      RETURN
  220 CONTINUE
      CALL XERMSG ('SLATEC', 'BINT4', 'IBCL IS NOT 1 OR 2', 2, 1)
      RETURN
  230 CONTINUE
      CALL XERMSG ('SLATEC', 'BINT4', 'IBCR IS NOT 1 OR 2', 2, 1)
      RETURN
  240 CONTINUE
      CALL XERMSG ('SLATEC', 'BINT4', 'KNTOPT IS NOT 1, 2, OR 3', 2, 1)
      RETURN
  250 CONTINUE
      CALL XERMSG ('SLATEC', 'BINT4',
     +   'KNOT INPUT THROUGH W ARRAY IS NOT ORDERED PROPERLY', 2, 1)
      RETURN
      END
