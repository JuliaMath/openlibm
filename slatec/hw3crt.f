*DECK HW3CRT
      SUBROUTINE HW3CRT (XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M,
     +   MBDCND, BDYS, BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, ELMBDA,
     +   LDIMF, MDIMF, F, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HW3CRT
C***PURPOSE  Solve the standard seven-point finite difference
C            approximation to the Helmholtz equation in Cartesian
C            coordinates.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HW3CRT-S)
C***KEYWORDS  CARTESIAN, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine HW3CRT solves the standard seven-point finite
C     difference approximation to the Helmholtz equation in Cartesian
C     coordinates:
C
C         (d/dX)(dU/dX) + (d/dY)(dU/dY) + (d/dZ)(dU/dZ)
C
C                    + LAMBDA*U = F(X,Y,Z) .
C
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C    * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C
C            * * * * * *   On Input    * * * * * *
C
C     XS,XF
C        The range of X, i.e. XS .LE. X .LE. XF .
C        XS must be less than XF.
C
C     L
C        The number of panels into which the interval (XS,XF) is
C        subdivided.  Hence, there will be L+1 grid points in the
C        X-direction given by X(I) = XS+(I-1)DX for I=1,2,...,L+1,
C        where DX = (XF-XS)/L is the panel width.  L must be at
C        least 5 .
C
C     LBDCND
C        Indicates the type of boundary conditions at X = XS and X = XF.
C
C        = 0  If the solution is periodic in X, i.e.
C             U(L+I,J,K) = U(I,J,K).
C        = 1  If the solution is specified at X = XS and X = XF.
C        = 2  If the solution is specified at X = XS and the derivative
C             of the solution with respect to X is specified at X = XF.
C        = 3  If the derivative of the solution with respect to X is
C             specified at X = XS and X = XF.
C        = 4  If the derivative of the solution with respect to X is
C             specified at X = XS and the solution is specified at X=XF.
C
C     BDXS
C        A two-dimensional array that specifies the values of the
C        derivative of the solution with respect to X at X = XS.
C        when LBDCND = 3 or 4,
C
C             BDXS(J,K) = (d/dX)U(XS,Y(J),Z(K)), J=1,2,...,M+1,
C                                                K=1,2,...,N+1.
C
C        When LBDCND has any other value, BDXS is a dummy variable.
C        BDXS must be dimensioned at least (M+1)*(N+1).
C
C     BDXF
C        A two-dimensional array that specifies the values of the
C        derivative of the solution with respect to X at X = XF.
C        When LBDCND = 2 or 3,
C
C             BDXF(J,K) = (d/dX)U(XF,Y(J),Z(K)), J=1,2,...,M+1,
C                                                K=1,2,...,N+1.
C
C        When LBDCND has any other value, BDXF is a dummy variable.
C        BDXF must be dimensioned at least (M+1)*(N+1).
C
C     YS,YF
C        The range of Y, i.e. YS .LE. Y .LE. YF.
C        YS must be less than YF.
C
C     M
C        The number of panels into which the interval (YS,YF) is
C        subdivided.  Hence, there will be M+1 grid points in the
C        Y-direction given by Y(J) = YS+(J-1)DY for J=1,2,...,M+1,
C        where DY = (YF-YS)/M is the panel width.  M must be at
C        least 5 .
C
C     MBDCND
C        Indicates the type of boundary conditions at Y = YS and Y = YF.
C
C        = 0  If the solution is periodic in Y, i.e.
C             U(I,M+J,K) = U(I,J,K).
C        = 1  If the solution is specified at Y = YS and Y = YF.
C        = 2  If the solution is specified at Y = YS and the derivative
C             of the solution with respect to Y is specified at Y = YF.
C        = 3  If the derivative of the solution with respect to Y is
C             specified at Y = YS and Y = YF.
C        = 4  If the derivative of the solution with respect to Y is
C             specified at Y = YS and the solution is specified at Y=YF.
C
C     BDYS
C        A two-dimensional array that specifies the values of the
C        derivative of the solution with respect to Y at Y = YS.
C        When MBDCND = 3 or 4,
C
C             BDYS(I,K) = (d/dY)U(X(I),YS,Z(K)), I=1,2,...,L+1,
C                                                K=1,2,...,N+1.
C
C        When MBDCND has any other value, BDYS is a dummy variable.
C        BDYS must be dimensioned at least (L+1)*(N+1).
C
C     BDYF
C        A two-dimensional array that specifies the values of the
C        derivative of the solution with respect to Y at Y = YF.
C        When MBDCND = 2 or 3,
C
C             BDYF(I,K) = (d/dY)U(X(I),YF,Z(K)), I=1,2,...,L+1,
C                                                K=1,2,...,N+1.
C
C        When MBDCND has any other value, BDYF is a dummy variable.
C        BDYF must be dimensioned at least (L+1)*(N+1).
C
C     ZS,ZF
C        The range of Z, i.e. ZS .LE. Z .LE. ZF.
C        ZS must be less than ZF.
C
C     N
C        The number of panels into which the interval (ZS,ZF) is
C        subdivided.  Hence, there will be N+1 grid points in the
C        Z-direction given by Z(K) = ZS+(K-1)DZ for K=1,2,...,N+1,
C        where DZ = (ZF-ZS)/N is the panel width.  N must be at least 5.
C
C     NBDCND
C        Indicates the type of boundary conditions at Z = ZS and Z = ZF.
C
C        = 0  If the solution is periodic in Z, i.e.
C             U(I,J,N+K) = U(I,J,K).
C        = 1  If the solution is specified at Z = ZS and Z = ZF.
C        = 2  If the solution is specified at Z = ZS and the derivative
C             of the solution with respect to Z is specified at Z = ZF.
C        = 3  If the derivative of the solution with respect to Z is
C             specified at Z = ZS and Z = ZF.
C        = 4  If the derivative of the solution with respect to Z is
C             specified at Z = ZS and the solution is specified at Z=ZF.
C
C     BDZS
C        A two-dimensional array that specifies the values of the
C        derivative of the solution with respect to Z at Z = ZS.
C        When NBDCND = 3 or 4,
C
C             BDZS(I,J) = (d/dZ)U(X(I),Y(J),ZS), I=1,2,...,L+1,
C                                                J=1,2,...,M+1.
C
C        When NBDCND has any other value, BDZS is a dummy variable.
C        BDZS must be dimensioned at least (L+1)*(M+1).
C
C     BDZF
C        A two-dimensional array that specifies the values of the
C        derivative of the solution with respect to Z at Z = ZF.
C        When NBDCND = 2 or 3,
C
C             BDZF(I,J) = (d/dZ)U(X(I),Y(J),ZF), I=1,2,...,L+1,
C                                                J=1,2,...,M+1.
C
C        When NBDCND has any other value, BDZF is a dummy variable.
C        BDZF must be dimensioned at least (L+1)*(M+1).
C
C     ELMBDA
C        The constant LAMBDA in the Helmholtz equation. If
C        LAMBDA .GT. 0, a solution may not exist.  However, HW3CRT will
C        attempt to find a solution.
C
C     F
C        A three-dimensional array that specifies the values of the
C        right side of the Helmholtz equation and boundary values (if
C        any).  For I=2,3,...,L, J=2,3,...,M, and K=2,3,...,N
C
C                   F(I,J,K) = F(X(I),Y(J),Z(K)).
C
C        On the boundaries F is defined by
C
C        LBDCND      F(1,J,K)         F(L+1,J,K)
C        ------   ---------------   ---------------
C
C          0      F(XS,Y(J),Z(K))   F(XS,Y(J),Z(K))
C          1      U(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
C          2      U(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))   J=1,2,...,M+1
C          3      F(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))   K=1,2,...,N+1
C          4      F(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
C
C        MBDCND      F(I,1,K)         F(I,M+1,K)
C        ------   ---------------   ---------------
C
C          0      F(X(I),YS,Z(K))   F(X(I),YS,Z(K))
C          1      U(X(I),YS,Z(K))   U(X(I),YF,Z(K))
C          2      U(X(I),YS,Z(K))   F(X(I),YF,Z(K))   I=1,2,...,L+1
C          3      F(X(I),YS,Z(K))   F(X(I),YF,Z(K))   K=1,2,...,N+1
C          4      F(X(I),YS,Z(K))   U(X(I),YF,Z(K))
C
C        NBDCND      F(I,J,1)         F(I,J,N+1)
C        ------   ---------------   ---------------
C
C          0      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZS)
C          1      U(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
C          2      U(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)   I=1,2,...,L+1
C          3      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)   J=1,2,...,M+1
C          4      F(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
C
C        F must be dimensioned at least (L+1)*(M+1)*(N+1).
C
C        NOTE:
C
C        If the table calls for both the solution U and the right side F
C        on a boundary, then the solution must be specified.
C
C     LDIMF
C        The row (or first) dimension of the arrays F,BDYS,BDYF,BDZS,
C        and BDZF as it appears in the program calling HW3CRT. this
C        parameter is used to specify the variable dimension of these
C        arrays.  LDIMF must be at least L+1.
C
C     MDIMF
C        The column (or second) dimension of the array F and the row (or
C        first) dimension of the arrays BDXS and BDXF as it appears in
C        the program calling HW3CRT.  This parameter is used to specify
C        the variable dimension of these arrays.
C        MDIMF must be at least M+1.
C
C     W
C        A one-dimensional array that must be provided by the user for
C        work space.  The length of W must be at least 30 + L + M + 5*N
C        + MAX(L,M,N) + 7*(INT((L+1)/2) + INT((M+1)/2))
C
C
C            * * * * * *   On Output   * * * * * *
C
C     F
C        Contains the solution U(I,J,K) of the finite difference
C        approximation for the grid point (X(I),Y(J),Z(K)) for
C        I=1,2,...,L+1, J=1,2,...,M+1, and K=1,2,...,N+1.
C
C     PERTRB
C        If a combination of periodic or derivative boundary conditions
C        is specified for a Poisson equation (LAMBDA = 0), a solution
C        may not exist.  PERTRB is a constant, calculated and subtracted
C        from F, which ensures that a solution exists.  PWSCRT then
C        computes this solution, which is a least squares solution to
C        the original approximation.  This solution is not unique and is
C        unnormalized.  The value of PERTRB should be small compared to
C        the right side F.  Otherwise, a solution is obtained to an
C        essentially different problem.  This comparison should always
C        be made to insure that a meaningful solution has been obtained.
C
C     IERROR
C        An error flag that indicates invalid input parameters.  Except
C        for numbers 0 and 12, a solution is not attempted.
C
C        =  0  No error
C        =  1  XS .GE. XF
C        =  2  L .LT. 5
C        =  3  LBDCND .LT. 0 .OR. LBDCND .GT. 4
C        =  4  YS .GE. YF
C        =  5  M .LT. 5
C        =  6  MBDCND .LT. 0 .OR. MBDCND .GT. 4
C        =  7  ZS .GE. ZF
C        =  8  N .LT. 5
C        =  9  NBDCND .LT. 0 .OR. NBDCND .GT. 4
C        = 10  LDIMF .LT. L+1
C        = 11  MDIMF .LT. M+1
C        = 12  LAMBDA .GT. 0
C
C        Since this is the only means of indicating a possibly incorrect
C        call to HW3CRT, the user should test IERROR after the call.
C
C *Long Description:
C
C    * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   BDXS(MDIMF,N+1),BDXF(MDIMF,N+1),BDYS(LDIMF,N+1),
C     Arguments      BDYF(LDIMF,N+1),BDZS(LDIMF,M+1),BDZF(LDIMF,M+1),
C                    F(LDIMF,MDIMF,N+1),W(see argument list)
C
C     Latest         December 1, 1978
C     Revision
C
C     Subprograms    HW3CRT,POIS3D,POS3D1,TRIDQ,RFFTI,RFFTF,RFFTF1,
C     Required       RFFTB,RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,
C                    COSQF1,COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,
C                    CFFTI1,CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,
C                    CFFTF,CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,
C                    PIMACH
C
C     Special        NONE
C     Conditions
C
C     Common         NONE
C     Blocks
C
C     I/O            NONE
C
C     Precision      Single
C
C     Specialist     Roland Sweet
C
C     Language       FORTRAN
C
C     History        Written by Roland Sweet at NCAR in July 1977
C
C     Algorithm      This subroutine defines the finite difference
C                    equations, incorporates boundary data, and
C                    adjusts the right side of singular systems and
C                    then calls POIS3D to solve the system.
C
C     Space          7862(decimal) = 17300(octal) locations on the
C     Required       NCAR Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine HW3CRT is roughly proportional
C                    to L*M*N*(log2(L)+log2(M)+5), but also depends on
C                    input parameters LBDCND and MBDCND.  Some typical
C                    values are listed in the table below.
C                       The solution process employed results in a loss
C                    of no more than three significant digits for L,M
C                    and N as large as 32.  More detailed information
C                    about accuracy can be found in the documentation
C                    for subroutine POIS3D which is the routine that
C                    actually solves the finite difference equations.
C
C
C                       L(=M=N)     LBDCND(=MBDCND=NBDCND)      T(MSECS)
C                       -------     ----------------------      --------
C
C                         16                  0                    300
C                         16                  1                    302
C                         16                  3                    348
C                         32                  0                   1925
C                         32                  1                   1929
C                         32                  3                   2109
C
C     Portability    American National Standards Institute FORTRAN.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C     Required       COS,SIN,ATAN
C     Resident
C     Routines
C
C     Reference      NONE
C
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  POIS3D
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  HW3CRT
C
C
      DIMENSION       BDXS(MDIMF,*)          ,BDXF(MDIMF,*)          ,
     1                BDYS(LDIMF,*)          ,BDYF(LDIMF,*)          ,
     2                BDZS(LDIMF,*)          ,BDZF(LDIMF,*)          ,
     3                F(LDIMF,MDIMF,*)       ,W(*)
C***FIRST EXECUTABLE STATEMENT  HW3CRT
      IERROR = 0
      IF (XF .LE. XS) IERROR = 1
      IF (L .LT. 5) IERROR = 2
      IF (LBDCND.LT.0 .OR. LBDCND.GT.4) IERROR = 3
      IF (YF .LE. YS) IERROR = 4
      IF (M .LT. 5) IERROR = 5
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) IERROR = 6
      IF (ZF .LE. ZS) IERROR = 7
      IF (N .LT. 5) IERROR = 8
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR = 9
      IF (LDIMF .LT. L+1) IERROR = 10
      IF (MDIMF .LT. M+1) IERROR = 11
      IF (IERROR .NE. 0) GO TO 188
      DY = (YF-YS)/M
      TWBYDY = 2./DY
      C2 = 1./(DY**2)
      MSTART = 1
      MSTOP = M
      MP1 = M+1
      MP = MBDCND+1
      GO TO (104,101,101,102,102),MP
  101 MSTART = 2
  102 GO TO (104,104,103,103,104),MP
  103 MSTOP = MP1
  104 MUNK = MSTOP-MSTART+1
      DZ = (ZF-ZS)/N
      TWBYDZ = 2./DZ
      NP = NBDCND+1
      C3 = 1./(DZ**2)
      NP1 = N+1
      NSTART = 1
      NSTOP = N
      GO TO (108,105,105,106,106),NP
  105 NSTART = 2
  106 GO TO (108,108,107,107,108),NP
  107 NSTOP = NP1
  108 NUNK = NSTOP-NSTART+1
      LP1 = L+1
      DX = (XF-XS)/L
      C1 = 1./(DX**2)
      TWBYDX = 2./DX
      LP = LBDCND+1
      LSTART = 1
      LSTOP = L
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      GO TO (122,109,109,112,112),LP
  109 LSTART = 2
      DO 111 J=MSTART,MSTOP
         DO 110 K=NSTART,NSTOP
            F(2,J,K) = F(2,J,K)-C1*F(1,J,K)
  110    CONTINUE
  111 CONTINUE
      GO TO 115
  112 DO 114 J=MSTART,MSTOP
         DO 113 K=NSTART,NSTOP
            F(1,J,K) = F(1,J,K)+TWBYDX*BDXS(J,K)
  113    CONTINUE
  114 CONTINUE
  115 GO TO (122,116,119,119,116),LP
  116 DO 118 J=MSTART,MSTOP
         DO 117 K=NSTART,NSTOP
            F(L,J,K) = F(L,J,K)-C1*F(LP1,J,K)
  117    CONTINUE
  118 CONTINUE
      GO TO 122
  119 LSTOP = LP1
      DO 121 J=MSTART,MSTOP
         DO 120 K=NSTART,NSTOP
            F(LP1,J,K) = F(LP1,J,K)-TWBYDX*BDXF(J,K)
  120    CONTINUE
  121 CONTINUE
  122 LUNK = LSTOP-LSTART+1
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (136,123,123,126,126),MP
  123 DO 125 I=LSTART,LSTOP
         DO 124 K=NSTART,NSTOP
            F(I,2,K) = F(I,2,K)-C2*F(I,1,K)
  124    CONTINUE
  125 CONTINUE
      GO TO 129
  126 DO 128 I=LSTART,LSTOP
         DO 127 K=NSTART,NSTOP
            F(I,1,K) = F(I,1,K)+TWBYDY*BDYS(I,K)
  127    CONTINUE
  128 CONTINUE
  129 GO TO (136,130,133,133,130),MP
  130 DO 132 I=LSTART,LSTOP
         DO 131 K=NSTART,NSTOP
            F(I,M,K) = F(I,M,K)-C2*F(I,MP1,K)
  131    CONTINUE
  132 CONTINUE
      GO TO 136
  133 DO 135 I=LSTART,LSTOP
         DO 134 K=NSTART,NSTOP
            F(I,MP1,K) = F(I,MP1,K)-TWBYDY*BDYF(I,K)
  134    CONTINUE
  135 CONTINUE
  136 CONTINUE
C
C     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
C
      GO TO (150,137,137,140,140),NP
  137 DO 139 I=LSTART,LSTOP
         DO 138 J=MSTART,MSTOP
            F(I,J,2) = F(I,J,2)-C3*F(I,J,1)
  138    CONTINUE
  139 CONTINUE
      GO TO 143
  140 DO 142 I=LSTART,LSTOP
         DO 141 J=MSTART,MSTOP
            F(I,J,1) = F(I,J,1)+TWBYDZ*BDZS(I,J)
  141    CONTINUE
  142 CONTINUE
  143 GO TO (150,144,147,147,144),NP
  144 DO 146 I=LSTART,LSTOP
         DO 145 J=MSTART,MSTOP
            F(I,J,N) = F(I,J,N)-C3*F(I,J,NP1)
  145    CONTINUE
  146 CONTINUE
      GO TO 150
  147 DO 149 I=LSTART,LSTOP
         DO 148 J=MSTART,MSTOP
            F(I,J,NP1) = F(I,J,NP1)-TWBYDZ*BDZF(I,J)
  148    CONTINUE
  149 CONTINUE
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
  150 CONTINUE
      IWB = NUNK+1
      IWC = IWB+NUNK
      IWW = IWC+NUNK
      DO 151 K=1,NUNK
         I = IWC+K-1
         W(K) = C3
         W(I) = C3
         I = IWB+K-1
         W(I) = -2.*C3+ELMBDA
  151 CONTINUE
      GO TO (155,155,153,152,152),NP
  152 W(IWC) = 2.*C3
  153 GO TO (155,155,154,154,155),NP
  154 W(IWB-1) = 2.*C3
  155 CONTINUE
      PERTRB = 0.
C
C     FOR SINGULAR PROBLEMS ADJUST DATA TO INSURE A SOLUTION WILL EXIST.
C
      GO TO (156,172,172,156,172),LP
  156 GO TO (157,172,172,157,172),MP
  157 GO TO (158,172,172,158,172),NP
  158 IF (ELMBDA) 172,160,159
  159 IERROR = 12
      GO TO 172
  160 CONTINUE
      MSTPM1 = MSTOP-1
      LSTPM1 = LSTOP-1
      NSTPM1 = NSTOP-1
      XLP = (2+LP)/3
      YLP = (2+MP)/3
      ZLP = (2+NP)/3
      S1 = 0.
      DO 164 K=2,NSTPM1
         DO 162 J=2,MSTPM1
            DO 161 I=2,LSTPM1
               S1 = S1+F(I,J,K)
  161       CONTINUE
            S1 = S1+(F(1,J,K)+F(LSTOP,J,K))/XLP
  162    CONTINUE
         S2 = 0.
         DO 163 I=2,LSTPM1
            S2 = S2+F(I,1,K)+F(I,MSTOP,K)
  163    CONTINUE
         S2 = (S2+(F(1,1,K)+F(1,MSTOP,K)+F(LSTOP,1,K)+F(LSTOP,MSTOP,K))/
     1                                                          XLP)/YLP
         S1 = S1+S2
  164 CONTINUE
      S = (F(1,1,1)+F(LSTOP,1,1)+F(1,1,NSTOP)+F(LSTOP,1,NSTOP)+
     1    F(1,MSTOP,1)+F(LSTOP,MSTOP,1)+F(1,MSTOP,NSTOP)+
     2                                   F(LSTOP,MSTOP,NSTOP))/(XLP*YLP)
      DO 166 J=2,MSTPM1
         DO 165 I=2,LSTPM1
            S = S+F(I,J,1)+F(I,J,NSTOP)
  165    CONTINUE
  166 CONTINUE
      S2 = 0.
      DO 167 I=2,LSTPM1
         S2 = S2+F(I,1,1)+F(I,1,NSTOP)+F(I,MSTOP,1)+F(I,MSTOP,NSTOP)
  167 CONTINUE
      S = S2/YLP+S
      S2 = 0.
      DO 168 J=2,MSTPM1
         S2 = S2+F(1,J,1)+F(1,J,NSTOP)+F(LSTOP,J,1)+F(LSTOP,J,NSTOP)
  168 CONTINUE
      S = S2/XLP+S
      PERTRB = (S/ZLP+S1)/((LUNK+1.-XLP)*(MUNK+1.-YLP)*
     1                                              (NUNK+1.-ZLP))
      DO 171 I=1,LUNK
         DO 170 J=1,MUNK
            DO 169 K=1,NUNK
               F(I,J,K) = F(I,J,K)-PERTRB
  169       CONTINUE
  170    CONTINUE
  171 CONTINUE
  172 CONTINUE
      NPEROD = 0
      IF (NBDCND .EQ. 0) GO TO 173
      NPEROD = 1
      W(1) = 0.
      W(IWW-1) = 0.
  173 CONTINUE
      CALL POIS3D (LBDCND,LUNK,C1,MBDCND,MUNK,C2,NPEROD,NUNK,W,W(IWB),
     1             W(IWC),LDIMF,MDIMF,F(LSTART,MSTART,NSTART),IR,W(IWW))
C
C     FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS.
C
      IF (LP .NE. 1) GO TO 180
      IF (MP .NE. 1) GO TO 175
      DO 174 K=NSTART,NSTOP
         F(1,MP1,K) = F(1,1,K)
  174 CONTINUE
      MSTOP = MP1
  175 IF (NP .NE. 1) GO TO 177
      DO 176 J=MSTART,MSTOP
         F(1,J,NP1) = F(1,J,1)
  176 CONTINUE
      NSTOP = NP1
  177 DO 179 J=MSTART,MSTOP
         DO 178 K=NSTART,NSTOP
            F(LP1,J,K) = F(1,J,K)
  178    CONTINUE
  179 CONTINUE
  180 CONTINUE
      IF (MP .NE. 1) GO TO 185
      IF (NP .NE. 1) GO TO 182
      DO 181 I=LSTART,LSTOP
         F(I,1,NP1) = F(I,1,1)
  181 CONTINUE
      NSTOP = NP1
  182 DO 184 I=LSTART,LSTOP
         DO 183 K=NSTART,NSTOP
            F(I,MP1,K) = F(I,1,K)
  183    CONTINUE
  184 CONTINUE
  185 CONTINUE
      IF (NP .NE. 1) GO TO 188
      DO 187 I=LSTART,LSTOP
         DO 186 J=MSTART,MSTOP
            F(I,J,NP1) = F(I,J,1)
  186    CONTINUE
  187 CONTINUE
  188 CONTINUE
      RETURN
      END
