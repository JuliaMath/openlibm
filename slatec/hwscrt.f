*DECK HWSCRT
      SUBROUTINE HWSCRT (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND,
     +   BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HWSCRT
C***PURPOSE  Solves the standard five-point finite difference
C            approximation to the Helmholtz equation in Cartesian
C            coordinates.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HWSCRT-S)
C***KEYWORDS  CARTESIAN, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine HWSCRT solves the standard five-point finite
C     difference approximation to the Helmholtz equation in Cartesian
C     coordinates:
C
C          (d/dX)(dU/dX) + (d/dY)(dU/dY) + LAMBDA*U = F(X,Y).
C
C
C
C     * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C             * * * * * *   On Input    * * * * * *
C
C     A,B
C       The range of X, i.e., A .LE. X .LE. B.  A must be less than B.
C
C     M
C       The number of panels into which the interval (A,B) is
C       subdivided.  Hence, there will be M+1 grid points in the
C       X-direction given by X(I) = A+(I-1)DX for I = 1,2,...,M+1,
C       where DX = (B-A)/M is the panel width. M must be greater than 3.
C
C     MBDCND
C       Indicates the type of boundary conditions at X = A and X = B.
C
C       = 0  If the solution is periodic in X, i.e., U(I,J) = U(M+I,J).
C       = 1  If the solution is specified at X = A and X = B.
C       = 2  If the solution is specified at X = A and the derivative of
C            the solution with respect to X is specified at X = B.
C       = 3  If the derivative of the solution with respect to X is
C            specified at X = A and X = B.
C       = 4  If the derivative of the solution with respect to X is
C            specified at X = A and the solution is specified at X = B.
C
C     BDA
C       A one-dimensional array of length N+1 that specifies the values
C       of the derivative of the solution with respect to X at X = A.
C       When MBDCND = 3 or 4,
C
C            BDA(J) = (d/dX)U(A,Y(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value, BDA is a dummy variable.
C
C     BDB
C       A one-dimensional array of length N+1 that specifies the values
C       of the derivative of the solution with respect to X at X = B.
C       When MBDCND = 2 or 3,
C
C            BDB(J) = (d/dX)U(B,Y(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value BDB is a dummy variable.
C
C     C,D
C       The range of Y, i.e., C .LE. Y .LE. D.  C must be less than D.
C
C     N
C       The number of panels into which the interval (C,D) is
C       subdivided.  Hence, there will be N+1 grid points in the
C       Y-direction given by Y(J) = C+(J-1)DY for J = 1,2,...,N+1, where
C       DY = (D-C)/N is the panel width.  N must be greater than 3.
C
C     NBDCND
C       Indicates the type of boundary conditions at Y = C and Y = D.
C
C       = 0  If the solution is periodic in Y, i.e., U(I,J) = U(I,N+J).
C       = 1  If the solution is specified at Y = C and Y = D.
C       = 2  If the solution is specified at Y = C and the derivative of
C            the solution with respect to Y is specified at Y = D.
C       = 3  If the derivative of the solution with respect to Y is
C            specified at Y = C and Y = D.
C       = 4  If the derivative of the solution with respect to Y is
C            specified at Y = C and the solution is specified at Y = D.
C
C     BDC
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to Y at Y = C.
C       When NBDCND = 3 or 4,
C
C            BDC(I) = (d/dY)U(X(I),C), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDC is a dummy variable.
C
C     BDD
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to Y at Y = D.
C       When NBDCND = 2 or 3,
C
C            BDD(I) = (d/dY)U(X(I),D), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDD is a dummy variable.
C
C     ELMBDA
C       The constant LAMBDA in the Helmholtz equation.  If
C       LAMBDA .GT. 0, a solution may not exist.  However, HWSCRT will
C       attempt to find a solution.
C
C     F
C       A two-dimensional array which specifies the values of the right
C       side of the Helmholtz equation and boundary values (if any).
C       For I = 2,3,...,M and J = 2,3,...,N
C
C            F(I,J) = F(X(I),Y(J)).
C
C       On the boundaries F is defined by
C
C            MBDCND     F(1,J)        F(M+1,J)
C            ------     ---------     --------
C
C              0        F(A,Y(J))     F(A,Y(J))
C              1        U(A,Y(J))     U(B,Y(J))
C              2        U(A,Y(J))     F(B,Y(J))     J = 1,2,...,N+1
C              3        F(A,Y(J))     F(B,Y(J))
C              4        F(A,Y(J))     U(B,Y(J))
C
C
C            NBDCND     F(I,1)        F(I,N+1)
C            ------     ---------     --------
C
C              0        F(X(I),C)     F(X(I),C)
C              1        U(X(I),C)     U(X(I),D)
C              2        U(X(I),C)     F(X(I),D)     I = 1,2,...,M+1
C              3        F(X(I),C)     F(X(I),D)
C              4        F(X(I),C)     U(X(I),D)
C
C       F must be dimensioned at least (M+1)*(N+1).
C
C       NOTE:
C
C       If the table calls for both the solution U and the right side F
C       at a corner then the solution must be specified.
C
C     IDIMF
C       The row (or first) dimension of the array F as it appears in the
C       program calling HWSCRT.  This parameter is used to specify the
C       variable dimension of F.  IDIMF must be at least M+1  .
C
C     W
C       A one-dimensional array that must be provided by the user for
C       work space.  W may require up to 4*(N+1) +
C       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of
C       locations used is computed by HWSCRT and is returned in location
C       W(1).
C
C
C             * * * * * *   On Output     * * * * * *
C
C     F
C       Contains the solution U(I,J) of the finite difference
C       approximation for the grid point (X(I),Y(J)), I = 1,2,...,M+1,
C       J = 1,2,...,N+1  .
C
C     PERTRB
C       If a combination of periodic or derivative boundary conditions
C       is specified for a Poisson equation (LAMBDA = 0), a solution may
C       not exist.  PERTRB is a constant, calculated and subtracted from
C       F, which ensures that a solution exists.  HWSCRT then computes
C       this solution, which is a least squares solution to the original
C       approximation.  This solution plus any constant is also a
C       solution.  Hence, the solution is not unique.  The value of
C       PERTRB should be small compared to the right side F.  Otherwise,
C       a solution is obtained to an essentially different problem.
C       This comparison should always be made to insure that a
C       meaningful solution has been obtained.
C
C     IERROR
C       An error flag that indicates invalid input parameters.  Except
C       for numbers 0 and 6, a solution is not attempted.
C
C       = 0  No error.
C       = 1  A .GE. B.
C       = 2  MBDCND .LT. 0 or MBDCND .GT. 4  .
C       = 3  C .GE. D.
C       = 4  N .LE. 3
C       = 5  NBDCND .LT. 0 or NBDCND .GT. 4  .
C       = 6  LAMBDA .GT. 0  .
C       = 7  IDIMF .LT. M+1  .
C       = 8  M .LE. 3
C
C       Since this is the only means of indicating a possibly incorrect
C       call to HWSCRT, the user should test IERROR after the call.
C
C     W
C       W(1) contains the required length of W.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C
C     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1),
C     Arguments      W(see argument list)
C
C     Latest         June 1, 1976
C     Revision
C
C     Subprograms    HWSCRT,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
C     Required       TRIX,TRI3,PIMACH
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
C     History        Standardized September 1, 1973
C                    Revised April 1, 1976
C
C     Algorithm      The routine defines the finite difference
C                    equations, incorporates boundary data, and adjusts
C                    the right side of singular systems and then calls
C                    GENBUN to solve the system.
C
C     Space          13110(octal) = 5704(decimal) locations on the NCAR
C     Required       Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine HWSCRT is roughly proportional
C                    to M*N*log2(N), but also depends on the input
C                    parameters NBDCND and MBDCND.  Some typical values
C                    are listed in the table below.
C                       The solution process employed results in a loss
C                    of no more than three significant digits for N and
C                    M as large as 64.  More detailed information about
C                    accuracy can be found in the documentation for
C                    subroutine GENBUN which is the routine that
C                    solves the finite difference equations.
C
C
C                       M(=N)    MBDCND    NBDCND    T(MSECS)
C                       -----    ------    ------    --------
C
C                        32        0         0          31
C                        32        1         1          23
C                        32        3         3          36
C                        64        0         0         128
C                        64        1         1          96
C                        64        3         3         142
C
C     Portability    American National Standards Institute FORTRAN.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C     Reference      Swarztrauber, P. and R. Sweet, 'Efficient FORTRAN
C                    Subprograms for The Solution Of Elliptic Equations'
C                    NCAR TN/IA-109, July, 1975, 138 pp.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
C                 subprograms for the solution of elliptic equations,
C                 NCAR TN/IA-109, July 1975, 138 pp.
C***ROUTINES CALLED  GENBUN
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  HWSCRT
C
C
      DIMENSION       F(IDIMF,*)
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
C***FIRST EXECUTABLE STATEMENT  HWSCRT
      IERROR = 0
      IF (A .GE. B) IERROR = 1
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) IERROR = 2
      IF (C .GE. D) IERROR = 3
      IF (N .LE. 3) IERROR = 4
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR = 5
      IF (IDIMF .LT. M+1) IERROR = 7
      IF (M .LE. 3) IERROR = 8
      IF (IERROR .NE. 0) RETURN
      NPEROD = NBDCND
      MPEROD = 0
      IF (MBDCND .GT. 0) MPEROD = 1
      DELTAX = (B-A)/M
      TWDELX = 2./DELTAX
      DELXSQ = 1./DELTAX**2
      DELTAY = (D-C)/N
      TWDELY = 2./DELTAY
      DELYSQ = 1./DELTAY**2
      NP = NBDCND+1
      NP1 = N+1
      MP = MBDCND+1
      MP1 = M+1
      NSTART = 1
      NSTOP = N
      NSKIP = 1
      GO TO (104,101,102,103,104),NP
  101 NSTART = 2
      GO TO 104
  102 NSTART = 2
  103 NSTOP = NP1
      NSKIP = 2
  104 NUNK = NSTOP-NSTART+1
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      MSTART = 1
      MSTOP = M
      MSKIP = 1
      GO TO (117,105,106,109,110),MP
  105 MSTART = 2
      GO TO 107
  106 MSTART = 2
      MSTOP = MP1
      MSKIP = 2
  107 DO 108 J=NSTART,NSTOP
         F(2,J) = F(2,J)-F(1,J)*DELXSQ
  108 CONTINUE
      GO TO 112
  109 MSTOP = MP1
      MSKIP = 2
  110 DO 111 J=NSTART,NSTOP
         F(1,J) = F(1,J)+BDA(J)*TWDELX
  111 CONTINUE
  112 GO TO (113,115),MSKIP
  113 DO 114 J=NSTART,NSTOP
         F(M,J) = F(M,J)-F(MP1,J)*DELXSQ
  114 CONTINUE
      GO TO 117
  115 DO 116 J=NSTART,NSTOP
         F(MP1,J) = F(MP1,J)-BDB(J)*TWDELX
  116 CONTINUE
  117 MUNK = MSTOP-MSTART+1
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (127,118,118,120,120),NP
  118 DO 119 I=MSTART,MSTOP
         F(I,2) = F(I,2)-F(I,1)*DELYSQ
  119 CONTINUE
      GO TO 122
  120 DO 121 I=MSTART,MSTOP
         F(I,1) = F(I,1)+BDC(I)*TWDELY
  121 CONTINUE
  122 GO TO (123,125),NSKIP
  123 DO 124 I=MSTART,MSTOP
         F(I,N) = F(I,N)-F(I,NP1)*DELYSQ
  124 CONTINUE
      GO TO 127
  125 DO 126 I=MSTART,MSTOP
         F(I,NP1) = F(I,NP1)-BDD(I)*TWDELY
  126 CONTINUE
C
C    MULTIPLY RIGHT SIDE BY DELTAY**2.
C
  127 DELYSQ = DELTAY*DELTAY
      DO 129 I=MSTART,MSTOP
         DO 128 J=NSTART,NSTOP
            F(I,J) = F(I,J)*DELYSQ
  128    CONTINUE
  129 CONTINUE
C
C     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = MUNK
      ID3 = ID2+MUNK
      ID4 = ID3+MUNK
      S = DELYSQ*DELXSQ
      ST2 = 2.*S
      DO 130 I=1,MUNK
         W(I) = S
         J = ID2+I
         W(J) = -ST2+ELMBDA*DELYSQ
         J = ID3+I
         W(J) = S
  130 CONTINUE
      IF (MP .EQ. 1) GO TO 131
      W(1) = 0.
      W(ID4) = 0.
  131 CONTINUE
      GO TO (135,135,132,133,134),MP
  132 W(ID2) = ST2
      GO TO 135
  133 W(ID2) = ST2
  134 W(ID3+1) = ST2
  135 CONTINUE
      PERTRB = 0.
      IF (ELMBDA) 144,137,136
  136 IERROR = 6
      GO TO 144
  137 IF ((NBDCND.EQ.0 .OR. NBDCND.EQ.3) .AND.
     1    (MBDCND.EQ.0 .OR. MBDCND.EQ.3)) GO TO 138
      GO TO 144
C
C     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
C     WILL EXIST.
C
  138 A1 = 1.
      A2 = 1.
      IF (NBDCND .EQ. 3) A2 = 2.
      IF (MBDCND .EQ. 3) A1 = 2.
      S1 = 0.
      MSP1 = MSTART+1
      MSTM1 = MSTOP-1
      NSP1 = NSTART+1
      NSTM1 = NSTOP-1
      DO 140 J=NSP1,NSTM1
         S = 0.
         DO 139 I=MSP1,MSTM1
            S = S+F(I,J)
  139    CONTINUE
         S1 = S1+S*A1+F(MSTART,J)+F(MSTOP,J)
  140 CONTINUE
      S1 = A2*S1
      S = 0.
      DO 141 I=MSP1,MSTM1
         S = S+F(I,NSTART)+F(I,NSTOP)
  141 CONTINUE
      S1 = S1+S*A1+F(MSTART,NSTART)+F(MSTART,NSTOP)+F(MSTOP,NSTART)+
     1     F(MSTOP,NSTOP)
      S = (2.+(NUNK-2)*A2)*(2.+(MUNK-2)*A1)
      PERTRB = S1/S
      DO 143 J=NSTART,NSTOP
         DO 142 I=MSTART,MSTOP
            F(I,J) = F(I,J)-PERTRB
  142    CONTINUE
  143 CONTINUE
      PERTRB = PERTRB/DELYSQ
C
C     SOLVE THE EQUATION.
C
  144 CALL GENBUN (NPEROD,NUNK,MPEROD,MUNK,W(1),W(ID2+1),W(ID3+1),
     1             IDIMF,F(MSTART,NSTART),IERR1,W(ID4+1))
      W(1) = W(ID4+1)+3*MUNK
C
C     FILL IN IDENTICAL VALUES WHEN HAVE PERIODIC BOUNDARY CONDITIONS.
C
      IF (NBDCND .NE. 0) GO TO 146
      DO 145 I=MSTART,MSTOP
         F(I,NP1) = F(I,1)
  145 CONTINUE
  146 IF (MBDCND .NE. 0) GO TO 148
      DO 147 J=NSTART,NSTOP
         F(MP1,J) = F(1,J)
  147 CONTINUE
      IF (NBDCND .EQ. 0) F(MP1,NP1) = F(1,NP1)
  148 CONTINUE
      RETURN
      END
