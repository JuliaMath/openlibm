*DECK HSTCRT
      SUBROUTINE HSTCRT (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND,
     +   BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HSTCRT
C***PURPOSE  Solve the standard five-point finite difference
C            approximation on a staggered grid to the Helmholtz equation
C            in Cartesian coordinates.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HSTCRT-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C      HSTCRT solves the standard five-point finite difference
C      approximation on a staggered grid to the Helmholtz equation in
C      Cartesian coordinates
C
C      (d/dX)(dU/dX) + (d/dY)(dU/dY) + LAMBDA*U = F(X,Y)
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C             * * * * * *   On Input    * * * * * *
C
C    A,B
C      The range of X, i.e. A .LE. X .LE. B.  A must be less than B.
C
C    M
C      The number of grid points in the interval (A,B).  The grid points
C      in the X-direction are given by X(I) = A + (I-0.5)dX for
C      I=1,2,...,M where dX =(B-A)/M.  M must be greater than 2.
C
C    MBDCND
C      Indicates the type of boundary conditions at X = A and X = B.
C
C      = 0  If the solution is periodic in X,
C           U(M+I,J) = U(I,J).
C
C      = 1  If the solution is specified at X = A and X = B.
C
C      = 2  If the solution is specified at X = A and the derivative
C           of the solution with respect to X is specified at X = B.
C
C      = 3  If the derivative of the solution with respect to X is
C           specified at X = A  and X = B.
C
C      = 4  If the derivative of the solution with respect to X is
C           specified at X = A  and the solution is specified at X = B.
C
C    BDA
C      A one-dimensional array of length N that specifies the boundary
C      values (if any) of the solution at X = A.  When MBDCND = 1 or 2,
C
C               BDA(J) = U(A,Y(J)) ,          J=1,2,...,N.
C
C      When MBDCND = 3 or 4,
C
C               BDA(J) = (d/dX)U(A,Y(J)) ,    J=1,2,...,N.
C
C    BDB
C      A one-dimensional array of length N that specifies the boundary
C      values of the solution at X = B.  When MBDCND = 1 or 4
C
C               BDB(J) = U(B,Y(J)) ,          J=1,2,...,N.
C
C      When MBDCND = 2 or 3
C
C               BDB(J) = (d/dX)U(B,Y(J)) ,    J=1,2,...,N.
C
C    C,D
C      The range of Y, i.e. C .LE. Y .LE. D.  C must be less
C      than D.
C
C    N
C      The number of unknowns in the interval (C,D).  The unknowns in
C      the Y-direction are given by Y(J) = C + (J-0.5)DY,
C      J=1,2,...,N, where DY = (D-C)/N.  N must be greater than 2.
C
C    NBDCND
C      Indicates the type of boundary conditions at Y = C
C      and Y = D.
C
C      = 0  If the solution is periodic in Y, i.e.
C           U(I,J) = U(I,N+J).
C
C      = 1  If the solution is specified at Y = C and Y = D.
C
C      = 2  If the solution is specified at Y = C and the derivative
C           of the solution with respect to Y is specified at Y = D.
C
C      = 3  If the derivative of the solution with respect to Y is
C           specified at Y = C and Y = D.
C
C      = 4  If the derivative of the solution with respect to Y is
C           specified at Y = C and the solution is specified at Y = D.
C
C    BDC
C      A one dimensional array of length M that specifies the boundary
C      values of the solution at Y = C.   When NBDCND = 1 or 2,
C
C               BDC(I) = U(X(I),C) ,              I=1,2,...,M.
C
C      When NBDCND = 3 or 4,
C
C               BDC(I) = (d/dY)U(X(I),C),     I=1,2,...,M.
C
C      When NBDCND = 0, BDC is a dummy variable.
C
C    BDD
C      A one-dimensional array of length M that specifies the boundary
C      values of the solution at Y = D.  When NBDCND = 1 or 4,
C
C               BDD(I) = U(X(I),D) ,              I=1,2,...,M.
C
C      When NBDCND = 2 or 3,
C
C               BDD(I) = (d/dY)U(X(I),D) ,    I=1,2,...,M.
C
C      When NBDCND = 0, BDD is a dummy variable.
C
C    ELMBDA
C      The constant LAMBDA in the Helmholtz equation.  If LAMBDA is
C      greater than 0, a solution may not exist.  However, HSTCRT will
C      attempt to find a solution.
C
C    F
C      A two-dimensional array that specifies the values of the right
C      side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N
C
C               F(I,J) = F(X(I),Y(J)) .
C
C      F must be dimensioned at least M X N.
C
C    IDIMF
C      The row (or first) dimension of the array F as it appears in the
C      program calling HSTCRT.  This parameter is used to specify the
C      variable dimension of F.  IDIMF must be at least M.
C
C    W
C      A one-dimensional array that must be provided by the user for
C      work space.  W may require up to 13M + 4N + M*INT(log2(N))
C      locations.  The actual number of locations used is computed by
C      HSTCRT and is returned in the location W(1).
C
C
C             * * * * * *   On Output   * * * * * *
C
C    F
C      Contains the solution U(I,J) of the finite difference
C      approximation for the grid point (X(I),Y(J)) for
C      I=1,2,...,M, J=1,2,...,N.
C
C    PERTRB
C      If a combination of periodic or derivative boundary conditions is
C      specified for a Poisson equation (LAMBDA = 0), a solution may not
C      exist.  PERTRB is a constant, calculated and subtracted from F,
C      which ensures that a solution exists.  HSTCRT then computes this
C      solution, which is a least squares solution to the original
C      approximation.  This solution plus any constant is also a
C      solution; hence, the solution is not unique.  The value of PERTRB
C      should be small compared to the right side F.  Otherwise, a
C      solution is obtained to an essentially different problem.  This
C      comparison should always be made to insure that a meaningful
C      solution has been obtained.
C
C    IERROR
C      An error flag that indicates invalid input parameters.
C       Except for numbers 0 and  6, a solution is not attempted.
C
C      =  0  No error
C
C      =  1  A .GE. B
C
C      =  2  MBDCND .LT. 0 or MBDCND .GT. 4
C
C      =  3  C .GE. D
C
C      =  4  N .LE. 2
C
C      =  5  NBDCND .LT. 0 or NBDCND .GT. 4
C
C      =  6  LAMBDA .GT. 0
C
C      =  7  IDIMF .LT. M
C
C      =  8  M .LE. 2
C
C      Since this is the only means of indicating a possibly
C      incorrect call to HSTCRT, the user should test IERROR after
C      the call.
C
C    W
C      W(1) contains the required length of W.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
C     Arguments      W(See argument list)
C
C     Latest         June 1, 1977
C     Revision
C
C     Subprograms    HSTCRT,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
C     Required       COSGEN,MERGE,TRIX,TRI3,PIMACH
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
C     History        Written by Roland Sweet at NCAR in January , 1977
C
C     Algorithm      This subroutine defines the finite-difference
C                    equations, incorporates boundary data, adjusts the
C                    right side when the system is singular and calls
C                    either POISTG or GENBUN which solves the linear
C                    system of equations.
C
C     Space          8131(decimal) = 17703(octal) locations on the
C     Required       NCAR Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine HSTCRT is roughly proportional
C                    to M*N*log2(N).  Some typical values are listed in
C                    the table below.
C                       The solution process employed results in a loss
C                    of no more than FOUR significant digits for N and M
C                    as large as 64.  More detailed information about
C                    accuracy can be found in the documentation for
C                    subroutine POISTG which is the routine that
C                    actually solves the finite difference equations.
C
C
C                       M(=N)    MBDCND    NBDCND    T(MSECS)
C                       -----    ------    ------    --------
C
C                        32       1-4       1-4         56
C                        64       1-4       1-4        230
C
C     Portability    American National Standards Institute Fortran.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C     Required       COS
C     Resident
C     Routines
C
C     Reference      Schumann, U. and R. Sweet,'A Direct Method For
C                    The Solution Of Poisson's Equation With Neumann
C                    Boundary Conditions On A Staggered Grid Of
C                    Arbitrary Size,' J. COMP. PHYS. 20(1976),
C                    PP. 171-182.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  U. Schumann and R. Sweet, A direct method for the
C                 solution of Poisson's equation with Neumann boundary
C                 conditions on a staggered grid of arbitrary size,
C                 Journal of Computational Physics 20, (1976),
C                 pp. 171-182.
C***ROUTINES CALLED  GENBUN, POISTG
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  HSTCRT
C
C
      DIMENSION       F(IDIMF,*) ,BDA(*)     ,BDB(*)     ,BDC(*)     ,
     1                BDD(*)     ,W(*)
C***FIRST EXECUTABLE STATEMENT  HSTCRT
      IERROR = 0
      IF (A .GE. B) IERROR = 1
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) IERROR = 2
      IF (C .GE. D) IERROR = 3
      IF (N .LE. 2) IERROR = 4
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR = 5
      IF (IDIMF .LT. M) IERROR = 7
      IF (M .LE. 2) IERROR = 8
      IF (IERROR .NE. 0) RETURN
      NPEROD = NBDCND
      MPEROD = 0
      IF (MBDCND .GT. 0) MPEROD = 1
      DELTAX = (B-A)/M
      TWDELX = 1./DELTAX
      DELXSQ = 2./DELTAX**2
      DELTAY = (D-C)/N
      TWDELY = 1./DELTAY
      DELYSQ = DELTAY**2
      TWDYSQ = 2./DELYSQ
      NP = NBDCND+1
      MP = MBDCND+1
C
C     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = M
      ID3 = ID2+M
      ID4 = ID3+M
      S = (DELTAY/DELTAX)**2
      ST2 = 2.*S
      DO 101 I=1,M
         W(I) = S
         J = ID2+I
         W(J) = -ST2+ELMBDA*DELYSQ
         J = ID3+I
         W(J) = S
  101 CONTINUE
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      GO TO (111,102,102,104,104),MP
  102 DO 103 J=1,N
         F(1,J) = F(1,J)-BDA(J)*DELXSQ
  103 CONTINUE
      W(ID2+1) = W(ID2+1)-W(1)
      GO TO 106
  104 DO 105 J=1,N
         F(1,J) = F(1,J)+BDA(J)*TWDELX
  105 CONTINUE
      W(ID2+1) = W(ID2+1)+W(1)
  106 GO TO (111,107,109,109,107),MP
  107 DO 108 J=1,N
         F(M,J) = F(M,J)-BDB(J)*DELXSQ
  108 CONTINUE
      W(ID3) = W(ID3)-W(1)
      GO TO 111
  109 DO 110 J=1,N
         F(M,J) = F(M,J)-BDB(J)*TWDELX
  110 CONTINUE
      W(ID3) = W(ID3)+W(1)
  111 CONTINUE
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (121,112,112,114,114),NP
  112 DO 113 I=1,M
         F(I,1) = F(I,1)-BDC(I)*TWDYSQ
  113 CONTINUE
      GO TO 116
  114 DO 115 I=1,M
         F(I,1) = F(I,1)+BDC(I)*TWDELY
  115 CONTINUE
  116 GO TO (121,117,119,119,117),NP
  117 DO 118 I=1,M
         F(I,N) = F(I,N)-BDD(I)*TWDYSQ
  118 CONTINUE
      GO TO 121
  119 DO 120 I=1,M
         F(I,N) = F(I,N)-BDD(I)*TWDELY
  120 CONTINUE
  121 CONTINUE
      DO 123 I=1,M
         DO 122 J=1,N
            F(I,J) = F(I,J)*DELYSQ
  122    CONTINUE
  123 CONTINUE
      IF (MPEROD .EQ. 0) GO TO 124
      W(1) = 0.
      W(ID4) = 0.
  124 CONTINUE
      PERTRB = 0.
      IF (ELMBDA) 133,126,125
  125 IERROR = 6
      GO TO 133
  126 GO TO (127,133,133,127,133),MP
  127 GO TO (128,133,133,128,133),NP
C
C     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
C     WILL EXIST.
C
  128 CONTINUE
      S = 0.
      DO 130 J=1,N
         DO 129 I=1,M
            S = S+F(I,J)
  129    CONTINUE
  130 CONTINUE
      PERTRB = S/(M*N)
      DO 132 J=1,N
         DO 131 I=1,M
            F(I,J) = F(I,J)-PERTRB
  131    CONTINUE
  132 CONTINUE
      PERTRB = PERTRB/DELYSQ
C
C     SOLVE THE EQUATION.
C
  133 CONTINUE
      IF (NPEROD .EQ. 0) GO TO 134
      CALL POISTG (NPEROD,N,MPEROD,M,W(1),W(ID2+1),W(ID3+1),IDIMF,F,
     1             IERR1,W(ID4+1))
      GO TO 135
  134 CONTINUE
      CALL GENBUN (NPEROD,N,MPEROD,M,W(1),W(ID2+1),W(ID3+1),IDIMF,F,
     1             IERR1,W(ID4+1))
  135 CONTINUE
      W(1) = W(ID4+1)+3*M
      RETURN
      END
