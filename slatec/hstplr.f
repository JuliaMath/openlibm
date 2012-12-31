*DECK HSTPLR
      SUBROUTINE HSTPLR (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND,
     +   BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HSTPLR
C***PURPOSE  Solve the standard five-point finite difference
C            approximation on a staggered grid to the Helmholtz equation
C            in polar coordinates.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HSTPLR-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, POLAR
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C      HSTPLR solves the standard five-point finite difference
C      approximation on a staggered grid to the Helmholtz equation in
C      polar coordinates
C
C      (1/R)(d/DR)(R(dU/DR)) + (1/R**2)(d/dTHETA)(dU/dTHETA)
C
C                      + LAMBDA*U = F(R,THETA)
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C             * * * * * *   On Input    * * * * * *
C
C    A,B
C      The range of R, i.e. A .LE. R .LE. B.  A must be less than B and
C      A must be non-negative.
C
C    M
C      The number of grid points in the interval (A,B).  The grid points
C      in the R-direction are given by R(I) = A + (I-0.5)DR for
C      I=1,2,...,M where DR =(B-A)/M.  M must be greater than 2.
C
C    MBDCND
C      Indicates the type of boundary conditions at R = A and R = B.
C
C      = 1  If the solution is specified at R = A and R = B.
C
C      = 2  If the solution is specified at R = A and the derivative
C           of the solution with respect to R is specified at R = B.
C           (see note 1 below)
C
C      = 3  If the derivative of the solution with respect to R is
C           specified at R = A (see note 2 below) and R = B.
C
C      = 4  If the derivative of the solution with respect to R is
C           specified at R = A (see note 2 below) and the solution is
C           specified at R = B.
C
C      = 5  If the solution is unspecified at R = A = 0 and the solution
C           is specified at R = B.
C
C      = 6  If the solution is unspecified at R = A = 0 and the
C           derivative of the solution with respect to R is specified at
C           R = B.
C
C      NOTE 1:  If A = 0, MBDCND = 2, and NBDCND = 0 or 3, the system of
C               equations to be solved is singular.  The unique solution
C               is determined by extrapolation to the specification of
C               U(0,THETA(1)).  But in this case the right side of the
C               system will be perturbed by the constant PERTRB.
C
C      NOTE 2:  If A = 0, do not use MBDCND = 3 or 4, but instead use
C               MBDCND = 1,2,5, or 6.
C
C    BDA
C      A one-dimensional array of length N that specifies the boundary
C      values (if any) of the solution at R = A.  When MBDCND = 1 or 2,
C
C               BDA(J) = U(A,THETA(J)) ,          J=1,2,...,N.
C
C      When MBDCND = 3 or 4,
C
C               BDA(J) = (d/dR)U(A,THETA(J)) ,    J=1,2,...,N.
C
C      When MBDCND = 5 or 6, BDA is a dummy variable.
C
C    BDB
C      A one-dimensional array of length N that specifies the boundary
C      values of the solution at R = B.  When MBDCND = 1,4, or 5,
C
C               BDB(J) = U(B,THETA(J)) ,          J=1,2,...,N.
C
C      When MBDCND = 2,3, or 6,
C
C               BDB(J) = (d/dR)U(B,THETA(J)) ,    J=1,2,...,N.
C
C    C,D
C      The range of THETA, i.e. C .LE. THETA .LE. D.  C must be less
C      than D.
C
C    N
C      The number of unknowns in the interval (C,D).  The unknowns in
C      the THETA-direction are given by THETA(J) = C + (J-0.5)DT,
C      J=1,2,...,N, where DT = (D-C)/N.  N must be greater than 2.
C
C    NBDCND
C      Indicates the type of boundary conditions at THETA = C
C      and THETA = D.
C
C      = 0  If the solution is periodic in THETA, i.e.
C           U(I,J) = U(I,N+J).
C
C      = 1  If the solution is specified at THETA = C and THETA = D
C           (see note below).
C
C      = 2  If the solution is specified at THETA = C and the derivative
C           of the solution with respect to THETA is specified at
C           THETA = D (see note below).
C
C      = 3  If the derivative of the solution with respect to THETA is
C           specified at THETA = C and THETA = D.
C
C      = 4  If the derivative of the solution with respect to THETA is
C           specified at THETA = C and the solution is specified at
C           THETA = d (see note below).
C
C      NOTE:  When NBDCND = 1, 2, or 4, do not use MBDCND = 5 or 6 (the
C      former indicates that the solution is specified at R =  0; the
C      latter indicates the solution is unspecified at R = 0).  Use
C      instead MBDCND = 1 or 2.
C
C    BDC
C      A one dimensional array of length M that specifies the boundary
C      values of the solution at THETA = C.   When NBDCND = 1 or 2,
C
C               BDC(I) = U(R(I),C) ,              I=1,2,...,M.
C
C      When NBDCND = 3 or 4,
C
C               BDC(I) = (d/dTHETA)U(R(I),C),     I=1,2,...,M.
C
C      When NBDCND = 0, BDC is a dummy variable.
C
C    BDD
C      A one-dimensional array of length M that specifies the boundary
C      values of the solution at THETA = D.  When NBDCND = 1 or 4,
C
C               BDD(I) = U(R(I),D) ,              I=1,2,...,M.
C
C      When NBDCND = 2 or 3,
C
C               BDD(I) = (d/dTHETA)U(R(I),D) ,    I=1,2,...,M.
C
C      When NBDCND = 0, BDD is a dummy variable.
C
C    ELMBDA
C      The constant LAMBDA in the Helmholtz equation.  If LAMBDA is
C      greater than 0, a solution may not exist.  However, HSTPLR will
C      attempt to find a solution.
C
C    F
C      A two-dimensional array that specifies the values of the right
C      side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N
C
C               F(I,J) = F(R(I),THETA(J)) .
C
C      F must be dimensioned at least M X N.
C
C    IDIMF
C      The row (or first) dimension of the array F as it appears in the
C      program calling HSTPLR.  This parameter is used to specify the
C      variable dimension of F.  IDIMF must be at least M.
C
C    W
C      A one-dimensional array that must be provided by the user for
C      work space.  W may require up to 13M + 4N + M*INT(log2(N))
C      locations.  The actual number of locations used is computed by
C      HSTPLR and is returned in the location W(1).
C
C
C             * * * * * *   On Output   * * * * * *
C
C    F
C      Contains the solution U(I,J) of the finite difference
C      approximation for the grid point (R(I),THETA(J)) for
C      I=1,2,...,M, J=1,2,...,N.
C
C    PERTRB
C      If a combination of periodic, derivative, or unspecified
C      boundary conditions is specified for a Poisson equation
C      (LAMBDA = 0), a solution may not exist.  PERTRB is a con-
C      stant, calculated and subtracted from F, which ensures
C      that a solution exists.  HSTPLR then computes this
C      solution, which is a least squares solution to the
C      original approximation.  This solution plus any constant is also
C      a solution; hence, the solution is not unique.  The value of
C      PERTRB should be small compared to the right side F.
C      Otherwise, a solution is obtained to an essentially different
C      problem.  This comparison should always be made to insure that
C      a meaningful solution has been obtained.
C
C    IERROR
C      An error flag that indicates invalid input parameters.
C      Except for numbers 0 and 11, a solution is not attempted.
C
C      =  0  No error
C
C      =  1  A .LT. 0
C
C      =  2  A .GE. B
C
C      =  3  MBDCND .LT. 1 or MBDCND .GT. 6
C
C      =  4  C .GE. D
C
C      =  5  N .LE. 2
C
C      =  6  NBDCND .LT. 0 or NBDCND .GT. 4
C
C      =  7  A = 0 and MBDCND = 3 or 4
C
C      =  8  A .GT. 0 and MBDCND .GE. 5
C
C      =  9  MBDCND .GE. 5 and NBDCND .NE. 0 or 3
C
C      = 10  IDIMF .LT. M
C
C      = 11  LAMBDA .GT. 0
C
C      = 12  M .LE. 2
C
C      Since this is the only means of indicating a possibly
C      incorrect call to HSTPLR, the user should test IERROR after
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
C     Arguments      W(see ARGUMENT LIST)
C
C     Latest         June 1, 1977
C     Revision
C
C     Subprograms    HSTPLR,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
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
C     History        Written by Roland Sweet at NCAR in February, 1977
C
C     Algorithm      This subroutine defines the finite-difference
C                    equations, incorporates boundary data, adjusts the
C                    right side when the system is singular and calls
C                    either POISTG or GENBUN which solves the linear
C                    system of equations.
C
C     Space          8265(decimal) = 20111(octal) LOCATIONS ON THE
C     Required       NCAR Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine HSTPLR is roughly proportional
C                    to M*N*log2(N).  Some typical values are listed in
C                    the table below.
C                       The solution process employed results in a loss
C                    of no more than four significant digits for N and M
C                    as large as 64.  More detailed information about
C                    accuracy can be found in the documentation for
C                    subroutine POISTG which is the routine that
C                    actually solves the finite difference equations.
C
C
C                       M(=N)    MBDCND    NBDCND    T(MSECS)
C                       -----    ------    ------    --------
C
C                        32       1-6       1-4         56
C                        64       1-6       1-4        230
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
C                    Boundary Conditions On A Staggered Grid of
C                    Arbitrary Size,' J. Comp. Phys. 20(1976),
C                    pp. 171-182.
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
C***END PROLOGUE  HSTPLR
C
C
      DIMENSION       F(IDIMF,*)
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
C***FIRST EXECUTABLE STATEMENT  HSTPLR
      IERROR = 0
      IF (A .LT. 0.) IERROR = 1
      IF (A .GE. B) IERROR = 2
      IF (MBDCND.LE.0 .OR. MBDCND.GE.7) IERROR = 3
      IF (C .GE. D) IERROR = 4
      IF (N .LE. 2) IERROR = 5
      IF (NBDCND.LT.0 .OR. NBDCND.GE.5) IERROR = 6
      IF (A.EQ.0. .AND. (MBDCND.EQ.3 .OR. MBDCND.EQ.4)) IERROR = 7
      IF (A.GT.0. .AND. MBDCND.GE.5) IERROR = 8
      IF (MBDCND.GE.5 .AND. NBDCND.NE.0 .AND. NBDCND.NE.3) IERROR = 9
      IF (IDIMF .LT. M) IERROR = 10
      IF (M .LE. 2) IERROR = 12
      IF (IERROR .NE. 0) RETURN
      DELTAR = (B-A)/M
      DLRSQ = DELTAR**2
      DELTHT = (D-C)/N
      DLTHSQ = DELTHT**2
      NP = NBDCND+1
      ISW = 1
      MB = MBDCND
      IF (A.EQ.0. .AND. MBDCND.EQ.2) MB = 6
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      IWB = M
      IWC = IWB+M
      IWR = IWC+M
      DO 101 I=1,M
         J = IWR+I
         W(J) = A+(I-0.5)*DELTAR
         W(I) = (A+(I-1)*DELTAR)/DLRSQ
         K = IWC+I
         W(K) = (A+I*DELTAR)/DLRSQ
         K = IWB+I
         W(K) = (ELMBDA-2./DLRSQ)*W(J)
  101 CONTINUE
      DO 103 I=1,M
         J = IWR+I
         A1 = W(J)
         DO 102 J=1,N
            F(I,J) = A1*F(I,J)
  102    CONTINUE
  103 CONTINUE
C
C     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
C
      GO TO (104,104,106,106,108,108),MB
  104 A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1)-W(1)
      DO 105 J=1,N
         F(1,J) = F(1,J)-A1*BDA(J)
  105 CONTINUE
      GO TO 108
  106 A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1)+W(1)
      DO 107 J=1,N
         F(1,J) = F(1,J)+A1*BDA(J)
  107 CONTINUE
  108 GO TO (109,111,111,109,109,111),MB
  109 A1 = 2.*W(IWR)
      W(IWC) = W(IWC)-W(IWR)
      DO 110 J=1,N
         F(M,J) = F(M,J)-A1*BDB(J)
  110 CONTINUE
      GO TO 113
  111 A1 = DELTAR*W(IWR)
      W(IWC) = W(IWC)+W(IWR)
      DO 112 J=1,N
         F(M,J) = F(M,J)-A1*BDB(J)
  112 CONTINUE
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  113 A1 = 2./DLTHSQ
      GO TO (123,114,114,116,116),NP
  114 DO 115 I=1,M
         J = IWR+I
         F(I,1) = F(I,1)-A1*BDC(I)/W(J)
  115 CONTINUE
      GO TO 118
  116 A1 = 1./DELTHT
      DO 117 I=1,M
         J = IWR+I
         F(I,1) = F(I,1)+A1*BDC(I)/W(J)
  117 CONTINUE
  118 A1 = 2./DLTHSQ
      GO TO (123,119,121,121,119),NP
  119 DO 120 I=1,M
         J = IWR+I
         F(I,N) = F(I,N)-A1*BDD(I)/W(J)
  120 CONTINUE
      GO TO 123
  121 A1 = 1./DELTHT
      DO 122 I=1,M
         J = IWR+I
         F(I,N) = F(I,N)-A1*BDD(I)/W(J)
  122 CONTINUE
  123 CONTINUE
C
C     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
C     SOLUTION.
C
      PERTRB = 0.
      IF (ELMBDA) 133,125,124
  124 IERROR = 11
      GO TO 133
  125 GO TO (133,133,126,133,133,126),MB
  126 GO TO (127,133,133,127,133),NP
  127 CONTINUE
      ISW = 2
      DO 129 J=1,N
         DO 128 I=1,M
            PERTRB = PERTRB+F(I,J)
  128    CONTINUE
  129 CONTINUE
      PERTRB = PERTRB/(M*N*0.5*(A+B))
      DO 131 I=1,M
         J = IWR+I
         A1 = PERTRB*W(J)
         DO 130 J=1,N
            F(I,J) = F(I,J)-A1
  130    CONTINUE
  131 CONTINUE
      A2 = 0.
      DO 132 J=1,N
         A2 = A2+F(1,J)
  132 CONTINUE
      A2 = A2/W(IWR+1)
  133 CONTINUE
C
C     MULTIPLY I-TH EQUATION THROUGH BY  R(I)*DELTHT**2
C
      DO 135 I=1,M
         J = IWR+I
         A1 = DLTHSQ*W(J)
         W(I) = A1*W(I)
         J = IWC+I
         W(J) = A1*W(J)
         J = IWB+I
         W(J) = A1*W(J)
         DO 134 J=1,N
            F(I,J) = A1*F(I,J)
  134    CONTINUE
  135 CONTINUE
      LP = NBDCND
      W(1) = 0.
      W(IWR) = 0.
C
C     CALL POISTG OR GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
C
      IF (LP .EQ. 0) GO TO 136
      CALL POISTG (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
      GO TO 137
  136 CALL GENBUN (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  137 CONTINUE
      W(1) = W(IWR+1)+3*M
      IF (A.NE.0. .OR. MBDCND.NE.2 .OR. ISW.NE.2) GO TO 141
      A1 = 0.
      DO 138 J=1,N
         A1 = A1+F(1,J)
  138 CONTINUE
      A1 = (A1-DLRSQ*A2/16.)/N
      IF (NBDCND .EQ. 3) A1 = A1+(BDD(1)-BDC(1))/(D-C)
      A1 = BDA(1)-A1
      DO 140 I=1,M
         DO 139 J=1,N
            F(I,J) = F(I,J)+A1
  139    CONTINUE
  140 CONTINUE
  141 CONTINUE
      RETURN
      END
