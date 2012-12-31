*DECK HSTCYL
      SUBROUTINE HSTCYL (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND,
     +   BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HSTCYL
C***PURPOSE  Solve the standard five-point finite difference
C            approximation on a staggered grid to the modified
C            Helmholtz equation in cylindrical coordinates.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HSTCYL-S)
C***KEYWORDS  CYLINDRICAL, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C      HSTCYL solves the standard five-point finite difference
C      approximation on a staggered grid to the modified Helmholtz
C      equation in cylindrical coordinates
C
C          (1/R)(d/dR)(R(dU/dR)) + (d/dZ)(dU/dZ)C
C                      + LAMBDA*(1/R**2)*U = F(R,Z)
C
C      This two-dimensional modified Helmholtz equation results
C      from the Fourier transform of a three-dimensional Poisson
C      equation.
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
C      = 1  If the solution is specified at R = A (see note below) and
C           R = B.
C
C      = 2  If the solution is specified at R = A (see note below) and
C           the derivative of the solution with respect to R is
C           specified at R = B.
C
C      = 3  If the derivative of the solution with respect to R is
C           specified at R = A (see note below) and R = B.
C
C      = 4  If the derivative of the solution with respect to R is
C           specified at R = A (see note below) and the solution is
C           specified at R = B.
C
C      = 5  If the solution is unspecified at R = A = 0 and the solution
C           is specified at R = B.
C
C      = 6  If the solution is unspecified at R = A = 0 and the
C           derivative of the solution with respect to R is specified at
C           R = B.
C
C      NOTE:  If A = 0, do not use MBDCND = 1,2,3, or 4, but instead
C             use MBDCND = 5 or 6.  The resulting approximation gives
C             the only meaningful boundary condition, i.e. dU/dR = 0.
C             (see D. Greenspan, 'Introductory Numerical Analysis Of
C             Elliptic Boundary Value Problems,' Harper and Row, 1965,
C             Chapter 5.)
C
C    BDA
C      A one-dimensional array of length N that specifies the boundary
C      values (if any) of the solution at R = A.  When MBDCND = 1 or 2,
C
C               BDA(J) = U(A,Z(J)) ,          J=1,2,...,N.
C
C      When MBDCND = 3 or 4,
C
C               BDA(J) = (d/dR)U(A,Z(J)) ,    J=1,2,...,N.
C
C      When MBDCND = 5 or 6, BDA is a dummy variable.
C
C    BDB
C      A one-dimensional array of length N that specifies the boundary
C      values of the solution at R = B.  When MBDCND = 1,4, or 5,
C
C               BDB(J) = U(B,Z(J)) ,          J=1,2,...,N.
C
C      When MBDCND = 2,3, or 6,
C
C               BDB(J) = (d/dR)U(B,Z(J)) ,    J=1,2,...,N.
C
C    C,D
C      The range of Z, i.e. C .LE. Z .LE. D.  C must be less
C      than D.
C
C    N
C      The number of unknowns in the interval (C,D).  The unknowns in
C      the Z-direction are given by Z(J) = C + (J-0.5)DZ,
C      J=1,2,...,N, where DZ = (D-C)/N.  N must be greater than 2.
C
C    NBDCND
C      Indicates the type of boundary conditions at Z = C
C      and Z = D.
C
C      = 0  If the solution is periodic in Z, i.e.
C           U(I,J) = U(I,N+J).
C
C      = 1  If the solution is specified at Z = C and Z = D.
C
C      = 2  If the solution is specified at Z = C and the derivative
C           of the solution with respect to Z is specified at
C           Z = D.
C
C      = 3  If the derivative of the solution with respect to Z is
C           specified at Z = C and Z = D.
C
C      = 4  If the derivative of the solution with respect to Z is
C           specified at Z = C and the solution is specified at
C           Z = D.
C
C    BDC
C      A one dimensional array of length M that specifies the boundary
C      values of the solution at Z = C.   When NBDCND = 1 or 2,
C
C               BDC(I) = U(R(I),C) ,              I=1,2,...,M.
C
C      When NBDCND = 3 or 4,
C
C               BDC(I) = (d/dZ)U(R(I),C),         I=1,2,...,M.
C
C      When NBDCND = 0, BDC is a dummy variable.
C
C    BDD
C      A one-dimensional array of length M that specifies the boundary
C      values of the solution at Z = D.  when NBDCND = 1 or 4,
C
C               BDD(I) = U(R(I),D) ,              I=1,2,...,M.
C
C      When NBDCND = 2 or 3,
C
C               BDD(I) = (d/dZ)U(R(I),D) ,        I=1,2,...,M.
C
C      When NBDCND = 0, BDD is a dummy variable.
C
C    ELMBDA
C      The constant LAMBDA in the modified Helmholtz equation.  If
C      LAMBDA is greater than 0, a solution may not exist.  However,
C      HSTCYL will attempt to find a solution.  LAMBDA must be zero
C      when MBDCND = 5 or 6.
C
C    F
C      A two-dimensional array that specifies the values of the right
C      side of the modified Helmholtz equation.  For I=1,2,...,M
C      and J=1,2,...,N
C
C               F(I,J) = F(R(I),Z(J)) .
C
C      F must be dimensioned at least M X N.
C
C    IDIMF
C      The row (or first) dimension of the array F as it appears in the
C      program calling HSTCYL.  This parameter is used to specify the
C      variable dimension of F.  IDIMF must be at least M.
C
C    W
C      A one-dimensional array that must be provided by the user for
C      work space.  W may require up to 13M + 4N + M*INT(log2(N))
C      locations.  The actual number of locations used is computed by
C      HSTCYL and is returned in the location W(1).
C
C
C             * * * * * *   On Output   * * * * * *
C
C    F
C      Contains the solution U(I,J) of the finite difference
C      approximation for the grid point (R(I),Z(J)) for
C      I=1,2,...,M, J=1,2,...,N.
C
C    PERTRB
C      If a combination of periodic, derivative, or unspecified
C      boundary conditions is specified for a Poisson equation
C      (LAMBDA = 0), a solution may not exist.  PERTRB is a con-
C      stant, calculated and subtracted from F, which ensures
C      that a solution exists.  HSTCYL then computes this
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
C      =  7  A = 0 and MBDCND = 1,2,3, or 4
C
C      =  8  A .GT. 0 and MBDCND .GE. 5
C
C      =  9  M .LE. 2
C
C      = 10  IDIMF .LT. M
C
C      = 11  LAMBDA .GT. 0
C
C      = 12  A=0, MBDCND .GE. 5, ELMBDA .NE. 0
C
C      Since this is the only means of indicating a possibly
C      incorrect call to HSTCYL, the user should test IERROR after
C      the call.
C
C    W
C      W(1) contains the required length of W.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension OF   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
C     Arguments      W(see argument list)
C
C     Latest         June 1, 1977
C     Revision
C
C     Subprograms    HSTCYL,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
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
C     History        Written by Roland Sweet at NCAR in March, 1977
C
C     Algorithm      This subroutine defines the finite-difference
C                    equations, incorporates boundary data, adjusts the
C                    right side when the system is singular and calls
C                    either POISTG or GENBUN which solves the linear
C                    system of equations.
C
C     Space          8228(decimal) = 20044(octal) locations on the
C     Required       NCAR Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine HSTCYL is roughly proportional
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
C                    The Solution of Poisson's Equation With Neumann
C                    Boundary Conditions On A Staggered Grid Of
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
C***END PROLOGUE  HSTCYL
C
C
      DIMENSION       F(IDIMF,*) ,BDA(*)     ,BDB(*)     ,BDC(*)     ,
     1                BDD(*)     ,W(*)
C***FIRST EXECUTABLE STATEMENT  HSTCYL
      IERROR = 0
      IF (A .LT. 0.) IERROR = 1
      IF (A .GE. B) IERROR = 2
      IF (MBDCND.LE.0 .OR. MBDCND.GE.7) IERROR = 3
      IF (C .GE. D) IERROR = 4
      IF (N .LE. 2) IERROR = 5
      IF (NBDCND.LT.0 .OR. NBDCND.GE.5) IERROR = 6
      IF (A.EQ.0. .AND. MBDCND.NE.5 .AND. MBDCND.NE.6) IERROR = 7
      IF (A.GT.0. .AND. MBDCND.GE.5) IERROR = 8
      IF (IDIMF .LT. M) IERROR = 10
      IF (M .LE. 2) IERROR = 9
      IF (A.EQ.0. .AND. MBDCND.GE.5 .AND. ELMBDA.NE.0.) IERROR = 12
      IF (IERROR .NE. 0) RETURN
      DELTAR = (B-A)/M
      DLRSQ = DELTAR**2
      DELTHT = (D-C)/N
      DLTHSQ = DELTHT**2
      NP = NBDCND+1
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      IWB = M
      IWC = IWB+M
      IWR = IWC+M
      DO 101 I=1,M
         J = IWR+I
         W(J) = A+(I-0.5)*DELTAR
         W(I) = (A+(I-1)*DELTAR)/(DLRSQ*W(J))
         K = IWC+I
         W(K) = (A+I*DELTAR)/(DLRSQ*W(J))
         K = IWB+I
         W(K) = ELMBDA/W(J)**2-2./DLRSQ
  101 CONTINUE
C
C     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
C
      GO TO (102,102,104,104,106,106),MBDCND
  102 A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1)-W(1)
      DO 103 J=1,N
         F(1,J) = F(1,J)-A1*BDA(J)
  103 CONTINUE
      GO TO 106
  104 A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1)+W(1)
      DO 105 J=1,N
         F(1,J) = F(1,J)+A1*BDA(J)
  105 CONTINUE
  106 CONTINUE
      GO TO (107,109,109,107,107,109),MBDCND
  107 W(IWC) = W(IWC)-W(IWR)
      A1 = 2.*W(IWR)
      DO 108 J=1,N
         F(M,J) = F(M,J)-A1*BDB(J)
  108 CONTINUE
      GO TO 111
  109 W(IWC) = W(IWC)+W(IWR)
      A1 = DELTAR*W(IWR)
      DO 110 J=1,N
         F(M,J) = F(M,J)-A1*BDB(J)
  110 CONTINUE
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  111 A1 = 2./DLTHSQ
      GO TO (121,112,112,114,114),NP
  112 DO 113 I=1,M
         F(I,1) = F(I,1)-A1*BDC(I)
  113 CONTINUE
      GO TO 116
  114 A1 = 1./DELTHT
      DO 115 I=1,M
         F(I,1) = F(I,1)+A1*BDC(I)
  115 CONTINUE
  116 A1 = 2./DLTHSQ
      GO TO (121,117,119,119,117),NP
  117 DO 118 I=1,M
         F(I,N) = F(I,N)-A1*BDD(I)
  118 CONTINUE
      GO TO 121
  119 A1 = 1./DELTHT
      DO 120 I=1,M
         F(I,N) = F(I,N)-A1*BDD(I)
  120 CONTINUE
  121 CONTINUE
C
C     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
C     SOLUTION.
C
      PERTRB = 0.
      IF (ELMBDA) 130,123,122
  122 IERROR = 11
      GO TO 130
  123 GO TO (130,130,124,130,130,124),MBDCND
  124 GO TO (125,130,130,125,130),NP
  125 CONTINUE
      DO 127 I=1,M
         A1 = 0.
         DO 126 J=1,N
            A1 = A1+F(I,J)
  126    CONTINUE
         J = IWR+I
         PERTRB = PERTRB+A1*W(J)
  127 CONTINUE
      PERTRB = PERTRB/(M*N*0.5*(A+B))
      DO 129 I=1,M
         DO 128 J=1,N
            F(I,J) = F(I,J)-PERTRB
  128    CONTINUE
  129 CONTINUE
  130 CONTINUE
C
C     MULTIPLY I-TH EQUATION THROUGH BY  DELTHT**2
C
      DO 132 I=1,M
         W(I) = W(I)*DLTHSQ
         J = IWC+I
         W(J) = W(J)*DLTHSQ
         J = IWB+I
         W(J) = W(J)*DLTHSQ
         DO 131 J=1,N
            F(I,J) = F(I,J)*DLTHSQ
  131    CONTINUE
  132 CONTINUE
      LP = NBDCND
      W(1) = 0.
      W(IWR) = 0.
C
C     CALL GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
C
      IF (NBDCND .EQ. 0) GO TO 133
      CALL POISTG (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
      GO TO 134
  133 CALL GENBUN (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  134 CONTINUE
      W(1) = W(IWR+1)+3*M
      RETURN
      END
