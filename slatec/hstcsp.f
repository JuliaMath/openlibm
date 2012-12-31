*DECK HSTCSP
      SUBROUTINE HSTCSP (INTL, A, B, M, MBDCND, BDA, BDB, C, D, N,
     +   NBDCND, BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HSTCSP
C***PURPOSE  Solve the standard five-point finite difference
C            approximation on a staggered grid to the modified Helmholtz
C            equation in spherical coordinates assuming axisymmetry
C            (no dependence on longitude).
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HSTCSP-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     HSTCSP solves the standard five-point finite difference
C     approximation on a staggered grid to the modified Helmholtz
C     equation spherical coordinates assuming axisymmetry (no dependence
C     on longitude).
C
C                  (1/R**2)(d/dR)(R**2(dU/dR)) +
C
C       1/(R**2*SIN(THETA))(d/dTHETA)(SIN(THETA)(dU/dTHETA)) +
C
C            (LAMBDA/(R*SIN(THETA))**2)U  =  F(THETA,R)
C
C     where THETA is colatitude and R is the radial coordinate.
C     This two-dimensional modified Helmholtz equation results from
C     the Fourier transform of the three-dimensional Poisson equation.
C
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C    * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C
C            * * * * * *   On Input    * * * * * *
C
C   INTL
C     = 0  On initial entry to HSTCSP or if any of the arguments
C          C, D, N, or NBDCND are changed from a previous call.
C
C     = 1  If C, D, N, and NBDCND are all unchanged from previous
C          call to HSTCSP.
C
C     NOTE:  A call with INTL = 0 takes approximately 1.5 times as much
C            time as a call with INTL = 1.  Once a call with INTL = 0
C            has been made then subsequent solutions corresponding to
C            different F, BDA, BDB, BDC, and BDD can be obtained
C            faster with INTL = 1 since initialization is not repeated.
C
C   A,B
C     The range of THETA (colatitude), i.e. A .LE. THETA .LE. B.  A
C     must be less than B and A must be non-negative.  A and B are in
C     radians.  A = 0 corresponds to the north pole and B = PI
C     corresponds to the south pole.
C
C                  * * *  IMPORTANT  * * *
C
C     If B is equal to PI, then B must be computed using the statement
C
C     B = PIMACH(DUM)
C
C     This insures that B in the user's program is equal to PI in this
C     program which permits several tests of the input parameters that
C     otherwise would not be possible.
C
C                  * * * * * * * * * * * *
C
C   M
C     The number of grid points in the interval (A,B).  The grid points
C     in the THETA-direction are given by THETA(I) = A + (I-0.5)DTHETA
C     for I=1,2,...,M where DTHETA =(B-A)/M.  M must be greater than 4.
C
C   MBDCND
C     Indicates the type of boundary conditions at THETA = A and
C     THETA = B.
C
C     = 1  If the solution is specified at THETA = A and THETA = B.
C          (See notes 1, 2 below)
C
C     = 2  If the solution is specified at THETA = A and the derivative
C          of the solution with respect to THETA is specified at
C          THETA = B (See notes 1, 2 below).
C
C     = 3  If the derivative of the solution with respect to THETA is
C          specified at THETA = A (See notes 1, 2 below) and THETA = B.
C
C     = 4  If the derivative of the solution with respect to THETA is
C          specified at THETA = A (See notes 1, 2 below) and the
C          solution is specified at THETA = B.
C
C     = 5  If the solution is unspecified at THETA = A = 0 and the
C          solution is specified at THETA = B. (See note 2 below)
C
C     = 6  If the solution is unspecified at THETA = A = 0 and the
C          derivative of the solution with respect to THETA is
C          specified at THETA = B (See note 2 below).
C
C     = 7  If the solution is specified at THETA = A and the
C          solution is unspecified at THETA = B = PI.
C
C     = 8  If the derivative of the solution with respect to
C          THETA is specified at THETA = A (See note 1 below)
C          and the solution is unspecified at THETA = B = PI.
C
C     = 9  If the solution is unspecified at THETA = A = 0 and
C          THETA = B = PI.
C
C     NOTES:  1.  If A = 0, do not use MBDCND = 1,2,3,4,7 or 8,
C                 but instead use MBDCND = 5, 6, or 9.
C
C             2.  if B = PI, do not use MBDCND = 1,2,3,4,5 or 6,
C                 but instead use MBDCND = 7, 8, or 9.
C
C             When A = 0  and/or B = PI the only meaningful boundary
C             condition is dU/dTHETA = 0.  (See D. Greenspan, 'Numerical
C             Analysis of Elliptic Boundary Value Problems,' Harper and
C             Row, 1965, Chapter 5.)
C
C   BDA
C     A one-dimensional array of length N that specifies the boundary
C     values (if any) of the solution at THETA = A.  When
C     MBDCND = 1, 2, or 7,
C
C              BDA(J) = U(A,R(J)) ,              J=1,2,...,N.
C
C     When MBDCND = 3, 4, or 8,
C
C              BDA(J) = (d/dTHETA)U(A,R(J)) ,    J=1,2,...,N.
C
C     When MBDCND has any other value, BDA is a dummy variable.
C
C   BDB
C     A one-dimensional array of length N that specifies the boundary
C     values of the solution at THETA = B.  When MBDCND = 1, 4, or 5,
C
C              BDB(J) = U(B,R(J)) ,              J=1,2,...,N.
C
C     When MBDCND = 2,3, or 6,
C
C              BDB(J) = (d/dTHETA)U(B,R(J)) ,    J=1,2,...,N.
C
C     When MBDCND has any other value, BDB is a dummy variable.
C
C   C,D
C     The range of R , i.e. C .LE. R .LE. D.
C     C must be less than D.  C must be non-negative.
C
C   N
C     The number of unknowns in the interval (C,D).  The unknowns in
C     the R-direction are given by R(J) = C + (J-0.5)DR,
C     J=1,2,...,N, where DR = (D-C)/N.  N must be greater than 4.
C
C   NBDCND
C     Indicates the type of boundary conditions at R = C
C     and R = D.
C
C     = 1  If the solution is specified at R = C and R = D.
C
C     = 2  If the solution is specified at R = C and the derivative
C          of the solution with respect to R is specified at
C          R = D. (See note 1 below)
C
C     = 3  If the derivative of the solution with respect to R is
C          specified at R = C and R = D.
C
C     = 4  If the derivative of the solution with respect to R is
C          specified at R = C and the solution is specified at
C          R = D.
C
C     = 5  If the solution is unspecified at R = C = 0 (See note 2
C          below) and the solution is specified at R = D.
C
C     = 6  If the solution is unspecified at R = C = 0 (See note 2
C          below) and the derivative of the solution with respect to R
C          is specified at R = D.
C
C     NOTE 1:  If C = 0 and MBDCND = 3,6,8 or 9, the system of equations
C              to be solved is singular.  The unique solution is
C              determined by extrapolation to the specification of
C              U(THETA(1),C).  But in these cases the right side of the
C              system will be perturbed by the constant PERTRB.
C
C     NOTE 2:  NBDCND = 5 or 6 cannot be used with MBDCND = 1, 2, 4, 5,
C              or 7 (the former indicates that the solution is
C              unspecified at R = 0; the latter indicates that the
C              solution is specified).  Use instead NBDCND = 1 or 2.
C
C   BDC
C     A one dimensional array of length M that specifies the boundary
C     values of the solution at R = C.   When NBDCND = 1 or 2,
C
C              BDC(I) = U(THETA(I),C) ,              I=1,2,...,M.
C
C     When NBDCND = 3 or 4,
C
C              BDC(I) = (d/dR)U(THETA(I),C),         I=1,2,...,M.
C
C     When NBDCND has any other value, BDC is a dummy variable.
C
C   BDD
C     A one-dimensional array of length M that specifies the boundary
C     values of the solution at R = D.  When NBDCND = 1 or 4,
C
C              BDD(I) = U(THETA(I),D) ,              I=1,2,...,M.
C
C     When NBDCND = 2 or 3,
C
C              BDD(I) = (d/dR)U(THETA(I),D) ,        I=1,2,...,M.
C
C     When NBDCND has any other value, BDD is a dummy variable.
C
C   ELMBDA
C     The constant LAMBDA in the modified Helmholtz equation.  If
C     LAMBDA is greater than 0, a solution may not exist.  However,
C     HSTCSP will attempt to find a solution.
C
C   F
C     A two-dimensional array that specifies the values of the right
C     side of the modified Helmholtz equation.  For I=1,2,...,M and
C     J=1,2,...,N
C
C              F(I,J) = F(THETA(I),R(J)) .
C
C     F must be dimensioned at least M X N.
C
C   IDIMF
C     The row (or first) dimension of the array F as it appears in the
C     program calling HSTCSP.  This parameter is used to specify the
C     variable dimension of F.  IDIMF must be at least M.
C
C   W
C     A one-dimensional array that must be provided by the user for
C     work space.  With K = INT(log2(N))+1 and L = 2**(K+1), W may
C     require up to (K-2)*L+K+MAX(2N,6M)+4(N+M)+5 locations.  The
C     actual number of locations used is computed by HSTCSP and is
C     returned in the location W(1).
C
C
C            * * * * * *   On Output   * * * * * *
C
C   F
C     Contains the solution U(I,J) of the finite difference
C     approximation for the grid point (THETA(I),R(J)) for
C     I=1,2,...,M, J=1,2,...,N.
C
C   PERTRB
C     If a combination of periodic, derivative, or unspecified
C     boundary conditions is specified for a Poisson equation
C     (LAMBDA = 0), a solution may not exist.  PERTRB is a con-
C     stant, calculated and subtracted from F, which ensures
C     that a solution exists.  HSTCSP then computes this
C     solution, which is a least squares solution to the
C     original approximation.  This solution plus any constant is also
C     a solution; hence, the solution is not unique.  The value of
C     PERTRB should be small compared to the right side F.
C     Otherwise, a solution is obtained to an essentially different
C     problem.  This comparison should always be made to insure that
C     a meaningful solution has been obtained.
C
C   IERROR
C     An error flag that indicates invalid input parameters.
C     Except for numbers 0 and 10, a solution is not attempted.
C
C     =  0  No error
C
C     =  1  A .LT. 0 or B .GT. PI
C
C     =  2  A .GE. B
C
C     =  3  MBDCND .LT. 1 or MBDCND .GT. 9
C
C     =  4  C .LT. 0
C
C     =  5  C .GE. D
C
C     =  6  NBDCND .LT. 1 or NBDCND .GT. 6
C
C     =  7  N .LT. 5
C
C     =  8  NBDCND = 5 or 6 and MBDCND = 1, 2, 4, 5, or 7
C
C     =  9  C .GT. 0 and NBDCND .GE. 5
C
C     = 10  ELMBDA .GT. 0
C
C     = 11  IDIMF .LT. M
C
C     = 12  M .LT. 5
C
C     = 13  A = 0 and MBDCND =1,2,3,4,7 or 8
C
C     = 14  B = PI and MBDCND .LE. 6
C
C     = 15  A .GT. 0 and MBDCND = 5, 6, or 9
C
C     = 16  B .LT. PI and MBDCND .GE. 7
C
C     = 17  LAMBDA .NE. 0 and NBDCND .GE. 5
C
C     Since this is the only means of indicating a possibly
C     incorrect call to HSTCSP, the user should test IERROR after
C     the call.
C
C   W
C     W(1) contains the required length of W.  Also  W contains
C     intermediate values that must not be destroyed if HSTCSP
C     will be called again with INTL = 1.
C
C *Long Description:
C
C    * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C    Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
C    Arguments      W(See argument list)
C
C    Latest         June 1979
C    Revision
C
C    Subprograms    HSTCSP,HSTCS1,BLKTRI,BLKTR1,INDXA,INDXB,INDXC,
C    Required       PROD,PRODP,CPROD,CPRODP,PPADD,PSGF,BSRH,PPSGF,
C                   PPSPF,COMPB,TEVLS,R1MACH
C
C    Special        NONE
C    Conditions
C
C    Common         CBLKT
C    Blocks
C
C    I/O            NONE
C
C    Precision      Single
C
C    Specialist     Roland Sweet
C
C    Language       FORTRAN
C
C    History        Written by Roland Sweet at NCAR in May, 1977
C
C    Algorithm      This subroutine defines the finite-difference
C                   equations, incorporates boundary data, adjusts the
C                   right side when the system is singular and calls
C                   BLKTRI which solves the linear system of equations.
C
C    Space          5269(decimal) = 12225(octal) locations on the
C    Required       NCAR Control Data 7600
C
C    Timing and        The execution time T on the NCAR Control Data
C    Accuracy       7600 for subroutine HSTCSP is roughly proportional
C                   to M*N*log2(N), but depends on the input parameter
C                   INTL.  Some values are listed in the table below.
C                      The solution process employed results in a loss
C                   of no more than FOUR significant digits for N and M
C                   as large as 64.  More detailed information about
C                   accuracy can be found in the documentation for
C                   subroutine BLKTRI which is the routine that
C                   actually solves the finite difference equations.
C
C
C                      M(=N)     INTL      MBDCND(=NBDCND)     T(MSECS)
C                      -----     ----      ---------------     --------
C
C                       32        0              1-6             132
C                       32        1              1-6              88
C                       64        0              1-6             546
C                       64        1              1-6             380
C
C    Portability    American National Standards Institute Fortran.
C                   The machine accuracy is set using function R1MACH.
C
C    Required       COS,SIN,ABS,SQRT
C    Resident
C    Routines
C
C    Reference      Swarztrauber, P.N., 'A Direct Method For The
C                   Discrete Solution Of Separable Elliptic Equations,'
C                   SIAM J. Numer. Anal. 11(1974), pp. 1136-1150.
C
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
C                 subprograms for the solution of elliptic equations,
C                 NCAR TN/IA-109, July 1975, 138 pp.
C               P. N. Swarztrauber, A direct method for the discrete
C                 solution of separable elliptic equations, SIAM Journal
C                 on Numerical Analysis 11, (1974), pp. 1136-1150.
C***ROUTINES CALLED  HSTCS1, PIMACH
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  HSTCSP
C
C
      DIMENSION       F(IDIMF,*) ,BDA(*)     ,BDB(*)     ,BDC(*)     ,
     1                BDD(*)     ,W(*)
C***FIRST EXECUTABLE STATEMENT  HSTCSP
      PI = PIMACH(DUM)
C
C     CHECK FOR INVALID INPUT PARAMETERS
C
      IERROR = 0
      IF (A.LT.0. .OR. B.GT.PI) IERROR = 1
      IF (A .GE. B) IERROR = 2
      IF (MBDCND.LT.1 .OR. MBDCND.GT.9) IERROR = 3
      IF (C .LT. 0.) IERROR = 4
      IF (C .GE. D) IERROR = 5
      IF (NBDCND.LT.1 .OR. NBDCND.GT.6) IERROR = 6
      IF (N .LT. 5) IERROR = 7
      IF ((NBDCND.EQ.5 .OR. NBDCND.EQ.6) .AND. (MBDCND.EQ.1 .OR.
     1    MBDCND.EQ.2 .OR. MBDCND.EQ.4 .OR. MBDCND.EQ.5 .OR.
     2                                                     MBDCND.EQ.7))
     3    IERROR = 8
      IF (C.GT.0. .AND. NBDCND.GE.5) IERROR = 9
      IF (IDIMF .LT. M) IERROR = 11
      IF (M .LT. 5) IERROR = 12
      IF (A.EQ.0. .AND. MBDCND.NE.5 .AND. MBDCND.NE.6 .AND. MBDCND.NE.9)
     1    IERROR = 13
      IF (B.EQ.PI .AND. MBDCND.LE.6) IERROR = 14
      IF (A.GT.0. .AND. (MBDCND.EQ.5 .OR. MBDCND.EQ.6 .OR. MBDCND.EQ.9))
     1    IERROR = 15
      IF (B.LT.PI .AND. MBDCND.GE.7) IERROR = 16
      IF (ELMBDA.NE.0. .AND. NBDCND.GE.5) IERROR = 17
      IF (IERROR .NE. 0) GO TO 101
      IWBM = M+1
      IWCM = IWBM+M
      IWAN = IWCM+M
      IWBN = IWAN+N
      IWCN = IWBN+N
      IWSNTH = IWCN+N
      IWRSQ = IWSNTH+M
      IWWRK = IWRSQ+N
      IERR1 = 0
      CALL HSTCS1 (INTL,A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     1             ELMBDA,F,IDIMF,PERTRB,IERR1,W,W(IWBM),W(IWCM),
     2             W(IWAN),W(IWBN),W(IWCN),W(IWSNTH),W(IWRSQ),W(IWWRK))
      W(1) = W(IWWRK)+IWWRK-1
      IERROR = IERR1
  101 CONTINUE
      RETURN
      END
