*DECK HSTSSP
      SUBROUTINE HSTSSP (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND,
     +   BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HSTSSP
C***PURPOSE  Solve the standard five-point finite difference
C            approximation on a staggered grid to the Helmholtz
C            equation in spherical coordinates and on the surface of
C            the unit sphere (radius of 1).
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HSTSSP-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     HSTSSP solves the standard five-point finite difference
C     approximation on a staggered grid to the Helmholtz equation in
C     spherical coordinates and on the surface of the unit sphere
C     (radius of 1)
C
C             (1/SIN(THETA))(d/dTHETA)(SIN(THETA)(dU/dTHETA)) +
C
C       (1/SIN(THETA)**2)(d/dPHI)(dU/dPHI) + LAMBDA*U = F(THETA,PHI)
C
C     where THETA is colatitude and PHI is longitude.
C
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C    * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C            * * * * * *   On Input    * * * * * *
C
C   A,B
C     The range of THETA (colatitude), i.e. A .LE. THETA .LE. B.  A
C     must be less than B and A must be non-negative.  A and B are in
C     radians.  A = 0 corresponds to the north pole and B = PI
C     corresponds to the south pole.
C
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
C
C
C   M
C     The number of grid points in the interval (A,B).  The grid points
C     in the THETA-direction are given by THETA(I) = A + (I-0.5)DTHETA
C     for I=1,2,...,M where DTHETA =(B-A)/M.  M must be greater than 2.
C
C   MBDCND
C     Indicates the type of boundary conditions at THETA = A and
C     THETA = B.
C
C     = 1  If the solution is specified at THETA = A and THETA = B.
C          (see note 3 below)
C
C     = 2  If the solution is specified at THETA = A and the derivative
C          of the solution with respect to THETA is specified at
C          THETA = B (see notes 2 and 3 below).
C
C     = 3  If the derivative of the solution with respect to THETA is
C          specified at THETA = A (see notes 1, 2 below) and THETA = B.
C
C     = 4  If the derivative of the solution with respect to THETA is
C          specified at THETA = A (see notes 1 and 2 below) and the
C          solution is specified at THETA = B.
C
C     = 5  If the solution is unspecified at THETA = A = 0 and the
C          solution is specified at THETA = B.  (see note 3 below)
C
C     = 6  If the solution is unspecified at THETA = A = 0 and the
C          derivative of the solution with respect to THETA is
C          specified at THETA = B (see note 2 below).
C
C     = 7  If the solution is specified at THETA = A and the
C          solution is unspecified at THETA = B = PI. (see note 3 below)
C
C     = 8  If the derivative of the solution with respect to
C          THETA is specified at THETA = A (see note 1 below)
C          and the solution is unspecified at THETA = B = PI.
C
C     = 9  If the solution is unspecified at THETA = A = 0 and
C          THETA = B = PI.
C
C     NOTES:  1.  If A = 0, do not use MBDCND = 3, 4, or 8,
C                 but instead use MBDCND = 5, 6, or 9.
C
C             2.  If B = PI, do not use MBDCND = 2, 3, or 6,
C                 but instead use MBDCND = 7, 8, or 9.
C
C             3.  When the solution is specified at THETA = 0 and/or
C                 THETA = PI and the other boundary conditions are
C                 combinations of unspecified, normal derivative, or
C                 periodicity a singular system results.  The unique
C                 solution is determined by extrapolation to the
C                 specification of the solution at either THETA = 0 or
C                 THETA = PI.  But in these cases the right side of the
C                 system will be perturbed by the constant PERTRB.
C
C   BDA
C     A one-dimensional array of length N that specifies the boundary
C     values (if any) of the solution at THETA = A.  When
C     MBDCND = 1, 2, or 7,
C
C              BDA(J) = U(A,PHI(J)) ,              J=1,2,...,N.
C
C     When MBDCND = 3, 4, or 8,
C
C              BDA(J) = (d/dTHETA)U(A,PHI(J)) ,    J=1,2,...,N.
C
C     When MBDCND has any other value, BDA is a dummy variable.
C
C   BDB
C     A one-dimensional array of length N that specifies the boundary
C     values of the solution at THETA = B.  When MBDCND = 1,4, or 5,
C
C              BDB(J) = U(B,PHI(J)) ,              J=1,2,...,N.
C
C     When MBDCND = 2,3, or 6,
C
C              BDB(J) = (d/dTHETA)U(B,PHI(J)) ,    J=1,2,...,N.
C
C     When MBDCND has any other value, BDB is a dummy variable.
C
C   C,D
C     The range of PHI (longitude), i.e. C .LE. PHI .LE. D.
C     C must be less than D.  If D-C = 2*PI, periodic boundary
C     conditions are usually prescribed.
C
C   N
C     The number of unknowns in the interval (C,D).  The unknowns in
C     the PHI-direction are given by PHI(J) = C + (J-0.5)DPHI,
C     J=1,2,...,N, where DPHI = (D-C)/N.  N must be greater than 2.
C
C   NBDCND
C     Indicates the type of boundary conditions at PHI = C
C     and PHI = D.
C
C     = 0  If the solution is periodic in PHI, i.e.
C          U(I,J) = U(I,N+J).
C
C     = 1  If the solution is specified at PHI = C and PHI = D
C          (see note below).
C
C     = 2  If the solution is specified at PHI = C and the derivative
C          of the solution with respect to PHI is specified at
C          PHI = D (see note below).
C
C     = 3  If the derivative of the solution with respect to PHI is
C          specified at PHI = C and PHI = D.
C
C     = 4  If the derivative of the solution with respect to PHI is
C          specified at PHI = C and the solution is specified at
C          PHI = D (see note below).
C
C     NOTE:  When NBDCND = 1, 2, or 4, do not use MBDCND = 5, 6, 7, 8,
C     or 9 (the former indicates that the solution is specified at
C     a pole; the latter indicates the solution is unspecified).  Use
C     instead MBDCND = 1 or 2.
C
C   BDC
C     A one dimensional array of length M that specifies the boundary
C     values of the solution at PHI = C.   When NBDCND = 1 or 2,
C
C              BDC(I) = U(THETA(I),C) ,              I=1,2,...,M.
C
C     When NBDCND = 3 or 4,
C
C              BDC(I) = (d/dPHI)U(THETA(I),C),       I=1,2,...,M.
C
C     When NBDCND = 0, BDC is a dummy variable.
C
C   BDD
C     A one-dimensional array of length M that specifies the boundary
C     values of the solution at PHI = D.  When NBDCND = 1 or 4,
C
C              BDD(I) = U(THETA(I),D) ,              I=1,2,...,M.
C
C     When NBDCND = 2 or 3,
C
C              BDD(I) = (d/dPHI)U(THETA(I),D) ,      I=1,2,...,M.
C
C     When NBDCND = 0, BDD is a dummy variable.
C
C   ELMBDA
C     The constant LAMBDA in the Helmholtz equation.  If LAMBDA is
C     greater than 0, a solution may not exist.  However, HSTSSP will
C     attempt to find a solution.
C
C   F
C     A two-dimensional array that specifies the values of the right
C     side of the Helmholtz equation.  For I=1,2,...,M and J=1,2,...,N
C
C              F(I,J) = F(THETA(I),PHI(J)) .
C
C     F must be dimensioned at least M X N.
C
C   IDIMF
C     The row (or first) dimension of the array F as it appears in the
C     program calling HSTSSP.  This parameter is used to specify the
C     variable dimension of F.  IDIMF must be at least M.
C
C   W
C     A one-dimensional array that must be provided by the user for
C     work space.  W may require up to 13M + 4N + M*INT(log2(N))
C     locations.  The actual number of locations used is computed by
C     HSTSSP and is returned in the location W(1).
C
C
C            * * * * * *   On Output   * * * * * *
C
C   F
C     Contains the solution U(I,J) of the finite difference
C     approximation for the grid point (THETA(I),PHI(J)) for
C     I=1,2,...,M, J=1,2,...,N.
C
C   PERTRB
C     If a combination of periodic, derivative, or unspecified
C     boundary conditions is specified for a Poisson equation
C     (LAMBDA = 0), a solution may not exist.  PERTRB is a con-
C     stant, calculated and subtracted from F, which ensures
C     that a solution exists.  HSTSSP then computes this
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
C      Except for numbers 0 and 14, a solution is not attempted.
C
C     =  0  No error
C
C     =  1  A .LT. 0 or B .GT. PI
C
C     =  2  A .GE. B
C
C     =  3  MBDCND .LT. 1 or MBDCND .GT. 9
C
C     =  4  C .GE. D
C
C     =  5  N .LE. 2
C
C     =  6  NBDCND .LT. 0 or NBDCND .GT. 4
C
C     =  7  A .GT. 0 and MBDCND = 5, 6, or 9
C
C     =  8  A = 0 and MBDCND = 3, 4, or 8
C
C     =  9  B .LT. PI and MBDCND .GE. 7
C
C     = 10  B = PI and MBDCND = 2,3, or 6
C
C     = 11  MBDCND .GE. 5 and NDBCND = 1, 2, or 4
C
C     = 12  IDIMF .LT. M
C
C     = 13  M .LE. 2
C
C     = 14  LAMBDA .GT. 0
C
C     Since this is the only means of indicating a possibly
C     incorrect call to HSTSSP, the user should test IERROR after
C     the call.
C
C   W
C     W(1) contains the required length of W.
C
C *Long Description:
C
C    * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C    Dimension of   BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
C    Arguments      W(see argument list)
C
C    Latest         June 1, 1977
C    Revision
C
C    Subprograms    HSTSSP,POISTG,POSTG2,GENBUN,POISD2,POISN2,POISP2,
C    Required       COSGEN,MERGE,TRIX,TRI3,PIMACH
C
C    Special        NONE
C    Conditions
C
C    Common         NONE
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
C    History        Written by Roland Sweet at NCAR in April, 1977
C
C    Algorithm      This subroutine defines the finite-difference
C                   equations, incorporates boundary data, adjusts the
C                   right side when the system is singular and calls
C                   either POISTG or GENBUN which solves the linear
C                   system of equations.
C
C    Space          8427(decimal) = 20353(octal) locations on the
C    Required       NCAR Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine HSTSSP is roughly proportional
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
C                        32       1-9       1-4         56
C                        64       1-9       1-4        230
C
C    Portability     American National Standards Institute FORTRAN.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C    Required       COS
C    Resident
C    Routines
C
C    Reference      Schumann, U. and R. Sweet,'A Direct Method For
C                   The Solution Of Poisson's Equation With Neumann
C                   Boundary Conditions On A Staggered Grid Of
C                   Arbitrary Size,' J. Comp. Phys. 20(1976),
C                   pp. 171-182.
C
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  U. Schumann and R. Sweet, A direct method for the
C                 solution of Poisson's equation with Neumann boundary
C                 conditions on a staggered grid of arbitrary size,
C                 Journal of Computational Physics 20, (1976),
C                 pp. 171-182.
C***ROUTINES CALLED  GENBUN, PIMACH, POISTG
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  HSTSSP
C
C
      DIMENSION       F(IDIMF,*) ,BDA(*)     ,BDB(*)     ,BDC(*)     ,
     1                BDD(*)     ,W(*)
C***FIRST EXECUTABLE STATEMENT  HSTSSP
      IERROR = 0
      PI = PIMACH(DUM)
      IF (A.LT.0. .OR. B.GT.PI) IERROR = 1
      IF (A .GE. B) IERROR = 2
      IF (MBDCND.LE.0 .OR. MBDCND.GT.9) IERROR = 3
      IF (C .GE. D) IERROR = 4
      IF (N .LE. 2) IERROR = 5
      IF (NBDCND.LT.0 .OR. NBDCND.GE.5) IERROR = 6
      IF (A.GT.0. .AND. (MBDCND.EQ.5 .OR. MBDCND.EQ.6 .OR. MBDCND.EQ.9))
     1    IERROR = 7
      IF (A.EQ.0. .AND. (MBDCND.EQ.3 .OR. MBDCND.EQ.4 .OR. MBDCND.EQ.8))
     1    IERROR = 8
      IF (B.LT.PI .AND. MBDCND.GE.7) IERROR = 9
      IF (B.EQ.PI .AND. (MBDCND.EQ.2 .OR. MBDCND.EQ.3 .OR. MBDCND.EQ.6))
     1    IERROR = 10
      IF (MBDCND.GE.5 .AND.
     1    (NBDCND.EQ.1 .OR. NBDCND.EQ.2 .OR. NBDCND.EQ.4)) IERROR = 11
      IF (IDIMF .LT. M) IERROR = 12
      IF (M .LE. 2) IERROR = 13
      IF (IERROR .NE. 0) RETURN
      DELTAR = (B-A)/M
      DLRSQ = DELTAR**2
      DELTHT = (D-C)/N
      DLTHSQ = DELTHT**2
      NP = NBDCND+1
      ISW = 1
      JSW = 1
      MB = MBDCND
      IF (ELMBDA .NE. 0.) GO TO 105
      GO TO (101,102,105,103,101,105,101,105,105),MBDCND
  101 IF (A.NE.0. .OR. B.NE.PI) GO TO 105
      MB = 9
      GO TO 104
  102 IF (A .NE. 0.) GO TO 105
      MB = 6
      GO TO 104
  103 IF (B .NE. PI) GO TO 105
      MB = 8
  104 JSW = 2
  105 CONTINUE
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      IWB = M
      IWC = IWB+M
      IWR = IWC+M
      IWS = IWR+M
      DO 106 I=1,M
         J = IWR+I
         W(J) = SIN(A+(I-0.5)*DELTAR)
         W(I) = SIN((A+(I-1)*DELTAR))/DLRSQ
  106 CONTINUE
      MM1 = M-1
      DO 107 I=1,MM1
         K = IWC+I
         W(K) = W(I+1)
         J = IWR+I
         K = IWB+I
         W(K) = ELMBDA*W(J)-(W(I)+W(I+1))
  107 CONTINUE
      W(IWR) = SIN(B)/DLRSQ
      W(IWC) = ELMBDA*W(IWS)-(W(M)+W(IWR))
      DO 109 I=1,M
         J = IWR+I
         A1 = W(J)
         DO 108 J=1,N
            F(I,J) = A1*F(I,J)
  108    CONTINUE
  109 CONTINUE
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
      GO TO (110,110,112,112,114,114,110,112,114),MB
  110 A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1)-W(1)
      DO 111 J=1,N
         F(1,J) = F(1,J)-A1*BDA(J)
  111 CONTINUE
      GO TO 114
  112 A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1)+W(1)
      DO 113 J=1,N
         F(1,J) = F(1,J)+A1*BDA(J)
  113 CONTINUE
  114 GO TO (115,117,117,115,115,117,119,119,119),MB
  115 A1 = 2.*W(IWR)
      W(IWC) = W(IWC)-W(IWR)
      DO 116 J=1,N
         F(M,J) = F(M,J)-A1*BDB(J)
  116 CONTINUE
      GO TO 119
  117 A1 = DELTAR*W(IWR)
      W(IWC) = W(IWC)+W(IWR)
      DO 118 J=1,N
         F(M,J) = F(M,J)-A1*BDB(J)
  118 CONTINUE
C
C     ENTER BOUNDARY DATA FOR PHI-BOUNDARIES.
C
  119 A1 = 2./DLTHSQ
      GO TO (129,120,120,122,122),NP
  120 DO 121 I=1,M
         J = IWR+I
         F(I,1) = F(I,1)-A1*BDC(I)/W(J)
  121 CONTINUE
      GO TO 124
  122 A1 = 1./DELTHT
      DO 123 I=1,M
         J = IWR+I
         F(I,1) = F(I,1)+A1*BDC(I)/W(J)
  123 CONTINUE
  124 A1 = 2./DLTHSQ
      GO TO (129,125,127,127,125),NP
  125 DO 126 I=1,M
         J = IWR+I
         F(I,N) = F(I,N)-A1*BDD(I)/W(J)
  126 CONTINUE
      GO TO 129
  127 A1 = 1./DELTHT
      DO 128 I=1,M
         J = IWR+I
         F(I,N) = F(I,N)-A1*BDD(I)/W(J)
  128 CONTINUE
  129 CONTINUE
C
C     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
C     SOLUTION.
C
      PERTRB = 0.
      IF (ELMBDA) 139,131,130
  130 IERROR = 14
      GO TO 139
  131 GO TO (139,139,132,139,139,132,139,132,132),MB
  132 GO TO (133,139,139,133,139),NP
  133 CONTINUE
      ISW = 2
      DO 135 J=1,N
         DO 134 I=1,M
            PERTRB = PERTRB+F(I,J)
  134    CONTINUE
  135 CONTINUE
      A1 = N*(COS(A)-COS(B))/(2.*SIN(0.5*DELTAR))
      PERTRB = PERTRB/A1
      DO 137 I=1,M
         J = IWR+I
         A1 = PERTRB*W(J)
         DO 136 J=1,N
            F(I,J) = F(I,J)-A1
  136    CONTINUE
  137 CONTINUE
      A2 = 0.
      A3 = 0.
      DO 138 J=1,N
         A2 = A2+F(1,J)
         A3 = A3+F(M,J)
  138 CONTINUE
      A2 = A2/W(IWR+1)
      A3 = A3/W(IWS)
  139 CONTINUE
C
C     MULTIPLY I-TH EQUATION THROUGH BY  R(I)*DELTHT**2
C
      DO 141 I=1,M
         J = IWR+I
         A1 = DLTHSQ*W(J)
         W(I) = A1*W(I)
         J = IWC+I
         W(J) = A1*W(J)
         J = IWB+I
         W(J) = A1*W(J)
         DO 140 J=1,N
            F(I,J) = A1*F(I,J)
  140    CONTINUE
  141 CONTINUE
      LP = NBDCND
      W(1) = 0.
      W(IWR) = 0.
C
C     CALL POISTG OR GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
C
      IF (NBDCND .EQ. 0) GO TO 142
      CALL POISTG (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
      GO TO 143
  142 CALL GENBUN (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  143 CONTINUE
      W(1) = W(IWR+1)+3*M
      IF (ISW.NE.2 .OR. JSW.NE.2) GO TO 150
      IF (MB .NE. 8) GO TO 145
      A1 = 0.
      DO 144 J=1,N
         A1 = A1+F(M,J)
  144 CONTINUE
      A1 = (A1-DLRSQ*A3/16.)/N
      IF (NBDCND .EQ. 3) A1 = A1+(BDD(M)-BDC(M))/(D-C)
      A1 = BDB(1)-A1
      GO TO 147
  145 A1 = 0.
      DO 146 J=1,N
         A1 = A1+F(1,J)
  146 CONTINUE
      A1 = (A1-DLRSQ*A2/16.)/N
      IF (NBDCND .EQ. 3) A1 = A1+(BDD(1)-BDC(1))/(D-C)
      A1 = BDA(1)-A1
  147 DO 149 I=1,M
         DO 148 J=1,N
            F(I,J) = F(I,J)+A1
  148    CONTINUE
  149 CONTINUE
  150 CONTINUE
      RETURN
      END
