*DECK HWSSSP
      SUBROUTINE HWSSSP (TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N,
     +   NBDCND, BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HWSSSP
C***PURPOSE  Solve a finite difference approximation to the Helmholtz
C            equation in spherical coordinates and on the surface of the
C            unit sphere (radius of 1).
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HWSSSP-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine HWSSSP solves a finite difference approximation to the
C     Helmholtz equation in spherical coordinates and on the surface of
C     the unit sphere (radius of 1):
C
C          (1/SIN(THETA))(d/dTHETA)(SIN(THETA)(dU/dTHETA))
C
C             + (1/SIN(THETA)**2)(d/dPHI)(dU/dPHI)
C
C             + LAMBDA*U = F(THETA,PHI)
C
C     Where THETA is colatitude and PHI is longitude.
C
C     * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C             * * * * * *   On Input    * * * * * *
C
C     TS,TF
C       The range of THETA (colatitude), i.e., TS .LE. THETA .LE. TF.
C       TS must be less than TF.  TS and TF are in radians.  A TS of
C       zero corresponds to the north pole and a TF of PI corresponds to
C       the south pole.
C
C     * * * * * * * * * * * * * * IMPORTANT * * * * * * * * * * * * * *
C
C     If TF is equal to PI then it must be computed using the statement
C     TF = PIMACH(DUM). This insures that TF in the users program is
C     equal to PI in this program which permits several tests of the
C     input parameters that otherwise would not be possible.
C
C
C     M
C       The number of panels into which the interval (TS,TF) is
C       subdivided.  Hence, there will be M+1 grid points in the
C       THETA-direction given by THETA(I) = (I-1)DTHETA+TS for
C       I = 1,2,...,M+1, where DTHETA = (TF-TS)/M is the panel width.
C       M must be greater than 5.
C
C     MBDCND
C       Indicates the type of boundary condition at THETA = TS and
C       THETA = TF.
C
C       = 1  If the solution is specified at THETA = TS and THETA = TF.
C       = 2  If the solution is specified at THETA = TS and the
C            derivative of the solution with respect to THETA is
C            specified at THETA = TF (see note 2 below).
C       = 3  If the derivative of the solution with respect to THETA is
C            specified at THETA = TS and THETA = TF (see notes 1,2
C            below).
C       = 4  If the derivative of the solution with respect to THETA is
C            specified at THETA = TS (see note 1 below) and the
C            solution is specified at THETA = TF.
C       = 5  If the solution is unspecified at THETA = TS = 0 and the
C            solution is specified at THETA = TF.
C       = 6  If the solution is unspecified at THETA = TS = 0 and the
C            derivative of the solution with respect to THETA is
C            specified at THETA = TF (see note 2 below).
C       = 7  If the solution is specified at THETA = TS and the
C            solution is unspecified at THETA = TF = PI.
C       = 8  If the derivative of the solution with respect to THETA is
C            specified at THETA = TS (see note 1 below) and the
C            solution is unspecified at THETA = TF = PI.
C       = 9  If the solution is unspecified at THETA = TS = 0 and
C            THETA = TF = PI.
C
C       NOTES:  1.  If TS = 0, do not use MBDCND = 3,4, or 8, but
C                   instead use MBDCND = 5,6, or 9  .
C               2.  If TF = PI, do not use MBDCND = 2,3, or 6, but
C                   instead use MBDCND = 7,8, or 9  .
C
C     BDTS
C       A one-dimensional array of length N+1 that specifies the values
C       of the derivative of the solution with respect to THETA at
C       THETA = TS.  When MBDCND = 3,4, or 8,
C
C            BDTS(J) = (d/dTHETA)U(TS,PHI(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value, BDTS is a dummy variable.
C
C     BDTF
C       A one-dimensional array of length N+1 that specifies the values
C       of the derivative of the solution with respect to THETA at
C       THETA = TF.  When MBDCND = 2,3, or 6,
C
C            BDTF(J) = (d/dTHETA)U(TF,PHI(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value, BDTF is a dummy variable.
C
C     PS,PF
C       The range of PHI (longitude), i.e., PS .LE. PHI .LE. PF.  PS
C       must be less than PF.  PS and PF are in radians.  If PS = 0 and
C       PF = 2*PI, periodic boundary conditions are usually prescribed.
C
C     * * * * * * * * * * * * * * IMPORTANT * * * * * * * * * * * * * *
C
C     If PF is equal to 2*PI then it must be computed using the
C     statement PF = 2.*PIMACH(DUM). This insures that PF in the users
C     program is equal to 2*PI in this program which permits tests of
C     the input parameters that otherwise would not be possible.
C
C
C     N
C       The number of panels into which the interval (PS,PF) is
C       subdivided.  Hence, there will be N+1 grid points in the
C       PHI-direction given by PHI(J) = (J-1)DPHI+PS  for
C       J = 1,2,...,N+1, where DPHI = (PF-PS)/N is the panel width.
C       N must be greater than 4.
C
C     NBDCND
C       Indicates the type of boundary condition at PHI = PS and
C       PHI = PF.
C
C       = 0  If the solution is periodic in PHI, i.e.,
C            U(I,J) = U(I,N+J).
C       = 1  If the solution is specified at PHI = PS and PHI = PF
C            (see note below).
C       = 2  If the solution is specified at PHI = PS (see note below)
C            and the derivative of the solution with respect to PHI is
C            specified at PHI = PF.
C       = 3  If the derivative of the solution with respect to PHI is
C            specified at PHI = PS and PHI = PF.
C       = 4  If the derivative of the solution with respect to PHI is
C            specified at PS and the solution is specified at PHI = PF
C            (see note below).
C
C       NOTE:  NBDCND = 1,2, or 4 cannot be used with
C              MBDCND = 5,6,7,8, or 9 (the former indicates that the
C                       solution is specified at a pole, the latter
C                       indicates that the solution is unspecified).
C                       Use instead
C              MBDCND = 1 or 2  .
C
C     BDPS
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to PHI at
C       PHI = PS.  When NBDCND = 3 or 4,
C
C            BDPS(I) = (d/dPHI)U(THETA(I),PS), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDPS is a dummy variable.
C
C     BDPF
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to PHI at
C       PHI = PF.  When NBDCND = 2 or 3,
C
C            BDPF(I) = (d/dPHI)U(THETA(I),PF), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDPF is a dummy variable.
C
C     ELMBDA
C       The constant LAMBDA in the Helmholtz equation.  If
C       LAMBDA .GT. 0, a solution may not exist.  However, HWSSSP will
C       attempt to find a solution.
C
C     F
C       A two-dimensional array that specifies the value of the right
C       side of the Helmholtz equation and boundary values (if any).
C       For I = 2,3,...,M  and  J = 2,3,...,N
C
C            F(I,J) = F(THETA(I),PHI(J)).
C
C       On the boundaries F is defined by
C
C            MBDCND   F(1,J)            F(M+1,J)
C            ------   ------------      ------------
C
C              1      U(TS,PHI(J))      U(TF,PHI(J))
C              2      U(TS,PHI(J))      F(TF,PHI(J))
C              3      F(TS,PHI(J))      F(TF,PHI(J))
C              4      F(TS,PHI(J))      U(TF,PHI(J))
C              5      F(0,PS)           U(TF,PHI(J))   J = 1,2,...,N+1
C              6      F(0,PS)           F(TF,PHI(J))
C              7      U(TS,PHI(J))      F(PI,PS)
C              8      F(TS,PHI(J))      F(PI,PS)
C              9      F(0,PS)           F(PI,PS)
C
C            NBDCND   F(I,1)            F(I,N+1)
C            ------   --------------    --------------
C
C              0      F(THETA(I),PS)    F(THETA(I),PS)
C              1      U(THETA(I),PS)    U(THETA(I),PF)
C              2      U(THETA(I),PS)    F(THETA(I),PF)   I = 1,2,...,M+1
C              3      F(THETA(I),PS)    F(THETA(I),PF)
C              4      F(THETA(I),PS)    U(THETA(I),PF)
C
C       F must be dimensioned at least (M+1)*(N+1).
C
C      *NOTE*
C
C       If the table calls for both the solution U and the right side F
C       at a corner then the solution must be specified.
C
C
C     IDIMF
C       The row (or first) dimension of the array F as it appears in the
C       program calling HWSSSP.  This parameter is used to specify the
C       variable dimension of F.  IDIMF must be at least M+1  .
C
C     W
C       A one-dimensional array that must be provided by the user for
C       work space. W may require up to 4*(N+1)+(16+INT(log2(N+1)))(M+1)
C       locations. The actual number of locations used is computed by
C       HWSSSP and is output in location W(1). INT( ) denotes the
C       FORTRAN integer function.
C
C
C     * * * * * * * * * *     On Output     * * * * * * * * * *
C
C     F
C       Contains the solution U(I,J) of the finite difference
C       approximation for the grid point (THETA(I),PHI(J)),
C       I = 1,2,...,M+1,   J = 1,2,...,N+1  .
C
C     PERTRB
C       If one specifies a combination of periodic, derivative or
C       unspecified boundary conditions for a Poisson equation
C       (LAMBDA = 0), a solution may not exist.  PERTRB is a constant,
C       calculated and subtracted from F, which ensures that a solution
C       exists.  HWSSSP then computes this solution, which is a least
C       squares solution to the original approximation.  This solution
C       is not unique and is unnormalized. The value of PERTRB should
C       be small compared to the right side F. Otherwise , a solution
C       is obtained to an essentially different problem. This comparison
C       should always be made to insure that a meaningful solution has
C       been obtained.
C
C     IERROR
C       An error flag that indicates invalid input parameters.  Except
C       for numbers 0 and 8, a solution is not attempted.
C
C       = 0  No error
C       = 1  TS.LT.0 or TF.GT.PI
C       = 2  TS.GE.TF
C       = 3  MBDCND.LT.1 or MBDCND.GT.9
C       = 4  PS.LT.0 or PS.GT.PI+PI
C       = 5  PS.GE.PF
C       = 6  N.LT.5
C       = 7  M.LT.5
C       = 8  NBDCND.LT.0 or NBDCND.GT.4
C       = 9  ELMBDA.GT.0
C       = 10 IDIMF.LT.M+1
C       = 11 NBDCND equals 1,2 or 4 and MBDCND.GE.5
C       = 12 TS.EQ.0 and MBDCND equals 3,4 or 8
C       = 13 TF.EQ.PI and MBDCND equals 2,3 or 6
C       = 14 MBDCND equals 5,6 or 9 and TS.NE.0
C       = 15 MBDCND.GE.7 and TF.NE.PI
C
C       Since this is the only means of indicating a possibly incorrect
C       call to HWSSSP, the user should test IERROR after a call.
C
C     W
C       Contains intermediate values that must not be destroyed if
C       HWSSSP will be called again with INTL = 1. W(1) contains the
C       required length of W .
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   BDTS(N+1),BDTF(N+1),BDPS(M+1),BDPF(M+1),
C     Arguments      F(IDIMF,N+1),W(see argument list)
C
C     Latest         January 1978
C     Revision
C
C
C     Subprograms    HWSSSP,HWSSS1,GENBUN,POISD2,POISN2,POISP2,COSGEN,ME
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
C     Specialist     Paul Swarztrauber
C
C     Language       FORTRAN
C
C     History        Version 1 - September 1973
C                    Version 2 - April     1976
C                    Version 3 - January   1978
C
C     Algorithm      The routine defines the finite difference
C                    equations, incorporates boundary data, and adjusts
C                    the right side of singular systems and then calls
C                    GENBUN to solve the system.
C
C     Space
C     Required       CONTROL DATA 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine HWSSSP is roughly proportional
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
C     Required       SIN,COS
C     Resident
C     Routines
C
C     References     P. N. Swarztrauber,'The Direct Solution Of The
C                    Discrete Poisson Equation On The Surface Of a
C                    Sphere, SIAM J. Numer. Anal.,15(1974), pp 212-215
C
C                    Swarztrauber,P. and R. Sweet, 'Efficient FORTRAN
C                    Subprograms for The Solution of Elliptic Equations'
C                    NCAR TN/IA-109, July, 1975, 138 pp.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
C                 subprograms for the solution of elliptic equations,
C                 NCAR TN/IA-109, July 1975, 138 pp.
C               P. N. Swarztrauber, The direct solution of the discrete
C                 Poisson equation on the surface of a sphere, SIAM
C                 Journal on Numerical Analysis 15 (1974), pp. 212-215.
C***ROUTINES CALLED  HWSSS1, PIMACH
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891009  Removed unreferenced variable.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  HWSSSP
C
      DIMENSION       F(IDIMF,*) ,BDTS(*)    ,BDTF(*)    ,BDPS(*)    ,
     1                BDPF(*)    ,W(*)
C***FIRST EXECUTABLE STATEMENT  HWSSSP
      PI = PIMACH(DUM)
      TPI = 2.*PI
      IERROR = 0
      IF (TS.LT.0. .OR. TF.GT.PI) IERROR = 1
      IF (TS .GE. TF) IERROR = 2
      IF (MBDCND.LT.1 .OR. MBDCND.GT.9) IERROR = 3
      IF (PS.LT.0. .OR. PF.GT.TPI) IERROR = 4
      IF (PS .GE. PF) IERROR = 5
      IF (N .LT. 5) IERROR = 6
      IF (M .LT. 5) IERROR = 7
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR = 8
      IF (ELMBDA .GT. 0.) IERROR = 9
      IF (IDIMF .LT. M+1) IERROR = 10
      IF ((NBDCND.EQ.1 .OR. NBDCND.EQ.2 .OR. NBDCND.EQ.4) .AND.
     1    MBDCND.GE.5) IERROR = 11
      IF (TS.EQ.0. .AND.
     1    (MBDCND.EQ.3 .OR. MBDCND.EQ.4 .OR. MBDCND.EQ.8)) IERROR = 12
      IF (TF.EQ.PI .AND.
     1    (MBDCND.EQ.2 .OR. MBDCND.EQ.3 .OR. MBDCND.EQ.6)) IERROR = 13
      IF ((MBDCND.EQ.5 .OR. MBDCND.EQ.6 .OR. MBDCND.EQ.9) .AND.
     1    TS.NE.0.) IERROR = 14
      IF (MBDCND.GE.7 .AND. TF.NE.PI) IERROR = 15
      IF (IERROR.NE.0 .AND. IERROR.NE.9) RETURN
      CALL HWSSS1 (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,BDPF,
     1             ELMBDA,F,IDIMF,PERTRB,W,W(M+2),W(2*M+3),W(3*M+4),
     2             W(4*M+5),W(5*M+6),W(6*M+7))
      W(1) = W(6*M+7)+6*(M+1)
      RETURN
      END
