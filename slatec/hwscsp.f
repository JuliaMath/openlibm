*DECK HWSCSP
      SUBROUTINE HWSCSP (INTL, TS, TF, M, MBDCND, BDTS, BDTF, RS, RF, N,
     +   NBDCND, BDRS, BDRF, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HWSCSP
C***PURPOSE  Solve a finite difference approximation to the modified
C            Helmholtz equation in spherical coordinates assuming
C            axisymmetry  (no dependence on longitude).
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HWSCSP-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine HWSCSP solves a finite difference approximation to the
C       modified Helmholtz equation in spherical coordinates assuming
C       axisymmetry  (no dependence on longitude)
C
C          (1/R**2)(d/dR)((R**2)(d/dR)U)
C
C             + (1/(R**2)SIN(THETA))(d/dTHETA)(SIN(THETA)(d/dTHETA)U)
C
C             + (LAMBDA/(RSIN(THETA))**2)U = F(THETA,R).
C
C     This two dimensional modified Helmholtz equation results from
C     the Fourier transform of the three dimensional Poisson equation
C
C     * * * * * * * * * *     On Input     * * * * * * * * * *
C
C     INTL
C       = 0  On initial entry to HWSCSP or if any of the arguments
C            RS, RF, N, NBDCND are changed from a previous call.
C       = 1  If RS, RF, N, NBDCND are all unchanged from previous call
C            to HWSCSP.
C
C       NOTE   A call with INTL=0 takes approximately 1.5 times as
C              much time as a call with INTL = 1.  Once a call with
C              INTL = 0 has been made then subsequent solutions
C              corresponding to different F, BDTS, BDTF, BDRS, BDRF can
C              be obtained faster with INTL = 1 since initialization is
C              not repeated.
C
C     TS,TF
C       The range of THETA (colatitude), i.e., TS .LE. THETA .LE. TF.
C       TS must be less than TF.  TS and TF are in radians.  A TS of
C       zero corresponds to the north pole and a TF of PI corresponds
C       to the south pole.
C
C     * * * * * * * * * * * * * * IMPORTANT * * * * * * * * * * * * * *
C
C     If TF is equal to PI then it must be computed using the statement
C     TF = PIMACH(DUM). This insures that TF in the users program is
C     equal to PI in this program which permits several tests of the
C     input parameters that otherwise would not be possible.
C
C     M
C       The number of panels into which the interval (TS,TF) is
C       subdivided.  Hence, there will be M+1 grid points in the
C       THETA-direction given by THETA(K) = (I-1)DTHETA+TS for
C       I = 1,2,...,M+1, where DTHETA = (TF-TS)/M is the panel width.
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
C            specified at THETA = TS (see note 1 below) and the solution
C            is unspecified at THETA = TF = PI.
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
C            BDTS(J) = (d/dTHETA)U(TS,R(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value, BDTS is a dummy variable.
C
C     BDTF
C       A one-dimensional array of length N+1 that specifies the values
C       of the derivative of the solution with respect to THETA at
C       THETA = TF.  When MBDCND = 2,3, or 6,
C
C            BDTF(J) = (d/dTHETA)U(TF,R(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value, BDTF is a dummy variable.
C
C     RS,RF
C       The range of R, i.e., RS .LE. R .LT. RF.  RS must be less than
C       RF.  RS must be non-negative.
C
C       N
C       The number of panels into which the interval (RS,RF) is
C       subdivided.  Hence, there will be N+1 grid points in the
C       R-direction given by R(J) = (J-1)DR+RS for J = 1,2,...,N+1,
C       where DR = (RF-RS)/N is the panel width.
C       N must be greater than 2
C
C     NBDCND
C       Indicates the type of boundary condition at R = RS and R = RF.
C
C       = 1  If the solution is specified at R = RS and R = RF.
C       = 2  If the solution is specified at R = RS and the derivative
C            of the solution with respect to R is specified at R = RF.
C       = 3  If the derivative of the solution with respect to R is
C            specified at R = RS and R = RF.
C       = 4  If the derivative of the solution with respect to R is
C            specified at RS and the solution is specified at R = RF.
C       = 5  If the solution is unspecified at R = RS = 0 (see note
C            below) and the solution is specified at R = RF.
C       = 6  If the solution is unspecified at R = RS = 0 (see note
C            below) and the derivative of the solution with respect to
C            R is specified at R = RF.
C
C       NOTE:  NBDCND = 5 or 6 cannot be used with
C              MBDCND = 1,2,4,5, or 7 (the former indicates that the
C                       solution is unspecified at R = 0, the latter
C                       indicates that the solution is specified).
C                       Use instead
C              NBDCND = 1 or 2  .
C
C     BDRS
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to R at R = RS.
C       When NBDCND = 3 or 4,
C
C            BDRS(I) = (d/dR)U(THETA(I),RS), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDRS is a dummy variable.
C
C     BDRF
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to R at R = RF.
C       When NBDCND = 2,3, or 6,
C
C            BDRF(I) = (d/dR)U(THETA(I),RF), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDRF is a dummy variable.
C
C     ELMBDA
C       The constant LAMBDA in the Helmholtz equation.  If
C       LAMBDA .GT. 0, a solution may not exist.  However, HWSCSP will
C       attempt to find a solution.  If NBDCND = 5 or 6 or
C       MBDCND = 5,6,7,8, or 9, ELMBDA must be zero.
C
C     F
C       A two-dimensional array that specifies the value of the right
C       side of the Helmholtz equation and boundary values (if any).
C       for I = 2,3,...,M and J = 2,3,...,N
C
C            F(I,J) = F(THETA(I),R(J)).
C
C       On the boundaries F is defined by
C
C            MBDCND   F(1,J)            F(M+1,J)
C            ------   ----------        ----------
C
C              1      U(TS,R(J))        U(TF,R(J))
C              2      U(TS,R(J))        F(TF,R(J))
C              3      F(TS,R(J))        F(TF,R(J))
C              4      F(TS,R(J))        U(TF,R(J))
C              5      F(0,R(J))         U(TF,R(J))   J = 1,2,...,N+1
C              6      F(0,R(J))         F(TF,R(J))
C              7      U(TS,R(J))        F(PI,R(J))
C              8      F(TS,R(J))        F(PI,R(J))
C              9      F(0,R(J))         F(PI,R(J))
C
C            NBDCND   F(I,1)            F(I,N+1)
C            ------   --------------    --------------
C
C              1      U(THETA(I),RS)    U(THETA(I),RF)
C              2      U(THETA(I),RS)    F(THETA(I),RF)
C              3      F(THETA(I),RS)    F(THETA(I),RF)
C              4      F(THETA(I),RS)    U(THETA(I),RF)   I = 1,2,...,M+1
C              5      F(TS,0)           U(THETA(I),RF)
C              6      F(TS,0)           F(THETA(I),RF)
C
C       F must be dimensioned at least (M+1)*(N+1).
C
C       NOTE
C
C       If the table calls for both the solution U and the right side F
C       at a corner then the solution must be specified.
C
C     IDIMF
C       The row (or first) dimension of the array F as it appears in the
C       program calling HWSCSP.  This parameter is used to specify the
C       variable dimension of F.  IDIMF must be at least M+1  .
C
C     W
C       A one-dimensional array that must be provided by the user for
C       work space. Its length can be computed from the formula below
C       which depends on the value of NBDCND.
C
C       If NBDCND=2,4 or 6 define NUNK=N
C       If NBDCND=1 or 5   define NUNK=N-1
C       If NBDCND=3        define NUNK=N+1
C
C       Now set K=INT(log2(NUNK))+1 and L=2**(K+1) then W must be
C       dimensioned at least (K-2)*L+K+5*(M+N)+MAX(2*N,6*M)+23
C
C       **IMPORTANT** For purposes of checking, the required length
C                     of W is computed by HWSCSP and stored in W(1)
C                     in floating point format.
C
C
C     * * * * * * * * * *     On Output     * * * * * * * * * *
C
C     F
C       Contains the solution U(I,J) of the finite difference
C       approximation for the grid point (THETA(I),R(J)),
C       I = 1,2,...,M+1,   J = 1,2,...,N+1  .
C
C     PERTRB
C       If a combination of periodic or derivative boundary conditions
C       is specified for a Poisson equation (LAMBDA = 0), a solution may
C       not exist.  PERTRB is a constant, calculated and subtracted from
C       F, which ensures that a solution exists.  HWSCSP then computes
C       this solution, which is a least squares solution to the original
C       approximation. This solution is not unique and is unnormalized.
C       The value of PERTRB should be small compared to the right side
C       F. Otherwise , a solution is obtained to an essentially
C       different problem. This comparison should always be made to
C       insure that a meaningful solution has been obtained.
C
C     IERROR
C       An error flag that indicates invalid input parameters.  Except
C       for numbers 0 and 10, a solution is not attempted.
C
C       = 1  TS.LT.0. or TF.GT.PI
C       = 2  TS.GE.TF
C       = 3  M.LT.5
C       = 4  MBDCND.LT.1 or MBDCND.GT.9
C       = 5  RS.LT.0
C       = 6  RS.GE.RF
C       = 7  N.LT.5
C       = 8  NBDCND.LT.1 or NBDCND.GT.6
C       = 9  ELMBDA.GT.0
C       = 10 IDIMF.LT.M+1
C       = 11 ELMBDA.NE.0 and MBDCND.GE.5
C       = 12 ELMBDA.NE.0 and NBDCND equals 5 or 6
C       = 13 MBDCND equals 5,6 or 9 and TS.NE.0
C       = 14 MBDCND.GE.7 and TF.NE.PI
C       = 15 TS.EQ.0 and MBDCND equals 3,4 or 8
C       = 16 TF.EQ.PI and MBDCND equals 2,3 or 6
C       = 17 NBDCND.GE.5 and RS.NE.0
C       = 18 NBDCND.GE.5 and MBDCND equals 1,2,4,5 or 7
C
C       Since this is the only means of indicating a possibly incorrect
C       call to HWSCSP, the user should test IERROR after a call.
C
C     W
C       Contains intermediate values that must not be destroyed if
C       HWSCSP will be called again with INTL = 1.  W(1) contains the
C       number of locations which W must have.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   BDTS(N+1),BDTF(N+1),BDRS(M+1),BDRF(M+1),
C     Arguments      F(IDIMF,N+1),W(see argument list)
C
C     Latest         June 1979
C     Revision
C
C     Subprograms    HWSCSP,HWSCS1,BLKTRI,BLKTR1,PROD,PRODP,CPROD,CPRODP
C     Required       ,COMBP,PPADD,PSGF,BSRH,PPSGF,PPSPF,TEVLS,INDXA,
C                    ,INDXB,INDXC,R1MACH
C
C     Special
C     Conditions
C
C     Common         CBLKT
C     Blocks
C
C     I/O            NONE
C
C     Precision      Single
C
C     Specialist     Paul N Swarztrauber
C
C     Language       FORTRAN
C
C     History        Version 1 September 1973
C                    Version 2 April     1976
C                    Version 3 June      1979
C
C     Algorithm      The routine defines the finite difference
C                    equations, incorporates boundary data, and adjusts
C                    the right side of singular systems and then calls
C                    BLKTRI to solve the system.
C
C     Space
C     Required
C
C     Portability    American National Standards Institute FORTRAN.
C                    The machine accuracy is set using function R1MACH.
C
C     Required       NONE
C     Resident
C     Routines
C
C     Reference      Swarztrauber,P. and R. Sweet, 'Efficient FORTRAN
C                    Subprograms for The Solution Of Elliptic Equations'
C                    NCAR TN/IA-109, July, 1975, 138 pp.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
C                 subprograms for the solution of elliptic equations,
C                 NCAR TN/IA-109, July 1975, 138 pp.
C***ROUTINES CALLED  HWSCS1, PIMACH
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  HWSCSP
C
      DIMENSION       F(IDIMF,*) ,BDTS(*)    ,BDTF(*)    ,BDRS(*)    ,
     1                BDRF(*)    ,W(*)
C***FIRST EXECUTABLE STATEMENT  HWSCSP
      PI = PIMACH(DUM)
      IERROR = 0
      IF (TS.LT.0. .OR. TF.GT.PI) IERROR = 1
      IF (TS .GE. TF) IERROR = 2
      IF (M .LT. 5) IERROR = 3
      IF (MBDCND.LT.1 .OR. MBDCND.GT.9) IERROR = 4
      IF (RS .LT. 0.) IERROR = 5
      IF (RS .GE. RF) IERROR = 6
      IF (N .LT. 5) IERROR = 7
      IF (NBDCND.LT.1 .OR. NBDCND.GT.6) IERROR = 8
      IF (ELMBDA .GT. 0.) IERROR = 9
      IF (IDIMF .LT. M+1) IERROR = 10
      IF (ELMBDA.NE.0. .AND. MBDCND.GE.5) IERROR = 11
      IF (ELMBDA.NE.0. .AND. (NBDCND.EQ.5 .OR. NBDCND.EQ.6)) IERROR = 12
      IF ((MBDCND.EQ.5 .OR. MBDCND.EQ.6 .OR. MBDCND.EQ.9) .AND.
     1    TS.NE.0.) IERROR = 13
      IF (MBDCND.GE.7 .AND. TF.NE.PI) IERROR = 14
      IF (TS.EQ.0. .AND.
     1    (MBDCND.EQ.4 .OR. MBDCND.EQ.8 .OR. MBDCND.EQ.3)) IERROR = 15
      IF (TF.EQ.PI .AND.
     1    (MBDCND.EQ.2 .OR. MBDCND.EQ.3 .OR. MBDCND.EQ.6)) IERROR = 16
      IF (NBDCND.GE.5 .AND. RS.NE.0.) IERROR = 17
      IF (NBDCND.GE.5 .AND. (MBDCND.EQ.1 .OR. MBDCND.EQ.2 .OR.
     1                                    MBDCND.EQ.5 .OR. MBDCND.EQ.7))
     2    IERROR = 18
      IF (IERROR.NE.0 .AND. IERROR.NE.9) RETURN
      NCK = N
      GO TO (101,103,102,103,101,103),NBDCND
  101 NCK = NCK-1
      GO TO 103
  102 NCK = NCK+1
  103 L = 2
      K = 1
  104 L = L+L
      K = K+1
      IF (NCK-L) 105,105,104
  105 L = L+L
      NP1 = N+1
      MP1 = M+1
      I1 = (K-2)*L+K+MAX(2*N,6*M)+13
      I2 = I1+NP1
      I3 = I2+NP1
      I4 = I3+NP1
      I5 = I4+NP1
      I6 = I5+NP1
      I7 = I6+MP1
      I8 = I7+MP1
      I9 = I8+MP1
      I10 = I9+MP1
      W(1) = I10+M
      CALL HWSCS1 (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     1             BDRF,ELMBDA,F,IDIMF,PERTRB,W(2),W(I1),W(I2),W(I3),
     2             W(I4),W(I5),W(I6),W(I7),W(I8),W(I9),W(I10))
      RETURN
      END
