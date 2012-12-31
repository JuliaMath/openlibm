*DECK SEPELI
      SUBROUTINE SEPELI (INTL, IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB,
     +   BETA, C, D, N, NBDCND, BDC, GAMA, BDD, XNU, COFX, COFY, GRHS,
     +   USOL, IDMN, W, PERTRB, IERROR)
C***BEGIN PROLOGUE  SEPELI
C***PURPOSE  Discretize and solve a second and, optionally, a fourth
C            order finite difference approximation on a uniform grid to
C            the general separable elliptic partial differential
C            equation on a rectangle with any combination of periodic or
C            mixed boundary conditions.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A2
C***TYPE      SINGLE PRECISION (SEPELI-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SEPARABLE
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C Dimension of           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
C Arguments              USOL(IDMN,N+1), GRHS(IDMN,N+1),
C                        W (see argument list)
C
C Latest Revision        March 1977
C
C Purpose                SEPELI solves for either the second-order
C                        finite difference approximation or a
C                        fourth-order approximation to a separable
C                        elliptic equation.
C
C                                    2    2
C                             AF(X)*d U/dX + BF(X)*dU/dX  + CF(X)*U +
C                                    2    2
C                             DF(Y)*d U/dY  + EF(Y)*dU/dY + FF(Y)*U
C
C                             = G(X,Y)
C
C                        on a rectangle (X greater than or equal to A
C                        and less than or equal to B; Y greater than
C                        or equal to C and less than or equal to D).
C                        Any combination of periodic or mixed boundary
C                        conditions is allowed.
C
C Purpose                The possible boundary conditions are:
C                        in the X-direction:
C                         (0) Periodic, U(X+B-A,Y)=U(X,Y) for all Y,X
C                         (1) U(A,Y), U(B,Y) are specified for all Y
C                         (2) U(A,Y), dU(B,Y)/dX+BETA*U(B,Y) are
C                             specified for all Y
C                         (3) dU(A,Y)/dX+ALPHA*U(A,Y),dU(B,Y)/dX+
C                             BETA*U(B,Y) are specified for all Y
C                         (4) dU(A,Y)/dX+ALPHA*U(A,Y),U(B,Y) are
C                             specified for all Y
C
C                        in the Y-direction:
C                         (0) Periodic, U(X,Y+D-C)=U(X,Y) for all X,Y
C                         (1) U(X,C),U(X,D) are specified for all X
C                         (2) U(X,C),dU(X,D)/dY+XNU*U(X,D) are specified
C                             for all X
C                         (3) dU(X,C)/dY+GAMA*U(X,C),dU(X,D)/dY+
C                             XNU*U(X,D) are specified for all X
C                         (4) dU(X,C)/dY+GAMA*U(X,C),U(X,D) are
C                             specified for all X
C
C Arguments
C
C On Input               INTL
C                          = 0 On initial entry to SEPELI or if any of
C                              the arguments C, D, N, NBDCND, COFY are
C                              changed from a previous call
C                          = 1 If C, D, N, NBDCND, COFY are unchanged
C                              from the previous call.
C
C                        IORDER
C                          = 2 If a second-order approximation is sought
C                          = 4 If a fourth-order approximation is sought
C
C                        A,B
C                          The range of the X-independent variable;
C                          i.e., X is greater than or equal to A and
C                          less than or equal to B.  A must be less than
C                          B.
C
C                        M
C                          The number of panels into which the interval
C                          [A,B] is subdivided.  Hence, there will be
C                          M+1 grid points in the X-direction given by
C                          XI=A+(I-1)*DLX for I=1,2,...,M+1 where
C                          DLX=(B-A)/M is the panel width.  M must be
C                          less than IDMN and greater than 5.
C
C                        MBDCND
C                          Indicates the type of boundary condition at
C                          X=A and X=B
C                          = 0 If the solution is periodic in X; i.e.,
C                              U(X+B-A,Y)=U(X,Y) for all Y,X
C                          = 1 If the solution is specified at X=A and
C                              X=B; i.e., U(A,Y) and U(B,Y) are
C                              specified for all Y
C                          = 2 If the solution is specified at X=A and
C                              the boundary condition is mixed at X=B;
C                              i.e., U(A,Y) and dU(B,Y)/dX+BETA*U(B,Y)
C                              are specified for all Y
C                          = 3 If the boundary conditions at X=A and X=B
C                              are mixed; i.e., dU(A,Y)/dX+ALPHA*U(A,Y)
C                              and dU(B,Y)/dX+BETA*U(B,Y) are specified
C                              for all Y
C                          = 4 If the boundary condition at X=A is mixed
C                              and the solution is specified at X=B;
C                              i.e., dU(A,Y)/dX+ALPHA*U(A,Y) and U(B,Y)
C                              are specified for all Y
C
C                        BDA
C                          A one-dimensional array of length N+1 that
C                          specifies the values of dU(A,Y)/dX+
C                          ALPHA*U(A,Y) at X=A, when MBDCND=3 or 4.
C                               BDA(J) = dU(A,YJ)/dX+ALPHA*U(A,YJ);
C                               J=1,2,...,N+1
C                          when MBDCND has any other value, BDA is a
C                          dummy parameter.
C
C On Input               ALPHA
C                          The scalar multiplying the solution in case
C                          of a mixed boundary condition at X=A (see
C                          argument BDA).  If MBDCND = 3,4 then ALPHA is
C                          a dummy parameter.
C
C                        BDB
C                          A one-dimensional array of length N+1 that
C                          specifies the values of dU(B,Y)/dX+
C                          BETA*U(B,Y) at X=B.  When MBDCND=2 or 3
C                               BDB(J) = dU(B,YJ)/dX+BETA*U(B,YJ);
C                               J=1,2,...,N+1
C                          When MBDCND has any other value, BDB is a
C                          dummy parameter.
C
C                        BETA
C                          The scalar multiplying the solution in case
C                          of a mixed boundary condition at X=B (see
C                          argument BDB).  If MBDCND=2,3 then BETA is a
C                          dummy parameter.
C
C                        C,D
C                          The range of the Y-independent variable;
C                          i.e., Y is greater than or equal to C and
C                          less than or equal to D.  C must be less than
C                          D.
C
C                        N
C                          The number of panels into which the interval
C                          [C,D] is subdivided.  Hence, there will be
C                          N+1 grid points in the Y-direction given by
C                          YJ=C+(J-1)*DLY for J=1,2,...,N+1 where
C                          DLY=(D-C)/N is the panel width.  In addition,
C                          N must be greater than 4.
C
C                        NBDCND
C                          Indicates the types of boundary conditions at
C                          Y=C and Y=D
C                          = 0 If the solution is periodic in Y; i.e.,
C                              U(X,Y+D-C)=U(X,Y) for all X,Y
C                          = 1 If the solution is specified at Y=C and
C                              Y = D, i.e., U(X,C) and U(X,D) are
C                              specified for all X
C                          = 2 If the solution is specified at Y=C and
C                              the boundary condition is mixed at Y=D;
C                              i.e., U(X,C) and dU(X,D)/dY+XNU*U(X,D)
C                              are specified for all X
C                          = 3 If the boundary conditions are mixed at
C                              Y=C and Y=D; i.e., dU(X,D)/dY+GAMA*U(X,C)
C                              and dU(X,D)/dY+XNU*U(X,D) are specified
C                              for all X
C                          = 4 If the boundary condition is mixed at Y=C
C                              and the solution is specified at Y=D;
C                              i.e. dU(X,C)/dY+GAMA*U(X,C) and U(X,D)
C                              are specified for all X
C
C                        BDC
C                          A one-dimensional array of length M+1 that
C                          specifies the value of dU(X,C)/dY+GAMA*U(X,C)
C                          at Y=C.  When NBDCND=3 or 4
C                             BDC(I) = dU(XI,C)/dY + GAMA*U(XI,C);
C                             I=1,2,...,M+1.
C                          When NBDCND has any other value, BDC is a
C                          dummy parameter.
C
C                        GAMA
C                          The scalar multiplying the solution in case
C                          of a mixed boundary condition at Y=C (see
C                          argument BDC).  If NBDCND=3,4 then GAMA is a
C                          dummy parameter.
C
C                        BDD
C                          A one-dimensional array of length M+1 that
C                          specifies the value of dU(X,D)/dY +
C                          XNU*U(X,D) at Y=C.  When NBDCND=2 or 3
C                            BDD(I) = dU(XI,D)/dY + XNU*U(XI,D);
C                            I=1,2,...,M+1.
C                          When NBDCND has any other value, BDD is a
C                          dummy parameter.
C
C                        XNU
C                          The scalar multiplying the solution in case
C                          of a mixed boundary condition at Y=D (see
C                          argument BDD).  If NBDCND=2 or 3 then XNU is
C                          a dummy parameter.
C
C                        COFX
C                          A user-supplied subprogram with
C                          parameters X, AFUN, BFUN, CFUN which
C                          returns the values of the X-dependent
C                          coefficients AF(X), BF(X), CF(X) in
C                          the elliptic equation at X.
C
C                        COFY
C                          A user-supplied subprogram with
C                          parameters Y, DFUN, EFUN, FFUN which
C                          returns the values of the Y-dependent
C                          coefficients DF(Y), EF(Y), FF(Y) in
C                          the elliptic equation at Y.
C
C                        NOTE:  COFX and COFY must be declared external
C                        in the calling routine.  The values returned in
C                        AFUN and DFUN must satisfy AFUN*DFUN greater
C                        than 0 for A less than X less than B,
C                        C less than Y less than D (see IERROR=10).
C                        The coefficients provided may lead to a matrix
C                        equation which is not diagonally dominant in
C                        which case solution may fail (see IERROR=4).
C
C                        GRHS
C                          A two-dimensional array that specifies the
C                          values of the right-hand side of the elliptic
C                          equation; i.e., GRHS(I,J)=G(XI,YI), for
C                          I=2,...,M; J=2,...,N.  At the boundaries,
C                          GRHS is defined by
C
C                          MBDCND   GRHS(1,J)   GRHS(M+1,J)
C                          ------   ---------   -----------
C                            0      G(A,YJ)     G(B,YJ)
C                            1         *           *
C                            2         *        G(B,YJ)  J=1,2,...,N+1
C                            3      G(A,YJ)     G(B,YJ)
C                            4      G(A,YJ)        *
C
C                          NBDCND   GRHS(I,1)   GRHS(I,N+1)
C                          ------   ---------   -----------
C                            0      G(XI,C)     G(XI,D)
C                            1         *           *
C                            2         *        G(XI,D)  I=1,2,...,M+1
C                            3      G(XI,C)     G(XI,D)
C                            4      G(XI,C)        *
C
C                          where * means these quantities are not used.
C                          GRHS should be dimensioned IDMN by at least
C                          N+1 in the calling routine.
C
C                        USOL
C                          A two-dimensional array that specifies the
C                          values of the solution along the boundaries.
C                          At the boundaries, USOL is defined by
C
C                          MBDCND   USOL(1,J)   USOL(M+1,J)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(A,YJ)     U(B,YJ)
C                            2      U(A,YJ)        *     J=1,2,...,N+1
C                            3         *           *
C                            4         *        U(B,YJ)
C
C                          NBDCND   USOL(I,1)   USOL(I,N+1)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(XI,C)     U(XI,D)
C                            2      U(XI,C)        *     I=1,2,...,M+1
C                            3         *           *
C                            4         *        U(XI,D)
C
C                          where * means the quantities are not used in
C                          the solution.
C
C                          If IORDER=2, the user may equivalence GRHS
C                          and USOL to save space.  Note that in this
C                          case the tables specifying the boundaries of
C                          the GRHS and USOL arrays determine the
C                          boundaries uniquely except at the corners.
C                          If the tables call for both G(X,Y) and
C                          U(X,Y) at a corner then the solution must be
C                          chosen.  For example, if MBDCND=2 and
C                          NBDCND=4, then U(A,C), U(A,D), U(B,D) must be
C                          chosen at the corners in addition to G(B,C).
C
C                          If IORDER=4, then the two arrays, USOL and
C                          GRHS, must be distinct.
C
C                          USOL should be dimensioned IDMN by at least
C                          N+1 in the calling routine.
C
C                        IDMN
C                          The row (or first) dimension of the arrays
C                          GRHS and USOL as it appears in the program
C                          calling SEPELI.  This parameter is used to
C                          specify the variable dimension of GRHS and
C                          USOL.  IDMN must be at least 7 and greater
C                          than or equal to M+1.
C
C                        W
C                          A one-dimensional array that must be provided
C                          by the user for work space.  Let
C                          K=INT(log2(N+1))+1 and set  L=2**(K+1).
C                          then (K-2)*L+K+10*N+12*M+27 will suffice
C                          as a length of W.  THE actual length of W in
C                          the calling routine must be set in W(1) (see
C                          IERROR=11).
C
C On Output              USOL
C                          Contains the approximate solution to the
C                          elliptic equation.  USOL(I,J) is the
C                          approximation to U(XI,YJ) for I=1,2...,M+1
C                          and J=1,2,...,N+1.  The approximation has
C                          error O(DLX**2+DLY**2) if called with
C                          IORDER=2 and O(DLX**4+DLY**4) if called with
C                          IORDER=4.
C
C                        W
C                          Contains intermediate values that must not be
C                          destroyed if SEPELI is called again with
C                          INTL=1.  In addition W(1) contains the exact
C                          minimal length (in floating point) required
C                          for the work space (see IERROR=11).
C
C                        PERTRB
C                          If a combination of periodic or derivative
C                          boundary conditions (i.e., ALPHA=BETA=0 if
C                          MBDCND=3; GAMA=XNU=0 if NBDCND=3) is
C                          specified and if the coefficients of U(X,Y)
C                          in the separable elliptic equation are zero
C                          (i.e., CF(X)=0 for X greater than or equal to
C                          A and less than or equal to B; FF(Y)=0 for
C                          Y greater than or equal to C and less than
C                          or equal to D) then a solution may not exist.
C                          PERTRB is a constant calculated and
C                          subtracted from the right-hand side of the
C                          matrix equations generated by SEPELI which
C                          insures that a solution exists.  SEPELI then
C                          computes this solution which is a weighted
C                          minimal least squares solution to the
C                          original problem.
C
C                        IERROR
C                          An error flag that indicates invalid input
C                          parameters or failure to find a solution
C                          = 0 No error
C                          = 1 If A greater than B or C greater than D
C                          = 2 If MBDCND less than 0 or MBDCND greater
C                              than 4
C                          = 3 If NBDCND less than 0 or NBDCND greater
C                              than 4
C                          = 4 If attempt to find a solution fails.
C                              (the linear system generated is not
C                              diagonally dominant.)
C                          = 5 If IDMN is too small (see discussion of
C                              IDMN)
C                          = 6 If M is too small or too large (see
C                              discussion of M)
C                          = 7 If N is too small (see discussion of N)
C                          = 8 If IORDER is not 2 or 4
C                          = 9 If INTL is not 0 or 1
C                          = 10 If AFUN*DFUN less than or equal to 0 for
C                               some interior mesh point (XI,YJ)
C                          = 11 If the work space length input in W(1)
C                               is less than the exact minimal work
C                               space length required output in W(1).
C
C                          NOTE (concerning IERROR=4):  for the
C                          coefficients input through COFX, COFY, the
C                          discretization may lead to a block
C                          tridiagonal linear system which is not
C                          diagonally dominant (for example, this
C                          happens if CFUN=0 and BFUN/(2.*DLX) greater
C                          than AFUN/DLX**2).  In this case solution may
C                          fail.  This cannot happen in the limit as
C                          DLX, DLY approach zero.  Hence, the condition
C                          may be remedied by taking larger values for M
C                          or N.
C
C Entry Points           SEPELI, SPELIP, CHKPRM, CHKSNG, ORTHOG, MINSOL,
C                        TRISP, DEFER, DX, DY, BLKTRI, BLKTR1, INDXB,
C                        INDXA, INDXC, PROD, PRODP, CPROD, CPRODP,
C                        PPADD, PSGF, BSRH, PPSGF, PPSPF, COMPB,
C                        TRUN1, STOR1, TQLRAT
C
C Special Conditions     NONE
C
C Common Blocks          SPLP, CBLKT
C
C I/O                    NONE
C
C Precision              Single
C
C Specialist             John C. Adams, NCAR, Boulder, Colorado  80307
C
C Language               FORTRAN
C
C History                Developed at NCAR during 1975-76.
C
C Algorithm              SEPELI automatically discretizes the separable
C                        elliptic equation which is then solved by a
C                        generalized cyclic reduction algorithm in the
C                        subroutine, BLKTRI.  The fourth-order solution
C                        is obtained using 'Deferred Corrections' which
C                        is described and referenced in sections,
C                        references and method.
C
C Space Required         14654 (octal) = 6572 (decimal)
C
C Accuracy and Timing    The following computational results were
C                        obtained by solving the sample problem at the
C                        end of this write-up on the Control Data 7600.
C                        The op count is proportional to M*N*log2(N).
C                        In contrast to the other routines in this
C                        chapter, accuracy is tested by computing and
C                        tabulating second- and fourth-order
C                        discretization errors.  Below is a table
C                        containing computational results.  The times
C                        given do not include initialization (i.e.,
C                        times are for INTL=1).  Note that the
C                        fourth-order accuracy is not realized until the
C                        mesh is sufficiently refined.
C
C              Second-order    Fourth-order   Second-order  Fourth-order
C    M    N   Execution Time  Execution Time    Error         Error
C               (M SEC)         (M SEC)
C     6    6         6              14          6.8E-1        1.2E0
C    14   14        23              58          1.4E-1        1.8E-1
C    30   30       100             247          3.2E-2        9.7E-3
C    62   62       445           1,091          7.5E-3        3.0E-4
C   126  126     2,002           4,772          1.8E-3        3.5E-6
C
C Portability            There are no machine-dependent constants.
C
C Required Resident      SQRT, ABS, LOG
C Routines
C
C References             Keller, H.B., 'Numerical Methods for Two-point
C                          Boundary-value Problems', Blaisdel (1968),
C                          Waltham, Mass.
C
C                        Swarztrauber, P., and R. Sweet (1975):
C                          'Efficient FORTRAN Subprograms for The
C                          Solution of Elliptic Partial Differential
C                          Equations'.  NCAR Technical Note
C                          NCAR-TN/IA-109, pp. 135-137.
C
C***REFERENCES  H. B. Keller, Numerical Methods for Two-point
C                 Boundary-value Problems, Blaisdel, Waltham, Mass.,
C                 1968.
C               P. N. Swarztrauber and R. Sweet, Efficient Fortran
C                 subprograms for the solution of elliptic equations,
C                 NCAR TN/IA-109, July 1975, 138 pp.
C***ROUTINES CALLED  CHKPRM, SPELIP
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SEPELI
C
      DIMENSION       GRHS(IDMN,*)           ,USOL(IDMN,*)
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
      EXTERNAL        COFX       ,COFY
C***FIRST EXECUTABLE STATEMENT  SEPELI
      CALL CHKPRM (INTL,IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,COFY,
     1             IDMN,IERROR)
      IF (IERROR .NE. 0) RETURN
C
C     COMPUTE MINIMUM WORK SPACE AND CHECK WORK SPACE LENGTH INPUT
C
      L = N+1
      IF (NBDCND .EQ. 0) L = N
      LOGB2N = INT(LOG(L+0.5)/LOG(2.0))+1
      LL = 2**(LOGB2N+1)
      K = M+1
      L = N+1
      LENGTH = (LOGB2N-2)*LL+LOGB2N+MAX(2*L,6*K)+5
      IF (NBDCND .EQ. 0) LENGTH = LENGTH+2*L
      IERROR = 11
      LINPUT = INT(W(1)+0.5)
      LOUTPT = LENGTH+6*(K+L)+1
      W(1) = LOUTPT
      IF (LOUTPT .GT. LINPUT) RETURN
      IERROR = 0
C
C     SET WORK SPACE INDICES
C
      I1 = LENGTH+2
      I2 = I1+L
      I3 = I2+L
      I4 = I3+L
      I5 = I4+L
      I6 = I5+L
      I7 = I6+L
      I8 = I7+K
      I9 = I8+K
      I10 = I9+K
      I11 = I10+K
      I12 = I11+K
      I13 = 2
      CALL SPELIP (INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     1             NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,W(I1),W(I2),W(I3),
     2             W(I4),W(I5),W(I6),W(I7),W(I8),W(I9),W(I10),W(I11),
     3             W(I12),GRHS,USOL,IDMN,W(I13),PERTRB,IERROR)
      RETURN
      END
