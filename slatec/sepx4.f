*DECK SEPX4
      SUBROUTINE SEPX4 (IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA,
     +   C, D, N, NBDCND, BDC, BDD, COFX, GRHS, USOL, IDMN, W, PERTRB,
     +   IERROR)
C***BEGIN PROLOGUE  SEPX4
C***PURPOSE  Solve for either the second or fourth order finite
C            difference approximation to the solution of a separable
C            elliptic partial differential equation on a rectangle.
C            Any combination of periodic or mixed boundary conditions is
C            allowed.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A2
C***TYPE      SINGLE PRECISION (SEPX4-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SEPARABLE
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C Purpose                SEPX4 solves for either the second-order
C                        finite difference approximation or a
C                        fourth-order approximation  to the
C                        solution of a separable elliptic equation
C                             AF(X)*UXX+BF(X)*UX+CF(X)*U+UYY = G(X,Y)
C
C                        on a rectangle (X greater than or equal to A
C                        and less than or equal to B; Y greater than
C                        or equal to C and less than or equal to D).
C                        Any combination of periodic or mixed boundary
C                        conditions is allowed.
C                        If boundary conditions in the X direction
C                        are periodic (see MBDCND=0 below) then the
C                        coefficients must satisfy
C                        AF(X)=C1,BF(X)=0,CF(X)=C2 for all X.
C                        Here C1,C2 are constants, C1.GT.0.
C
C                        The possible boundary conditions are
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
C                        In the Y-direction:
C                         (0) Periodic, U(X,Y+D-C)=U(X,Y) for all X,Y
C                         (1) U(X,C),U(X,D) are specified for all X
C                         (2) U(X,C),dU(X,D)/dY are specified for all X
C                         (3) dU(X,C)/DY,dU(X,D)/dY are specified for
C                            all X
C                        (4) dU(X,C)/DY,U(X,D) are specified for all X
C
C Usage                  Call SEPX4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,
C                                  BETA,C,D,N,NBDCND,BDC,BDD,COFX,
C                                  GRHS,USOL,IDMN,W,PERTRB,IERROR)
C
C Arguments
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
C                          When MBDCND has any other value, BDA is a
C                          dummy parameter.
C
C On Input               ALPHA
C                          The scalar multiplying the solution in case
C                          of a mixed boundary condition AT X=A (see
C                          argument BDA).  If MBDCND = 3,4 then ALPHA is
C                          a dummy parameter.
C
C                        BDB
C                          A one-dimensional array of length N+1 that
C                          specifies the values of dU(B,Y)/dX+
C                          BETA*U(B,Y) at X=B.  when MBDCND=2 or 3
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
C                              i.e., dU(X,C)/dY and U(X,D)
C                              are specified for all X
C                          = 3 If the boundary conditions are mixed at
C                              Y= C and Y=D i.e., dU(X,D)/DY
C                              and dU(X,D)/dY are specified
C                              for all X
C                          = 4 If the boundary condition is mixed at Y=C
C                              and the solution is specified at Y=D;
C                              i.e. dU(X,C)/dY+GAMA*U(X,C) and U(X,D)
C                              are specified for all X
C
C                        BDC
C                          A one-dimensional array of length M+1 that
C                          specifies the value dU(X,C)/DY
C                          at Y=C.  When NBDCND=3 or 4
C                            BDC(I) = dU(XI,C)/DY
C                             I=1,2,...,M+1.
C                          When NBDCND has any other value, BDC is a
C                          dummy parameter.
C
C
C                        BDD
C                          A one-dimensional array of length M+1 that
C                          specifies the value of dU(X,D)/DY
C                          at Y=D.  When NBDCND=2 or 3
C                            BDD(I)=dU(XI,D)/DY
C                             I=1,2,...,M+1.
C                          When NBDCND has any other value, BDD is a
C                          dummy parameter.
C
C
C                        COFX
C                          A user-supplied subprogram with
C                          parameters X, AFUN, BFUN, CFUN which
C                          returns the values of the X-dependent
C                          coefficients AF(X), BF(X), CF(X) in
C                          the elliptic equation at X.
C                          If boundary conditions in the X direction
C                          are periodic then the coefficients
C                          must satisfy AF(X)=C1,BF(X)=0,CF(X)=C2 for
C                          all X.  Here C1.GT.0 and C2 are constants.
C
C                          Note that COFX must be declared external
C                          in the calling routine.
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
C                          calling SEPX4.  This parameter is used to
C                          specify the variable dimension of GRHS and
C                          USOL.  IDMN must be at least 7 and greater
C                          than or equal to M+1.
C
C                        W
C                          A one-dimensional array that must be provided
C                          by the user for work space.
C                          10*N+(16+INT(log2(N)))*(M+1)+23 will suffice
C                          as a length for W.  The actual length of
C                          W in the calling routine must be set in W(1)
C                          (see IERROR=11).
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
C                          W(1) contains the exact minimal length (in
C                          floating point) required for the work space
C                          (see IERROR=11).
C
C                        PERTRB
C                          If a combination of periodic or derivative
C                          boundary conditions (i.e., ALPHA=BETA=0 if
C                          MBDCND=3) is specified and if CF(X)=0 for all
C                          X, then a solution to the discretized matrix
C                          equation may not exist (reflecting the non-
C                          uniqueness of solutions to the PDE).  PERTRB
C                          is a constant calculated and subtracted from
C                          the right hand side of the matrix equation
C                          insuring the existence of a solution.
C                          SEPX4 computes this solution which is a
C                          weighted minimal least squares solution to
C                          the original problem.  If singularity is
C                          not detected PERTRB=0.0 is returned by
C                          SEPX4.
C
C                        IERROR
C                          An error flag that indicates invalid input
C                          parameters or failure to find a solution
C                          = 0  No error
C                          = 1  If A greater than B or C greater than D
C                          = 2  If MBDCND less than 0 or MBDCND greater
C                               than 4
C                          = 3  If NBDCND less than 0 or NBDCND greater
C                               than 4
C                          = 4  If attempt to find a solution fails.
C                               (the linear system generated is not
C                               diagonally dominant.)
C                          = 5  If IDMN is too small (see discussion of
C                               IDMN)
C                          = 6  If M is too small or too large (see
C                               discussion of M)
C                          = 7  If N is too small (see discussion of N)
C                          = 8  If IORDER is not 2 or 4
C                          = 10 If AFUN is less than or equal to zero
C                               for some interior mesh point XI
C                          = 11 If the work space length input in W(1)
C                               is less than the exact minimal work
C                               space length required output in W(1).
C                          = 12 If MBDCND=0 and AF(X)=CF(X)=constant
C                               or BF(X)=0 for all X is not true.
C
C *Long Description:
C
C Dimension of           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
C Arguments              USOL(IDMN,N+1), GRHS(IDMN,N+1),
C                        W (see argument list)
C
C Latest Revision        October 1980
C
C Special Conditions     NONE
C
C Common Blocks          SPL4
C
C I/O                    NONE
C
C Precision              Single
C
C Required Library       NONE
C Files
C
C Specialist             John C. Adams, NCAR, Boulder, Colorado  80307
C
C Language               FORTRAN
C
C
C Entry Points           SEPX4,SPELI4,CHKPR4,CHKSN4,ORTHO4,MINSO4,TRIS4,
C                        DEFE4,DX4,DY4
C
C History                SEPX4 was developed by modifying the ULIB
C                        routine SEPELI during October 1978.
C                        It should be used instead of SEPELI whenever
C                        possible.  The increase in speed is at least
C                        a factor of three.
C
C Algorithm              SEPX4 automatically discretizes the separable
C                        elliptic equation which is then solved by a
C                        generalized cyclic reduction algorithm in the
C                        subroutine POIS.  The fourth order solution
C                        is obtained using the technique of
C                        deferred corrections referenced below.
C
C
C References             Keller, H.B., 'Numerical Methods for Two-point
C                          Boundary-value Problems', Blaisdel (1968),
C                          Waltham, Mass.
C
C                        Swarztrauber, P., and R. Sweet (1975):
C                          'Efficient FORTRAN Subprograms For The
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
C***ROUTINES CALLED  CHKPR4, SPELI4
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920122  Minor corrections and modifications to prologue.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SEPX4
C
      DIMENSION       GRHS(IDMN,*)           ,USOL(IDMN,*)
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
      EXTERNAL COFX
C***FIRST EXECUTABLE STATEMENT  SEPX4
      CALL CHKPR4(IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,IDMN,IERROR)
      IF (IERROR .NE. 0) RETURN
C
C     COMPUTE MINIMUM WORK SPACE AND CHECK WORK SPACE LENGTH INPUT
C
      L = N+1
      IF (NBDCND .EQ. 0) L = N
      K = M+1
      L = N+1
C     ESTIMATE LOG BASE 2 OF N
      LOG2N=INT(LOG(REAL(N+1))/LOG(2.0)+0.5)
      LENGTH=4*(N+1)+(10+LOG2N)*(M+1)
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
      CALL SPELI4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     1NBDCND,BDC,BDD,COFX,W(I1),W(I2),W(I3),
     2             W(I4),W(I5),W(I6),W(I7),W(I8),W(I9),W(I10),W(I11),
     3             W(I12),GRHS,USOL,IDMN,W(I13),PERTRB,IERROR)
      RETURN
      END
