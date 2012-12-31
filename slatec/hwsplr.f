*DECK HWSPLR
      SUBROUTINE HWSPLR (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND,
     +   BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HWSPLR
C***PURPOSE  Solve a finite difference approximation to the Helmholtz
C            equation in polar coordinates.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HWSPLR-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, POLAR
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine HWSPLR solves a finite difference approximation to the
C     Helmholtz equation in polar coordinates:
C
C          (1/R)(d/dR)(R(dU/dR)) + (1/R**2)(d/dTHETA)(dU/dTHETA)
C
C                                + LAMBDA*U = F(R,THETA).
C
C
C
C
C     * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C             * * * * * *   On Input    * * * * * *
C
C     A,B
C       The range of R, i.e., A .LE. R .LE. B.  A must be less than B
C       and A must be non-negative.
C
C     M
C       The number of panels into which the interval (A,B) is
C       subdivided.  Hence, there will be M+1 grid points in the
C       R-direction given by R(I) = A+(I-1)DR, for I = 1,2,...,M+1,
C       where DR = (B-A)/M is the panel width. M must be greater than 3.
C
C     MBDCND
C       Indicates the type of boundary condition at R = A and R = B.
C
C       = 1  If the solution is specified at R = A and R = B.
C       = 2  If the solution is specified at R = A and the derivative of
C            the solution with respect to R is specified at R = B.
C       = 3  If the derivative of the solution with respect to R is
C            specified at R = A (see note below) and R = B.
C       = 4  If the derivative of the solution with respect to R is
C            specified at R = A (see note below) and the solution is
C            specified at R = B.
C       = 5  If the solution is unspecified at R = A = 0 and the
C            solution is specified at R = B.
C       = 6  If the solution is unspecified at R = A = 0 and the
C            derivative of the solution with respect to R is specified
C            at R = B.
C
C       NOTE:  If A = 0, do not use MBDCND = 3 or 4, but instead use
C              MBDCND = 1,2,5, or 6  .
C
C     BDA
C       A one-dimensional array of length N+1 that specifies the values
C       of the derivative of the solution with respect to R at R = A.
C       When MBDCND = 3 or 4,
C
C            BDA(J) = (d/dR)U(A,THETA(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value, BDA is a dummy variable.
C
C     BDB
C       A one-dimensional array of length N+1 that specifies the values
C       of the derivative of the solution with respect to R at R = B.
C       When MBDCND = 2,3, or 6,
C
C            BDB(J) = (d/dR)U(B,THETA(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value, BDB is a dummy variable.
C
C     C,D
C       The range of THETA, i.e., C .LE. THETA .LE. D.  C must be less
C       than D.
C
C     N
C       The number of panels into which the interval (C,D) is
C       subdivided.  Hence, there will be N+1 grid points in the
C       THETA-direction given by THETA(J) = C+(J-1)DTHETA for
C       J = 1,2,...,N+1, where DTHETA = (D-C)/N is the panel width.  N
C       must be greater than 3.
C
C     NBDCND
C       Indicates the type of boundary conditions at THETA = C and
C       at THETA = D.
C
C       = 0  If the solution is periodic in THETA, i.e.,
C            U(I,J) = U(I,N+J).
C       = 1  If the solution is specified at THETA = C and THETA = D
C            (see note below).
C       = 2  If the solution is specified at THETA = C and the
C            derivative of the solution with respect to THETA is
C            specified at THETA = D (see note below).
C       = 4  If the derivative of the solution with respect to THETA is
C            specified at THETA = C and the solution is specified at
C            THETA = D (see note below).
C
C       NOTE:  When NBDCND = 1,2, or 4, do not use MBDCND = 5 or 6
C              (the former indicates that the solution is specified at
C              R = 0, the latter indicates the solution is unspecified
C              at R = 0).  Use instead MBDCND = 1 or 2  .
C
C     BDC
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to THETA at
C       THETA = C.  When NBDCND = 3 or 4,
C
C            BDC(I) = (d/dTHETA)U(R(I),C), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDC is a dummy variable.
C
C     BDD
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to THETA at
C       THETA = D.  When NBDCND = 2 or 3,
C
C            BDD(I) = (d/dTHETA)U(R(I),D), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDD is a dummy variable.
C
C     ELMBDA
C       The constant LAMBDA in the Helmholtz equation.  If
C       LAMBDA .LT. 0, a solution may not exist.  However, HWSPLR will
C       attempt to find a solution.
C
C     F
C       A two-dimensional array that specifies the values of the right
C       side of the Helmholtz equation and boundary values (if any).
C       For I = 2,3,...,M and J = 2,3,...,N
C
C            F(I,J) = F(R(I),THETA(J)).
C
C       On the boundaries F is defined by
C
C            MBDCND   F(1,J)            F(M+1,J)
C            ------   -------------     -------------
C
C              1      U(A,THETA(J))     U(B,THETA(J))
C              2      U(A,THETA(J))     F(B,THETA(J))
C              3      F(A,THETA(J))     F(B,THETA(J))
C              4      F(A,THETA(J))     U(B,THETA(J))   J = 1,2,...,N+1
C              5      F(0,0)            U(B,THETA(J))
C              6      F(0,0)            F(B,THETA(J))
C
C            NBDCND   F(I,1)            F(I,N+1)
C            ------   ---------         ---------
C
C              0      F(R(I),C)         F(R(I),C)
C              1      U(R(I),C)         U(R(I),D)
C              2      U(R(I),C)         F(R(I),D)   I = 1,2,...,M+1
C              3      F(R(I),C)         F(R(I),D)
C              4      F(R(I),C)         U(R(I),D)
C
C       F must be dimensioned at least (M+1)*(N+1).
C
C       NOTE
C
C       If the table calls for both the solution U and the right side F
C       at a corner then the solution must be specified.
C
C
C     IDIMF
C       The row (or first) dimension of the array F as it appears in the
C       program calling HWSPLR.  This parameter is used to specify the
C       variable dimension of F.  IDIMF must be at least M+1  .
C
C     W
C       A one-dimensional array that must be provided by the user for
C       work space.  W may require up to 4*(N+1) +
C       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of
C       locations used is computed by HWSPLR and is returned in location
C       W(1).
C
C
C             * * * * * *   On Output     * * * * * *
C
C     F
C       Contains the solution U(I,J) of the finite difference
C       approximation for the grid point (R(I),THETA(J)),
C       I = 1,2,...,M+1, J = 1,2,...,N+1  .
C
C     PERTRB
C       If a combination of periodic, derivative, or unspecified
C       boundary conditions is specified for a Poisson equation
C       (LAMBDA = 0), a solution may not exist.  PERTRB is a constant,
C       calculated and subtracted from F, which ensures that a solution
C       exists.  HWSPLR then computes this solution, which is a least
C       squares solution to the original approximation.  This solution
C       plus any constant is also a solution.  Hence, the solution is
C       not unique.  PERTRB should be small compared to the right side.
C       Otherwise, a solution is obtained to an essentially different
C       problem.  This comparison should always be made to insure that a
C       meaningful solution has been obtained.
C
C     IERROR
C       An error flag that indicates invalid input parameters.  Except
C       for numbers 0 and 11, a solution is not attempted.
C
C       =  0  No error.
C       =  1  A .LT. 0  .
C       =  2  A .GE. B.
C       =  3  MBDCND .LT. 1 or MBDCND .GT. 6  .
C       =  4  C .GE. D.
C       =  5  N .LE. 3
C       =  6  NBDCND .LT. 0 or .GT. 4  .
C       =  7  A = 0, MBDCND = 3 or 4  .
C       =  8  A .GT. 0, MBDCND .GE. 5  .
C       =  9  MBDCND .GE. 5, NBDCND .NE. 0 and NBDCND .NE. 3  .
C       = 10  IDIMF .LT. M+1  .
C       = 11  LAMBDA .GT. 0  .
C       = 12  M .LE. 3
C
C       Since this is the only means of indicating a possibly incorrect
C       call to HWSPLR, the user should test IERROR after the call.
C
C     W
C       W(1) contains the required length of W.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1),
C     Arguments      W(see argument list)
C
C     Latest         June 1, 1976
C     Revision
C
C     Subprograms    HWSPLR,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
C     Required       TRIX,TRI3,PIMACH
C
C     Special        None
C     Conditions
C
C     Common         NONE
C     Blocks
C
C     I/O
C
C     Precision      Single
C
C     Specialist     Roland Sweet
C
C     Language       FORTRAN
C
C     History        Standardized April 1, 1973
C                    Revised January 1, 1976
C
C     Algorithm      The routine defines the finite difference
C                    equations, incorporates boundary data, and adjusts
C                    the right side of singular systems and then calls
C                    GENBUN to solve the system.
C
C     Space          13430(octal) = 5912(decimal)  locations on the NCAR
C     Required       Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine HWSPLR is roughly proportional
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
C                        32        1         0          31
C                        32        1         1          23
C                        32        3         3          36
C                        64        1         0         128
C                        64        1         1          96
C                        64        3         3         142
C
C     Portability    American National Standards Institute FORTRAN.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C     Required       COS
C     Resident
C     Routines
C
C     Reference      Swarztrauber, P. and R. Sweet, 'Efficient FORTRAN
C                    Subprograms For The Solution Of Elliptic Equations'
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
C***END PROLOGUE  HWSPLR
C
C
      DIMENSION       F(IDIMF,*)
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
C***FIRST EXECUTABLE STATEMENT  HWSPLR
      IERROR = 0
      IF (A .LT. 0.) IERROR = 1
      IF (A .GE. B) IERROR = 2
      IF (MBDCND.LE.0 .OR. MBDCND.GE.7) IERROR = 3
      IF (C .GE. D) IERROR = 4
      IF (N .LE. 3) IERROR = 5
      IF (NBDCND.LE.-1 .OR. NBDCND.GE.5) IERROR = 6
      IF (A.EQ.0. .AND. (MBDCND.EQ.3 .OR. MBDCND.EQ.4)) IERROR = 7
      IF (A.GT.0. .AND. MBDCND.GE.5) IERROR = 8
      IF (MBDCND.GE.5 .AND. NBDCND.NE.0 .AND. NBDCND.NE.3) IERROR = 9
      IF (IDIMF .LT. M+1) IERROR = 10
      IF (M .LE. 3) IERROR = 12
      IF (IERROR .NE. 0) RETURN
      MP1 = M+1
      DELTAR = (B-A)/M
      DLRBY2 = DELTAR/2.
      DLRSQ = DELTAR**2
      NP1 = N+1
      DELTHT = (D-C)/N
      DLTHSQ = DELTHT**2
      NP = NBDCND+1
C
C     DEFINE RANGE OF INDICES I AND J FOR UNKNOWNS U(I,J).
C
      MSTART = 2
      MSTOP = MP1
      GO TO (101,105,102,103,104,105),MBDCND
  101 MSTOP = M
      GO TO 105
  102 MSTART = 1
      GO TO 105
  103 MSTART = 1
  104 MSTOP = M
  105 MUNK = MSTOP-MSTART+1
      NSTART = 1
      NSTOP = N
      GO TO (109,106,107,108,109),NP
  106 NSTART = 2
      GO TO 109
  107 NSTART = 2
  108 NSTOP = NP1
  109 NUNK = NSTOP-NSTART+1
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = MUNK
      ID3 = ID2+MUNK
      ID4 = ID3+MUNK
      ID5 = ID4+MUNK
      ID6 = ID5+MUNK
      A1 = 2./DLRSQ
      IJ = 0
      IF (MBDCND.EQ.3 .OR. MBDCND.EQ.4) IJ = 1
      DO 110 I=1,MUNK
         R = A+(I-IJ)*DELTAR
         J = ID5+I
         W(J) = R
         J = ID6+I
         W(J) = 1./R**2
         W(I) = (R-DLRBY2)/(R*DLRSQ)
         J = ID3+I
         W(J) = (R+DLRBY2)/(R*DLRSQ)
         J = ID2+I
         W(J) = -A1+ELMBDA
  110 CONTINUE
      GO TO (114,111,112,113,114,111),MBDCND
  111 W(ID2) = A1
      GO TO 114
  112 W(ID2) = A1
  113 W(ID3+1) = A1
  114 CONTINUE
C
C     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
C
      GO TO (115,115,117,117,119,119),MBDCND
  115 A1 = W(1)
      DO 116 J=NSTART,NSTOP
         F(2,J) = F(2,J)-A1*F(1,J)
  116 CONTINUE
      GO TO 119
  117 A1 = 2.*DELTAR*W(1)
      DO 118 J=NSTART,NSTOP
         F(1,J) = F(1,J)+A1*BDA(J)
  118 CONTINUE
  119 GO TO (120,122,122,120,120,122),MBDCND
  120 A1 = W(ID4)
      DO 121 J=NSTART,NSTOP
         F(M,J) = F(M,J)-A1*F(MP1,J)
  121 CONTINUE
      GO TO 124
  122 A1 = 2.*DELTAR*W(ID4)
      DO 123 J=NSTART,NSTOP
         F(MP1,J) = F(MP1,J)-A1*BDB(J)
  123 CONTINUE
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  124 A1 = 1./DLTHSQ
      L = ID5-MSTART+1
      LP = ID6-MSTART+1
      GO TO (134,125,125,127,127),NP
  125 DO 126 I=MSTART,MSTOP
         J = I+LP
         F(I,2) = F(I,2)-A1*W(J)*F(I,1)
  126 CONTINUE
      GO TO 129
  127 A1 = 2./DELTHT
      DO 128 I=MSTART,MSTOP
         J = I+LP
         F(I,1) = F(I,1)+A1*W(J)*BDC(I)
  128 CONTINUE
  129 A1 = 1./DLTHSQ
      GO TO (134,130,132,132,130),NP
  130 DO 131 I=MSTART,MSTOP
         J = I+LP
         F(I,N) = F(I,N)-A1*W(J)*F(I,NP1)
  131 CONTINUE
      GO TO 134
  132 A1 = 2./DELTHT
      DO 133 I=MSTART,MSTOP
         J = I+LP
         F(I,NP1) = F(I,NP1)-A1*W(J)*BDD(I)
  133 CONTINUE
  134 CONTINUE
C
C     ADJUST RIGHT SIDE OF EQUATION FOR UNKNOWN AT POLE WHEN HAVE
C     DERIVATIVE SPECIFIED BOUNDARY CONDITIONS.
C
      IF (MBDCND.GE.5 .AND. NBDCND.EQ.3)
     1    F(1,1) = F(1,1)-(BDD(2)-BDC(2))*4./(N*DELTHT*DLRSQ)
C
C     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
C     SOLUTION.
C
      PERTRB = 0.
      IF (ELMBDA) 144,136,135
  135 IERROR = 11
      GO TO 144
  136 IF (NBDCND.NE.0 .AND. NBDCND.NE.3) GO TO 144
      S2 = 0.
      GO TO (144,144,137,144,144,138),MBDCND
  137 W(ID5+1) = .5*(W(ID5+2)-DLRBY2)
      S2 = .25*DELTAR
  138 A2 = 2.
      IF (NBDCND .EQ. 0) A2 = 1.
      J = ID5+MUNK
      W(J) = .5*(W(J-1)+DLRBY2)
      S = 0.
      DO 140 I=MSTART,MSTOP
         S1 = 0.
         IJ = NSTART+1
         K = NSTOP-1
         DO 139 J=IJ,K
            S1 = S1+F(I,J)
  139    CONTINUE
         J = I+L
         S = S+(A2*S1+F(I,NSTART)+F(I,NSTOP))*W(J)
  140 CONTINUE
      S2 = M*A+DELTAR*((M-1)*(M+1)*.5+.25)+S2
      S1 = (2.+A2*(NUNK-2))*S2
      IF (MBDCND .EQ. 3) GO TO 141
      S2 = N*A2*DELTAR/8.
      S = S+F(1,1)*S2
      S1 = S1+S2
  141 CONTINUE
      PERTRB = S/S1
      DO 143 I=MSTART,MSTOP
         DO 142 J=NSTART,NSTOP
            F(I,J) = F(I,J)-PERTRB
  142    CONTINUE
  143 CONTINUE
  144 CONTINUE
C
C     MULTIPLY I-TH EQUATION THROUGH BY (R(I)*DELTHT)**2.
C
      DO 146 I=MSTART,MSTOP
         K = I-MSTART+1
         J = I+LP
         A1 = DLTHSQ/W(J)
         W(K) = A1*W(K)
         J = ID2+K
         W(J) = A1*W(J)
         J = ID3+K
         W(J) = A1*W(J)
         DO 145 J=NSTART,NSTOP
            F(I,J) = A1*F(I,J)
  145    CONTINUE
  146 CONTINUE
      W(1) = 0.
      W(ID4) = 0.
C
C     CALL GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
C
      CALL GENBUN (NBDCND,NUNK,1,MUNK,W(1),W(ID2+1),W(ID3+1),IDIMF,
     1             F(MSTART,NSTART),IERR1,W(ID4+1))
      IWSTOR = W(ID4+1)+3*MUNK
      GO TO (157,157,157,157,148,147),MBDCND
C
C     ADJUST THE SOLUTION AS NECESSARY FOR THE PROBLEMS WHERE A = 0.
C
  147 IF (ELMBDA .NE. 0.) GO TO 148
      YPOLE = 0.
      GO TO 155
  148 CONTINUE
      J = ID5+MUNK
      W(J) = W(ID2)/W(ID3)
      DO 149 IP=3,MUNK
         I = MUNK-IP+2
         J = ID5+I
         LP = ID2+I
         K = ID3+I
         W(J) = W(I)/(W(LP)-W(K)*W(J+1))
  149 CONTINUE
      W(ID5+1) = -.5*DLTHSQ/(W(ID2+1)-W(ID3+1)*W(ID5+2))
      DO 150 I=2,MUNK
         J = ID5+I
         W(J) = -W(J)*W(J-1)
  150 CONTINUE
      S = 0.
      DO 151 J=NSTART,NSTOP
         S = S+F(2,J)
  151 CONTINUE
      A2 = NUNK
      IF (NBDCND .EQ. 0) GO TO 152
      S = S-.5*(F(2,NSTART)+F(2,NSTOP))
      A2 = A2-1.
  152 YPOLE = (.25*DLRSQ*F(1,1)-S/A2)/(W(ID5+1)-1.+ELMBDA*DLRSQ*.25)
      DO 154 I=MSTART,MSTOP
         K = L+I
         DO 153 J=NSTART,NSTOP
            F(I,J) = F(I,J)+YPOLE*W(K)
  153    CONTINUE
  154 CONTINUE
  155 DO 156 J=1,NP1
         F(1,J) = YPOLE
  156 CONTINUE
  157 CONTINUE
      IF (NBDCND .NE. 0) GO TO 159
      DO 158 I=MSTART,MSTOP
         F(I,NP1) = F(I,1)
  158 CONTINUE
  159 CONTINUE
      W(1) = IWSTOR
      RETURN
      END
