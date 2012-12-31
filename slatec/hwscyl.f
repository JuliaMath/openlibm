*DECK HWSCYL
      SUBROUTINE HWSCYL (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND,
     +   BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C***BEGIN PROLOGUE  HWSCYL
C***PURPOSE  Solve a standard finite difference approximation
C            to the Helmholtz equation in cylindrical coordinates.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B1A1A
C***TYPE      SINGLE PRECISION (HWSCYL-S)
C***KEYWORDS  CYLINDRICAL, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine HWSCYL solves a finite difference approximation to the
C     Helmholtz equation in cylindrical coordinates:
C
C          (1/R)(d/dR)(R(dU/dR)) + (d/dZ)(dU/dZ)
C
C                                + (LAMBDA/R**2)U = F(R,Z)
C
C     This modified Helmholtz equation results from the Fourier
C     transform of the three-dimensional Poisson equation.
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
C       Indicates the type of boundary conditions at R = A and R = B.
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
C            BDA(J) = (d/dR)U(A,Z(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value, BDA is a dummy variable.
C
C     BDB
C       A one-dimensional array of length N+1 that specifies the values
C       of the derivative of the solution with respect to R at R = B.
C       When MBDCND = 2,3, or 6,
C
C            BDB(J) = (d/dR)U(B,Z(J)), J = 1,2,...,N+1  .
C
C       When MBDCND has any other value, BDB is a dummy variable.
C
C     C,D
C       The range of Z, i.e., C .LE. Z .LE. D.  C must be less than D.
C
C     N
C       The number of panels into which the interval (C,D) is
C       subdivided.  Hence, there will be N+1 grid points in the
C       Z-direction given by Z(J) = C+(J-1)DZ, for J = 1,2,...,N+1,
C       where DZ = (D-C)/N is the panel width. N must be greater than 3.
C
C     NBDCND
C       Indicates the type of boundary conditions at Z = C and Z = D.
C
C       = 0  If the solution is periodic in Z, i.e., U(I,1) = U(I,N+1).
C       = 1  If the solution is specified at Z = C and Z = D.
C       = 2  If the solution is specified at Z = C and the derivative of
C            the solution with respect to Z is specified at Z = D.
C       = 3  If the derivative of the solution with respect to Z is
C            specified at Z = C and Z = D.
C       = 4  If the derivative of the solution with respect to Z is
C            specified at Z = C and the solution is specified at Z = D.
C
C     BDC
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to Z at Z = C.
C       When NBDCND = 3 or 4,
C
C            BDC(I) = (d/dZ)U(R(I),C), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDC is a dummy variable.
C
C     BDD
C       A one-dimensional array of length M+1 that specifies the values
C       of the derivative of the solution with respect to Z at Z = D.
C       When NBDCND = 2 or 3,
C
C            BDD(I) = (d/dZ)U(R(I),D), I = 1,2,...,M+1  .
C
C       When NBDCND has any other value, BDD is a dummy variable.
C
C     ELMBDA
C       The constant LAMBDA in the Helmholtz equation.  If
C       LAMBDA .GT. 0, a solution may not exist.  However, HWSCYL will
C       attempt to find a solution.  LAMBDA must be zero when
C       MBDCND = 5 or 6  .
C
C     F
C       A two-dimensional array that specifies the values of the right
C       side of the Helmholtz equation and boundary data (if any).  For
C       I = 2,3,...,M and J = 2,3,...,N
C
C            F(I,J) = F(R(I),Z(J)).
C
C       On the boundaries F is defined by
C
C            MBDCND   F(1,J)            F(M+1,J)
C            ------   ---------         ---------
C
C              1      U(A,Z(J))         U(B,Z(J))
C              2      U(A,Z(J))         F(B,Z(J))
C              3      F(A,Z(J))         F(B,Z(J))   J = 1,2,...,N+1
C              4      F(A,Z(J))         U(B,Z(J))
C              5      F(0,Z(J))         U(B,Z(J))
C              6      F(0,Z(J))         F(B,Z(J))
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
C     IDIMF
C       The row (or first) dimension of the array F as it appears in the
C       program calling HWSCYL.  This parameter is used to specify the
C       variable dimension of F.  IDIMF must be at least M+1  .
C
C     W
C       A one-dimensional array that must be provided by the user for
C       work space.  W may require up to 4*(N+1) +
C       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of
C       locations used is computed by HWSCYL and is returned in location
C       W(1).
C
C
C             * * * * * *   On Output     * * * * * *
C
C     F
C       Contains the solution U(I,J) of the finite difference
C       approximation for the grid point (R(I),Z(J)), I = 1,2,...,M+1,
C       J = 1,2,...,N+1  .
C
C     PERTRB
C       If one specifies a combination of periodic, derivative, and
C       unspecified boundary conditions for a Poisson equation
C       (LAMBDA = 0), a solution may not exist.  PERTRB is a constant,
C       calculated and subtracted from F, which ensures that a solution
C       exists.  HWSCYL then computes this solution, which is a least
C       squares solution to the original approximation.  This solution
C       plus any constant is also a solution.  Hence, the solution is
C       not unique.  The value of PERTRB should be small compared to the
C       right side F.  Otherwise, a solution is obtained to an
C       essentially different problem.  This comparison should always
C       be made to insure that a meaningful solution has been obtained.
C
C     IERROR
C       An error flag which indicates invalid input parameters.  Except
C       for numbers 0 and 11, a solution is not attempted.
C
C       =  0  No error.
C       =  1  A .LT. 0  .
C       =  2  A .GE. B.
C       =  3  MBDCND .LT. 1 or MBDCND .GT. 6  .
C       =  4  C .GE. D.
C       =  5  N .LE. 3
C       =  6  NBDCND .LT. 0 or NBDCND .GT. 4  .
C       =  7  A = 0, MBDCND = 3 or 4  .
C       =  8  A .GT. 0, MBDCND .GE. 5  .
C       =  9  A = 0, LAMBDA .NE. 0, MBDCND .GE. 5  .
C       = 10  IDIMF .LT. M+1  .
C       = 11  LAMBDA .GT. 0  .
C       = 12  M .LE. 3
C
C       Since this is the only means of indicating a possibly incorrect
C       call to HWSCYL, the user should test IERROR after the call.
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
C     Subprograms    HWSCYL,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
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
C     Specialist     Roland Sweet
C
C     Language       FORTRAN
C
C     History        Standardized September 1, 1973
C                    Revised April 1, 1976
C
C     Algorithm      The routine defines the finite difference
C                    equations, incorporates boundary data, and adjusts
C                    the right side of singular systems and then calls
C                    GENBUN to solve the system.
C
C     Space          5818(decimal) = 13272(octal) locations on the NCAR
C     Required       Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine HWSCYL is roughly proportional
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
C                    Subprograms for the Solution of Elliptic Equations'
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
C***END PROLOGUE  HWSCYL
C
C
      DIMENSION       F(IDIMF,*)
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
C***FIRST EXECUTABLE STATEMENT  HWSCYL
      IERROR = 0
      IF (A .LT. 0.) IERROR = 1
      IF (A .GE. B) IERROR = 2
      IF (MBDCND.LE.0 .OR. MBDCND.GE.7) IERROR = 3
      IF (C .GE. D) IERROR = 4
      IF (N .LE. 3) IERROR = 5
      IF (NBDCND.LE.-1 .OR. NBDCND.GE.5) IERROR = 6
      IF (A.EQ.0. .AND. (MBDCND.EQ.3 .OR. MBDCND.EQ.4)) IERROR = 7
      IF (A.GT.0. .AND. MBDCND.GE.5) IERROR = 8
      IF (A.EQ.0. .AND. ELMBDA.NE.0. .AND. MBDCND.GE.5) IERROR = 9
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
      MSTOP = M
      GO TO (104,103,102,101,101,102),MBDCND
  101 MSTART = 1
      GO TO 104
  102 MSTART = 1
  103 MSTOP = MP1
  104 MUNK = MSTOP-MSTART+1
      NSTART = 1
      NSTOP = N
      GO TO (108,105,106,107,108),NP
  105 NSTART = 2
      GO TO 108
  106 NSTART = 2
  107 NSTOP = NP1
  108 NUNK = NSTOP-NSTART+1
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = MUNK
      ID3 = ID2+MUNK
      ID4 = ID3+MUNK
      ID5 = ID4+MUNK
      ID6 = ID5+MUNK
      ISTART = 1
      A1 = 2./DLRSQ
      IJ = 0
      IF (MBDCND.EQ.3 .OR. MBDCND.EQ.4) IJ = 1
      IF (MBDCND .LE. 4) GO TO 109
      W(1) = 0.
      W(ID2+1) = -2.*A1
      W(ID3+1) = 2.*A1
      ISTART = 2
      IJ = 1
  109 DO 110 I=ISTART,MUNK
         R = A+(I-IJ)*DELTAR
         J = ID5+I
         W(J) = R
         J = ID6+I
         W(J) = 1./R**2
         W(I) = (R-DLRBY2)/(R*DLRSQ)
         J = ID3+I
         W(J) = (R+DLRBY2)/(R*DLRSQ)
         K = ID6+I
         J = ID2+I
         W(J) = -A1+ELMBDA*W(K)
  110 CONTINUE
      GO TO (114,111,112,113,114,112),MBDCND
  111 W(ID2) = A1
      GO TO 114
  112 W(ID2) = A1
  113 W(ID3+1) = A1*ISTART
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
C     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
C
  124 A1 = 1./DLTHSQ
      L = ID5-MSTART+1
      GO TO (134,125,125,127,127),NP
  125 DO 126 I=MSTART,MSTOP
         F(I,2) = F(I,2)-A1*F(I,1)
  126 CONTINUE
      GO TO 129
  127 A1 = 2./DELTHT
      DO 128 I=MSTART,MSTOP
         F(I,1) = F(I,1)+A1*BDC(I)
  128 CONTINUE
  129 A1 = 1./DLTHSQ
      GO TO (134,130,132,132,130),NP
  130 DO 131 I=MSTART,MSTOP
         F(I,N) = F(I,N)-A1*F(I,NP1)
  131 CONTINUE
      GO TO 134
  132 A1 = 2./DELTHT
      DO 133 I=MSTART,MSTOP
         F(I,NP1) = F(I,NP1)-A1*BDD(I)
  133 CONTINUE
  134 CONTINUE
C
C     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
C     SOLUTION.
C
      PERTRB = 0.
      IF (ELMBDA) 146,136,135
  135 IERROR = 11
      GO TO 146
  136 W(ID5+1) = .5*(W(ID5+2)-DLRBY2)
      GO TO (146,146,138,146,146,137),MBDCND
  137 W(ID5+1) = .5*W(ID5+1)
  138 GO TO (140,146,146,139,146),NP
  139 A2 = 2.
      GO TO 141
  140 A2 = 1.
  141 K = ID5+MUNK
      W(K) = .5*(W(K-1)+DLRBY2)
      S = 0.
      DO 143 I=MSTART,MSTOP
         S1 = 0.
         NSP1 = NSTART+1
         NSTM1 = NSTOP-1
         DO 142 J=NSP1,NSTM1
            S1 = S1+F(I,J)
  142    CONTINUE
         K = I+L
         S = S+(A2*S1+F(I,NSTART)+F(I,NSTOP))*W(K)
  143 CONTINUE
      S2 = M*A+(.75+(M-1)*(M+1))*DLRBY2
      IF (MBDCND .EQ. 3) S2 = S2+.25*DLRBY2
      S1 = (2.+A2*(NUNK-2))*S2
      PERTRB = S/S1
      DO 145 I=MSTART,MSTOP
         DO 144 J=NSTART,NSTOP
            F(I,J) = F(I,J)-PERTRB
  144    CONTINUE
  145 CONTINUE
  146 CONTINUE
C
C     MULTIPLY I-TH EQUATION THROUGH BY DELTHT**2 TO PUT EQUATION INTO
C     CORRECT FORM FOR SUBROUTINE GENBUN.
C
      DO 148 I=MSTART,MSTOP
         K = I-MSTART+1
         W(K) = W(K)*DLTHSQ
         J = ID2+K
         W(J) = W(J)*DLTHSQ
         J = ID3+K
         W(J) = W(J)*DLTHSQ
         DO 147 J=NSTART,NSTOP
            F(I,J) = F(I,J)*DLTHSQ
  147    CONTINUE
  148 CONTINUE
      W(1) = 0.
      W(ID4) = 0.
C
C     CALL GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
C
      CALL GENBUN (NBDCND,NUNK,1,MUNK,W(1),W(ID2+1),W(ID3+1),IDIMF,
     1             F(MSTART,NSTART),IERR1,W(ID4+1))
      W(1) = W(ID4+1)+3*MUNK
      IF (NBDCND .NE. 0) GO TO 150
      DO 149 I=MSTART,MSTOP
         F(I,NP1) = F(I,1)
  149 CONTINUE
  150 CONTINUE
      RETURN
      END
