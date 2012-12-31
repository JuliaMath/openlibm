*DECK POISTG
      SUBROUTINE POISTG (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y,
     +   IERROR, W)
C***BEGIN PROLOGUE  POISTG
C***PURPOSE  Solve a block tridiagonal system of linear equations
C            that results from a staggered grid finite difference
C            approximation to 2-D elliptic PDE's.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B4B
C***TYPE      SINGLE PRECISION (POISTG-S)
C***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, TRIDIAGONAL
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine POISTG solves the linear system of equations
C
C       A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C       + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
C
C       for I=1,2,...,M and J=1,2,...,N.
C
C     The indices I+1 and I-1 are evaluated modulo M, i.e.
C     X(0,J) = X(M,J) and X(M+1,J) = X(1,J), and X(I,0) may be equal to
C     X(I,1) or -X(I,1) and X(I,N+1) may be equal to X(I,N) or -X(I,N)
C     depending on an input parameter.
C
C
C     * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C             * * * * * *   On Input    * * * * * *
C
C   NPEROD
C     Indicates the values which X(I,0) and X(I,N+1) are assumed
C     to have.
C     = 1 If X(I,0) = -X(I,1) and X(I,N+1) = -X(I,N)
C     = 2 If X(I,0) = -X(I,1) and X(I,N+1) =  X(I,N)
C     = 3 If X(I,0) =  X(I,1) and X(I,N+1) =  X(I,N)
C     = 4 If X(I,0) =  X(I,1) and X(I,N+1) = -X(I,N)
C
C   N
C     The number of unknowns in the J-direction.  N must
C     be greater than 2.
C
C   MPEROD
C     = 0 If A(1) and C(M) are not zero
C     = 1 If A(1) = C(M) = 0
C
C   M
C     The number of unknowns in the I-direction.  M must
C     be greater than 2.
C
C   A,B,C
C     One-dimensional arrays of length M that specify the coefficients
C     in the linear equations given above.  If MPEROD = 0 the array
C     elements must not depend on the index I, but must be constant.
C     Specifically, the subroutine checks the following condition
C
C           A(I) = C(1)
C           B(I) = B(1)
C           C(I) = C(1)
C
C     for I = 1, 2, ..., M.
C
C   IDIMY
C     The row (or first) dimension of the two-dimensional array Y as
C     it appears in the program calling POISTG.  This parameter is
C     used to specify the variable dimension of Y.  IDIMY must be at
C     least M.
C
C   Y
C     A two-dimensional array that specifies the values of the
C     right side of the linear system of equations given above.
C     Y must be dimensioned at least M X N.
C
C   W
C     A one-dimensional work array that must be provided by the user
C     for work space.  W may require up to 9M + 4N + M(INT(log2(N)))
C     locations.  The actual number of locations used is computed by
C     POISTG and returned in location W(1).
C
C
C             * * * * * *   On Output     * * * * * *
C
C   Y
C     Contains the solution X.
C
C   IERROR
C     An error flag that indicates invalid input parameters.  Except
C     for number zero, a solution is not attempted.
C     = 0  No error
C     = 1  If M .LE. 2
C     = 2  If N .LE. 2
C     = 3  IDIMY .LT. M
C     = 4  If NPEROD .LT. 1 or NPEROD .GT. 4
C     = 5  If MPEROD .LT. 0 or MPEROD .GT. 1
C     = 6  If MPEROD = 0 and
C          A(I) .NE. C(1) or B(I) .NE. B(1) or C(I) .NE. C(1)
C          for some I = 1, 2, ..., M.
C       = 7 If MPEROD .EQ. 1 .AND. (A(1).NE.0 .OR. C(M).NE.0)
C
C   W
C     W(1) contains the required length of W.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   A(M),B(M),C(M),Y(IDIMY,N),
C     Arguments      W(see argument list)
C
C     Latest         June 1, 1977
C     Revision
C
C     Subprograms    POISTG,POSTG2,COSGEN,MERGE,TRIX,TRI3,PIMACH
C     Required
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
C     History        Written by Roland Sweet in 1973
C                    Revised by Roland Sweet in 1977
C
C
C     Space          3297(decimal) = 6341(octal) locations on the
C     Required       NCAR Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine POISTG is roughly proportional
C                    to M*N*log2(N).  Some typical values are listed
C                    in the table below.  More comprehensive timing
C                    charts may be found in the reference.
C                       To measure the accuracy of the algorithm a
C                    uniform random number generator was used to create
C                    a solution array X for the system given in the
C                    'PURPOSE ' with
C
C                       A(I) = C(I) = -0.5*B(I) = 1,       I=1,2,...,M
C
C                    and, when MPEROD = 1
C
C                       A(1) = C(M) = 0
C                       B(1) = B(M) =-1.
C
C                    The solution X was substituted into the given sys-
C                    tem and, using double precision, a right side Y was
C                    computed.  Using this array Y subroutine POISTG was
C                    called to produce an approximate solution Z.  Then
C                    the relative error, defined as
C
C                       E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J)))
C
C                    where the two maxima are taken over all I=1,2,...,M
C                    and J=1,2,...,N, was computed.  The value of E is
C                    given in the table below for some typical values of
C                    M and N.
C
C
C                       M (=N)    MPEROD    NPEROD    T(MSECS)    E
C                       ------    ------    ------    --------  ------
C
C                         31        0-1       1-4        45     9.E-13
C                         31        1         1          21     4.E-13
C                         31        1         3          41     3.E-13
C                         32        0-1       1-4        51     3.E-12
C                         32        1         1          32     3.E-13
C                         32        1         3          48     1.E-13
C                         33        0-1       1-4        42     1.E-12
C                         33        1         1          30     4.E-13
C                         33        1         3          34     1.E-13
C                         63        0-1       1-4       186     3.E-12
C                         63        1         1          91     1.E-12
C                         63        1         3         173     2.E-13
C                         64        0-1       1-4       209     4.E-12
C                         64        1         1         128     1.E-12
C                         64        1         3         199     6.E-13
C                         65        0-1       1-4       143     2.E-13
C                         65        1         1         160     1.E-11
C                         65        1         3         138     4.E-13
C
C     Portability    American National Standards Institute FORTRAN.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C     Required       COS
C     Resident
C     Routines
C
C     Reference      Schumann, U. and R. Sweet,'A Direct Method for
C                    the Solution of Poisson's Equation With Neumann
C                    Boundary Conditions on a Staggered Grid of
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
C***ROUTINES CALLED  POSTG2
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  POISTG
C
C
      DIMENSION       Y(IDIMY,*)
      DIMENSION       W(*)       ,B(*)       ,A(*)       ,C(*)
C***FIRST EXECUTABLE STATEMENT  POISTG
      IERROR = 0
      IF (M .LE. 2) IERROR = 1
      IF (N .LE. 2) IERROR = 2
      IF (IDIMY .LT. M) IERROR = 3
      IF (NPEROD.LT.1 .OR. NPEROD.GT.4) IERROR = 4
      IF (MPEROD.LT.0 .OR. MPEROD.GT.1) IERROR = 5
      IF (MPEROD .EQ. 1) GO TO 103
      DO 101 I=1,M
         IF (A(I) .NE. C(1)) GO TO 102
         IF (C(I) .NE. C(1)) GO TO 102
         IF (B(I) .NE. B(1)) GO TO 102
  101 CONTINUE
      GO TO 104
  102 IERROR = 6
      RETURN
  103 IF (A(1).NE.0. .OR. C(M).NE.0.) IERROR = 7
  104 IF (IERROR .NE. 0) RETURN
      IWBA = M+1
      IWBB = IWBA+M
      IWBC = IWBB+M
      IWB2 = IWBC+M
      IWB3 = IWB2+M
      IWW1 = IWB3+M
      IWW2 = IWW1+M
      IWW3 = IWW2+M
      IWD = IWW3+M
      IWTCOS = IWD+M
      IWP = IWTCOS+4*N
      DO 106 I=1,M
         K = IWBA+I-1
         W(K) = -A(I)
         K = IWBC+I-1
         W(K) = -C(I)
         K = IWBB+I-1
         W(K) = 2.-B(I)
         DO 105 J=1,N
            Y(I,J) = -Y(I,J)
  105    CONTINUE
  106 CONTINUE
      NP = NPEROD
      MP = MPEROD+1
      GO TO (110,107),MP
  107 CONTINUE
      GO TO (108,108,108,119),NPEROD
  108 CONTINUE
      CALL POSTG2 (NP,N,M,W(IWBA),W(IWBB),W(IWBC),IDIMY,Y,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD .EQ. 4) GO TO 120
  109 CONTINUE
      GO TO (123,129),MP
  110 CONTINUE
C
C     REORDER UNKNOWNS WHEN MP =0
C
      MH = (M+1)/2
      MHM1 = MH-1
      MODD = 1
      IF (MH*2 .EQ. M) MODD = 2
      DO 115 J=1,N
         DO 111 I=1,MHM1
            MHPI = MH+I
            MHMI = MH-I
            W(I) = Y(MHMI,J)-Y(MHPI,J)
            W(MHPI) = Y(MHMI,J)+Y(MHPI,J)
  111    CONTINUE
         W(MH) = 2.*Y(MH,J)
         GO TO (113,112),MODD
  112    W(M) = 2.*Y(M,J)
  113    CONTINUE
         DO 114 I=1,M
            Y(I,J) = W(I)
  114    CONTINUE
  115 CONTINUE
      K = IWBC+MHM1-1
      I = IWBA+MHM1
      W(K) = 0.
      W(I) = 0.
      W(K+1) = 2.*W(K+1)
      GO TO (116,117),MODD
  116 CONTINUE
      K = IWBB+MHM1-1
      W(K) = W(K)-W(I-1)
      W(IWBC-1) = W(IWBC-1)+W(IWBB-1)
      GO TO 118
  117 W(IWBB-1) = W(K+1)
  118 CONTINUE
      GO TO 107
  119 CONTINUE
C
C     REVERSE COLUMNS WHEN NPEROD = 4.
C
      IREV = 1
      NBY2 = N/2
      NP = 2
  120 DO 122 J=1,NBY2
         MSKIP = N+1-J
         DO 121 I=1,M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
  121    CONTINUE
  122 CONTINUE
      GO TO (108,109),IREV
  123 CONTINUE
      DO 128 J=1,N
         DO 124 I=1,MHM1
            MHMI = MH-I
            MHPI = MH+I
            W(MHMI) = .5*(Y(MHPI,J)+Y(I,J))
            W(MHPI) = .5*(Y(MHPI,J)-Y(I,J))
  124    CONTINUE
         W(MH) = .5*Y(MH,J)
         GO TO (126,125),MODD
  125    W(M) = .5*Y(M,J)
  126    CONTINUE
         DO 127 I=1,M
            Y(I,J) = W(I)
  127    CONTINUE
  128 CONTINUE
  129 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR W ARRAY.
C
      W(1) = IPSTOR+IWP-1
      RETURN
      END
