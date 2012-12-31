*DECK GENBUN
      SUBROUTINE GENBUN (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y,
     +   IERROR, W)
C***BEGIN PROLOGUE  GENBUN
C***PURPOSE  Solve by a cyclic reduction algorithm the linear system
C            of equations that results from a finite difference
C            approximation to certain 2-d elliptic PDE's on a centered
C            grid .
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B4B
C***TYPE      SINGLE PRECISION (GENBUN-S, CMGNBN-C)
C***KEYWORDS  ELLIPTIC, FISHPACK, PDE, TRIDIAGONAL
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine GENBUN solves the linear system of equations
C
C          A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C
C          + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
C
C               for I = 1,2,...,M  and  J = 1,2,...,N.
C
C     The indices I+1 and I-1 are evaluated modulo M, i.e.,
C     X(0,J) = X(M,J) and X(M+1,J) = X(1,J), and X(I,0) may be equal to
C     0, X(I,2), or X(I,N) and X(I,N+1) may be equal to 0, X(I,N-1), or
C     X(I,1) depending on an input parameter.
C
C
C     * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C             * * * * * *   On Input    * * * * * *
C
C     NPEROD
C       Indicates the values that X(I,0) and X(I,N+1) are assumed to
C       have.
C
C       = 0  If X(I,0) = X(I,N) and X(I,N+1) = X(I,1).
C       = 1  If X(I,0) = X(I,N+1) = 0  .
C       = 2  If X(I,0) = 0 and X(I,N+1) = X(I,N-1).
C       = 3  If X(I,0) = X(I,2) and X(I,N+1) = X(I,N-1).
C       = 4  If X(I,0) = X(I,2) and X(I,N+1) = 0.
C
C     N
C       The number of unknowns in the J-direction.  N must be greater
C       than 2.
C
C     MPEROD
C       = 0 if A(1) and C(M) are not zero.
C       = 1 if A(1) = C(M) = 0.
C
C     M
C       The number of unknowns in the I-direction.  M must be greater
C       than 2.
C
C     A,B,C
C       One-dimensional arrays of length M that specify the
C       coefficients in the linear equations given above.  If MPEROD = 0
C       the array elements must not depend upon the index I, but must be
C       constant.  Specifically, the subroutine checks the following
C       condition
C
C             A(I) = C(1)
C             C(I) = C(1)
C             B(I) = B(1)
C
C       for I=1,2,...,M.
C
C     IDIMY
C       The row (or first) dimension of the two-dimensional array Y as
C       it appears in the program calling GENBUN.  This parameter is
C       used to specify the variable dimension of Y.  IDIMY must be at
C       least M.
C
C     Y
C       A two-dimensional array that specifies the values of the right
C       side of the linear system of equations given above.  Y must be
C       dimensioned at least M*N.
C
C     W
C       A one-dimensional array that must be provided by the user for
C       work space.  W may require up to 4*N + (10 + INT(log2(N)))*M
C       locations.  The actual number of locations used is computed by
C       GENBUN and is returned in location W(1).
C
C
C             * * * * * *   On Output     * * * * * *
C
C     Y
C       Contains the solution X.
C
C     IERROR
C       An error flag that indicates invalid input parameters.  Except
C       for number zero, a solution is not attempted.
C
C       = 0  No error.
C       = 1  M .LE. 2
C       = 2  N .LE. 2
C       = 3  IDIMY .LT. M
C       = 4  NPEROD .LT. 0 or NPEROD .GT. 4
C       = 5  MPEROD .LT. 0 or MPEROD .GT. 1
C       = 6  A(I) .NE. C(1) or C(I) .NE. C(1) or B(I) .NE. B(1) for
C            some I=1,2,...,M.
C       = 7  A(1) .NE. 0 or C(M) .NE. 0 and MPEROD = 1
C
C     W
C       W(1) contains the required length of W.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   A(M),B(M),C(M),Y(IDIMY,N),W(see parameter list)
C     Arguments
C
C     Latest         June 1, 1976
C     Revision
C
C     Subprograms    GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,TRIX,TRI3,
C     Required       PIMACH
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
C     History        Standardized April 1, 1973
C                    Revised August 20,1973
C                    Revised January 1, 1976
C
C     Algorithm      The linear system is solved by a cyclic reduction
C                    algorithm described in the reference.
C
C     Space          4944(decimal) = 11520(octal) locations on the NCAR
C     Required       Control Data 7600.
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine GENBUN is roughly proportional
C                    to M*N*log2(N), but also depends on the input
C                    parameter NPEROD.  Some typical values are listed
C                    in the table below.  More comprehensive timing
C                    charts may be found in the reference.
C                       To measure the accuracy of the algorithm a
C                    uniform random number generator was used to create
C                    a solution array X for the system given in the
C                    'PURPOSE' with
C
C                       A(I) = C(I) = -0.5*B(I) = 1,       I=1,2,...,M
C
C                    and, when MPEROD = 1
C
C                       A(1) = C(M) = 0
C                       A(M) = C(1) = 2.
C
C                    The solution X was substituted into the given sys-
C                    tem and, using double precision, a right side Y was
C                    computed.  Using this array Y subroutine GENBUN was
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
C                         31        0         0          36     6.E-14
C                         31        1         1          21     4.E-13
C                         31        1         3          41     3.E-13
C                         32        0         0          29     9.E-14
C                         32        1         1          32     3.E-13
C                         32        1         3          48     1.E-13
C                         33        0         0          36     9.E-14
C                         33        1         1          30     4.E-13
C                         33        1         3          34     1.E-13
C                         63        0         0         150     1.E-13
C                         63        1         1          91     1.E-12
C                         63        1         3         173     2.E-13
C                         64        0         0         122     1.E-13
C                         64        1         1         128     1.E-12
C                         64        1         3         199     6.E-13
C                         65        0         0         143     2.E-13
C                         65        1         1         120     1.E-12
C                         65        1         3         138     4.E-13
C
C     Portability    American National Standards Institute Fortran.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C     Required       COS
C     Resident
C     Routines
C
C     Reference      Sweet, R., 'A Cyclic Reduction Algorithm For
C                    Solving Block Tridiagonal Systems Of Arbitrary
C                    Dimensions,' SIAM J. on Numer. Anal.,
C                    14(Sept., 1977), PP. 706-720.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  R. Sweet, A cyclic reduction algorithm for solving
C                 block tridiagonal systems of arbitrary dimensions,
C                 SIAM Journal on Numerical Analysis 14, (September
C                 1977), pp. 706-720.
C***ROUTINES CALLED  POISD2, POISN2, POISP2
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  GENBUN
C
C
      DIMENSION       Y(IDIMY,*)
      DIMENSION       W(*)       ,B(*)       ,A(*)       ,C(*)
C***FIRST EXECUTABLE STATEMENT  GENBUN
      IERROR = 0
      IF (M .LE. 2) IERROR = 1
      IF (N .LE. 2) IERROR = 2
      IF (IDIMY .LT. M) IERROR = 3
      IF (NPEROD.LT.0 .OR. NPEROD.GT.4) IERROR = 4
      IF (MPEROD.LT.0 .OR. MPEROD.GT.1) IERROR = 5
      IF (MPEROD .EQ. 1) GO TO 102
      DO 101 I=2,M
         IF (A(I) .NE. C(1)) GO TO 103
         IF (C(I) .NE. C(1)) GO TO 103
         IF (B(I) .NE. B(1)) GO TO 103
  101 CONTINUE
      GO TO 104
  102 IF (A(1).NE.0. .OR. C(M).NE.0.) IERROR = 7
      GO TO 104
  103 IERROR = 6
  104 IF (IERROR .NE. 0) RETURN
      MP1 = M+1
      IWBA = MP1
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
      MP = MPEROD+1
      NP = NPEROD+1
      GO TO (114,107),MP
  107 GO TO (108,109,110,111,123),NP
  108 CALL POISP2 (M,N,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      GO TO 112
  109 CALL POISD2 (M,N,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWW1),
     1             W(IWD),W(IWTCOS),W(IWP))
      GO TO 112
  110 CALL POISN2 (M,N,1,2,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      GO TO 112
  111 CALL POISN2 (M,N,1,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
  112 IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD .EQ. 4) GO TO 124
  113 GO TO (127,133),MP
  114 CONTINUE
C
C     REORDER UNKNOWNS WHEN MP =0
C
      MH = (M+1)/2
      MHM1 = MH-1
      MODD = 1
      IF (MH*2 .EQ. M) MODD = 2
      DO 119 J=1,N
         DO 115 I=1,MHM1
            MHPI = MH+I
            MHMI = MH-I
            W(I) = Y(MHMI,J)-Y(MHPI,J)
            W(MHPI) = Y(MHMI,J)+Y(MHPI,J)
  115    CONTINUE
         W(MH) = 2.*Y(MH,J)
         GO TO (117,116),MODD
  116    W(M) = 2.*Y(M,J)
  117    CONTINUE
         DO 118 I=1,M
            Y(I,J) = W(I)
  118    CONTINUE
  119 CONTINUE
      K = IWBC+MHM1-1
      I = IWBA+MHM1
      W(K) = 0.
      W(I) = 0.
      W(K+1) = 2.*W(K+1)
      GO TO (120,121),MODD
  120 CONTINUE
      K = IWBB+MHM1-1
      W(K) = W(K)-W(I-1)
      W(IWBC-1) = W(IWBC-1)+W(IWBB-1)
      GO TO 122
  121 W(IWBB-1) = W(K+1)
  122 CONTINUE
      GO TO 107
C
C     REVERSE COLUMNS WHEN NPEROD = 4.
C
  123 IREV = 1
      NBY2 = N/2
  124 DO 126 J=1,NBY2
         MSKIP = N+1-J
         DO 125 I=1,M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
  125    CONTINUE
  126 CONTINUE
      GO TO (110,113),IREV
  127 CONTINUE
      DO 132 J=1,N
         DO 128 I=1,MHM1
            MHMI = MH-I
            MHPI = MH+I
            W(MHMI) = .5*(Y(MHPI,J)+Y(I,J))
            W(MHPI) = .5*(Y(MHPI,J)-Y(I,J))
  128    CONTINUE
         W(MH) = .5*Y(MH,J)
         GO TO (130,129),MODD
  129    W(M) = .5*Y(M,J)
  130    CONTINUE
         DO 131 I=1,M
            Y(I,J) = W(I)
  131    CONTINUE
  132 CONTINUE
  133 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR W ARRAY.
C
      W(1) = IPSTOR+IWP-1
      RETURN
      END
