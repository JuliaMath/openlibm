*DECK POIS3D
      SUBROUTINE POIS3D (LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B,
     +   C, LDIMF, MDIMF, F, IERROR, W)
C***BEGIN PROLOGUE  POIS3D
C***PURPOSE  Solve a three-dimensional block tridiagonal linear system
C            which arises from a finite difference approximation to a
C            three-dimensional Poisson equation using the Fourier
C            transform package FFTPAK written by Paul Swarztrauber.
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B4B
C***TYPE      SINGLE PRECISION (POIS3D-S)
C***KEYWORDS  ELLIPTIC PDE, FISHPACK, HELMHOLTZ, POISSON
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine POIS3D solves the linear system of equations
C
C       C1*(X(I-1,J,K)-2.*X(I,J,K)+X(I+1,J,K))
C     + C2*(X(I,J-1,K)-2.*X(I,J,K)+X(I,J+1,K))
C     + A(K)*X(I,J,K-1)+B(K)*X(I,J,K)+C(K)*X(I,J,K+1) = F(I,J,K)
C
C     for  I=1,2,...,L , J=1,2,...,M , and K=1,2,...,N .
C
C     The indices K-1 and K+1 are evaluated modulo N, i.e.
C     X(I,J,0) = X(I,J,N) and X(I,J,N+1) = X(I,J,1). The unknowns
C     X(0,J,K), X(L+1,J,K), X(I,0,K), and X(I,M+1,K) are assumed to take
C     on certain prescribed values described below.
C
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C    * * * * * * * *    Parameter Description     * * * * * * * * * *
C
C
C            * * * * * *   On Input    * * * * * *
C
C     LPEROD   Indicates the values that X(0,J,K) and X(L+1,J,K) are
C              assumed to have.
C
C              = 0  If X(0,J,K) = X(L,J,K) and X(L+1,J,K) = X(1,J,K).
C              = 1  If X(0,J,K) = X(L+1,J,K) = 0.
C              = 2  If X(0,J,K) = 0  and X(L+1,J,K) = X(L-1,J,K).
C              = 3  If X(0,J,K) = X(2,J,K) and X(L+1,J,K) = X(L-1,J,K).
C              = 4  If X(0,J,K) = X(2,J,K) and X(L+1,J,K) = 0.
C
C     L        The number of unknowns in the I-direction. L must be at
C              least 3.
C
C     C1       The real constant that appears in the above equation.
C
C     MPEROD   Indicates the values that X(I,0,K) and X(I,M+1,K) are
C              assumed to have.
C
C              = 0  If X(I,0,K) = X(I,M,K) and X(I,M+1,K) = X(I,1,K).
C              = 1  If X(I,0,K) = X(I,M+1,K) = 0.
C              = 2  If X(I,0,K) = 0 and X(I,M+1,K) = X(I,M-1,K).
C              = 3  If X(I,0,K) = X(I,2,K) and X(I,M+1,K) = X(I,M-1,K).
C              = 4  If X(I,0,K) = X(I,2,K) and X(I,M+1,K) = 0.
C
C     M        The number of unknowns in the J-direction. M must be at
C              least 3.
C
C     C2       The real constant which appears in the above equation.
C
C     NPEROD   = 0  If A(1) and C(N) are not zero.
C              = 1  If A(1) = C(N) = 0.
C
C     N        The number of unknowns in the K-direction. N must be at
C              least 3.
C
C
C     A,B,C    One-dimensional arrays of length N that specify the
C              coefficients in the linear equations given above.
C
C              If NPEROD = 0 the array elements must not depend upon the
C              index K, but must be constant.  Specifically, the
C              subroutine checks the following condition
C
C                          A(K) = C(1)
C                          C(K) = C(1)
C                          B(K) = B(1)
C
C                  for K=1,2,...,N.
C
C     LDIMF    The row (or first) dimension of the three-dimensional
C              array F as it appears in the program calling POIS3D.
C              This parameter is used to specify the variable dimension
C              of F.  LDIMF must be at least L.
C
C     MDIMF    The column (or second) dimension of the three-dimensional
C              array F as it appears in the program calling POIS3D.
C              This parameter is used to specify the variable dimension
C              of F.  MDIMF must be at least M.
C
C     F        A three-dimensional array that specifies the values of
C              the right side of the linear system of equations given
C              above.  F must be dimensioned at least L x M x N.
C
C     W        A one-dimensional array that must be provided by the
C              user for work space.  The length of W must be at least
C              30 + L + M + 2*N + MAX(L,M,N) +
C              7*(INT((L+1)/2) + INT((M+1)/2)).
C
C
C            * * * * * *   On Output   * * * * * *
C
C     F        Contains the solution X.
C
C     IERROR   An error flag that indicates invalid input parameters.
C              Except for number zero, a solution is not attempted.
C              = 0  No error
C              = 1  If LPEROD .LT. 0 or .GT. 4
C              = 2  If L .LT. 3
C              = 3  If MPEROD .LT. 0 or .GT. 4
C              = 4  If M .LT. 3
C              = 5  If NPEROD .LT. 0 or .GT. 1
C              = 6  If N .LT. 3
C              = 7  If LDIMF .LT. L
C              = 8  If MDIMF .LT. M
C              = 9  If A(K) .NE. C(1) or C(K) .NE. C(1) or B(I) .NE.B(1)
C                      for some K=1,2,...,N.
C              = 10 If NPEROD = 1 and A(1) .NE. 0 or C(N) .NE. 0
C
C              Since this is the only means of indicating a possibly
C              incorrect call to POIS3D, the user should test IERROR
C              after the call.
C
C *Long Description:
C
C    * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   A(N),B(N),C(N),F(LDIMF,MDIMF,N),
C     Arguments      W(see argument list)
C
C     Latest         December 1, 1978
C     Revision
C
C     Subprograms    POIS3D,POS3D1,TRIDQ,RFFTI,RFFTF,RFFTF1,RFFTB,
C     Required       RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,COSQF1
C                    COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,CFFTI1,
C                    CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,CFFTF,
C                    CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,PIMACH,
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
C     History        Written by Roland Sweet at NCAR in July 1977
C
C     Algorithm      This subroutine solves three-dimensional block
C                    tridiagonal linear systems arising from finite
C                    difference approximations to three-dimensional
C                    Poisson equations using the Fourier transform
C                    package FFTPAK written by Paul Swarztrauber.
C
C     Space          6561(decimal) = 14641(octal) locations on the
C     Required       NCAR Control Data 7600
C
C     Timing and        The execution time T on the NCAR Control Data
C     Accuracy       7600 for subroutine POIS3D is roughly proportional
C                    to L*M*N*(log2(L)+log2(M)+5), but also depends on
C                    input parameters LPEROD and MPEROD.  Some typical
C                    values are listed in the table below when NPEROD=0.
C                       To measure the accuracy of the algorithm a
C                    uniform random number generator was used to create
C                    a solution array X for the system given in the
C                    'PURPOSE' with
C
C                       A(K) = C(K) = -0.5*B(K) = 1,       K=1,2,...,N
C
C                    and, when NPEROD = 1
C
C                       A(1) = C(N) = 0
C                       A(N) = C(1) = 2.
C
C                    The solution X was substituted into the given sys-
C                    tem and, using double precision, a right side Y was
C                    computed.  Using this array Y subroutine POIS3D was
C                    called to produce an approximate solution Z.  Then
C                    the relative error, defined as
C
C                    E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K)))
C
C                    where the two maxima are taken over I=1,2,...,L,
C                    J=1,2,...,M and K=1,2,...,N, was computed.  The
C                    value of E is given in the table below for some
C                    typical values of L,M and N.
C
C
C                       L(=M=N)   LPEROD    MPEROD    T(MSECS)    E
C                       ------    ------    ------    --------  ------
C
C                         16        0         0         272     1.E-13
C                         15        1         1         287     4.E-13
C                         17        3         3         338     2.E-13
C                         32        0         0        1755     2.E-13
C                         31        1         1        1894     2.E-12
C                         33        3         3        2042     7.E-13
C
C
C     Portability    American National Standards Institute FORTRAN.
C                    The machine dependent constant PI is defined in
C                    function PIMACH.
C
C     Required       COS,SIN,ATAN
C     Resident
C     Routines
C
C     Reference      NONE
C
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  POS3D1
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  POIS3D
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,
     1                F(LDIMF,MDIMF,*)       ,W(*)       ,SAVE(6)
C***FIRST EXECUTABLE STATEMENT  POIS3D
      LP = LPEROD+1
      MP = MPEROD+1
      NP = NPEROD+1
C
C     CHECK FOR INVALID INPUT.
C
      IERROR = 0
      IF (LP.LT.1 .OR. LP.GT.5) IERROR = 1
      IF (L .LT. 3) IERROR = 2
      IF (MP.LT.1 .OR. MP.GT.5) IERROR = 3
      IF (M .LT. 3) IERROR = 4
      IF (NP.LT.1 .OR. NP.GT.2) IERROR = 5
      IF (N .LT. 3) IERROR = 6
      IF (LDIMF .LT. L) IERROR = 7
      IF (MDIMF .LT. M) IERROR = 8
      IF (NP .NE. 1) GO TO 103
      DO 101 K=1,N
         IF (A(K) .NE. C(1)) GO TO 102
         IF (C(K) .NE. C(1)) GO TO 102
         IF (B(K) .NE. B(1)) GO TO 102
  101 CONTINUE
      GO TO 104
  102 IERROR = 9
  103 IF (NPEROD.EQ.1 .AND. (A(1).NE.0. .OR. C(N).NE.0.)) IERROR = 10
  104 IF (IERROR .NE. 0) GO TO 122
      IWYRT = L+1
      IWT = IWYRT+M
      IWD = IWT+MAX(L,M,N)+1
      IWBB = IWD+N
      IWX = IWBB+N
      IWY = IWX+7*((L+1)/2)+15
      GO TO (105,114),NP
C
C     REORDER UNKNOWNS WHEN NPEROD = 0.
C
  105 NH = (N+1)/2
      NHM1 = NH-1
      NODD = 1
      IF (2*NH .EQ. N) NODD = 2
      DO 111 I=1,L
         DO 110 J=1,M
            DO 106 K=1,NHM1
               NHPK = NH+K
               NHMK = NH-K
               W(K) = F(I,J,NHMK)-F(I,J,NHPK)
               W(NHPK) = F(I,J,NHMK)+F(I,J,NHPK)
  106       CONTINUE
            W(NH) = 2.*F(I,J,NH)
            GO TO (108,107),NODD
  107       W(N) = 2.*F(I,J,N)
  108       DO 109 K=1,N
               F(I,J,K) = W(K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      SAVE(1) = C(NHM1)
      SAVE(2) = A(NH)
      SAVE(3) = C(NH)
      SAVE(4) = B(NHM1)
      SAVE(5) = B(N)
      SAVE(6) = A(N)
      C(NHM1) = 0.
      A(NH) = 0.
      C(NH) = 2.*C(NH)
      GO TO (112,113),NODD
  112 B(NHM1) = B(NHM1)-A(NH-1)
      B(N) = B(N)+A(N)
      GO TO 114
  113 A(N) = C(NH)
  114 CONTINUE
      CALL POS3D1 (LP,L,MP,M,N,A,B,C,LDIMF,MDIMF,F,W,W(IWYRT),W(IWT),
     1             W(IWD),W(IWX),W(IWY),C1,C2,W(IWBB))
      GO TO (115,122),NP
  115 DO 121 I=1,L
         DO 120 J=1,M
            DO 116 K=1,NHM1
               NHMK = NH-K
               NHPK = NH+K
               W(NHMK) = .5*(F(I,J,NHPK)+F(I,J,K))
               W(NHPK) = .5*(F(I,J,NHPK)-F(I,J,K))
  116       CONTINUE
            W(NH) = .5*F(I,J,NH)
            GO TO (118,117),NODD
  117       W(N) = .5*F(I,J,N)
  118       DO 119 K=1,N
               F(I,J,K) = W(K)
  119       CONTINUE
  120    CONTINUE
  121 CONTINUE
      C(NHM1) = SAVE(1)
      A(NH) = SAVE(2)
      C(NH) = SAVE(3)
      B(NHM1) = SAVE(4)
      B(N) = SAVE(5)
      A(N) = SAVE(6)
  122 CONTINUE
      RETURN
      END
