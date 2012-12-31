*DECK DULSIA
      SUBROUTINE DULSIA (A, MDA, M, N, B, MDB, NB, RE, AE, KEY, MODE,
     +   NP, KRANK, KSURE, RNORM, W, LW, IWORK, LIW, INFO)
C***BEGIN PROLOGUE  DULSIA
C***PURPOSE  Solve an underdetermined linear system of equations by
C            performing an LQ factorization of the matrix using
C            Householder transformations.  Emphasis is put on detecting
C            possible rank deficiency.
C***LIBRARY   SLATEC
C***CATEGORY  D9
C***TYPE      DOUBLE PRECISION (ULSIA-S, DULSIA-D)
C***KEYWORDS  LINEAR LEAST SQUARES, LQ FACTORIZATION,
C             UNDERDETERMINED LINEAR SYSTEM
C***AUTHOR  Manteuffel, T. A., (LANL)
C***DESCRIPTION
C
C     DULSIA computes the minimal length solution(s) to the problem AX=B
C     where A is an M by N matrix with M.LE.N and B is the M by NB
C     matrix of right hand sides.  User input bounds on the uncertainty
C     in the elements of A are used to detect numerical rank deficiency.
C     The algorithm employs a row and column pivot strategy to
C     minimize the growth of uncertainty and round-off errors.
C
C     DULSIA requires (MDA+1)*N + (MDB+1)*NB + 6*M dimensioned space
C
C   ******************************************************************
C   *                                                                *
C   *         WARNING - All input arrays are changed on exit.        *
C   *                                                                *
C   ******************************************************************
C
C     Input.. All TYPE REAL variables are DOUBLE PRECISION
C
C     A(,)          Linear coefficient matrix of AX=B, with MDA the
C      MDA,M,N      actual first dimension of A in the calling program.
C                   M is the row dimension (no. of EQUATIONS of the
C                   problem) and N the col dimension (no. of UNKNOWNS).
C                   Must have MDA.GE.M and M.LE.N.
C
C     B(,)          Right hand side(s), with MDB the actual first
C      MDB,NB       dimension of B in the calling program. NB is the
C                   number of M by 1 right hand sides.  Since the
C                   solution is returned in B, must have MDB.GE.N.  If
C                   NB = 0, B is never accessed.
C
C   ******************************************************************
C   *                                                                *
C   *         Note - Use of RE and AE are what make this             *
C   *                code significantly different from               *
C   *                other linear least squares solvers.             *
C   *                However, the inexperienced user is              *
C   *                advised to set RE=0.,AE=0.,KEY=0.               *
C   *                                                                *
C   ******************************************************************
C
C     RE(),AE(),KEY
C     RE()          RE() is a vector of length N such that RE(I) is
C                   the maximum relative uncertainty in row I of
C                   the matrix A. The values of RE() must be between
C                   0 and 1. A minimum of 10*machine precision will
C                   be enforced.
C
C     AE()          AE() is a vector of length N such that AE(I) is
C                   the maximum absolute uncertainty in row I of
C                   the matrix A. The values of AE() must be greater
C                   than or equal to 0.
C
C     KEY           For ease of use, RE and AE may be input as either
C                   vectors or scalars. If a scalar is input, the algo-
C                   rithm will use that value for each column of A.
C                   The parameter KEY indicates whether scalars or
C                   vectors are being input.
C                        KEY=0     RE scalar  AE scalar
C                        KEY=1     RE vector  AE scalar
C                        KEY=2     RE scalar  AE vector
C                        KEY=3     RE vector  AE vector
C
C
C     MODE          The integer MODE indicates how the routine
C                   is to react if rank deficiency is detected.
C                   If MODE = 0 return immediately, no solution
C                             1 compute truncated solution
C                             2 compute minimal length least squares sol
C                   The inexperienced user is advised to set MODE=0
C
C     NP            The first NP rows of A will not be interchanged
C                   with other rows even though the pivot strategy
C                   would suggest otherwise.
C                   The inexperienced user is advised to set NP=0.
C
C     WORK()        A real work array dimensioned 5*M.  However, if
C                   RE or AE have been specified as vectors, dimension
C                   WORK 4*M. If both RE and AE have been specified
C                   as vectors, dimension WORK 3*M.
C
C     LW            Actual dimension of WORK
C
C     IWORK()       Integer work array dimensioned at least N+M.
C
C     LIW           Actual dimension of IWORK.
C
C
C     INFO          Is a flag which provides for the efficient
C                   solution of subsequent problems involving the
C                   same A but different B.
C                   If INFO = 0 original call
C                      INFO = 1 subsequent calls
C                   On subsequent calls, the user must supply A, KRANK,
C                   LW, IWORK, LIW, and the first 2*M locations of WORK
C                   as output by the original call to DULSIA. MODE must
C                   be equal to the value of MODE in the original call.
C                   If MODE.LT.2, only the first N locations of WORK
C                   are accessed. AE, RE, KEY, and NP are not accessed.
C
C
C
C
C     Output..All TYPE REAL variables are DOUBLE PRECISION
C
C     A(,)          Contains the lower triangular part of the reduced
C                   matrix and the transformation information. It togeth
C                   with the first M elements of WORK (see below)
C                   completely specify the LQ factorization of A.
C
C     B(,)          Contains the N by NB solution matrix for X.
C
C     KRANK,KSURE   The numerical rank of A,  based upon the relative
C                   and absolute bounds on uncertainty, is bounded
C                   above by KRANK and below by KSURE. The algorithm
C                   returns a solution based on KRANK. KSURE provides
C                   an indication of the precision of the rank.
C
C     RNORM()       Contains the Euclidean length of the NB residual
C                   vectors  B(I)-AX(I), I=1,NB. If the matrix A is of
C                   full rank, then RNORM=0.0.
C
C     WORK()        The first M locations of WORK contain values
C                   necessary to reproduce the Householder
C                   transformation.
C
C     IWORK()       The first N locations contain the order in
C                   which the columns of A were used. The next
C                   M locations contain the order in which the
C                   rows of A were used.
C
C     INFO          Flag to indicate status of computation on completion
C                  -1   Parameter error(s)
C                   0 - Rank deficient, no solution
C                   1 - Rank deficient, truncated solution
C                   2 - Rank deficient, minimal length least squares sol
C                   3 - Numerical rank 0, zero solution
C                   4 - Rank .LT. NP
C                   5 - Full rank
C
C***REFERENCES  T. Manteuffel, An interval analysis approach to rank
C                 determination in linear least squares problems,
C                 Report SAND80-0655, Sandia Laboratories, June 1980.
C***ROUTINES CALLED  D1MACH, DU11US, DU12US, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   810801  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891009  Removed unreferenced variable.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Fixed an error message.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DULSIA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION D1MACH
      DIMENSION A(MDA,*),B(MDB,*),RE(*),AE(*),RNORM(*),W(*)
      INTEGER IWORK(*)
C
C***FIRST EXECUTABLE STATEMENT  DULSIA
      IF(INFO.LT.0 .OR. INFO.GT.1) GO TO 514
      IT=INFO
      INFO=-1
      IF(NB.EQ.0 .AND. IT.EQ.1) GO TO 501
      IF(M.LT.1) GO TO 502
      IF(N.LT.1) GO TO 503
      IF(N.LT.M) GO TO 504
      IF(MDA.LT.M) GO TO 505
      IF(LIW.LT.M+N) GO TO 506
      IF(MODE.LT.0 .OR. MODE.GT.3) GO TO 515
      IF(NB.EQ.0) GO TO 4
      IF(NB.LT.0) GO TO 507
      IF(MDB.LT.N) GO TO 508
      IF(IT.EQ.0) GO TO 4
      GO TO 400
    4 IF(KEY.LT.0.OR.KEY.GT.3) GO TO 509
      IF(KEY.EQ.0 .AND. LW.LT.5*M) GO TO 510
      IF(KEY.EQ.1 .AND. LW.LT.4*M) GO TO 510
      IF(KEY.EQ.2 .AND. LW.LT.4*M) GO TO 510
      IF(KEY.EQ.3 .AND. LW.LT.3*M) GO TO 510
      IF(NP.LT.0 .OR. NP.GT.M) GO TO 516
C
      EPS=10.*D1MACH(3)
      M1=1
      M2=M1+M
      M3=M2+M
      M4=M3+M
      M5=M4+M
C
      IF(KEY.EQ.1) GO TO 100
      IF(KEY.EQ.2) GO TO 200
      IF(KEY.EQ.3) GO TO 300
C
      IF(RE(1).LT.0.D00) GO TO 511
      IF(RE(1).GT.1.0D0) GO TO 512
      IF(RE(1).LT.EPS) RE(1)=EPS
      IF(AE(1).LT.0.0D0) GO TO 513
      DO 20 I=1,M
      W(M4-1+I)=RE(1)
      W(M5-1+I)=AE(1)
   20 CONTINUE
      CALL DU11US(A,MDA,M,N,W(M4),W(M5),MODE,NP,KRANK,KSURE,
     1            W(M1),W(M2),W(M3),IWORK(M1),IWORK(M2))
      GO TO 400
C
  100 CONTINUE
      IF(AE(1).LT.0.0D0) GO TO 513
      DO 120 I=1,M
      IF(RE(I).LT.0.0D0) GO TO 511
      IF(RE(I).GT.1.0D0) GO TO 512
      IF(RE(I).LT.EPS) RE(I)=EPS
      W(M4-1+I)=AE(1)
  120 CONTINUE
      CALL DU11US(A,MDA,M,N,RE,W(M4),MODE,NP,KRANK,KSURE,
     1            W(M1),W(M2),W(M3),IWORK(M1),IWORK(M2))
      GO TO 400
C
  200 CONTINUE
      IF(RE(1).LT.0.0D0) GO TO 511
      IF(RE(1).GT.1.0D0) GO TO 512
      IF(RE(1).LT.EPS) RE(1)=EPS
      DO 220 I=1,M
      W(M4-1+I)=RE(1)
      IF(AE(I).LT.0.0D0) GO TO 513
  220 CONTINUE
      CALL DU11US(A,MDA,M,N,W(M4),AE,MODE,NP,KRANK,KSURE,
     1            W(M1),W(M2),W(M3),IWORK(M1),IWORK(M2))
      GO TO 400
C
  300 CONTINUE
      DO 320 I=1,M
      IF(RE(I).LT.0.0D0) GO TO 511
      IF(RE(I).GT.1.0D0) GO TO 512
      IF(RE(I).LT.EPS) RE(I)=EPS
      IF(AE(I).LT.0.0D0) GO TO 513
  320 CONTINUE
      CALL DU11US(A,MDA,M,N,RE,AE,MODE,NP,KRANK,KSURE,
     1            W(M1),W(M2),W(M3),IWORK(M1),IWORK(M2))
C
C     DETERMINE INFO
C
  400 IF(KRANK.NE.M) GO TO 402
          INFO=5
          GO TO 410
  402 IF(KRANK.NE.0) GO TO 404
          INFO=3
          GO TO 410
  404 IF(KRANK.GE.NP) GO TO 406
          INFO=4
          RETURN
  406 INFO=MODE
      IF(MODE.EQ.0) RETURN
  410 IF(NB.EQ.0) RETURN
C
C
C     SOLUTION PHASE
C
      M1=1
      M2=M1+M
      M3=M2+M
      IF(INFO.EQ.2) GO TO 420
      IF(LW.LT.M2-1) GO TO 510
      CALL DU12US(A,MDA,M,N,B,MDB,NB,MODE,KRANK,
     1            RNORM,W(M1),W(M1),IWORK(M1),IWORK(M2))
      RETURN
C
  420 IF(LW.LT.M3-1) GO TO 510
      CALL DU12US(A,MDA,M,N,B,MDB,NB,MODE,KRANK,
     1            RNORM,W(M1),W(M2),IWORK(M1),IWORK(M2))
      RETURN
C
C     ERROR MESSAGES
C
  501 CALL XERMSG ('SLATEC', 'DULSIA',
     +   'SOLUTION ONLY (INFO=1) BUT NO RIGHT HAND SIDE (NB=0)', 1, 0)
      RETURN
  502 CALL XERMSG ('SLATEC', 'DULSIA', 'M.LT.1', 2, 1)
      RETURN
  503 CALL XERMSG ('SLATEC', 'DULSIA', 'N.LT.1', 2, 1)
      RETURN
  504 CALL XERMSG ('SLATEC', 'DULSIA', 'N.LT.M', 2, 1)
      RETURN
  505 CALL XERMSG ('SLATEC', 'DULSIA', 'MDA.LT.M', 2, 1)
      RETURN
  506 CALL XERMSG ('SLATEC', 'DULSIA', 'LIW.LT.M+N', 2, 1)
      RETURN
  507 CALL XERMSG ('SLATEC', 'DULSIA', 'NB.LT.0', 2, 1)
      RETURN
  508 CALL XERMSG ('SLATEC', 'DULSIA', 'MDB.LT.N', 2, 1)
      RETURN
  509 CALL XERMSG ('SLATEC', 'DULSIA', 'KEY OUT OF RANGE', 2, 1)
      RETURN
  510 CALL XERMSG ('SLATEC', 'DULSIA', 'INSUFFICIENT WORK SPACE', 8, 1)
      INFO=-1
      RETURN
  511 CALL XERMSG ('SLATEC', 'DULSIA', 'RE(I) .LT. 0', 2, 1)
      RETURN
  512 CALL XERMSG ('SLATEC', 'DULSIA', 'RE(I) .GT. 1', 2, 1)
      RETURN
  513 CALL XERMSG ('SLATEC', 'DULSIA', 'AE(I) .LT. 0', 2, 1)
      RETURN
  514 CALL XERMSG ('SLATEC', 'DULSIA', 'INFO OUT OF RANGE', 2, 1)
      RETURN
  515 CALL XERMSG ('SLATEC', 'DULSIA', 'MODE OUT OF RANGE', 2, 1)
      RETURN
  516 CALL XERMSG ('SLATEC', 'DULSIA', 'NP OUT OF RANGE', 2, 1)
      RETURN
      END
