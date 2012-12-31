*DECK SGLSS
      SUBROUTINE SGLSS (A, MDA, M, N, B, MDB, NB, RNORM, WORK, LW,
     +   IWORK, LIW, INFO)
C***BEGIN PROLOGUE  SGLSS
C***PURPOSE  Solve a linear least squares problems by performing a QR
C            factorization of the matrix using Householder
C            transformations.  Emphasis is put on detecting possible
C            rank deficiency.
C***LIBRARY   SLATEC
C***CATEGORY  D9, D5
C***TYPE      SINGLE PRECISION (SGLSS-S, DGLSS-D)
C***KEYWORDS  LINEAR LEAST SQUARES, LQ FACTORIZATION, QR FACTORIZATION,
C             UNDERDETERMINED LINEAR SYSTEMS
C***AUTHOR  Manteuffel, T. A., (LANL)
C***DESCRIPTION
C
C     SGLSS solves both underdetermined and overdetermined
C     LINEAR systems AX = B, where A is an M by N matrix
C     and B is an M by NB matrix of right hand sides. If
C     M.GE.N, the least squares solution is computed by
C     decomposing the matrix A into the product of an
C     orthogonal matrix Q and an upper triangular matrix
C     R (QR factorization). If M.LT.N, the minimal
C     length solution is computed by factoring the
C     matrix A into the product of a lower triangular
C     matrix L and an orthogonal matrix Q (LQ factor-
C     ization). If the matrix A is determined to be rank
C     deficient, that is the rank of A is less than
C     MIN(M,N), then the minimal length least squares
C     solution is computed.
C
C     SGLSS assumes full machine precision in the data.
C     If more control over the uncertainty in the data
C     is desired, the codes LLSIA and ULSIA are
C     recommended.
C
C     SGLSS requires MDA*N + (MDB + 1)*NB + 5*MIN(M,N) dimensioned
C     real space and M+N dimensioned integer space.
C
C
C   ******************************************************************
C   *                                                                *
C   *         WARNING - All input arrays are changed on exit.        *
C   *                                                                *
C   ******************************************************************
C     SUBROUTINE SGLSS(A,MDA,M,N,B,MDB,NB,RNORM,WORK,LW,IWORK,LIW,INFO)
C
C     Input..
C
C     A(,)          Linear coefficient matrix of AX=B, with MDA the
C      MDA,M,N      actual first dimension of A in the calling program.
C                   M is the row dimension (no. of EQUATIONS of the
C                   problem) and N the col dimension (no. of UNKNOWNS).
C
C     B(,)          Right hand side(s), with MDB the actual first
C      MDB,NB       dimension of B in the calling program. NB is the
C                   number of M by 1 right hand sides. Must have
C                   MDB.GE.MAX(M,N). If NB = 0, B is never accessed.
C
C
C     RNORM()       Vector of length at least NB.  On input the contents
C                   of RNORM are unused.
C
C     WORK()        A real work array dimensioned 5*MIN(M,N).
C
C     LW            Actual dimension of WORK.
C
C     IWORK()       Integer work array dimensioned at least N+M.
C
C     LIW           Actual dimension of IWORK.
C
C
C     INFO          A flag which provides for the efficient
C                   solution of subsequent problems involving the
C                   same A but different B.
C                   If INFO = 0 original call
C                      INFO = 1 subsequent calls
C                   On subsequent calls, the user must supply A, INFO,
C                   LW, IWORK, LIW, and the first 2*MIN(M,N) locations
C                   of WORK as output by the original call to SGLSS.
C
C
C     Output..
C
C     A(,)          Contains the triangular part of the reduced matrix
C                   and the transformation information. It together with
C                   the first 2*MIN(M,N) elements of WORK (see below)
C                   completely specify the factorization of A.
C
C     B(,)          Contains the N by NB solution matrix X.
C
C
C     RNORM()       Contains the Euclidean length of the NB residual
C                   vectors  B(I)-AX(I), I=1,NB.
C
C     WORK()        The first 2*MIN(M,N) locations of WORK contain value
C                   necessary to reproduce the factorization of A.
C
C     IWORK()       The first M+N locations contain the order in
C                   which the rows and columns of A were used.
C                   If M.GE.N columns then rows. If M.LT.N rows
C                   then columns.
C
C     INFO          Flag to indicate status of computation on completion
C                  -1   Parameter error(s)
C                   0 - Full rank
C                   N.GT.0 - Reduced rank  rank=MIN(M,N)-INFO
C
C***REFERENCES  T. Manteuffel, An interval analysis approach to rank
C                 determination in linear least squares problems,
C                 Report SAND80-0655, Sandia Laboratories, June 1980.
C***ROUTINES CALLED  LLSIA, ULSIA
C***REVISION HISTORY  (YYMMDD)
C   810801  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SGLSS
      DIMENSION A(MDA,*),B(MDB,*),RNORM(*),WORK(*)
      INTEGER IWORK(*)
C
C***FIRST EXECUTABLE STATEMENT  SGLSS
      RE=0.
      AE=0.
      KEY=0
      MODE=2
      NP=0
C
C     IF M.GE.N CALL LLSIA
C     IF M.LT.N CALL ULSIA
C
      IF(M.LT.N) GO TO 10
      CALL LLSIA(A,MDA,M,N,B,MDB,NB,RE,AE,KEY,MODE,NP,
     1            KRANK,KSURE,RNORM,WORK,LW,IWORK,LIW,INFO)
      IF(INFO.EQ.-1) RETURN
      INFO=N-KRANK
      RETURN
   10 CALL ULSIA(A,MDA,M,N,B,MDB,NB,RE,AE,KEY,MODE,NP,
     1            KRANK,KSURE,RNORM,WORK,LW,IWORK,LIW,INFO)
      IF(INFO.EQ.-1) RETURN
      INFO=M-KRANK
      RETURN
      END
