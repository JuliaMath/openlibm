*DECK SCOV
      SUBROUTINE SCOV (FCN, IOPT, M, N, X, FVEC, R, LDR, INFO, WA1, WA2,
     +   WA3, WA4)
C***BEGIN PROLOGUE  SCOV
C***PURPOSE  Calculate the covariance matrix for a nonlinear data
C            fitting problem.  It is intended to be used after a
C            successful return from either SNLS1 or SNLS1E.
C***LIBRARY   SLATEC
C***CATEGORY  K1B1
C***TYPE      SINGLE PRECISION (SCOV-S, DCOV-D)
C***KEYWORDS  COVARIANCE MATRIX, NONLINEAR DATA FITTING,
C             NONLINEAR LEAST SQUARES
C***AUTHOR  Hiebert, K. L., (SNLA)
C***DESCRIPTION
C
C  1. Purpose.
C
C     SCOV calculates the covariance matrix for a nonlinear data
C     fitting problem.  It is intended to be used after a
C     successful return from either SNLS1 or SNLS1E. SCOV
C     and SNLS1 (and SNLS1E) have compatible parameters.  The
C     required external subroutine, FCN, is the same
C     for all three codes, SCOV, SNLS1, and SNLS1E.
C
C  2. Subroutine and Type Statements.
C
C     SUBROUTINE SCOV(FCN,IOPT,M,N,X,FVEC,R,LDR,INFO,
C                     WA1,WA2,WA3,WA4)
C     INTEGER IOPT,M,N,LDR,INFO
C     REAL X(N),FVEC(M),R(LDR,N),WA1(N),WA2(N),WA3(N),WA4(M)
C     EXTERNAL FCN
C
C  3. Parameters.
C
C       FCN is the name of the user-supplied subroutine which calculates
C         the functions.  If the user wants to supply the Jacobian
C         (IOPT=2 or 3), then FCN must be written to calculate the
C         Jacobian, as well as the functions.  See the explanation
C         of the IOPT argument below.  FCN must be declared in an
C         EXTERNAL statement in the calling program and should be
C         written as follows.
C
C         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C         INTEGER IFLAG,LDFJAC,M,N
C         REAL X(N),FVEC(M)
C         ----------
C         FJAC and LDFJAC may be ignored     , if IOPT=1.
C         REAL FJAC(LDFJAC,N)                , if IOPT=2.
C         REAL FJAC(N)                       , if IOPT=3.
C         ----------
C           IFLAG will never be zero when FCN is called by SCOV.
C         RETURN
C         ----------
C           If IFLAG=1, calculate the functions at X and return
C           this vector in FVEC.
C         RETURN
C         ----------
C           If IFLAG=2, calculate the full Jacobian at X and return
C           this matrix in FJAC.  Note that IFLAG will never be 2 unless
C           IOPT=2.  FVEC contains the function values at X and must
C           not be altered.  FJAC(I,J) must be set to the derivative
C           of FVEC(I) with respect to X(J).
C         RETURN
C         ----------
C           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian
C           and return this vector in FJAC.  Note that IFLAG will
C           never be 3 unless IOPT=3.  FJAC(J) must be set to
C           the derivative of FVEC(LDFJAC) with respect to X(J).
C         RETURN
C         ----------
C         END
C
C
C         The value of IFLAG should not be changed by FCN unless the
C         user wants to terminate execution of SCOV.  In this case, set
C         IFLAG to a negative integer.
C
C
C    IOPT is an input variable which specifies how the Jacobian will
C         be calculated.  If IOPT=2 or 3, then the user must supply the
C         Jacobian, as well as the function values, through the
C         subroutine FCN.  If IOPT=2, the user supplies the full
C         Jacobian with one call to FCN.  If IOPT=3, the user supplies
C         one row of the Jacobian with each call.  (In this manner,
C         storage can be saved because the full Jacobian is not stored.)
C         If IOPT=1, the code will approximate the Jacobian by forward
C         differencing.
C
C       M is a positive integer input variable set to the number of
C         functions.
C
C       N is a positive integer input variable set to the number of
C         variables.  N must not exceed M.
C
C       X is an array of length N.  On input X must contain the value
C         at which the covariance matrix is to be evaluated.  This is
C         usually the value for X returned from a successful run of
C         SNLS1 (or SNLS1E).  The value of X will not be changed.
C
C    FVEC is an output array of length M which contains the functions
C         evaluated at X.
C
C       R is an output array.  For IOPT=1 and 2, R is an M by N array.
C         For IOPT=3, R is an N by N array.  On output, if INFO=1,
C         the upper N by N submatrix of R contains the covariance
C         matrix evaluated at X.
C
C     LDR is a positive integer input variable which specifies
C         the leading dimension of the array R.  For IOPT=1 and 2,
C         LDR must not be less than M.  For IOPT=3, LDR must not
C         be less than N.
C
C    INFO is an integer output variable.  If the user has terminated
C         execution, INFO is set to the (negative) value of IFLAG.  See
C         description of FCN. Otherwise, INFO is set as follows.
C
C         INFO = 0 Improper input parameters (M.LE.0 or N.LE.0).
C
C         INFO = 1 Successful return.  The covariance matrix has been
C                  calculated and stored in the upper N by N
C                  submatrix of R.
C
C         INFO = 2 The Jacobian matrix is singular for the input value
C                  of X.  The covariance matrix cannot be calculated.
C                  The upper N by N submatrix of R contains the QR
C                  factorization of the Jacobian (probably not of
C                  interest to the user).
C
C     WA1 is a work array of length N.
C     WA2 is a work array of length N.
C     WA3 is a work array of length N.
C     WA4 is a work array of length M.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ENORM, FDJAC3, QRFAC, RWUPDT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   810522  DATE WRITTEN
C   890505  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Fixed an error message.  (RWC)
C***END PROLOGUE  SCOV
C
C     REVISED 820707-1100
C     REVISED YYMMDD HHMM
C
      INTEGER I,IDUM,IFLAG,INFO,IOPT,J,K,KP1,LDR,M,N,NM1,NROW
      REAL X(*),R(LDR,*),FVEC(*),WA1(*),WA2(*),WA3(*),WA4(*)
      EXTERNAL FCN
      REAL ONE,SIGMA,TEMP,ZERO
      LOGICAL SING
      SAVE ZERO, ONE
      DATA ZERO/0.E0/,ONE/1.E0/
C***FIRST EXECUTABLE STATEMENT  SCOV
      SING=.FALSE.
      IFLAG=0
      IF (M.LE.0 .OR. N.LE.0) GO TO 300
C
C     CALCULATE SIGMA = (SUM OF THE SQUARED RESIDUALS) / (M-N)
      IFLAG=1
      CALL FCN(IFLAG,M,N,X,FVEC,R,LDR)
      IF (IFLAG.LT.0) GO TO 300
      TEMP=ENORM(M,FVEC)
      SIGMA=ONE
      IF (M.NE.N) SIGMA=TEMP*TEMP/(M-N)
C
C     CALCULATE THE JACOBIAN
      IF (IOPT.EQ.3) GO TO 200
C
C     STORE THE FULL JACOBIAN USING M*N STORAGE
      IF (IOPT.EQ.1) GO TO 100
C
C     USER SUPPLIES THE JACOBIAN
      IFLAG=2
      CALL FCN(IFLAG,M,N,X,FVEC,R,LDR)
      GO TO 110
C
C     CODE APPROXIMATES THE JACOBIAN
100   CALL FDJAC3(FCN,M,N,X,FVEC,R,LDR,IFLAG,ZERO,WA4)
110   IF (IFLAG.LT.0) GO TO 300
C
C     COMPUTE THE QR DECOMPOSITION
      CALL QRFAC(M,N,R,LDR,.FALSE.,IDUM,1,WA1,WA1,WA1)
      DO 120 I=1,N
120   R(I,I)=WA1(I)
      GO TO 225
C
C     COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX CALCULATED ONE
C     ROW AT A TIME AND STORED IN THE UPPER TRIANGLE OF R.
C     ( (Q TRANSPOSE)*FVEC IS ALSO CALCULATED BUT NOT USED.)
200   CONTINUE
      DO 210 J=1,N
      WA2(J)=ZERO
      DO 205 I=1,N
      R(I,J)=ZERO
205   CONTINUE
210   CONTINUE
      IFLAG=3
      DO 220 I=1,M
      NROW = I
      CALL FCN(IFLAG,M,N,X,FVEC,WA1,NROW)
      IF (IFLAG.LT.0) GO TO 300
      TEMP=FVEC(I)
      CALL RWUPDT(N,R,LDR,WA1,WA2,TEMP,WA3,WA4)
220   CONTINUE
C
C     CHECK IF R IS SINGULAR.
225   CONTINUE
      DO 230 I=1,N
      IF (R(I,I).EQ.ZERO) SING=.TRUE.
230   CONTINUE
      IF (SING) GO TO 300
C
C     R IS UPPER TRIANGULAR.  CALCULATE (R TRANSPOSE) INVERSE AND STORE
C     IN THE UPPER TRIANGLE OF R.
      IF (N.EQ.1) GO TO 275
      NM1=N-1
      DO 270 K=1,NM1
C
C     INITIALIZE THE RIGHT-HAND SIDE (WA1(*)) AS THE K-TH COLUMN OF THE
C     IDENTITY MATRIX.
      DO 240 I=1,N
      WA1(I)=ZERO
240   CONTINUE
      WA1(K)=ONE
C
      R(K,K)=WA1(K)/R(K,K)
      KP1=K+1
      DO 260 I=KP1,N
C
C     SUBTRACT R(K,I-1)*R(I-1,*) FROM THE RIGHT-HAND SIDE, WA1(*).
      DO 250 J=I,N
      WA1(J)=WA1(J)-R(K,I-1)*R(I-1,J)
250   CONTINUE
      R(K,I)=WA1(I)/R(I,I)
260   CONTINUE
270   CONTINUE
275   R(N,N)=ONE/R(N,N)
C
C     CALCULATE R-INVERSE * (R TRANSPOSE) INVERSE AND STORE IN THE UPPER
C     TRIANGLE OF R.
      DO 290 I=1,N
      DO 290 J=I,N
      TEMP=ZERO
      DO 280 K=J,N
      TEMP=TEMP+R(I,K)*R(J,K)
280   CONTINUE
      R(I,J)=TEMP*SIGMA
290   CONTINUE
      INFO=1
C
300   CONTINUE
      IF (M.LE.0 .OR. N.LE.0) INFO=0
      IF (IFLAG.LT.0) INFO=IFLAG
      IF (SING) INFO=2
      IF (INFO .LT. 0) CALL XERMSG ('SLATEC', 'SCOV',
     +   'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'SCOV',
     +   'INVALID INPUT PARAMETER.', 2, 1)
      IF (INFO .EQ. 2) CALL XERMSG ('SLATEC', 'SCOV',
     +   'SINGULAR JACOBIAN MATRIX, COVARIANCE MATRIX CANNOT BE ' //
     +   'CALCULATED.', 1, 1)
      RETURN
      END
