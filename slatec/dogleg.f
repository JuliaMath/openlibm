*DECK DOGLEG
      SUBROUTINE DOGLEG (N, R, LR, DIAG, QTB, DELTA, X, WA1, WA2)
C***BEGIN PROLOGUE  DOGLEG
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SNSQ and SNSQE
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (DOGLEG-S, DDOGLG-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Given an M by N matrix A, an N by N nonsingular DIAGONAL
C     matrix D, an M-vector B, and a positive number DELTA, the
C     problem is to determine the convex combination X of the
C     Gauss-Newton and scaled gradient directions that minimizes
C     (A*X - B) in the least squares sense, subject to the
C     restriction that the Euclidean norm of D*X be at most DELTA.
C
C     This subroutine completes the solution of the problem
C     if it is provided with the necessary information from the
C     QR factorization of A. That is, if A = Q*R, where Q has
C     orthogonal columns and R is an upper triangular matrix,
C     then DOGLEG expects the full upper triangle of R and
C     the first N components of (Q TRANSPOSE)*B.
C
C     The subroutine statement is
C
C       SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
C
C     where
C
C       N is a positive integer input variable set to the order of R.
C
C       R is an input array of length LR which must contain the upper
C         triangular matrix R stored by rows.
C
C       LR is a positive integer input variable not less than
C         (N*(N+1))/2.
C
C       DIAG is an input array of length N which must contain the
C         diagonal elements of the matrix D.
C
C       QTB is an input array of length N which must contain the first
C         N elements of the vector (Q TRANSPOSE)*B.
C
C       DELTA is a positive input variable which specifies an upper
C         bound on the Euclidean norm of D*X.
C
C       X is an output array of length N which contains the desired
C         convex combination of the Gauss-Newton direction and the
C         scaled gradient direction.
C
C       WA1 and WA2 are work arrays of length N.
C
C***SEE ALSO  SNSQ, SNSQE
C***ROUTINES CALLED  ENORM, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DOGLEG
      INTEGER N,LR
      REAL DELTA
      REAL R(LR),DIAG(*),QTB(*),X(*),WA1(*),WA2(*)
      INTEGER I,J,JJ,JP1,K,L
      REAL ALPHA,BNORM,EPSMCH,GNORM,ONE,QNORM,SGNORM,SUM,TEMP,ZERO
      REAL R1MACH,ENORM
      SAVE ONE, ZERO
      DATA ONE,ZERO /1.0E0,0.0E0/
C***FIRST EXECUTABLE STATEMENT  DOGLEG
      EPSMCH = R1MACH(4)
C
C     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.
C
      JJ = (N*(N + 1))/2 + 1
      DO 50 K = 1, N
         J = N - K + 1
         JP1 = J + 1
         JJ = JJ - K
         L = JJ + 1
         SUM = ZERO
         IF (N .LT. JP1) GO TO 20
         DO 10 I = JP1, N
            SUM = SUM + R(L)*X(I)
            L = L + 1
   10       CONTINUE
   20    CONTINUE
         TEMP = R(JJ)
         IF (TEMP .NE. ZERO) GO TO 40
         L = J
         DO 30 I = 1, J
            TEMP = MAX(TEMP,ABS(R(L)))
            L = L + N - I
   30       CONTINUE
         TEMP = EPSMCH*TEMP
         IF (TEMP .EQ. ZERO) TEMP = EPSMCH
   40    CONTINUE
         X(J) = (QTB(J) - SUM)/TEMP
   50    CONTINUE
C
C     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.
C
      DO 60 J = 1, N
         WA1(J) = ZERO
         WA2(J) = DIAG(J)*X(J)
   60    CONTINUE
      QNORM = ENORM(N,WA2)
      IF (QNORM .LE. DELTA) GO TO 140
C
C     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
C     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
C
      L = 1
      DO 80 J = 1, N
         TEMP = QTB(J)
         DO 70 I = J, N
            WA1(I) = WA1(I) + R(L)*TEMP
            L = L + 1
   70       CONTINUE
         WA1(J) = WA1(J)/DIAG(J)
   80    CONTINUE
C
C     CALCULATE THE NORM OF THE SCALED GRADIENT DIRECTION,
C     NORMALIZE, AND RESCALE THE GRADIENT.
C
      GNORM = ENORM(N,WA1)
      SGNORM = ZERO
      ALPHA = DELTA/QNORM
      IF (GNORM .EQ. ZERO) GO TO 120
      DO 90 J = 1, N
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)
   90    CONTINUE
C
C     CALCULATE THE POINT ALONG THE SCALED GRADIENT
C     AT WHICH THE QUADRATIC IS MINIMIZED.
C
      L = 1
      DO 110 J = 1, N
         SUM = ZERO
         DO 100 I = J, N
            SUM = SUM + R(L)*WA1(I)
            L = L + 1
  100       CONTINUE
         WA2(J) = SUM
  110    CONTINUE
      TEMP = ENORM(N,WA2)
      SGNORM = (GNORM/TEMP)/TEMP
C
C     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
C
      ALPHA = ZERO
      IF (SGNORM .GE. DELTA) GO TO 120
C
C     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
C     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
C     AT WHICH THE QUADRATIC IS MINIMIZED.
C
      BNORM = ENORM(N,QTB)
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
      TEMP = TEMP - (DELTA/QNORM)*(SGNORM/DELTA)**2
     1       + SQRT((TEMP-(DELTA/QNORM))**2
     2              +(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))
      ALPHA = ((DELTA/QNORM)*(ONE - (SGNORM/DELTA)**2))/TEMP
  120 CONTINUE
C
C     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
C     DIRECTION AND THE SCALED GRADIENT DIRECTION.
C
      TEMP = (ONE - ALPHA)*MIN(SGNORM,DELTA)
      DO 130 J = 1, N
         X(J) = TEMP*WA1(J) + ALPHA*X(J)
  130    CONTINUE
  140 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE DOGLEG.
C
      END
