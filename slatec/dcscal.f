*DECK DCSCAL
      SUBROUTINE DCSCAL (A, NRDA, NROW, NCOL, COLS, COLSAV, ROWS,
     +   ROWSAV, ANORM, SCALES, ISCALE, IC)
C***BEGIN PROLOGUE  DCSCAL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP and DSUDS
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (CSCALE-S, DCSCAL-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     This routine scales the matrix A by columns when needed.
C
C***SEE ALSO  DBVSUP, DSUDS
C***ROUTINES CALLED  DDOT
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DCSCAL
      DOUBLE PRECISION DDOT
      INTEGER IC, IP, ISCALE, J, K, NCOL, NRDA, NROW
      DOUBLE PRECISION A(NRDA,*), ALOG2, ANORM, ASCALE, COLS(*),
     1     COLSAV(*), CS, P, ROWS(*), ROWSAV(*), S,
     2     SCALES(*), TEN20, TEN4
C
      SAVE TEN4, TEN20
      DATA TEN4,TEN20 /1.0D4,1.0D20/
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 130
C        BEGIN BLOCK PERMITTING ...EXITS TO 60
C***FIRST EXECUTABLE STATEMENT  DCSCAL
            IF (ISCALE .NE. (-1)) GO TO 40
C
               IF (IC .EQ. 0) GO TO 20
                  DO 10 K = 1, NCOL
                     COLS(K) = DDOT(NROW,A(1,K),1,A(1,K),1)
   10             CONTINUE
   20          CONTINUE
C
               ASCALE = ANORM/NCOL
               DO 30 K = 1, NCOL
                  CS = COLS(K)
C        .........EXIT
                  IF ((CS .GT. TEN4*ASCALE) .OR. (TEN4*CS .LT. ASCALE))
     1               GO TO 60
C        .........EXIT
                  IF ((CS .LT. 1.0D0/TEN20) .OR. (CS .GT. TEN20))
     1               GO TO 60
   30          CONTINUE
   40       CONTINUE
C
            DO 50 K = 1, NCOL
               SCALES(K) = 1.0D0
   50       CONTINUE
C     ......EXIT
            GO TO 130
   60    CONTINUE
C
         ALOG2 = LOG(2.0D0)
         ANORM = 0.0D0
         DO 110 K = 1, NCOL
            CS = COLS(K)
            IF (CS .NE. 0.0D0) GO TO 70
               SCALES(K) = 1.0D0
            GO TO 100
   70       CONTINUE
               P = LOG(CS)/ALOG2
               IP = -0.5D0*P
               S = 2.0D0**IP
               SCALES(K) = S
               IF (IC .EQ. 1) GO TO 80
                  COLS(K) = S*S*COLS(K)
                  ANORM = ANORM + COLS(K)
                  COLSAV(K) = COLS(K)
   80          CONTINUE
               DO 90 J = 1, NROW
                  A(J,K) = S*A(J,K)
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
C
C     ...EXIT
         IF (IC .EQ. 0) GO TO 130
C
         DO 120 K = 1, NROW
            ROWS(K) = DDOT(NCOL,A(K,1),NRDA,A(K,1),NRDA)
            ROWSAV(K) = ROWS(K)
            ANORM = ANORM + ROWS(K)
  120    CONTINUE
  130 CONTINUE
      RETURN
      END
