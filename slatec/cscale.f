*DECK CSCALE
      SUBROUTINE CSCALE (A, NRDA, NROW, NCOL, COLS, COLSAV, ROWS,
     +   ROWSAV, ANORM, SCALES, ISCALE, IC)
C***BEGIN PROLOGUE  CSCALE
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (CSCALE-S, DCSCAL-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     This routine scales the matrix A by columns when needed
C
C***SEE ALSO  BVSUP
C***ROUTINES CALLED  SDOT
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  CSCALE
      DIMENSION A(NRDA,*),COLS(*),COLSAV(*),SCALES(*),
     1          ROWS(*),ROWSAV(*)
C
      SAVE TEN4, TEN20
      DATA TEN4,TEN20/1.E+4,1.E+20/
C
C***FIRST EXECUTABLE STATEMENT  CSCALE
      IF (ISCALE .NE. (-1)) GO TO 25
C
      IF (IC .EQ. 0) GO TO 10
      DO 5 K=1,NCOL
    5    COLS(K)=SDOT(NROW,A(1,K),1,A(1,K),1)
C
   10 ASCALE=ANORM/NCOL
      DO 20 K=1,NCOL
         CS=COLS(K)
         IF ((CS .GT. TEN4*ASCALE) .OR. (TEN4*CS .LT. ASCALE)) GO TO 50
         IF ((CS .LT. 1./TEN20) .OR. (CS .GT. TEN20)) GO TO 50
   20 CONTINUE
C
   25 DO 30 K=1,NCOL
   30    SCALES(K)=1.
      RETURN
C
   50 ALOG2=LOG(2.)
      ANORM=0.
      DO 100 K=1,NCOL
         CS=COLS(K)
         IF (CS .NE. 0.) GO TO 60
         SCALES(K)=1.
         GO TO 100
   60    P=LOG(CS)/ALOG2
         IP=-0.5*P
         S=2.**IP
         SCALES(K)=S
         IF (IC .EQ. 1) GO TO 70
         COLS(K)=S*S*COLS(K)
         ANORM=ANORM+COLS(K)
         COLSAV(K)=COLS(K)
   70    DO 80 J=1,NROW
   80       A(J,K)=S*A(J,K)
  100 CONTINUE
C
      IF (IC .EQ. 0) RETURN
C
      DO 200 K=1,NROW
         ROWS(K)=SDOT(NCOL,A(K,1),NRDA,A(K,1),NRDA)
         ROWSAV(K)=ROWS(K)
  200    ANORM=ANORM+ROWS(K)
      RETURN
      END
