*DECK ORTHOL
      SUBROUTINE ORTHOL (A, M, N, NRDA, IFLAG, IRANK, ISCALE, DIAG,
     +   KPIVOT, SCALES, COLS, CS)
C***BEGIN PROLOGUE  ORTHOL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (ORTHOL-S)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   Reduction of the matrix A to upper triangular form by a sequence of
C   orthogonal HOUSEHOLDER transformations pre-multiplying A
C
C   Modeled after the ALGOL codes in the articles in the REFERENCES
C   section.
C
C **********************************************************************
C   INPUT
C **********************************************************************
C
C     A -- Contains the matrix to be decomposed, must be dimensioned
C           NRDA by N
C     M -- Number of rows in the matrix, M greater or equal to N
C     N -- Number of columns in the matrix, N greater or equal to 1
C     IFLAG -- Indicates the uncertainty in the matrix data
C             = 0 when the data is to be treated as exact
C             =-K when the data is assumed to be accurate to about
C                 K digits
C     ISCALE -- Scaling indicator
C               =-1 if the matrix A is to be pre-scaled by
C               columns when appropriate.
C               Otherwise no scaling will be attempted
C     NRDA -- Row dimension of A, NRDA greater or equal to M
C     DIAG,KPIVOT,COLS -- Arrays of length at least n used internally
C         ,CS,SCALES
C
C **********************************************************************
C   OUTPUT
C **********************************************************************
C
C     IFLAG - Status indicator
C            =1 for successful decomposition
C            =2 if improper input is detected
C            =3 if rank of the matrix is less than N
C     A -- Contains the reduced matrix in the strictly upper triangular
C          part and transformation information in the lower part
C     IRANK -- Contains the numerically determined matrix rank
C     DIAG -- Contains the diagonal elements of the reduced
C             triangular matrix
C     KPIVOT -- Contains the pivotal information, the column
C               interchanges performed on the original matrix are
C               recorded here.
C     SCALES -- Contains the column scaling parameters
C
C **********************************************************************
C
C***SEE ALSO  BVSUP
C***REFERENCES  G. Golub, Numerical methods for solving linear least
C                 squares problems, Numerische Mathematik 7, (1965),
C                 pp. 206-216.
C               P. Businger and G. Golub, Linear least squares
C                 solutions by Householder transformations, Numerische
C                 Mathematik  7, (1965), pp. 269-276.
C***ROUTINES CALLED  CSCALE, R1MACH, SDOT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900402  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  ORTHOL
      DIMENSION A(NRDA,*),DIAG(*),KPIVOT(*),COLS(*),CS(*),SCALES(*)
C
C **********************************************************************
C
C     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
C     BY THE FUNCTION R1MACH.
C
C***FIRST EXECUTABLE STATEMENT  ORTHOL
      URO = R1MACH(3)
C
C **********************************************************************
C
      IF (M .GE. N  .AND.  N .GE. 1  .AND.  NRDA .GE. M) GO TO 1
      IFLAG=2
      CALL XERMSG ('SLATEC', 'ORTHOL', 'INVALID INPUT PARAMETERS.', 2,
     +   1)
      RETURN
C
    1 ACC=10.*URO
      IF (IFLAG .LT. 0) ACC=MAX(ACC,10.**IFLAG)
      SRURO=SQRT(URO)
      IFLAG=1
      IRANK=N
C
C     COMPUTE NORM**2 OF JTH COLUMN AND A MATRIX NORM
C
      ANORM=0.
      DO 2 J=1,N
         KPIVOT(J)=J
         COLS(J)=SDOT(M,A(1,J),1,A(1,J),1)
         CS(J)=COLS(J)
         ANORM=ANORM+COLS(J)
    2 CONTINUE
C
C     PERFORM COLUMN SCALING ON A WHEN SPECIFIED
C
      CALL CSCALE(A,NRDA,M,N,COLS,CS,DUM,DUM,ANORM,SCALES,ISCALE,0)
C
      ANORM=SQRT(ANORM)
C
C
C     CONSTRUCTION OF UPPER TRIANGULAR MATRIX AND RECORDING OF
C     ORTHOGONAL TRANSFORMATIONS
C
C
      DO 50 K=1,N
         MK=M-K+1
         IF (K .EQ. N) GO TO 25
         KP=K+1
C
C        SEARCHING FOR PIVOTAL COLUMN
C
         DO 10 J=K,N
            IF (COLS(J) .GE. SRURO*CS(J)) GO TO 5
            COLS(J)=SDOT(MK,A(K,J),1,A(K,J),1)
            CS(J)=COLS(J)
    5       IF (J .EQ. K) GO TO 7
            IF (SIGMA .GE. 0.99*COLS(J)) GO TO 10
    7       SIGMA=COLS(J)
            JCOL=J
   10    CONTINUE
         IF (JCOL .EQ. K) GO TO 25
C
C        PERFORM COLUMN INTERCHANGE
C
         L=KPIVOT(K)
         KPIVOT(K)=KPIVOT(JCOL)
         KPIVOT(JCOL)=L
         COLS(JCOL)=COLS(K)
         COLS(K)=SIGMA
         CSS=CS(K)
         CS(K)=CS(JCOL)
         CS(JCOL)=CSS
         SC=SCALES(K)
         SCALES(K)=SCALES(JCOL)
         SCALES(JCOL)=SC
         DO 20 L=1,M
            ASAVE=A(L,K)
            A(L,K)=A(L,JCOL)
   20       A(L,JCOL)=ASAVE
C
C        CHECK RANK OF THE MATRIX
C
   25    SIG=SDOT(MK,A(K,K),1,A(K,K),1)
         DIAGK=SQRT(SIG)
         IF (DIAGK .GT. ACC*ANORM) GO TO 30
C
C        RANK DEFICIENT PROBLEM
         IFLAG=3
         IRANK=K-1
         CALL XERMSG ('SLATEC', 'ORTHOL',
     +      'RANK OF MATRIX IS LESS THAN THE NUMBER OF COLUMNS.', 1, 1)
         RETURN
C
C        CONSTRUCT AND APPLY TRANSFORMATION TO MATRIX A
C
   30    AKK=A(K,K)
         IF (AKK .GT. 0.) DIAGK=-DIAGK
         DIAG(K)=DIAGK
         A(K,K)=AKK-DIAGK
         IF (K .EQ. N) GO TO 50
         SAD=DIAGK*AKK-SIG
         DO 40 J=KP,N
            AS=SDOT(MK,A(K,K),1,A(K,J),1)/SAD
            DO 35 L=K,M
   35          A(L,J)=A(L,J)+AS*A(L,K)
   40       COLS(J)=COLS(J)-A(K,J)**2
   50 CONTINUE
C
C
      RETURN
      END
