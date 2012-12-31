*DECK DORTHR
      SUBROUTINE DORTHR (A, N, M, NRDA, IFLAG, IRANK, ISCALE, DIAG,
     +   KPIVOT, SCALES, ROWS, RS)
C***BEGIN PROLOGUE  DORTHR
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP and DSUDS
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (ORTHOR-S, DORTHR-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   Reduction of the matrix A to lower triangular form by a sequence of
C   orthogonal HOUSEHOLDER transformations post-multiplying A.
C
C *********************************************************************
C   INPUT
C *********************************************************************
C
C     A -- Contains the matrix to be decomposed, must be dimensioned
C           NRDA by N.
C     N -- Number of rows in the matrix, N greater or equal to 1.
C     M -- Number of columns in the matrix, M greater or equal to N.
C     IFLAG -- Indicates the uncertainty in the matrix data.
C             = 0 when the data is to be treated as exact.
C             =-K when the data is assumed to be accurate to about
C                 K digits.
C     ISCALE -- Scaling indicator.
C               =-1 if the matrix is to be pre-scaled by
C               columns when appropriate.
C               Otherwise no scaling will be attempted.
C     NRDA -- Row dimension of A, NRDA greater or equal to N.
C     DIAG,KPIVOT,ROWS, -- Arrays of length at least N used internally
C          RS,SCALES         (except for SCALES which is M).
C
C *********************************************************************
C   OUTPUT
C *********************************************************************
C
C     IFLAG - Status indicator
C            =1 for successful decomposition.
C            =2 if improper input is detected.
C            =3 if rank of the matrix is less than N.
C     A -- Contains the reduced matrix in the strictly lower triangular
C          part and transformation information.
C     IRANK -- Contains the numerically determined matrix rank.
C     DIAG -- Contains the diagonal elements of the reduced
C             triangular matrix.
C     KPIVOT -- Contains the pivotal information, the column
C               interchanges performed on the original matrix are
C               recorded here.
C     SCALES -- Contains the column scaling parameters.
C
C *********************************************************************
C
C***SEE ALSO  DBVSUP, DSUDS
C***REFERENCES  G. Golub, Numerical methods for solving linear least
C                 squares problems, Numerische Mathematik 7, (1965),
C                 pp. 206-216.
C               P. Businger and G. Golub, Linear least squares
C                 solutions by Householder transformations, Numerische
C                 Mathematik  7, (1965), pp. 269-276.
C***ROUTINES CALLED  D1MACH, DCSCAL, DDOT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DORTHR
      DOUBLE PRECISION DDOT, D1MACH
      INTEGER IFLAG, IRANK, ISCALE, J, JROW, K, KP, KPIVOT(*), L, M,
     1     MK, N, NRDA
      DOUBLE PRECISION A(NRDA,*), ACC, AKK, ANORM, AS, ASAVE, DIAG(*),
     1     DIAGK, DUM, ROWS(*), RS(*), RSS, SAD, SCALES(*), SIG, SIGMA,
     2     SRURO, URO
C
C     ******************************************************************
C
C          MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
C          BY THE FUNCTION D1MACH.
C
C     ******************************************************************
C
C***FIRST EXECUTABLE STATEMENT  DORTHR
      URO = D1MACH(4)
      IF (M .GE. N .AND. N .GE. 1 .AND. NRDA .GE. N) GO TO 10
         IFLAG = 2
         CALL XERMSG ('SLATEC', 'DORTHR', 'INVALID INPUT PARAMETERS.',
     +      2, 1)
      GO TO 150
   10 CONTINUE
C
         ACC = 10.0D0*URO
         IF (IFLAG .LT. 0) ACC = MAX(ACC,10.0D0**IFLAG)
         SRURO = SQRT(URO)
         IFLAG = 1
         IRANK = N
C
C        COMPUTE NORM**2 OF JTH ROW AND A MATRIX NORM
C
         ANORM = 0.0D0
         DO 20 J = 1, N
            KPIVOT(J) = J
            ROWS(J) = DDOT(M,A(J,1),NRDA,A(J,1),NRDA)
            RS(J) = ROWS(J)
            ANORM = ANORM + ROWS(J)
   20    CONTINUE
C
C        PERFORM COLUMN SCALING ON A WHEN SPECIFIED
C
         CALL DCSCAL(A,NRDA,N,M,SCALES,DUM,ROWS,RS,ANORM,SCALES,ISCALE,
     1               1)
C
         ANORM = SQRT(ANORM)
C
C
C        CONSTRUCTION OF LOWER TRIANGULAR MATRIX AND RECORDING OF
C        ORTHOGONAL TRANSFORMATIONS
C
C
         DO 130 K = 1, N
C           BEGIN BLOCK PERMITTING ...EXITS TO 80
               MK = M - K + 1
C           ...EXIT
               IF (K .EQ. N) GO TO 80
               KP = K + 1
C
C              SEARCHING FOR PIVOTAL ROW
C
               DO 60 J = K, N
C                 BEGIN BLOCK PERMITTING ...EXITS TO 50
                     IF (ROWS(J) .GE. SRURO*RS(J)) GO TO 30
                        ROWS(J) = DDOT(MK,A(J,K),NRDA,A(J,K),NRDA)
                        RS(J) = ROWS(J)
   30                CONTINUE
                     IF (J .EQ. K) GO TO 40
C                 ......EXIT
                        IF (SIGMA .GE. 0.99D0*ROWS(J)) GO TO 50
   40                CONTINUE
                     SIGMA = ROWS(J)
                     JROW = J
   50             CONTINUE
   60          CONTINUE
C           ...EXIT
               IF (JROW .EQ. K) GO TO 80
C
C              PERFORM ROW INTERCHANGE
C
               L = KPIVOT(K)
               KPIVOT(K) = KPIVOT(JROW)
               KPIVOT(JROW) = L
               ROWS(JROW) = ROWS(K)
               ROWS(K) = SIGMA
               RSS = RS(K)
               RS(K) = RS(JROW)
               RS(JROW) = RSS
               DO 70 L = 1, M
                  ASAVE = A(K,L)
                  A(K,L) = A(JROW,L)
                  A(JROW,L) = ASAVE
   70          CONTINUE
   80       CONTINUE
C
C           CHECK RANK OF THE MATRIX
C
            SIG = DDOT(MK,A(K,K),NRDA,A(K,K),NRDA)
            DIAGK = SQRT(SIG)
            IF (DIAGK .GT. ACC*ANORM) GO TO 90
C
C              RANK DEFICIENT PROBLEM
               IFLAG = 3
               IRANK = K - 1
               CALL XERMSG ('SLATEC', 'DORTHR',
     +            'RANK OF MATRIX IS LESS THAN THE NUMBER OF ROWS.', 1,
     +            1)
C        ......EXIT
               GO TO 140
   90       CONTINUE
C
C           CONSTRUCT AND APPLY TRANSFORMATION TO MATRIX A
C
            AKK = A(K,K)
            IF (AKK .GT. 0.0D0) DIAGK = -DIAGK
            DIAG(K) = DIAGK
            A(K,K) = AKK - DIAGK
            IF (K .EQ. N) GO TO 120
               SAD = DIAGK*AKK - SIG
               DO 110 J = KP, N
                  AS = DDOT(MK,A(K,K),NRDA,A(J,K),NRDA)/SAD
                  DO 100 L = K, M
                     A(J,L) = A(J,L) + AS*A(K,L)
  100             CONTINUE
                  ROWS(J) = ROWS(J) - A(J,K)**2
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
C
C
      RETURN
      END
