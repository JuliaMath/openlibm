*DECK CTRMM
      SUBROUTINE CTRMM (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $   B, LDB)
C***BEGIN PROLOGUE  CTRMM
C***PURPOSE  Multiply a complex general matrix by a complex triangular
C            matrix.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B6
C***TYPE      COMPLEX (STRMM-S, DTRMM-D, CTRMM-C)
C***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S. (NAG)
C***DESCRIPTION
C
C  CTRMM  performs one of the matrix-matrix operations
C
C     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
C
C  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
C  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
C
C     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C           On entry,  SIDE specifies whether  op( A ) multiplies B from
C           the left or right as follows:
C
C              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
C
C              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
C
C           Unchanged on exit.
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix A is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n'   op( A ) = A.
C
C              TRANSA = 'T' or 't'   op( A ) = A'.
C
C              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit triangular
C           as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of B. M must be at
C           least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of B.  N must be
C           at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
C           zero then  A is not referenced and  B need not be set before
C           entry.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, k ), where k is m
C           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
C           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
C           upper triangular part of the array  A must contain the upper
C           triangular matrix  and the strictly lower triangular part of
C           A is not referenced.
C           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
C           lower triangular part of the array  A must contain the lower
C           triangular matrix  and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
C           A  are not referenced either,  but are assumed to be  unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
C           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
C           then LDA must be at least max( 1, n ).
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, n ).
C           Before entry,  the leading  m by n part of the array  B must
C           contain the matrix  B,  and  on exit  is overwritten  by the
C           transformed matrix.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in  the  calling  (sub)  program.   LDB  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C***REFERENCES  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
C                 A set of level 3 basic linear algebra subprograms.
C                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
C***ROUTINES CALLED  LSAME, XERBLA
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910605  Modified to meet SLATEC prologue standards.  Only comment
C           lines were modified.  (BKS)
C***END PROLOGUE  CTRMM
C     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX            ALPHA
C     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
C     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX            TEMP
C     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
C***FIRST EXECUTABLE STATEMENT  CTRMM
C
C     Test the input parameters.
C
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
C
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRMM ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
C     And when  alpha.eq.zero.
C
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
C
C     Start the operations.
C
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
C
C           Form  B := alpha*A*B.
C
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
C
            IF( UPPER )THEN
               DO 120, J = 1, N
                  DO 110, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOCONJ )THEN
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( I, I )
                        DO 90, K = 1, I - 1
                           TEMP = TEMP + A( K, I )*B( K, J )
   90                   CONTINUE
                     ELSE
                        IF( NOUNIT )
     $                     TEMP = TEMP*CONJG( A( I, I ) )
                        DO 100, K = 1, I - 1
                           TEMP = TEMP + CONJG( A( K, I ) )*B( K, J )
  100                   CONTINUE
                     END IF
                     B( I, J ) = ALPHA*TEMP
  110             CONTINUE
  120          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = 1, M
                     TEMP = B( I, J )
                     IF( NOCONJ )THEN
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( I, I )
                        DO 130, K = I + 1, M
                           TEMP = TEMP + A( K, I )*B( K, J )
  130                   CONTINUE
                     ELSE
                        IF( NOUNIT )
     $                     TEMP = TEMP*CONJG( A( I, I ) )
                        DO 140, K = I + 1, M
                           TEMP = TEMP + CONJG( A( K, I ) )*B( K, J )
  140                   CONTINUE
                     END IF
                     B( I, J ) = ALPHA*TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
C
C           Form  B := alpha*B*A.
C
            IF( UPPER )THEN
               DO 200, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 170, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  170             CONTINUE
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
  200          CONTINUE
            ELSE
               DO 240, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 210, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  210             CONTINUE
                  DO 230, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 220, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  220                   CONTINUE
                     END IF
  230             CONTINUE
  240          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
C
            IF( UPPER )THEN
               DO 280, K = 1, N
                  DO 260, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = ALPHA*A( J, K )
                        ELSE
                           TEMP = ALPHA*CONJG( A( J, K ) )
                        END IF
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = TEMP*A( K, K )
                     ELSE
                        TEMP = TEMP*CONJG( A( K, K ) )
                     END IF
                  END IF
                  IF( TEMP.NE.ONE )THEN
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
  280          CONTINUE
            ELSE
               DO 320, K = N, 1, -1
                  DO 300, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = ALPHA*A( J, K )
                        ELSE
                           TEMP = ALPHA*CONJG( A( J, K ) )
                        END IF
                        DO 290, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  290                   CONTINUE
                     END IF
  300             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = TEMP*A( K, K )
                     ELSE
                        TEMP = TEMP*CONJG( A( K, K ) )
                     END IF
                  END IF
                  IF( TEMP.NE.ONE )THEN
                     DO 310, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  310                CONTINUE
                  END IF
  320          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of CTRMM .
C
      END
