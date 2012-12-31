*DECK DSYR2K
      SUBROUTINE DSYR2K (UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA,
     $   C, LDC)
C***BEGIN PROLOGUE  DSYR2K
C***PURPOSE  Perform one of the symmetric rank 2k operations.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B6
C***TYPE      DOUBLE PRECISION (SSYR2-S, DSYR2-D, CSYR2-C, DSYR2K-D)
C***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S. (NAG)
C***DESCRIPTION
C
C  DSYR2K  performs one of the symmetric rank 2k operations
C
C     C := alpha*A*B' + alpha*B*A' + beta*C,
C
C  or
C
C     C := alpha*A'*B + alpha*B'*A + beta*C,
C
C  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
C  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
C  matrices in the second case.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On  entry,   UPLO  specifies  whether  the  upper  or  lower
C           triangular  part  of the  array  C  is to be  referenced  as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry,  TRANS  specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
C                                        beta*C.
C
C              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
C                                        beta*C.
C
C              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
C                                        beta*C.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry,  N specifies the order of the matrix C.  N must be
C           at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
C           of  columns  of the  matrices  A and B,  and on  entry  with
C           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
C           of rows of the matrices  A and B.  K must be at least  zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
C           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
C           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
C           part of the array  A  must contain the matrix  A,  otherwise
C           the leading  k by n  part of the array  A  must contain  the
C           matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
C           then  LDA must be at least  max( 1, n ), otherwise  LDA must
C           be at least  max( 1, k ).
C           Unchanged on exit.
C
C  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
C           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
C           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
C           part of the array  B  must contain the matrix  B,  otherwise
C           the leading  k by n  part of the array  B  must contain  the
C           matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
C           then  LDB must be at least  max( 1, n ), otherwise  LDB must
C           be at least  max( 1, k ).
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta.
C           Unchanged on exit.
C
C  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
C           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
C           upper triangular part of the array C must contain the upper
C           triangular part  of the  symmetric matrix  and the strictly
C           lower triangular part of C is not referenced.  On exit, the
C           upper triangular part of the array  C is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
C           lower triangular part of the array C must contain the lower
C           triangular part  of the  symmetric matrix  and the strictly
C           upper triangular part of C is not referenced.  On exit, the
C           lower triangular part of the array  C is overwritten by the
C           lower triangular part of the updated matrix.
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
C           max( 1, n ).
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
C***END PROLOGUE  DSYR2K
C     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
C
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      DOUBLE PRECISION   TEMP1, TEMP2
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C***FIRST EXECUTABLE STATEMENT  DSYR2K
C
C     Test the input parameters.
C
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
C
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR2K', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     And when  alpha.eq.zero.
C
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF( LSAME( TRANS, 'N' ) )THEN
C
C        Form  C := alpha*A*B' + alpha*B*A' + C.
C
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) +
     $                              A( I, L )*TEMP1 + B( I, L )*TEMP2
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) +
     $                              A( I, L )*TEMP1 + B( I, L )*TEMP2
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
C
C        Form  C := alpha*A'*B + alpha*B'*A + C.
C
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 190, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 220, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of DSYR2K.
C
      END
