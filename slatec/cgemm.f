*DECK CGEMM
      SUBROUTINE CGEMM (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $   BETA, C, LDC)
C***BEGIN PROLOGUE  CGEMM
C***PURPOSE  Multiply a complex general matrix by a complex general
C            matrix.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B6
C***TYPE      COMPLEX (SGEMM-S, DGEMM-D, CGEMM-C)
C***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S. (NAG)
C***DESCRIPTION
C
C  CGEMM  performs one of the matrix-matrix operations
C
C     C := alpha*op( A )*op( B ) + beta*C,
C
C  where  op( X ) is one of
C
C     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
C
C  alpha and beta are scalars, and A, B and C are matrices, with op( A )
C  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n',  op( A ) = A.
C
C              TRANSA = 'T' or 't',  op( A ) = A'.
C
C              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
C
C           Unchanged on exit.
C
C  TRANSB - CHARACTER*1.
C           On entry, TRANSB specifies the form of op( B ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSB = 'N' or 'n',  op( B ) = B.
C
C              TRANSB = 'T' or 't',  op( B ) = B'.
C
C              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry,  M  specifies  the number  of rows  of the  matrix
C           op( A )  and of the  matrix  C.  M  must  be at least  zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry,  N  specifies the number  of columns of the matrix
C           op( B ) and the number of columns of the matrix C. N must be
C           at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry,  K  specifies  the number of columns of the matrix
C           op( A ) and the number of rows of the matrix op( B ). K must
C           be at least  zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
C           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
C           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
C           part of the array  A  must contain the matrix  A,  otherwise
C           the leading  k by m  part of the array  A  must contain  the
C           matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
C           LDA must be at least  max( 1, m ), otherwise  LDA must be at
C           least  max( 1, k ).
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
C           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
C           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
C           part of the array  B  must contain the matrix  B,  otherwise
C           the leading  n by k  part of the array  B  must contain  the
C           matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
C           LDB must be at least  max( 1, k ), otherwise  LDB must be at
C           least  max( 1, n ).
C           Unchanged on exit.
C
C  BETA   - COMPLEX         .
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then C need not be set on input.
C           Unchanged on exit.
C
C  C      - COMPLEX          array of DIMENSION ( LDC, n ).
C           Before entry, the leading  m by n  part of the array  C must
C           contain the matrix  C,  except when  beta  is zero, in which
C           case C need not be set on entry.
C           On exit, the array  C  is overwritten by the  m by n  matrix
C           ( alpha*op( A )*op( B ) + beta*C ).
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
C***END PROLOGUE  CGEMM
C     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX            ALPHA, BETA
C     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
C     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX            TEMP
C     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
C***FIRST EXECUTABLE STATEMENT  CGEMM
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
C     B  respectively are to be  transposed but  not conjugated  and set
C     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
C     and the number of rows of  B  respectively.
C
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
C
C     Test the input parameters.
C
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.CONJA                ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.CONJB                ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMM ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     And when  alpha.eq.zero.
C
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF( NOTB )THEN
         IF( NOTA )THEN
C
C           Form  C := alpha*A*B + beta*C.
C
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
C
C           Form  C := alpha*conjg( A' )*B + beta*C.
C
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B + beta*C
C
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
C
C           Form  C := alpha*A*conjg( B' ) + beta*C.
C
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*CONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
C
C           Form  C := alpha*A*B'          + beta*C
C
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
C
C           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
C
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*CONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
C
C           Form  C := alpha*conjg( A' )*B' + beta*C
C
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
C
C           Form  C := alpha*A'*conjg( B' ) + beta*C
C
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*CONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B' + beta*C
C
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of CGEMM .
C
      END
