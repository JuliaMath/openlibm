*DECK DSYR
      SUBROUTINE DSYR (UPLO, N, ALPHA, X, INCX, A, LDA)
C***BEGIN PROLOGUE  DSYR
C***PURPOSE  Perform the symmetric rank 1 operation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B4
C***TYPE      DOUBLE PRECISION (DSYR-D)
C***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J. J., (ANL)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  DSYR   performs the symmetric rank 1 operation
C
C     A := alpha*x*x' + A,
C
C  where alpha is a real scalar, x is an n element vector and A is an
C  n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1)*abs( INCX)).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced. On exit, the
C           upper triangular part of the array A is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced. On exit, the
C           lower triangular part of the array A is overwritten by the
C           lower triangular part of the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and
C                 Hanson, R. J.  An extended set of Fortran basic linear
C                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1,
C                 pp. 1-17, March 1988.
C***ROUTINES CALLED  LSAME, XERBLA
C***REVISION HISTORY  (YYMMDD)
C   861022  DATE WRITTEN
C   910605  Modified to meet SLATEC prologue standards.  Only comment
C           lines were modified.  (BKS)
C***END PROLOGUE  DSYR
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C***FIRST EXECUTABLE STATEMENT  DSYR
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR  ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Set the start point in X if the increment is not unity.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
      IF( LSAME( UPLO, 'U' ) )THEN
C
C        Form  A  when A is stored in upper triangle.
C
         IF( INCX.EQ.1 )THEN
            DO 20, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            JX = KX
            DO 40, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = KX
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
C
C        Form  A  when A is stored in lower triangle.
C
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of DSYR  .
C
      END
