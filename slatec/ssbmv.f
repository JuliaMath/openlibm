*DECK SSBMV
      SUBROUTINE SSBMV (UPLO, N, K, ALPHA, A, LDA, X, INCX, BETA, Y,
     $   INCY)
C***BEGIN PROLOGUE  SSBMV
C***PURPOSE  Multiply a real vector by a real symmetric band matrix.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B4
C***TYPE      SINGLE PRECISION (SSBMV-S, DSBMV-D, CSBMV-C)
C***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J. J., (ANL)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  SSBMV  performs the matrix-vector  operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n symmetric band matrix, with k super-diagonals.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the band matrix A is being supplied as
C           follows:
C
C              UPLO = 'U' or 'u'   The upper triangular part of A is
C                                  being supplied.
C
C              UPLO = 'L' or 'l'   The lower triangular part of A is
C                                  being supplied.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry, K specifies the number of super-diagonals of the
C           matrix A. K must satisfy  0 .le. K.
C           Unchanged on exit.
C
C  ALPHA  - REAL            .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n ).
C           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
C           by n part of the array A must contain the upper triangular
C           band part of the symmetric matrix, supplied column by
C           column, with the leading diagonal of the matrix in row
C           ( k + 1 ) of the array, the first super-diagonal starting at
C           position 2 in row k, and so on. The top left k by k triangle
C           of the array A is not referenced.
C           The following program segment will transfer the upper
C           triangular part of a symmetric band matrix from conventional
C           full matrix storage to band storage:
C
C                 DO 20, J = 1, N
C                    M = K + 1 - J
C                    DO 10, I = MAX( 1, J - K ), J
C                       A( M + I, J ) = matrix( I, J )
C              10    CONTINUE
C              20 CONTINUE
C
C           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
C           by n part of the array A must contain the lower triangular
C           band part of the symmetric matrix, supplied column by
C           column, with the leading diagonal of the matrix in row 1 of
C           the array, the first sub-diagonal starting at position 1 in
C           row 2, and so on. The bottom right k by k triangle of the
C           array A is not referenced.
C           The following program segment will transfer the lower
C           triangular part of a symmetric band matrix from conventional
C           full matrix storage to band storage:
C
C                 DO 20, J = 1, N
C                    M = 1 - J
C                    DO 10, I = J, MIN( N, J + K )
C                       A( M + I, J ) = matrix( I, J )
C              10    CONTINUE
C              20 CONTINUE
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           ( k + 1 ).
C           Unchanged on exit.
C
C  X      - REAL             array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - REAL            .
C           On entry, BETA specifies the scalar beta.
C           Unchanged on exit.
C
C  Y      - REAL             array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the
C           vector y. On exit, Y is overwritten by the updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
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
C***END PROLOGUE  SSBMV
C     .. Scalar Arguments ..
      REAL               ALPHA, BETA
      INTEGER            INCX, INCY, K, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      REAL               A( LDA, * ), X( * ), Y( * )
C     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
C     .. Local Scalars ..
      REAL               TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KPLUS1, KX, KY, L
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C***FIRST EXECUTABLE STATEMENT  SSBMV
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( K.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.( K + 1 ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSBMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set up the start points in  X  and  Y.
C
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of the array A
C     are accessed sequentially with one pass through A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( UPLO, 'U' ) )THEN
C
C        Form  y  when upper triangle of A is stored.
C
         KPLUS1 = K + 1
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               L     = KPLUS1 - J
               DO 50, I = MAX( 1, J - K ), J - 1
                  Y( I ) = Y( I ) + TEMP1*A( L + I, J )
                  TEMP2  = TEMP2  + A( L + I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( KPLUS1, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               L     = KPLUS1 - J
               DO 70, I = MAX( 1, J - K ), J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( L + I, J )
                  TEMP2   = TEMP2   + A( L + I, J )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( KPLUS1, J ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
               IF( J.GT.K )THEN
                  KX = KX + INCX
                  KY = KY + INCY
               END IF
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y  when lower triangle of A is stored.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*A( 1, J )
               L      = 1            - J
               DO 90, I = J + 1, MIN( N, J + K )
                  Y( I ) = Y( I ) + TEMP1*A( L + I, J )
                  TEMP2  = TEMP2  + A( L + I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*A( 1, J )
               L       = 1             - J
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, MIN( N, J + K )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( L + I, J )
                  TEMP2   = TEMP2   + A( L + I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of SSBMV .
C
      END
