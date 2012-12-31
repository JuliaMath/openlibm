*DECK CGBMV
      SUBROUTINE CGBMV (TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
     $   BETA, Y, INCY)
C***BEGIN PROLOGUE  CGBMV
C***PURPOSE  Multiply a complex vector by a complex general band matrix.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B4
C***TYPE      COMPLEX (SGBMV-S, DGBMV-D, CGBMV-C)
C***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J. J., (ANL)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  CGBMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
C
C     y := alpha*conjg( A' )*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  KL     - INTEGER.
C           On entry, KL specifies the number of sub-diagonals of the
C           matrix A. KL must satisfy  0 .le. KL.
C           Unchanged on exit.
C
C  KU     - INTEGER.
C           On entry, KU specifies the number of super-diagonals of the
C           matrix A. KU must satisfy  0 .le. KU.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, n ).
C           Before entry, the leading ( kl + ku + 1 ) by n part of the
C           array A must contain the matrix of coefficients, supplied
C           column by column, with the leading diagonal of the matrix in
C           row ( ku + 1 ) of the array, the first super-diagonal
C           starting at position 2 in row ku, the first sub-diagonal
C           starting at position 1 in row ( ku + 2 ), and so on.
C           Elements in the array A that do not correspond to elements
C           in the band matrix (such as the top left ku by ku triangle)
C           are not referenced.
C           The following program segment will transfer a band matrix
C           from conventional full matrix storage to band storage:
C
C                 DO 20, J = 1, N
C                    K = KU + 1 - J
C                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
C                       A( K + I, J ) = matrix( I, J )
C              10    CONTINUE
C              20 CONTINUE
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           ( kl + ku + 1 ).
C           Unchanged on exit.
C
C  X      - COMPLEX          array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - COMPLEX         .
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - COMPLEX          array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry, the incremented array Y must contain the
C           vector y. On exit, Y is overwritten by the updated vector y.
C
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
C***END PROLOGUE  CGBMV
C     .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
      INTEGER            INCX, INCY, KL, KU, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
C     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
C     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KUP1, KX, KY,
     $                   LENX, LENY
      LOGICAL            NOCONJ
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, MIN
C***FIRST EXECUTABLE STATEMENT  CGBMV
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( KL.LT.0 )THEN
         INFO = 4
      ELSE IF( KU.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.( KL + KU + 1 ) )THEN
         INFO = 8
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 10
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGBMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
      NOCONJ = LSAME( TRANS, 'T' )
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the band part of A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      KUP1 = KU + 1
      IF( LSAME( TRANS, 'N' ) )THEN
C
C        Form  y := alpha*A*x + y.
C
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  K    = KUP1 - J
                  DO 50, I = MAX( 1, J - KU ), MIN( M, J + KL )
                     Y( I ) = Y( I ) + TEMP*A( K + I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  K    = KUP1 - J
                  DO 70, I = MAX( 1, J - KU ), MIN( M, J + KL )
                     Y( IY ) = Y( IY ) + TEMP*A( K + I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
               IF( J.GT.KU )
     $            KY = KY + INCY
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
C
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               K    = KUP1 - J
               IF( NOCONJ )THEN
                  DO 90, I = MAX( 1, J - KU ), MIN( M, J + KL )
                     TEMP = TEMP + A( K + I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = MAX( 1, J - KU ), MIN( M, J + KL )
                     TEMP = TEMP + CONJG( A( K + I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               K    = KUP1 - J
               IF( NOCONJ )THEN
                  DO 120, I = MAX( 1, J - KU ), MIN( M, J + KL )
                     TEMP = TEMP + A( K + I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = MAX( 1, J - KU ), MIN( M, J + KL )
                     TEMP = TEMP + CONJG( A( K + I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
               IF( J.GT.KU )
     $            KX = KX + INCX
  140       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of CGBMV .
C
      END
