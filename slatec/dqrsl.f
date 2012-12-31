*DECK DQRSL
      SUBROUTINE DQRSL (X, LDX, N, K, QRAUX, Y, QY, QTY, B, RSD, XB,
     +   JOB, INFO)
C***BEGIN PROLOGUE  DQRSL
C***PURPOSE  Apply the output of DQRDC to compute coordinate transfor-
C            mations, projections, and least squares solutions.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D9, D2A1
C***TYPE      DOUBLE PRECISION (SQRSL-S, DQRSL-D, CQRSL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR,
C             SOLVE
C***AUTHOR  Stewart, G. W., (U. of Maryland)
C***DESCRIPTION
C
C     DQRSL applies the output of DQRDC to compute coordinate
C     transformations, projections, and least squares solutions.
C     For K .LE. MIN(N,P), let XK be the matrix
C
C            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
C
C     formed from columns JPVT(1), ... ,JPVT(K) of the original
C     N X P matrix X that was input to DQRDC (if no pivoting was
C     done, XK consists of the first K columns of X in their
C     original order).  DQRDC produces a factored orthogonal matrix Q
C     and an upper triangular matrix R such that
C
C              XK = Q * (R)
C                       (0)
C
C     This information is contained in coded form in the arrays
C     X and QRAUX.
C
C     On Entry
C
C        X      DOUBLE PRECISION(LDX,P).
C               X contains the output of DQRDC.
C
C        LDX    INTEGER.
C               LDX is the leading dimension of the array X.
C
C        N      INTEGER.
C               N is the number of rows of the matrix XK.  It must
C               have the same value as N in DQRDC.
C
C        K      INTEGER.
C               K is the number of columns of the matrix XK.  K
C               must not be greater than MIN(N,P), where P is the
C               same as in the calling sequence to DQRDC.
C
C        QRAUX  DOUBLE PRECISION(P).
C               QRAUX contains the auxiliary output from DQRDC.
C
C        Y      DOUBLE PRECISION(N)
C               Y contains an N-vector that is to be manipulated
C               by DQRSL.
C
C        JOB    INTEGER.
C               JOB specifies what is to be computed.  JOB has
C               the decimal expansion ABCDE, with the following
C               meaning.
C
C                    If A .NE. 0, compute QY.
C                    If B,C,D, or E .NE. 0, compute QTY.
C                    If C .NE. 0, compute B.
C                    If D .NE. 0, compute RSD.
C                    If E .NE. 0, compute XB.
C
C               Note that a request to compute B, RSD, or XB
C               automatically triggers the computation of QTY, for
C               which an array must be provided in the calling
C               sequence.
C
C     On Return
C
C        QY     DOUBLE PRECISION(N).
C               QY contains Q*Y, if its computation has been
C               requested.
C
C        QTY    DOUBLE PRECISION(N).
C               QTY contains TRANS(Q)*Y, if its computation has
C               been requested.  Here TRANS(Q) is the
C               transpose of the matrix Q.
C
C        B      DOUBLE PRECISION(K)
C               B contains the solution of the least squares problem
C
C                    minimize norm2(Y - XK*B),
C
C               if its computation has been requested.  (Note that
C               if pivoting was requested in DQRDC, the J-th
C               component of B will be associated with column JPVT(J)
C               of the original matrix X that was input into DQRDC.)
C
C        RSD    DOUBLE PRECISION(N).
C               RSD contains the least squares residual Y - XK*B,
C               if its computation has been requested.  RSD is
C               also the orthogonal projection of Y onto the
C               orthogonal complement of the column space of XK.
C
C        XB     DOUBLE PRECISION(N).
C               XB contains the least squares approximation XK*B,
C               if its computation has been requested.  XB is also
C               the orthogonal projection of Y onto the column space
C               of X.
C
C        INFO   INTEGER.
C               INFO is zero unless the computation of B has
C               been requested and R is exactly singular.  In
C               this case, INFO is the index of the first zero
C               diagonal element of R and B is left unaltered.
C
C     The parameters QY, QTY, B, RSD, and XB are not referenced
C     if their computation is not requested and in this case
C     can be replaced by dummy variables in the calling program.
C     To save storage, the user may in some cases use the same
C     array for different parameters in the calling sequence.  A
C     frequently occurring example is when one wishes to compute
C     any of B, RSD, or XB and does not need Y or QTY.  In this
C     case one may identify Y, QTY, and one of B, RSD, or XB, while
C     providing separate arrays for anything else that is to be
C     computed.  Thus the calling sequence
C
C          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
C
C     will result in the computation of B and RSD, with RSD
C     overwriting Y.  More generally, each item in the following
C     list contains groups of permissible identifications for
C     a single calling sequence.
C
C          1. (Y,QTY,B) (RSD) (XB) (QY)
C
C          2. (Y,QTY,RSD) (B) (XB) (QY)
C
C          3. (Y,QTY,XB) (B) (RSD) (QY)
C
C          4. (Y,QY) (QTY,B) (RSD) (XB)
C
C          5. (Y,QY) (QTY,RSD) (B) (XB)
C
C          6. (Y,QY) (QTY,XB) (B) (RSD)
C
C     In any group the value returned in the array allocated to
C     the group corresponds to the last member of the group.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DCOPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DQRSL
      INTEGER LDX,N,K,JOB,INFO
      DOUBLE PRECISION X(LDX,*),QRAUX(*),Y(*),QY(*),QTY(*),B(*),RSD(*),
     1                 XB(*)
C
      INTEGER I,J,JJ,JU,KP1
      DOUBLE PRECISION DDOT,T,TEMP
      LOGICAL CB,CQY,CQTY,CR,CXB
C***FIRST EXECUTABLE STATEMENT  DQRSL
C
C     SET INFO FLAG.
C
      INFO = 0
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      CQY = JOB/10000 .NE. 0
      CQTY = MOD(JOB,10000) .NE. 0
      CB = MOD(JOB,1000)/100 .NE. 0
      CR = MOD(JOB,100)/10 .NE. 0
      CXB = MOD(JOB,10) .NE. 0
      JU = MIN(K,N-1)
C
C     SPECIAL ACTION WHEN N=1.
C
      IF (JU .NE. 0) GO TO 40
         IF (CQY) QY(1) = Y(1)
         IF (CQTY) QTY(1) = Y(1)
         IF (CXB) XB(1) = Y(1)
         IF (.NOT.CB) GO TO 30
            IF (X(1,1) .NE. 0.0D0) GO TO 10
               INFO = 1
            GO TO 20
   10       CONTINUE
               B(1) = Y(1)/X(1,1)
   20       CONTINUE
   30    CONTINUE
         IF (CR) RSD(1) = 0.0D0
      GO TO 250
   40 CONTINUE
C
C        SET UP TO COMPUTE QY OR QTY.
C
         IF (CQY) CALL DCOPY(N,Y,1,QY,1)
         IF (CQTY) CALL DCOPY(N,Y,1,QTY,1)
         IF (.NOT.CQY) GO TO 70
C
C           COMPUTE QY.
C
            DO 60 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 50
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QY(J),1)
                  X(J,J) = TEMP
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
         IF (.NOT.CQTY) GO TO 100
C
C           COMPUTE TRANS(Q)*Y.
C
            DO 90 J = 1, JU
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 80
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QTY(J),1)
                  X(J,J) = TEMP
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        SET UP TO COMPUTE B, RSD, OR XB.
C
         IF (CB) CALL DCOPY(K,QTY,1,B,1)
         KP1 = K + 1
         IF (CXB) CALL DCOPY(K,QTY,1,XB,1)
         IF (CR .AND. K .LT. N) CALL DCOPY(N-K,QTY(KP1),1,RSD(KP1),1)
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120
            DO 110 I = KP1, N
               XB(I) = 0.0D0
  110       CONTINUE
  120    CONTINUE
         IF (.NOT.CR) GO TO 140
            DO 130 I = 1, K
               RSD(I) = 0.0D0
  130       CONTINUE
  140    CONTINUE
         IF (.NOT.CB) GO TO 190
C
C           COMPUTE B.
C
            DO 170 JJ = 1, K
               J = K - JJ + 1
               IF (X(J,J) .NE. 0.0D0) GO TO 150
                  INFO = J
                  GO TO 180
  150          CONTINUE
               B(J) = B(J)/X(J,J)
               IF (J .EQ. 1) GO TO 160
                  T = -B(J)
                  CALL DAXPY(J-1,T,X(1,J),1,B,1)
  160          CONTINUE
  170       CONTINUE
  180       CONTINUE
  190    CONTINUE
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240
C
C           COMPUTE RSD OR XB AS REQUIRED.
C
            DO 230 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 220
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  IF (.NOT.CR) GO TO 200
                     T = -DDOT(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,RSD(J),1)
  200             CONTINUE
                  IF (.NOT.CXB) GO TO 210
                     T = -DDOT(N-J+1,X(J,J),1,XB(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,XB(J),1)
  210             CONTINUE
                  X(J,J) = TEMP
  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RETURN
      END
