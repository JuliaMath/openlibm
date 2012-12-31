*DECK CTRDI
      SUBROUTINE CTRDI (T, LDT, N, DET, JOB, INFO)
C***BEGIN PROLOGUE  CTRDI
C***PURPOSE  Compute the determinant and inverse of a triangular matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C3, D3C3
C***TYPE      COMPLEX (STRDI-S, DTRDI-D, CTRDI-C)
C***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
C             TRIANGULAR MATRIX
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CTRDI computes the determinant and inverse of a complex
C     triangular matrix.
C
C     On Entry
C
C        T       COMPLEX(LDT,N)
C                T contains the triangular matrix.  The zero
C                elements of the matrix are not referenced, and
C                the corresponding elements of the array can be
C                used to store other information.
C
C        LDT     INTEGER
C                LDT is the leading dimension of the array T.
C
C        N       INTEGER
C                N is the order of the system.
C
C        JOB     INTEGER
C                = 010       no det, inverse of lower triangular.
C                = 011       no det, inverse of upper triangular.
C                = 100       det, no inverse.
C                = 110       det, inverse of lower triangular.
C                = 111       det, inverse of upper triangular.
C
C     On Return
C
C        T       inverse of original matrix if requested.
C                Otherwise unchanged.
C
C        DET     COMPLEX(2)
C                determinant of original matrix if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. CABS1(DET(1)) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C        INFO    INTEGER
C                INFO contains zero if the system is nonsingular
C                and the inverse is requested.
C                Otherwise INFO contains the index of
C                a zero diagonal element of T.
C
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CSCAL
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CTRDI
      INTEGER LDT,N,JOB,INFO
      COMPLEX T(LDT,*),DET(2)
C
      COMPLEX TEMP
      REAL TEN
      INTEGER I,J,K,KB,KM1,KP1
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C***FIRST EXECUTABLE STATEMENT  CTRDI
C
C        COMPUTE DETERMINANT
C
         IF (JOB/100 .EQ. 0) GO TO 70
            DET(1) = (1.0E0,0.0E0)
            DET(2) = (0.0E0,0.0E0)
            TEN = 10.0E0
            DO 50 I = 1, N
               DET(1) = T(I,I)*DET(1)
               IF (CABS1(DET(1)) .EQ. 0.0E0) GO TO 60
   10          IF (CABS1(DET(1)) .GE. 1.0E0) GO TO 20
                  DET(1) = CMPLX(TEN,0.0E0)*DET(1)
                  DET(2) = DET(2) - (1.0E0,0.0E0)
               GO TO 10
   20          CONTINUE
   30          IF (CABS1(DET(1)) .LT. TEN) GO TO 40
                  DET(1) = DET(1)/CMPLX(TEN,0.0E0)
                  DET(2) = DET(2) + (1.0E0,0.0E0)
               GO TO 30
   40          CONTINUE
   50       CONTINUE
   60       CONTINUE
   70    CONTINUE
C
C        COMPUTE INVERSE OF UPPER TRIANGULAR
C
         IF (MOD(JOB/10,10) .EQ. 0) GO TO 170
            IF (MOD(JOB,10) .EQ. 0) GO TO 120
                  DO 100 K = 1, N
                     INFO = K
                     IF (CABS1(T(K,K)) .EQ. 0.0E0) GO TO 110
                     T(K,K) = (1.0E0,0.0E0)/T(K,K)
                     TEMP = -T(K,K)
                     CALL CSCAL(K-1,TEMP,T(1,K),1)
                     KP1 = K + 1
                     IF (N .LT. KP1) GO TO 90
                     DO 80 J = KP1, N
                        TEMP = T(K,J)
                        T(K,J) = (0.0E0,0.0E0)
                        CALL CAXPY(K,TEMP,T(1,K),1,T(1,J),1)
   80                CONTINUE
   90                CONTINUE
  100             CONTINUE
                  INFO = 0
  110          CONTINUE
            GO TO 160
  120       CONTINUE
C
C              COMPUTE INVERSE OF LOWER TRIANGULAR
C
               DO 150 KB = 1, N
                  K = N + 1 - KB
                  INFO = K
                  IF (CABS1(T(K,K)) .EQ. 0.0E0) GO TO 180
                  T(K,K) = (1.0E0,0.0E0)/T(K,K)
                  TEMP = -T(K,K)
                  IF (K .NE. N) CALL CSCAL(N-K,TEMP,T(K+1,K),1)
                  KM1 = K - 1
                  IF (KM1 .LT. 1) GO TO 140
                  DO 130 J = 1, KM1
                     TEMP = T(K,J)
                     T(K,J) = (0.0E0,0.0E0)
                     CALL CAXPY(N-K+1,TEMP,T(K,K),1,T(K,J),1)
  130             CONTINUE
  140             CONTINUE
  150          CONTINUE
               INFO = 0
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
      RETURN
      END
