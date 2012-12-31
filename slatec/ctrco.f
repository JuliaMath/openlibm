*DECK CTRCO
      SUBROUTINE CTRCO (T, LDT, N, RCOND, Z, JOB)
C***BEGIN PROLOGUE  CTRCO
C***PURPOSE  Estimate the condition number of a triangular matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C3
C***TYPE      COMPLEX (STRCO-S, DTRCO-D, CTRCO-C)
C***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
C             TRIANGULAR MATRIX
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CTRCO estimates the condition of a complex triangular matrix.
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
C                = 0         T  is lower triangular.
C                = nonzero   T  is upper triangular.
C
C     On Return
C
C        RCOND   REAL
C                an estimate of the reciprocal condition of  T .
C                For the system  T*X = B , relative perturbations
C                in  T  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  T  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       COMPLEX(N)
C                a work vector whose contents are usually unimportant.
C                If  T  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CSSCAL, SCASUM
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CTRCO
      INTEGER LDT,N,JOB
      COMPLEX T(LDT,*),Z(*)
      REAL RCOND
C
      COMPLEX W,WK,WKM,EK
      REAL TNORM,YNORM,S,SM,SCASUM
      INTEGER I1,J,J1,J2,K,KK,L
      LOGICAL LOWER
      COMPLEX ZDUM,ZDUM1,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))
C
C***FIRST EXECUTABLE STATEMENT  CTRCO
      LOWER = JOB .EQ. 0
C
C     COMPUTE 1-NORM OF T
C
      TNORM = 0.0E0
      DO 10 J = 1, N
         L = J
         IF (LOWER) L = N + 1 - J
         I1 = 1
         IF (LOWER) I1 = J
         TNORM = MAX(TNORM,SCASUM(L,T(I1,J),1))
   10 CONTINUE
C
C     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  CTRANS(T)*Y = E .
C     CTRANS(T)  IS THE CONJUGATE TRANSPOSE OF T .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF Y .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE CTRANS(T)*Y = E
C
      EK = (1.0E0,0.0E0)
      DO 20 J = 1, N
         Z(J) = (0.0E0,0.0E0)
   20 CONTINUE
      DO 100 KK = 1, N
         K = KK
         IF (LOWER) K = N + 1 - KK
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))
         IF (CABS1(EK-Z(K)) .LE. CABS1(T(K,K))) GO TO 30
            S = CABS1(T(K,K))/CABS1(EK-Z(K))
            CALL CSSCAL(N,S,Z,1)
            EK = CMPLX(S,0.0E0)*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = CABS1(WK)
         SM = CABS1(WKM)
         IF (CABS1(T(K,K)) .EQ. 0.0E0) GO TO 40
            WK = WK/CONJG(T(K,K))
            WKM = WKM/CONJG(T(K,K))
         GO TO 50
   40    CONTINUE
            WK = (1.0E0,0.0E0)
            WKM = (1.0E0,0.0E0)
   50    CONTINUE
         IF (KK .EQ. N) GO TO 90
            J1 = K + 1
            IF (LOWER) J1 = 1
            J2 = N
            IF (LOWER) J2 = K - 1
            DO 60 J = J1, J2
               SM = SM + CABS1(Z(J)+WKM*CONJG(T(K,J)))
               Z(J) = Z(J) + WK*CONJG(T(K,J))
               S = S + CABS1(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               W = WKM - WK
               WK = WKM
               DO 70 J = J1, J2
                  Z(J) = Z(J) + W*CONJG(T(K,J))
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE T*Z = Y
C
      DO 130 KK = 1, N
         K = N + 1 - KK
         IF (LOWER) K = KK
         IF (CABS1(Z(K)) .LE. CABS1(T(K,K))) GO TO 110
            S = CABS1(T(K,K))/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  110    CONTINUE
         IF (CABS1(T(K,K)) .NE. 0.0E0) Z(K) = Z(K)/T(K,K)
         IF (CABS1(T(K,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         I1 = 1
         IF (LOWER) I1 = K + 1
         IF (KK .GE. N) GO TO 120
            W = -Z(K)
            CALL CAXPY(N-KK,W,T(I1,K),1,Z(I1),1)
  120    CONTINUE
  130 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (TNORM .NE. 0.0E0) RCOND = YNORM/TNORM
      IF (TNORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
