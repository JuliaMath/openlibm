*DECK CGEDI
      SUBROUTINE CGEDI (A, LDA, N, IPVT, DET, WORK, JOB)
C***BEGIN PROLOGUE  CGEDI
C***PURPOSE  Compute the determinant and inverse of a matrix using the
C            factors computed by CGECO or CGEFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C1, D3C1
C***TYPE      COMPLEX (SGEDI-S, DGEDI-D, CGEDI-C)
C***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CGEDI computes the determinant and inverse of a matrix
C     using the factors computed by CGECO or CGEFA.
C
C     On Entry
C
C        A       COMPLEX(LDA, N)
C                the output from CGECO or CGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from CGECO or CGEFA.
C
C        WORK    COMPLEX(N)
C                work vector.  Contents destroyed.
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        A       inverse of original matrix if requested.
C                Otherwise unchanged.
C
C        DET     COMPLEX(2)
C                determinant of original matrix if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. CABS1(DET(1)) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if CGECO has set RCOND .GT. 0.0 or CGEFA has set
C        INFO .EQ. 0 .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CSCAL, CSWAP
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CGEDI
      INTEGER LDA,N,IPVT(*),JOB
      COMPLEX A(LDA,*),DET(2),WORK(*)
C
      COMPLEX T
      REAL TEN
      INTEGER I,J,K,KB,KP1,L,NM1
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C***FIRST EXECUTABLE STATEMENT  CGEDI
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = (1.0E0,0.0E0)
         DET(2) = (0.0E0,0.0E0)
         TEN = 10.0E0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
            IF (CABS1(DET(1)) .EQ. 0.0E0) GO TO 60
   10       IF (CABS1(DET(1)) .GE. 1.0E0) GO TO 20
               DET(1) = CMPLX(TEN,0.0E0)*DET(1)
               DET(2) = DET(2) - (1.0E0,0.0E0)
            GO TO 10
   20       CONTINUE
   30       IF (CABS1(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/CMPLX(TEN,0.0E0)
               DET(2) = DET(2) + (1.0E0,0.0E0)
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = (1.0E0,0.0E0)/A(K,K)
            T = -A(K,K)
            CALL CSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = (0.0E0,0.0E0)
               CALL CAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = (0.0E0,0.0E0)
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL CAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL CSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
