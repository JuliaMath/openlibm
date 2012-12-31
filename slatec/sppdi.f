*DECK SPPDI
      SUBROUTINE SPPDI (AP, N, DET, JOB)
C***BEGIN PROLOGUE  SPPDI
C***PURPOSE  Compute the determinant and inverse of a real symmetric
C            positive definite matrix using factors from SPPCO or SPPFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2B1B, D3B1B
C***TYPE      SINGLE PRECISION (SPPDI-S, DPPDI-D, CPPDI-C)
C***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
C             PACKED, POSITIVE DEFINITE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     SPPDI computes the determinant and inverse
C     of a real symmetric positive definite matrix
C     using the factors computed by SPPCO or SPPFA .
C
C     On Entry
C
C        AP      REAL (N*(N+1)/2)
C                the output from SPPCO or SPPFA.
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        AP      the upper triangular half of the inverse .
C                The strict lower triangle is unaltered.
C
C        DET     REAL(2)
C                determinant of original matrix if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. DET(1) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if SPOCO or SPOFA has set INFO .EQ. 0 .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  SAXPY, SSCAL
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SPPDI
      INTEGER N,JOB
      REAL AP(*)
      REAL DET(2)
C
      REAL T
      REAL S
      INTEGER I,II,J,JJ,JM1,J1,K,KJ,KK,KP1,K1
C***FIRST EXECUTABLE STATEMENT  SPPDI
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0E0
         DET(2) = 0.0E0
         S = 10.0E0
         II = 0
         DO 50 I = 1, N
            II = II + I
            DET(1) = AP(II)**2*DET(1)
            IF (DET(1) .EQ. 0.0E0) GO TO 60
   10       IF (DET(1) .GE. 1.0E0) GO TO 20
               DET(1) = S*DET(1)
               DET(2) = DET(2) - 1.0E0
            GO TO 10
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40
               DET(1) = DET(1)/S
               DET(2) = DET(2) + 1.0E0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(R)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 140
         KK = 0
         DO 100 K = 1, N
            K1 = KK + 1
            KK = KK + K
            AP(KK) = 1.0E0/AP(KK)
            T = -AP(KK)
            CALL SSCAL(K-1,T,AP(K1),1)
            KP1 = K + 1
            J1 = KK + 1
            KJ = KK + K
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = AP(KJ)
               AP(KJ) = 0.0E0
               CALL SAXPY(K,T,AP(K1),1,AP(J1),1)
               J1 = J1 + J
               KJ = KJ + J
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM  INVERSE(R) * TRANS(INVERSE(R))
C
         JJ = 0
         DO 130 J = 1, N
            J1 = JJ + 1
            JJ = JJ + J
            JM1 = J - 1
            K1 = 1
            KJ = J1
            IF (JM1 .LT. 1) GO TO 120
            DO 110 K = 1, JM1
               T = AP(KJ)
               CALL SAXPY(K,T,AP(J1),1,AP(K1),1)
               K1 = K1 + K
               KJ = KJ + 1
  110       CONTINUE
  120       CONTINUE
            T = AP(JJ)
            CALL SSCAL(J,T,AP(J1),1)
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
