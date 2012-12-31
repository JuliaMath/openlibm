*DECK QZVAL
      SUBROUTINE QZVAL (NM, N, A, B, ALFR, ALFI, BETA, MATZ, Z)
C***BEGIN PROLOGUE  QZVAL
C***PURPOSE  The third step of the QZ algorithm for generalized
C            eigenproblems.  Accepts a pair of real matrices, one in
C            quasi-triangular form and the other in upper triangular
C            form and computes the eigenvalues of the associated
C            eigenproblem.  Usually preceded by QZHES, QZIT, and
C            followed by QZVEC.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C2C
C***TYPE      SINGLE PRECISION (QZVAL-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is the third step of the QZ algorithm
C     for solving generalized matrix eigenvalue problems,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART.
C
C     This subroutine accepts a pair of REAL matrices, one of them
C     in quasi-triangular form and the other in upper triangular form.
C     It reduces the quasi-triangular matrix further, so that any
C     remaining 2-by-2 blocks correspond to pairs of complex
C     eigenvalues, and returns quantities whose ratios give the
C     generalized eigenvalues.  It is usually preceded by  QZHES
C     and  QZIT  and may be followed by  QZVEC.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A, B, and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrices A and B.  N is an INTEGER
C          variable.  N must be less than or equal to NM.
C
C        A contains a real upper quasi-triangular matrix.  A is a two-
C          dimensional REAL array, dimensioned A(NM,N).
C
C        B contains a real upper triangular matrix.  In addition,
C          location B(N,1) contains the tolerance quantity (EPSB)
C          computed and saved in  QZIT.  B is a two-dimensional REAL
C          array, dimensioned B(NM,N).
C
C        MATZ should be set to .TRUE. if the right hand transformations
C          are to be accumulated for later use in computing
C          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL
C          variable.
C
C        Z contains, if MATZ has been set to .TRUE., the transformation
C          matrix produced in the reductions by  QZHES  and  QZIT,  if
C          performed, or else the identity matrix.  If MATZ has been set
C          to .FALSE., Z is not referenced.  Z is a two-dimensional REAL
C          array, dimensioned Z(NM,N).
C
C     On Output
C
C        A has been reduced further to a quasi-triangular matrix in
C          which all nonzero subdiagonal elements correspond to pairs
C          of complex eigenvalues.
C
C        B is still in upper triangular form, although its elements
C          have been altered.  B(N,1) is unaltered.
C
C        ALFR and ALFI contain the real and imaginary parts of the
C          diagonal elements of the triangular matrix that would be
C          obtained if A were reduced completely to triangular form
C          by unitary transformations.  Non-zero values of ALFI occur
C          in pairs, the first member positive and the second negative.
C          ALFR and ALFI are one-dimensional REAL arrays, dimensioned
C          ALFR(N) and ALFI(N).
C
C        BETA contains the diagonal elements of the corresponding B,
C          normalized to be real and non-negative.  The generalized
C          eigenvalues are then the ratios ((ALFR+I*ALFI)/BETA).
C          BETA is a one-dimensional REAL array, dimensioned BETA(N).
C
C        Z contains the product of the right hand transformations
C          (for all three steps) if MATZ has been set to .TRUE.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  QZVAL
C
      INTEGER I,J,N,EN,NA,NM,NN,ISW
      REAL A(NM,*),B(NM,*),ALFR(*),ALFI(*),BETA(*),Z(NM,*)
      REAL C,D,E,R,S,T,AN,A1,A2,BN,CQ,CZ,DI,DR,EI,TI,TR
      REAL U1,U2,V1,V2,A1I,A11,A12,A2I,A21,A22,B11,B12,B22
      REAL SQI,SQR,SSI,SSR,SZI,SZR,A11I,A11R,A12I,A12R
      REAL A22I,A22R,EPSB
      LOGICAL MATZ
C
C***FIRST EXECUTABLE STATEMENT  QZVAL
      EPSB = B(N,1)
      ISW = 1
C     .......... FIND EIGENVALUES OF QUASI-TRIANGULAR MATRICES.
C                FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 510 NN = 1, N
         EN = N + 1 - NN
         NA = EN - 1
         IF (ISW .EQ. 2) GO TO 505
         IF (EN .EQ. 1) GO TO 410
         IF (A(EN,NA) .NE. 0.0E0) GO TO 420
C     .......... 1-BY-1 BLOCK, ONE REAL ROOT ..........
  410    ALFR(EN) = A(EN,EN)
         IF (B(EN,EN) .LT. 0.0E0) ALFR(EN) = -ALFR(EN)
         BETA(EN) = ABS(B(EN,EN))
         ALFI(EN) = 0.0E0
         GO TO 510
C     .......... 2-BY-2 BLOCK ..........
  420    IF (ABS(B(NA,NA)) .LE. EPSB) GO TO 455
         IF (ABS(B(EN,EN)) .GT. EPSB) GO TO 430
         A1 = A(EN,EN)
         A2 = A(EN,NA)
         BN = 0.0E0
         GO TO 435
  430    AN = ABS(A(NA,NA)) + ABS(A(NA,EN)) + ABS(A(EN,NA))
     1      + ABS(A(EN,EN))
         BN = ABS(B(NA,NA)) + ABS(B(NA,EN)) + ABS(B(EN,EN))
         A11 = A(NA,NA) / AN
         A12 = A(NA,EN) / AN
         A21 = A(EN,NA) / AN
         A22 = A(EN,EN) / AN
         B11 = B(NA,NA) / BN
         B12 = B(NA,EN) / BN
         B22 = B(EN,EN) / BN
         E = A11 / B11
         EI = A22 / B22
         S = A21 / (B11 * B22)
         T = (A22 - E * B22) / B22
         IF (ABS(E) .LE. ABS(EI)) GO TO 431
         E = EI
         T = (A11 - E * B11) / B11
  431    C = 0.5E0 * (T - S * B12)
         D = C * C + S * (A12 - E * B12)
         IF (D .LT. 0.0E0) GO TO 480
C     .......... TWO REAL ROOTS.
C                ZERO BOTH A(EN,NA) AND B(EN,NA) ..........
         E = E + (C + SIGN(SQRT(D),C))
         A11 = A11 - E * B11
         A12 = A12 - E * B12
         A22 = A22 - E * B22
         IF (ABS(A11) + ABS(A12) .LT.
     1       ABS(A21) + ABS(A22)) GO TO 432
         A1 = A12
         A2 = A11
         GO TO 435
  432    A1 = A22
         A2 = A21
C     .......... CHOOSE AND APPLY REAL Z ..........
  435    S = ABS(A1) + ABS(A2)
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 440 I = 1, EN
            T = A(I,EN) + U2 * A(I,NA)
            A(I,EN) = A(I,EN) + T * V1
            A(I,NA) = A(I,NA) + T * V2
            T = B(I,EN) + U2 * B(I,NA)
            B(I,EN) = B(I,EN) + T * V1
            B(I,NA) = B(I,NA) + T * V2
  440    CONTINUE
C
         IF (.NOT. MATZ) GO TO 450
C
         DO 445 I = 1, N
            T = Z(I,EN) + U2 * Z(I,NA)
            Z(I,EN) = Z(I,EN) + T * V1
            Z(I,NA) = Z(I,NA) + T * V2
  445    CONTINUE
C
  450    IF (BN .EQ. 0.0E0) GO TO 475
         IF (AN .LT. ABS(E) * BN) GO TO 455
         A1 = B(NA,NA)
         A2 = B(EN,NA)
         GO TO 460
  455    A1 = A(NA,NA)
         A2 = A(EN,NA)
C     .......... CHOOSE AND APPLY REAL Q ..........
  460    S = ABS(A1) + ABS(A2)
         IF (S .EQ. 0.0E0) GO TO 475
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 470 J = NA, N
            T = A(NA,J) + U2 * A(EN,J)
            A(NA,J) = A(NA,J) + T * V1
            A(EN,J) = A(EN,J) + T * V2
            T = B(NA,J) + U2 * B(EN,J)
            B(NA,J) = B(NA,J) + T * V1
            B(EN,J) = B(EN,J) + T * V2
  470    CONTINUE
C
  475    A(EN,NA) = 0.0E0
         B(EN,NA) = 0.0E0
         ALFR(NA) = A(NA,NA)
         ALFR(EN) = A(EN,EN)
         IF (B(NA,NA) .LT. 0.0E0) ALFR(NA) = -ALFR(NA)
         IF (B(EN,EN) .LT. 0.0E0) ALFR(EN) = -ALFR(EN)
         BETA(NA) = ABS(B(NA,NA))
         BETA(EN) = ABS(B(EN,EN))
         ALFI(EN) = 0.0E0
         ALFI(NA) = 0.0E0
         GO TO 505
C     .......... TWO COMPLEX ROOTS ..........
  480    E = E + C
         EI = SQRT(-D)
         A11R = A11 - E * B11
         A11I = EI * B11
         A12R = A12 - E * B12
         A12I = EI * B12
         A22R = A22 - E * B22
         A22I = EI * B22
         IF (ABS(A11R) + ABS(A11I) + ABS(A12R) + ABS(A12I) .LT.
     1       ABS(A21) + ABS(A22R) + ABS(A22I)) GO TO 482
         A1 = A12R
         A1I = A12I
         A2 = -A11R
         A2I = -A11I
         GO TO 485
  482    A1 = A22R
         A1I = A22I
         A2 = -A21
         A2I = 0.0E0
C     .......... CHOOSE COMPLEX Z ..........
  485    CZ = SQRT(A1*A1+A1I*A1I)
         IF (CZ .EQ. 0.0E0) GO TO 487
         SZR = (A1 * A2 + A1I * A2I) / CZ
         SZI = (A1 * A2I - A1I * A2) / CZ
         R = SQRT(CZ*CZ+SZR*SZR+SZI*SZI)
         CZ = CZ / R
         SZR = SZR / R
         SZI = SZI / R
         GO TO 490
  487    SZR = 1.0E0
         SZI = 0.0E0
  490    IF (AN .LT. (ABS(E) + EI) * BN) GO TO 492
         A1 = CZ * B11 + SZR * B12
         A1I = SZI * B12
         A2 = SZR * B22
         A2I = SZI * B22
         GO TO 495
  492    A1 = CZ * A11 + SZR * A12
         A1I = SZI * A12
         A2 = CZ * A21 + SZR * A22
         A2I = SZI * A22
C     .......... CHOOSE COMPLEX Q ..........
  495    CQ = SQRT(A1*A1+A1I*A1I)
         IF (CQ .EQ. 0.0E0) GO TO 497
         SQR = (A1 * A2 + A1I * A2I) / CQ
         SQI = (A1 * A2I - A1I * A2) / CQ
         R = SQRT(CQ*CQ+SQR*SQR+SQI*SQI)
         CQ = CQ / R
         SQR = SQR / R
         SQI = SQI / R
         GO TO 500
  497    SQR = 1.0E0
         SQI = 0.0E0
C     .......... COMPUTE DIAGONAL ELEMENTS THAT WOULD RESULT
C                IF TRANSFORMATIONS WERE APPLIED ..........
  500    SSR = SQR * SZR + SQI * SZI
         SSI = SQR * SZI - SQI * SZR
         I = 1
         TR = CQ * CZ * A11 + CQ * SZR * A12 + SQR * CZ * A21
     1      + SSR * A22
         TI = CQ * SZI * A12 - SQI * CZ * A21 + SSI * A22
         DR = CQ * CZ * B11 + CQ * SZR * B12 + SSR * B22
         DI = CQ * SZI * B12 + SSI * B22
         GO TO 503
  502    I = 2
         TR = SSR * A11 - SQR * CZ * A12 - CQ * SZR * A21
     1      + CQ * CZ * A22
         TI = -SSI * A11 - SQI * CZ * A12 + CQ * SZI * A21
         DR = SSR * B11 - SQR * CZ * B12 + CQ * CZ * B22
         DI = -SSI * B11 - SQI * CZ * B12
  503    T = TI * DR - TR * DI
         J = NA
         IF (T .LT. 0.0E0) J = EN
         R = SQRT(DR*DR+DI*DI)
         BETA(J) = BN * R
         ALFR(J) = AN * (TR * DR + TI * DI) / R
         ALFI(J) = AN * T / R
         IF (I .EQ. 1) GO TO 502
  505    ISW = 3 - ISW
  510 CONTINUE
C
      RETURN
      END
