*DECK QZIT
      SUBROUTINE QZIT (NM, N, A, B, EPS1, MATZ, Z, IERR)
C***BEGIN PROLOGUE  QZIT
C***PURPOSE  The second step of the QZ algorithm for generalized
C            eigenproblems.  Accepts an upper Hessenberg and an upper
C            triangular matrix and reduces the former to
C            quasi-triangular form while preserving the form of the
C            latter.  Usually preceded by QZHES and followed by QZVAL
C            and QZVEC.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1B3
C***TYPE      SINGLE PRECISION (QZIT-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is the second step of the QZ algorithm
C     for solving generalized matrix eigenvalue problems,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART,
C     as modified in technical note NASA TN D-7305(1973) by WARD.
C
C     This subroutine accepts a pair of REAL matrices, one of them
C     in upper Hessenberg form and the other in upper triangular form.
C     It reduces the Hessenberg matrix to quasi-triangular form using
C     orthogonal transformations while maintaining the triangular form
C     of the other matrix.  It is usually preceded by  QZHES  and
C     followed by  QZVAL  and, possibly,  QZVEC.
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
C        A contains a real upper Hessenberg matrix.  A is a two-
C          dimensional REAL array, dimensioned A(NM,N).
C
C        B contains a real upper triangular matrix.  B is a two-
C          dimensional REAL array, dimensioned B(NM,N).
C
C        EPS1 is a tolerance used to determine negligible elements.
C          EPS1 = 0.0 (or negative) may be input, in which case an
C          element will be neglected only if it is less than roundoff
C          error times the norm of its matrix.  If the input EPS1 is
C          positive, then an element will be considered negligible
C          if it is less than EPS1 times the norm of its matrix.  A
C          positive value of EPS1 may result in faster execution,
C          but less accurate results.  EPS1 is a REAL variable.
C
C        MATZ should be set to .TRUE. if the right hand transformations
C          are to be accumulated for later use in computing
C          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL
C          variable.
C
C        Z contains, if MATZ has been set to .TRUE., the transformation
C          matrix produced in the reduction by  QZHES, if performed, or
C          else the identity matrix.  If MATZ has been set to .FALSE.,
C          Z is not referenced.  Z is a two-dimensional REAL array,
C          dimensioned Z(NM,N).
C
C     On Output
C
C        A has been reduced to quasi-triangular form.  The elements
C          below the first subdiagonal are still zero, and no two
C          consecutive subdiagonal elements are nonzero.
C
C        B is still in upper triangular form, although its elements
C          have been altered.  The location B(N,1) is used to store
C          EPS1 times the norm of B for later use by  QZVAL  and  QZVEC.
C
C        Z contains the product of the right hand transformations
C          (for both steps) if MATZ has been set to .TRUE.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if neither A(J,J-1) nor A(J-1,J-2) has become
C                     zero after a total of 30*N iterations.
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
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  QZIT
C
      INTEGER I,J,K,L,N,EN,K1,K2,LD,LL,L1,NA,NM,ISH,ITN,ITS,KM1,LM1
      INTEGER ENM2,IERR,LOR1,ENORN
      REAL A(NM,*),B(NM,*),Z(NM,*)
      REAL R,S,T,A1,A2,A3,EP,SH,U1,U2,U3,V1,V2,V3,ANI
      REAL A11,A12,A21,A22,A33,A34,A43,A44,BNI,B11
      REAL B12,B22,B33,B34,B44,EPSA,EPSB,EPS1,ANORM,BNORM
      LOGICAL MATZ,NOTLAS
C
C***FIRST EXECUTABLE STATEMENT  QZIT
      IERR = 0
C     .......... COMPUTE EPSA,EPSB ..........
      ANORM = 0.0E0
      BNORM = 0.0E0
C
      DO 30 I = 1, N
         ANI = 0.0E0
         IF (I .NE. 1) ANI = ABS(A(I,I-1))
         BNI = 0.0E0
C
         DO 20 J = I, N
            ANI = ANI + ABS(A(I,J))
            BNI = BNI + ABS(B(I,J))
   20    CONTINUE
C
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   30 CONTINUE
C
      IF (ANORM .EQ. 0.0E0) ANORM = 1.0E0
      IF (BNORM .EQ. 0.0E0) BNORM = 1.0E0
      EP = EPS1
      IF (EP .GT. 0.0E0) GO TO 50
C     .......... COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO ..........
      EP = 1.0E0
   40 EP = EP / 2.0E0
      IF (1.0E0 + EP .GT. 1.0E0) GO TO 40
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
C     .......... REDUCE A TO QUASI-TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR ..........
      LOR1 = 1
      ENORN = N
      EN = N
      ITN = 30*N
C     .......... BEGIN QZ STEP ..........
   60 IF (EN .LE. 2) GO TO 1001
      IF (.NOT. MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
   70 ISH = 2
C     .......... CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- ..........
      DO 80 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L .EQ. 1) GO TO 95
         IF (ABS(A(L,LM1)) .LE. EPSA) GO TO 90
   80 CONTINUE
C
   90 A(L,LM1) = 0.0E0
      IF (L .LT. NA) GO TO 95
C     .......... 1-BY-1 OR 2-BY-2 BLOCK ISOLATED ..........
      EN = LM1
      GO TO 60
C     .......... CHECK FOR SMALL TOP OF B ..........
   95 LD = L
  100 L1 = L + 1
      B11 = B(L,L)
      IF (ABS(B11) .GT. EPSB) GO TO 120
      B(L,L) = 0.0E0
      S = ABS(A(L,L)) + ABS(A(L1,L))
      U1 = A(L,L) / S
      U2 = A(L1,L) / S
      R = SIGN(SQRT(U1*U1+U2*U2),U1)
      V1 = -(U1 + R) / R
      V2 = -U2 / R
      U2 = V2 / V1
C
      DO 110 J = L, ENORN
         T = A(L,J) + U2 * A(L1,J)
         A(L,J) = A(L,J) + T * V1
         A(L1,J) = A(L1,J) + T * V2
         T = B(L,J) + U2 * B(L1,J)
         B(L,J) = B(L,J) + T * V1
         B(L1,J) = B(L1,J) + T * V2
  110 CONTINUE
C
      IF (L .NE. 1) A(L,LM1) = -A(L,LM1)
      LM1 = L
      L = L1
      GO TO 90
  120 A11 = A(L,L) / B11
      A21 = A(L1,L) / B11
      IF (ISH .EQ. 1) GO TO 140
C     .......... ITERATION STRATEGY ..........
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10) GO TO 155
C     .......... DETERMINE TYPE OF SHIFT ..........
      B22 = B(L1,L1)
      IF (ABS(B22) .LT. EPSB) B22 = EPSB
      B33 = B(NA,NA)
      IF (ABS(B33) .LT. EPSB) B33 = EPSB
      B44 = B(EN,EN)
      IF (ABS(B44) .LT. EPSB) B44 = EPSB
      A33 = A(NA,NA) / B33
      A34 = A(NA,EN) / B44
      A43 = A(EN,NA) / B33
      A44 = A(EN,EN) / B44
      B34 = B(NA,EN) / B44
      T = 0.5E0 * (A43 * B34 - A33 - A44)
      R = T * T + A34 * A43 - A33 * A44
      IF (R .LT. 0.0E0) GO TO 150
C     .......... DETERMINE SINGLE SHIFT ZEROTH COLUMN OF A ..........
      ISH = 1
      R = SQRT(R)
      SH = -T + R
      S = -T - R
      IF (ABS(S-A44) .LT. ABS(SH-A44)) SH = S
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS OF A.
C                FOR L=EN-2 STEP -1 UNTIL LD DO -- ..........
      DO 130 LL = LD, ENM2
         L = ENM2 + LD - LL
         IF (L .EQ. LD) GO TO 140
         LM1 = L - 1
         L1 = L + 1
         T = A(L,L)
         IF (ABS(B(L,L)) .GT. EPSB) T = T - SH * B(L,L)
         IF (ABS(A(L,LM1)) .LE. ABS(T/A(L1,L)) * EPSA) GO TO 100
  130 CONTINUE
C
  140 A1 = A11 - SH
      A2 = A21
      IF (L .NE. LD) A(L,LM1) = -A(L,LM1)
      GO TO 160
C     .......... DETERMINE DOUBLE SHIFT ZEROTH COLUMN OF A ..........
  150 A12 = A(L,L1) / B22
      A22 = A(L1,L1) / B22
      B12 = B(L,L1) / B22
      A1 = ((A33 - A11) * (A44 - A11) - A34 * A43 + A43 * B34 * A11)
     1     / A21 + A12 - A11 * B12
      A2 = (A22 - A11) - A21 * B12 - (A33 - A11) - (A44 - A11)
     1     + A43 * B34
      A3 = A(L1+1,L1) / B22
      GO TO 160
C     .......... AD HOC SHIFT ..........
  155 A1 = 0.0E0
      A2 = 1.0E0
      A3 = 1.1605E0
  160 ITS = ITS + 1
      ITN = ITN - 1
      IF (.NOT. MATZ) LOR1 = LD
C     .......... MAIN LOOP ..........
      DO 260 K = L, NA
         NOTLAS = K .NE. NA .AND. ISH .EQ. 2
         K1 = K + 1
         K2 = K + 2
         KM1 = MAX(K-1,L)
         LL = MIN(EN,K1+ISH)
         IF (NOTLAS) GO TO 190
C     .......... ZERO A(K+1,K-1) ..........
         IF (K .EQ. L) GO TO 170
         A1 = A(K,KM1)
         A2 = A(K1,KM1)
  170    S = ABS(A1) + ABS(A2)
         IF (S .EQ. 0.0E0) GO TO 70
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 180 J = KM1, ENORN
            T = A(K,J) + U2 * A(K1,J)
            A(K,J) = A(K,J) + T * V1
            A(K1,J) = A(K1,J) + T * V2
            T = B(K,J) + U2 * B(K1,J)
            B(K,J) = B(K,J) + T * V1
            B(K1,J) = B(K1,J) + T * V2
  180    CONTINUE
C
         IF (K .NE. L) A(K1,KM1) = 0.0E0
         GO TO 240
C     .......... ZERO A(K+1,K-1) AND A(K+2,K-1) ..........
  190    IF (K .EQ. L) GO TO 200
         A1 = A(K,KM1)
         A2 = A(K1,KM1)
         A3 = A(K2,KM1)
  200    S = ABS(A1) + ABS(A2) + ABS(A3)
         IF (S .EQ. 0.0E0) GO TO 260
         U1 = A1 / S
         U2 = A2 / S
         U3 = A3 / S
         R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         V3 = -U3 / R
         U2 = V2 / V1
         U3 = V3 / V1
C
         DO 210 J = KM1, ENORN
            T = A(K,J) + U2 * A(K1,J) + U3 * A(K2,J)
            A(K,J) = A(K,J) + T * V1
            A(K1,J) = A(K1,J) + T * V2
            A(K2,J) = A(K2,J) + T * V3
            T = B(K,J) + U2 * B(K1,J) + U3 * B(K2,J)
            B(K,J) = B(K,J) + T * V1
            B(K1,J) = B(K1,J) + T * V2
            B(K2,J) = B(K2,J) + T * V3
  210    CONTINUE
C
         IF (K .EQ. L) GO TO 220
         A(K1,KM1) = 0.0E0
         A(K2,KM1) = 0.0E0
C     .......... ZERO B(K+2,K+1) AND B(K+2,K) ..........
  220    S = ABS(B(K2,K2)) + ABS(B(K2,K1)) + ABS(B(K2,K))
         IF (S .EQ. 0.0E0) GO TO 240
         U1 = B(K2,K2) / S
         U2 = B(K2,K1) / S
         U3 = B(K2,K) / S
         R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         V3 = -U3 / R
         U2 = V2 / V1
         U3 = V3 / V1
C
         DO 230 I = LOR1, LL
            T = A(I,K2) + U2 * A(I,K1) + U3 * A(I,K)
            A(I,K2) = A(I,K2) + T * V1
            A(I,K1) = A(I,K1) + T * V2
            A(I,K) = A(I,K) + T * V3
            T = B(I,K2) + U2 * B(I,K1) + U3 * B(I,K)
            B(I,K2) = B(I,K2) + T * V1
            B(I,K1) = B(I,K1) + T * V2
            B(I,K) = B(I,K) + T * V3
  230    CONTINUE
C
         B(K2,K) = 0.0E0
         B(K2,K1) = 0.0E0
         IF (.NOT. MATZ) GO TO 240
C
         DO 235 I = 1, N
            T = Z(I,K2) + U2 * Z(I,K1) + U3 * Z(I,K)
            Z(I,K2) = Z(I,K2) + T * V1
            Z(I,K1) = Z(I,K1) + T * V2
            Z(I,K) = Z(I,K) + T * V3
  235    CONTINUE
C     .......... ZERO B(K+1,K) ..........
  240    S = ABS(B(K1,K1)) + ABS(B(K1,K))
         IF (S .EQ. 0.0E0) GO TO 260
         U1 = B(K1,K1) / S
         U2 = B(K1,K) / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 250 I = LOR1, LL
            T = A(I,K1) + U2 * A(I,K)
            A(I,K1) = A(I,K1) + T * V1
            A(I,K) = A(I,K) + T * V2
            T = B(I,K1) + U2 * B(I,K)
            B(I,K1) = B(I,K1) + T * V1
            B(I,K) = B(I,K) + T * V2
  250    CONTINUE
C
         B(K1,K) = 0.0E0
         IF (.NOT. MATZ) GO TO 260
C
         DO 255 I = 1, N
            T = Z(I,K1) + U2 * Z(I,K)
            Z(I,K1) = Z(I,K1) + T * V1
            Z(I,K) = Z(I,K) + T * V2
  255    CONTINUE
C
  260 CONTINUE
C     .......... END QZ STEP ..........
      GO TO 70
C     .......... SET ERROR -- NEITHER BOTTOM SUBDIAGONAL ELEMENT
C                HAS BECOME NEGLIGIBLE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
C     .......... SAVE EPSB FOR USE BY QZVAL AND QZVEC ..........
 1001 IF (N .GT. 1) B(N,1) = EPSB
      RETURN
      END
