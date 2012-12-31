*DECK QZHES
      SUBROUTINE QZHES (NM, N, A, B, MATZ, Z)
C***BEGIN PROLOGUE  QZHES
C***PURPOSE  The first step of the QZ algorithm for solving generalized
C            matrix eigenproblems.  Accepts a pair of real general
C            matrices and reduces one of them to upper Hessenberg
C            and the other to upper triangular form using orthogonal
C            transformations. Usually followed by QZIT, QZVAL, QZVEC.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1B3
C***TYPE      SINGLE PRECISION (QZHES-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is the first step of the QZ algorithm
C     for solving generalized matrix eigenvalue problems,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART.
C
C     This subroutine accepts a pair of REAL GENERAL matrices and
C     reduces one of them to upper Hessenberg form and the other
C     to upper triangular form using orthogonal transformations.
C     It is usually followed by  QZIT,  QZVAL  and, possibly,  QZVEC.
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
C        A contains a real general matrix.  A is a two-dimensional
C          REAL array, dimensioned A(NM,N).
C
C        B contains a real general matrix.  B is a two-dimensional
C          REAL array, dimensioned B(NM,N).
C
C        MATZ should be set to .TRUE. if the right hand transformations
C          are to be accumulated for later use in computing
C          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL
C          variable.
C
C     On Output
C
C        A has been reduced to upper Hessenberg form.  The elements
C          below the first subdiagonal have been set to zero.
C
C        B has been reduced to upper triangular form.  The elements
C          below the main diagonal have been set to zero.
C
C        Z contains the product of the right hand transformations if
C          MATZ has been set to .TRUE.  Otherwise, Z is not referenced.
C          Z is a two-dimensional REAL array, dimensioned Z(NM,N).
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
C***END PROLOGUE  QZHES
C
      INTEGER I,J,K,L,N,LB,L1,NM,NK1,NM1,NM2
      REAL A(NM,*),B(NM,*),Z(NM,*)
      REAL R,S,T,U1,U2,V1,V2,RHO
      LOGICAL MATZ
C
C     .......... INITIALIZE Z ..........
C***FIRST EXECUTABLE STATEMENT  QZHES
      IF (.NOT. MATZ) GO TO 10
C
      DO 3 I = 1, N
C
         DO 2 J = 1, N
            Z(I,J) = 0.0E0
    2    CONTINUE
C
         Z(I,I) = 1.0E0
    3 CONTINUE
C     .......... REDUCE B TO UPPER TRIANGULAR FORM ..........
   10 IF (N .LE. 1) GO TO 170
      NM1 = N - 1
C
      DO 100 L = 1, NM1
         L1 = L + 1
         S = 0.0E0
C
         DO 20 I = L1, N
            S = S + ABS(B(I,L))
   20    CONTINUE
C
         IF (S .EQ. 0.0E0) GO TO 100
         S = S + ABS(B(L,L))
         R = 0.0E0
C
         DO 25 I = L, N
            B(I,L) = B(I,L) / S
            R = R + B(I,L)**2
   25    CONTINUE
C
         R = SIGN(SQRT(R),B(L,L))
         B(L,L) = B(L,L) + R
         RHO = R * B(L,L)
C
         DO 50 J = L1, N
            T = 0.0E0
C
            DO 30 I = L, N
               T = T + B(I,L) * B(I,J)
   30       CONTINUE
C
            T = -T / RHO
C
            DO 40 I = L, N
               B(I,J) = B(I,J) + T * B(I,L)
   40       CONTINUE
C
   50    CONTINUE
C
         DO 80 J = 1, N
            T = 0.0E0
C
            DO 60 I = L, N
               T = T + B(I,L) * A(I,J)
   60       CONTINUE
C
            T = -T / RHO
C
            DO 70 I = L, N
               A(I,J) = A(I,J) + T * B(I,L)
   70       CONTINUE
C
   80    CONTINUE
C
         B(L,L) = -S * R
C
         DO 90 I = L1, N
            B(I,L) = 0.0E0
   90    CONTINUE
C
  100 CONTINUE
C     .......... REDUCE A TO UPPER HESSENBERG FORM, WHILE
C                KEEPING B TRIANGULAR ..........
      IF (N .EQ. 2) GO TO 170
      NM2 = N - 2
C
      DO 160 K = 1, NM2
         NK1 = NM1 - K
C     .......... FOR L=N-1 STEP -1 UNTIL K+1 DO -- ..........
         DO 150 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C     .......... ZERO A(L+1,K) ..........
            S = ABS(A(L,K)) + ABS(A(L1,K))
            IF (S .EQ. 0.0E0) GO TO 150
            U1 = A(L,K) / S
            U2 = A(L1,K) / S
            R = SIGN(SQRT(U1*U1+U2*U2),U1)
            V1 =  -(U1 + R) / R
            V2 = -U2 / R
            U2 = V2 / V1
C
            DO 110 J = K, N
               T = A(L,J) + U2 * A(L1,J)
               A(L,J) = A(L,J) + T * V1
               A(L1,J) = A(L1,J) + T * V2
  110       CONTINUE
C
            A(L1,K) = 0.0E0
C
            DO 120 J = L, N
               T = B(L,J) + U2 * B(L1,J)
               B(L,J) = B(L,J) + T * V1
               B(L1,J) = B(L1,J) + T * V2
  120       CONTINUE
C     .......... ZERO B(L+1,L) ..........
            S = ABS(B(L1,L1)) + ABS(B(L1,L))
            IF (S .EQ. 0.0E0) GO TO 150
            U1 = B(L1,L1) / S
            U2 = B(L1,L) / S
            R = SIGN(SQRT(U1*U1+U2*U2),U1)
            V1 =  -(U1 + R) / R
            V2 = -U2 / R
            U2 = V2 / V1
C
            DO 130 I = 1, L1
               T = B(I,L1) + U2 * B(I,L)
               B(I,L1) = B(I,L1) + T * V1
               B(I,L) = B(I,L) + T * V2
  130       CONTINUE
C
            B(L1,L) = 0.0E0
C
            DO 140 I = 1, N
               T = A(I,L1) + U2 * A(I,L)
               A(I,L1) = A(I,L1) + T * V1
               A(I,L) = A(I,L) + T * V2
  140       CONTINUE
C
            IF (.NOT. MATZ) GO TO 150
C
            DO 145 I = 1, N
               T = Z(I,L1) + U2 * Z(I,L)
               Z(I,L1) = Z(I,L1) + T * V1
               Z(I,L) = Z(I,L) + T * V2
  145       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C
  170 RETURN
      END
