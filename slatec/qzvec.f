*DECK QZVEC
      SUBROUTINE QZVEC (NM, N, A, B, ALFR, ALFI, BETA, Z)
C***BEGIN PROLOGUE  QZVEC
C***PURPOSE  The optional fourth step of the QZ algorithm for
C            generalized eigenproblems.  Accepts a matrix in
C            quasi-triangular form and another in upper triangular
C            and computes the eigenvectors of the triangular problem
C            and transforms them back to the original coordinates
C            Usually preceded by QZHES, QZIT, and QZVAL.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C3
C***TYPE      SINGLE PRECISION (QZVEC-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is the optional fourth step of the QZ algorithm
C     for solving generalized matrix eigenvalue problems,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART.
C
C     This subroutine accepts a pair of REAL matrices, one of them in
C     quasi-triangular form (in which each 2-by-2 block corresponds to
C     a pair of complex eigenvalues) and the other in upper triangular
C     form.  It computes the eigenvectors of the triangular problem and
C     transforms the results back to the original coordinate system.
C     It is usually preceded by  QZHES,  QZIT, and  QZVAL.
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
C        ALFR, ALFI, and BETA are one-dimensional REAL arrays with
C          components whose ratios ((ALFR+I*ALFI)/BETA) are the
C          generalized eigenvalues.  They are usually obtained from
C          QZVAL.  They are dimensioned ALFR(N), ALFI(N), and BETA(N).
C
C        Z contains the transformation matrix produced in the reductions
C          by  QZHES,  QZIT, and  QZVAL,  if performed.  If the
C          eigenvectors of the triangular problem are desired, Z must
C          contain the identity matrix.  Z is a two-dimensional REAL
C          array, dimensioned Z(NM,N).
C
C     On Output
C
C        A is unaltered.  Its subdiagonal elements provide information
C           about the storage of the complex eigenvectors.
C
C        B has been destroyed.
C
C        ALFR, ALFI, and BETA are unaltered.
C
C        Z contains the real and imaginary parts of the eigenvectors.
C          If ALFI(J) .EQ. 0.0, the J-th eigenvalue is real and
C            the J-th column of Z contains its eigenvector.
C          If ALFI(J) .NE. 0.0, the J-th eigenvalue is complex.
C            If ALFI(J) .GT. 0.0, the eigenvalue is the first of
C              a complex pair and the J-th and (J+1)-th columns
C              of Z contain its eigenvector.
C            If ALFI(J) .LT. 0.0, the eigenvalue is the second of
C              a complex pair and the (J-1)-th and J-th columns
C              of Z contain the conjugate of its eigenvector.
C          Each eigenvector is normalized so that the modulus
C          of its largest component is 1.0 .
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
C***END PROLOGUE  QZVEC
C
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN,ISW,ENM2
      REAL A(NM,*),B(NM,*),ALFR(*),ALFI(*),BETA(*),Z(NM,*)
      REAL D,Q,R,S,T,W,X,Y,DI,DR,RA,RR,SA,TI,TR,T1,T2
      REAL W1,X1,ZZ,Z1,ALFM,ALMI,ALMR,BETM,EPSB
C
C***FIRST EXECUTABLE STATEMENT  QZVEC
      EPSB = B(N,1)
      ISW = 1
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 800 NN = 1, N
         EN = N + 1 - NN
         NA = EN - 1
         IF (ISW .EQ. 2) GO TO 795
         IF (ALFI(EN) .NE. 0.0E0) GO TO 710
C     .......... REAL VECTOR ..........
         M = EN
         B(EN,EN) = 1.0E0
         IF (NA .EQ. 0) GO TO 800
         ALFM = ALFR(M)
         BETM = BETA(M)
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 700 II = 1, NA
            I = EN - II
            W = BETM * A(I,I) - ALFM * B(I,I)
            R = 0.0E0
C
            DO 610 J = M, EN
  610       R = R + (BETM * A(I,J) - ALFM * B(I,J)) * B(J,EN)
C
            IF (I .EQ. 1 .OR. ISW .EQ. 2) GO TO 630
            IF (BETM * A(I,I-1) .EQ. 0.0E0) GO TO 630
            ZZ = W
            S = R
            GO TO 690
  630       M = I
            IF (ISW .EQ. 2) GO TO 640
C     .......... REAL 1-BY-1 BLOCK ..........
            T = W
            IF (W .EQ. 0.0E0) T = EPSB
            B(I,EN) = -R / T
            GO TO 700
C     .......... REAL 2-BY-2 BLOCK ..........
  640       X = BETM * A(I,I+1) - ALFM * B(I,I+1)
            Y = BETM * A(I+1,I)
            Q = W * ZZ - X * Y
            T = (X * S - ZZ * R) / Q
            B(I,EN) = T
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            B(I+1,EN) = (-R - W * T) / X
            GO TO 690
  650       B(I+1,EN) = (-S - Y * T) / ZZ
  690       ISW = 3 - ISW
  700    CONTINUE
C     .......... END REAL VECTOR ..........
         GO TO 800
C     .......... COMPLEX VECTOR ..........
  710    M = NA
         ALMR = ALFR(M)
         ALMI = ALFI(M)
         BETM = BETA(M)
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
         Y = BETM * A(EN,NA)
         B(NA,NA) = -ALMI * B(EN,EN) / Y
         B(NA,EN) = (ALMR * B(EN,EN) - BETM * A(EN,EN)) / Y
         B(EN,NA) = 0.0E0
         B(EN,EN) = 1.0E0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 795
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
         DO 790 II = 1, ENM2
            I = NA - II
            W = BETM * A(I,I) - ALMR * B(I,I)
            W1 = -ALMI * B(I,I)
            RA = 0.0E0
            SA = 0.0E0
C
            DO 760 J = M, EN
               X = BETM * A(I,J) - ALMR * B(I,J)
               X1 = -ALMI * B(I,J)
               RA = RA + X * B(J,NA) - X1 * B(J,EN)
               SA = SA + X * B(J,EN) + X1 * B(J,NA)
  760       CONTINUE
C
            IF (I .EQ. 1 .OR. ISW .EQ. 2) GO TO 770
            IF (BETM * A(I,I-1) .EQ. 0.0E0) GO TO 770
            ZZ = W
            Z1 = W1
            R = RA
            S = SA
            ISW = 2
            GO TO 790
  770       M = I
            IF (ISW .EQ. 2) GO TO 780
C     .......... COMPLEX 1-BY-1 BLOCK ..........
            TR = -RA
            TI = -SA
  773       DR = W
            DI = W1
C     .......... COMPLEX DIVIDE (T1,T2) = (TR,TI) / (DR,DI) ..........
  775       IF (ABS(DI) .GT. ABS(DR)) GO TO 777
            RR = DI / DR
            D = DR + DI * RR
            T1 = (TR + TI * RR) / D
            T2 = (TI - TR * RR) / D
            GO TO (787,782), ISW
  777       RR = DR / DI
            D = DR * RR + DI
            T1 = (TR * RR + TI) / D
            T2 = (TI * RR - TR) / D
            GO TO (787,782), ISW
C     .......... COMPLEX 2-BY-2 BLOCK ..........
  780       X = BETM * A(I,I+1) - ALMR * B(I,I+1)
            X1 = -ALMI * B(I,I+1)
            Y = BETM * A(I+1,I)
            TR = Y * RA - W * R + W1 * S
            TI = Y * SA - W * S - W1 * R
            DR = W * ZZ - W1 * Z1 - X * Y
            DI = W * Z1 + W1 * ZZ - X1 * Y
            IF (DR .EQ. 0.0E0 .AND. DI .EQ. 0.0E0) DR = EPSB
            GO TO 775
  782       B(I+1,NA) = T1
            B(I+1,EN) = T2
            ISW = 1
            IF (ABS(Y) .GT. ABS(W) + ABS(W1)) GO TO 785
            TR = -RA - X * B(I+1,NA) + X1 * B(I+1,EN)
            TI = -SA - X * B(I+1,EN) - X1 * B(I+1,NA)
            GO TO 773
  785       T1 = (-R - ZZ * B(I+1,NA) + Z1 * B(I+1,EN)) / Y
            T2 = (-S - ZZ * B(I+1,EN) - Z1 * B(I+1,NA)) / Y
  787       B(I,NA) = T1
            B(I,EN) = T2
  790    CONTINUE
C     .......... END COMPLEX VECTOR ..........
  795    ISW = 3 - ISW
  800 CONTINUE
C     .......... END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 1 DO -- ..........
      DO 880 JJ = 1, N
         J = N + 1 - JJ
C
         DO 880 I = 1, N
            ZZ = 0.0E0
C
            DO 860 K = 1, J
  860       ZZ = ZZ + Z(I,K) * B(K,J)
C
            Z(I,J) = ZZ
  880 CONTINUE
C     .......... NORMALIZE SO THAT MODULUS OF LARGEST
C                COMPONENT OF EACH VECTOR IS 1.
C                (ISW IS 1 INITIALLY FROM BEFORE) ..........
      DO 950 J = 1, N
         D = 0.0E0
         IF (ISW .EQ. 2) GO TO 920
         IF (ALFI(J) .NE. 0.0E0) GO TO 945
C
         DO 890 I = 1, N
            IF (ABS(Z(I,J)) .GT. D) D = ABS(Z(I,J))
  890    CONTINUE
C
         DO 900 I = 1, N
  900    Z(I,J) = Z(I,J) / D
C
         GO TO 950
C
  920    DO 930 I = 1, N
            R = ABS(Z(I,J-1)) + ABS(Z(I,J))
            IF (R .NE. 0.0E0) R = R * SQRT((Z(I,J-1)/R)**2
     1                                     +(Z(I,J)/R)**2)
            IF (R .GT. D) D = R
  930    CONTINUE
C
         DO 940 I = 1, N
            Z(I,J-1) = Z(I,J-1) / D
            Z(I,J) = Z(I,J) / D
  940    CONTINUE
C
  945    ISW = 3 - ISW
  950 CONTINUE
C
      RETURN
      END
