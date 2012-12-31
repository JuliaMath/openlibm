*DECK HQR2
      SUBROUTINE HQR2 (NM, N, LOW, IGH, H, WR, WI, Z, IERR)
C***BEGIN PROLOGUE  HQR2
C***PURPOSE  Compute the eigenvalues and eigenvectors of a real upper
C            Hessenberg matrix using QR method.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C2B
C***TYPE      SINGLE PRECISION (HQR2-S, COMQR2-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure HQR2,
C     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a REAL UPPER Hessenberg matrix by the QR method.  The
C     eigenvectors of a REAL GENERAL matrix can also be found
C     if  ELMHES  and  ELTRAN  or  ORTHES  and  ORTRAN  have
C     been used to reduce this general matrix to Hessenberg form
C     and to accumulate the similarity transformations.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, H and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix H.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        LOW and IGH are two INTEGER variables determined by the
C          balancing subroutine  BALANC.  If  BALANC  has not been
C          used, set LOW=1 and IGH equal to the order of the matrix, N.
C
C        H contains the upper Hessenberg matrix.  H is a two-dimensional
C          REAL array, dimensioned H(NM,N).
C
C        Z contains the transformation matrix produced by  ELTRAN
C          after the reduction by  ELMHES, or by  ORTRAN  after the
C          reduction by  ORTHES, if performed.  If the eigenvectors
C          of the Hessenberg matrix are desired, Z must contain the
C          identity matrix.  Z is a two-dimensional REAL array,
C          dimensioned Z(NM,M).
C
C     On OUTPUT
C
C        H has been destroyed.
C
C        WR and WI contain the real and imaginary parts, respectively,
C          of the eigenvalues.  The eigenvalues are unordered except
C          that complex conjugate pairs of values appear consecutively
C          with the eigenvalue having the positive imaginary part first.
C          If an error exit is made, the eigenvalues should be correct
C          for indices IERR+1, IERR+2, ..., N.  WR and WI are one-
C          dimensional REAL arrays, dimensioned WR(N) and WI(N).
C
C        Z contains the real and imaginary parts of the eigenvectors.
C          If the J-th eigenvalue is real, the J-th column of Z
C          contains its eigenvector.  If the J-th eigenvalue is complex
C          with positive imaginary part, the J-th and (J+1)-th
C          columns of Z contain the real and imaginary parts of its
C          eigenvector.  The eigenvectors are unnormalized.  If an
C          error exit is made, none of the eigenvectors has been found.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after a total of 30*N iterations.
C                     The eigenvalues should be correct for indices
C                     IERR+1, IERR+2, ..., N, but no eigenvectors are
C                     computed.
C
C     Calls CDIV for complex division.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  CDIV
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  HQR2
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN
      INTEGER IGH,ITN,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,*),WR(*),WI(*),Z(NM,*)
      REAL P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,S1,S2
      LOGICAL NOTLAS
C
C***FIRST EXECUTABLE STATEMENT  HQR2
      IERR = 0
      NORM = 0.0E0
      K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
      DO 50 I = 1, N
C
         DO 40 J = K, N
   40    NORM = NORM + ABS(H(I,J))
C
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0E0
   50 CONTINUE
C
      EN = IGH
      T = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
   60 IF (EN .LT. LOW) GO TO 340
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S .EQ. 0.0E0) S = NORM
         S2 = S + ABS(H(L,L-1))
         IF (S2 .EQ. S) GO TO 100
   80 CONTINUE
C     .......... FORM SHIFT ..........
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     .......... FORM EXCEPTIONAL SHIFT ..........
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75E0 * S
      Y = X
      W = -0.4375E0 * S * S
  130 ITS = ITS + 1
      ITN = ITN - 1
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         S1 = ABS(P) * (ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
         S2 = S1 + ABS(H(M,M-1)) * (ABS(Q) + ABS(R))
         IF (S2 .EQ. S1) GO TO 150
  140 CONTINUE
C
  150 MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0E0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0E0
  160 CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0E0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0E0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
C     .......... ROW MODIFICATION ..........
         DO 210 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            IF (.NOT. NOTLAS) GO TO 200
            P = P + R * H(K+2,J)
            H(K+2,J) = H(K+2,J) - P * ZZ
  200       H(K+1,J) = H(K+1,J) - P * Y
            H(K,J) = H(K,J) - P * X
  210    CONTINUE
C
         J = MIN(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 230 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            IF (.NOT. NOTLAS) GO TO 220
            P = P + ZZ * H(I,K+2)
            H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
            H(I,K) = H(I,K) - P
  230    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
         DO 250 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1)
            IF (.NOT. NOTLAS) GO TO 240
            P = P + ZZ * Z(I,K+2)
            Z(I,K+2) = Z(I,K+2) - P * R
  240       Z(I,K+1) = Z(I,K+1) - P * Q
            Z(I,K) = Z(I,K) - P
  250    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     .......... ONE ROOT FOUND ..........
  270 H(EN,EN) = X + T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0E0
      EN = NA
      GO TO 60
C     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0E0
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0E0) GO TO 320
C     .......... REAL PAIR ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0E0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0E0
      WI(EN) = 0.0E0
      X = H(EN,NA)
      S = ABS(X) + ABS(ZZ)
      P = X / S
      Q = ZZ / S
      R = SQRT(P*P+Q*Q)
      P = P / R
      Q = Q / R
C     .......... ROW MODIFICATION ..........
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
C     .......... COLUMN MODIFICATION ..........
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
      DO 310 I = LOW, IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q * ZZ + P * Z(I,EN)
         Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
C
      GO TO 330
C     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  340 IF (NORM .EQ. 0.0E0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 800 NN = 1, N
         EN = N + 1 - NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q) 710, 600, 800
C     .......... REAL VECTOR ..........
  600    M = EN
         H(EN,EN) = 1.0E0
         IF (NA .EQ. 0) GO TO 800
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 700 II = 1, NA
            I = EN - II
            W = H(I,I) - P
            R = H(I,EN)
            IF (M .GT. NA) GO TO 620
C
            DO 610 J = M, NA
  610       R = R + H(I,J) * H(J,EN)
C
  620       IF (WI(I) .GE. 0.0E0) GO TO 630
            ZZ = W
            S = R
            GO TO 700
  630       M = I
            IF (WI(I) .NE. 0.0E0) GO TO 640
            T = W
            IF (T .NE. 0.0E0) GO TO 635
            T = NORM
  632       T = 0.5E0*T
            IF (NORM + T .GT. NORM) GO TO 632
            T = 2.0E0*T
  635       H(I,EN) = -R / T
            GO TO 700
C     .......... SOLVE REAL EQUATIONS ..........
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
            T = (X * S - ZZ * R) / Q
            H(I,EN) = T
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            H(I+1,EN) = (-R - W * T) / X
            GO TO 700
  650       H(I+1,EN) = (-S - Y * T) / ZZ
  700    CONTINUE
C     .......... END REAL VECTOR ..........
         GO TO 800
C     .......... COMPLEX VECTOR ..........
  710    M = NA
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
         IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GO TO 720
         H(NA,NA) = Q / H(EN,NA)
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
         GO TO 730
  720    CALL CDIV(0.0E0,-H(NA,EN),H(NA,NA)-P,Q,H(NA,NA),H(NA,EN))
  730    H(EN,NA) = 0.0E0
         H(EN,EN) = 1.0E0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 800
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
         DO 790 II = 1, ENM2
            I = NA - II
            W = H(I,I) - P
            RA = 0.0E0
            SA = H(I,EN)
C
            DO 760 J = M, NA
               RA = RA + H(I,J) * H(J,NA)
               SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
C
            IF (WI(I) .GE. 0.0E0) GO TO 770
            ZZ = W
            R = RA
            S = SA
            GO TO 790
  770       M = I
            IF (WI(I) .NE. 0.0E0) GO TO 780
            CALL CDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
            GO TO 790
C     .......... SOLVE COMPLEX EQUATIONS ..........
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
            VI = (WR(I) - P) * 2.0E0 * Q
            IF (VR .NE. 0.0E0 .OR. VI .NE. 0.0E0) GO TO 783
            S1 = NORM * (ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(ZZ))
            VR = S1
  782       VR = 0.5E0*VR
            IF (S1 + VR .GT. S1) GO TO 782
            VR = 2.0E0*VR
  783       CALL CDIV(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI,
     1                H(I,NA),H(I,EN))
            IF (ABS(X) .LE. ABS(ZZ) + ABS(Q)) GO TO 785
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
            GO TO 790
  785       CALL CDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q,
     1                H(I+1,NA),H(I+1,EN))
  790    CONTINUE
C     .......... END COMPLEX VECTOR ..........
  800 CONTINUE
C     .......... END BACK SUBSTITUTION.
C                VECTORS OF ISOLATED ROOTS ..........
      DO 840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
C
         DO 820 J = I, N
  820    Z(I,J) = H(I,J)
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW DO -- ..........
      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZ = 0.0E0
C
            DO 860 K = LOW, M
  860       ZZ = ZZ + Z(I,K) * H(K,J)
C
            Z(I,J) = ZZ
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
