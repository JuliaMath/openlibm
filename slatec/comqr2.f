*DECK COMQR2
      SUBROUTINE COMQR2 (NM, N, LOW, IGH, ORTR, ORTI, HR, HI, WR, WI,
     +   ZR, ZI, IERR)
C***BEGIN PROLOGUE  COMQR2
C***PURPOSE  Compute the eigenvalues and eigenvectors of a complex upper
C            Hessenberg matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C2B
C***TYPE      COMPLEX (HQR2-S, COMQR2-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of a unitary analogue of the
C     ALGOL procedure  COMLR2, NUM. MATH. 16, 181-204(1970) by Peters
C     and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C     The unitary analogue substitutes the QR algorithm of Francis
C     (COMP. JOUR. 4, 332-345(1962)) for the LR algorithm.
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a COMPLEX UPPER Hessenberg matrix by the QR
C     method.  The eigenvectors of a COMPLEX GENERAL matrix
C     can also be found if  CORTH  has been used to reduce
C     this general matrix to Hessenberg form.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, HR, HI, ZR, and ZI, as declared in the
C          calling program dimension statement.  NM is an INTEGER
C          variable.
C
C        N is the order of the matrix H=(HR,HI).  N is an INTEGER
C          variable.  N must be less than or equal to NM.
C
C        LOW and IGH are two INTEGER variables determined by the
C          balancing subroutine  CBAL.  If  CBAL  has not been used,
C          set LOW=1 and IGH equal to the order of the matrix, N.
C
C        ORTR and ORTI contain information about the unitary trans-
C          formations used in the reduction by  CORTH, if performed.
C          Only elements LOW through IGH are used.  If the eigenvectors
C          of the Hessenberg matrix are desired, set ORTR(J) and
C          ORTI(J) to 0.0E0 for these elements.  ORTR and ORTI are
C          one-dimensional REAL arrays, dimensioned ORTR(IGH) and
C          ORTI(IGH).
C
C        HR and HI contain the real and imaginary parts, respectively,
C          of the complex upper Hessenberg matrix.  Their lower
C          triangles below the subdiagonal contain information about
C          the unitary transformations used in the reduction by  CORTH,
C          if performed.  If the eigenvectors of the Hessenberg matrix
C          are desired, these elements may be arbitrary.  HR and HI
C          are two-dimensional REAL arrays, dimensioned HR(NM,N) and
C          HI(NM,N).
C
C     On OUTPUT
C
C        ORTR, ORTI, and the upper Hessenberg portions of HR and HI
C          have been destroyed.
C
C        WR and WI contain the real and imaginary parts, respectively,
C          of the eigenvalues of the upper Hessenberg matrix.  If an
C          error exit is made, the eigenvalues should be correct for
C          indices IERR+1, IERR+2, ..., N.  WR and WI are one-
C          dimensional REAL arrays, dimensioned WR(N) and WI(N).
C
C        ZR and ZI contain the real and imaginary parts, respectively,
C          of the eigenvectors.  The eigenvectors are unnormalized.
C          If an error exit is made, none of the eigenvectors has been
C          found.  ZR and ZI are two-dimensional REAL arrays,
C          dimensioned ZR(NM,N) and ZI(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after a total of 30*N iterations.
C                     The eigenvalues should be correct for indices
C                     IERR+1, IERR+2, ..., N, but no eigenvectors are
C                     computed.
C
C     Calls CSROOT for complex square root.
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
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
C***ROUTINES CALLED  CDIV, CSROOT, PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  COMQR2
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,NM,NN,IGH,IP1
      INTEGER ITN,ITS,LOW,LP1,ENM1,IEND,IERR
      REAL HR(NM,*),HI(NM,*),WR(*),WI(*),ZR(NM,*),ZI(NM,*)
      REAL ORTR(*),ORTI(*)
      REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,S1,S2
      REAL PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  COMQR2
      IERR = 0
C     .......... INITIALIZE EIGENVECTOR MATRIX ..........
      DO 100 I = 1, N
C
         DO 100 J = 1, N
            ZR(I,J) = 0.0E0
            ZI(I,J) = 0.0E0
            IF (I .EQ. J) ZR(I,J) = 1.0E0
  100 CONTINUE
C     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
C                FROM THE INFORMATION LEFT BY CORTH ..........
      IEND = IGH - LOW - 1
      IF (IEND) 180, 150, 105
C     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  105 DO 140 II = 1, IEND
         I = IGH - II
         IF (ORTR(I) .EQ. 0.0E0 .AND. ORTI(I) .EQ. 0.0E0) GO TO 140
         IF (HR(I,I-1) .EQ. 0.0E0 .AND. HI(I,I-1) .EQ. 0.0E0) GO TO 140
C     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
         NORM = HR(I,I-1) * ORTR(I) + HI(I,I-1) * ORTI(I)
         IP1 = I + 1
C
         DO 110 K = IP1, IGH
            ORTR(K) = HR(K,I-1)
            ORTI(K) = HI(K,I-1)
  110    CONTINUE
C
         DO 130 J = I, IGH
            SR = 0.0E0
            SI = 0.0E0
C
            DO 115 K = I, IGH
               SR = SR + ORTR(K) * ZR(K,J) + ORTI(K) * ZI(K,J)
               SI = SI + ORTR(K) * ZI(K,J) - ORTI(K) * ZR(K,J)
  115       CONTINUE
C
            SR = SR / NORM
            SI = SI / NORM
C
            DO 120 K = I, IGH
               ZR(K,J) = ZR(K,J) + SR * ORTR(K) - SI * ORTI(K)
               ZI(K,J) = ZI(K,J) + SR * ORTI(K) + SI * ORTR(K)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  150 L = LOW + 1
C
      DO 170 I = L, IGH
         LL = MIN(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0E0) GO TO 170
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0E0
C
         DO 155 J = I, N
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  155    CONTINUE
C
         DO 160 J = 1, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  160    CONTINUE
C
         DO 165 J = LOW, IGH
            SI = YR * ZI(J,I) + YI * ZR(J,I)
            ZR(J,I) = YR * ZR(J,I) - YI * ZI(J,I)
            ZI(J,I) = SI
  165    CONTINUE
C
  170 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0E0
      TI = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 680
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         S1 = ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1))
     1             + ABS(HR(L,L)) +ABS(HI(L,L))
         S2 = S1 + ABS(HR(L,L-1))
         IF (S2 .EQ. S1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0E0
      YI = (HI(ENM1,ENM1) - SI) / 2.0E0
      CALL CSROOT(YR**2-YI**2+XR,2.0E0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0E0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = 0.0E0
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
      DO 500 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0E0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0E0
         HI(I,I-1) = SR / NORM
C
         DO 490 J = I, N
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
C
  500 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0E0) GO TO 540
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0E0
      IF (EN .EQ. N) GO TO 540
      IP1 = EN + 1
C
      DO 520 J = IP1, N
         YR = HR(EN,J)
         YI = HI(EN,J)
         HR(EN,J) = SR * YR + SI * YI
         HI(EN,J) = SR * YI - SI * YR
  520 CONTINUE
C     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 580 I = 1, J
            YR = HR(I,J-1)
            YI = 0.0E0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GO TO 560
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
C
         DO 590 I = LOW, IGH
            YR = ZR(I,J-1)
            YI = ZI(I,J-1)
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
            ZR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            ZI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
            ZR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            ZI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  590    CONTINUE
C
  600 CONTINUE
C
      IF (SI .EQ. 0.0E0) GO TO 240
C
      DO 630 I = 1, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
C
      DO 640 I = LOW, IGH
         YR = ZR(I,EN)
         YI = ZI(I,EN)
         ZR(I,EN) = SR * YR - SI * YI
         ZI(I,EN) = SR * YI + SI * YR
  640 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GO TO 220
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  680 NORM = 0.0E0
C
      DO 720 I = 1, N
C
         DO 720 J = I, N
            NORM = NORM + ABS(HR(I,J)) + ABS(HI(I,J))
  720 CONTINUE
C
      IF (N .EQ. 1 .OR. NORM .EQ. 0.0E0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
      DO 800 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         ENM1 = EN - 1
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 780 II = 1, ENM1
            I = EN - II
            ZZR = HR(I,EN)
            ZZI = HI(I,EN)
            IF (I .EQ. ENM1) GO TO 760
            IP1 = I + 1
C
            DO 740 J = IP1, ENM1
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
C
  760       YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR .NE. 0.0E0 .OR. YI .NE. 0.0E0) GO TO 775
            YR = NORM
  770       YR = 0.5E0*YR
            IF (NORM + YR .GT. NORM) GO TO 770
            YR = 2.0E0*YR
  775       CALL CDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN))
  780    CONTINUE
C
  800 CONTINUE
C     .......... END BACKSUBSTITUTION ..........
      ENM1 = N - 1
C     .......... VECTORS OF ISOLATED ROOTS ..........
      DO  840 I = 1, ENM1
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
         IP1 = I + 1
C
         DO 820 J = IP1, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  820    CONTINUE
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW+1 DO -- ..........
      DO 880 JJ = LOW, ENM1
         J = N + LOW - JJ
         M = MIN(J-1,IGH)
C
         DO 880 I = LOW, IGH
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
C
            DO 860 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
C
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
