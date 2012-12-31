*DECK COMLR2
      SUBROUTINE COMLR2 (NM, N, LOW, IGH, INT, HR, HI, WR, WI, ZR, ZI,
     +   IERR)
C***BEGIN PROLOGUE  COMLR2
C***PURPOSE  Compute the eigenvalues and eigenvectors of a complex upper
C            Hessenberg matrix using the modified LR method.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C2B
C***TYPE      COMPLEX (COMLR2-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK, LR METHOD
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure COMLR2,
C     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a COMPLEX UPPER Hessenberg matrix by the modified LR
C     method.  The eigenvectors of a COMPLEX GENERAL matrix
C     can also be found if  COMHES  has been used to reduce
C     this general matrix to Hessenberg form.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, HR, HI, ZR and ZI, as declared in the
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
C        INT contains information on the rows and columns
C          interchanged in the reduction by  COMHES, if performed.
C          Only elements LOW through IGH are used.  If you want the
C          eigenvectors of a complex general matrix, leave INT as it
C          came from  COMHES.  If the eigenvectors of the Hessenberg
C          matrix are desired, set INT(J)=J for these elements.  INT
C          is a one-dimensional INTEGER array, dimensioned INT(IGH).
C
C        HR and HI contain the real and imaginary parts, respectively,
C          of the complex upper Hessenberg matrix.  Their lower
C          triangles below the subdiagonal contain the multipliers
C          which were used in the reduction by  COMHES, if performed.
C          If the eigenvectors of a complex general matrix are
C          desired, leave these multipliers in the lower triangles.
C          If the eigenvectors of the Hessenberg matrix are desired,
C          these elements must be set to zero.  HR and HI are
C          two-dimensional REAL arrays, dimensioned HR(NM,N) and
C          HI(NM,N).
C
C     On OUTPUT
C
C        The upper Hessenberg portions of HR and HI have been
C          destroyed, but the location HR(1,1) contains the norm
C          of the triangularized matrix.
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
C***ROUTINES CALLED  CDIV, CSROOT
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  COMLR2
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NM,NN,IGH,IM1,IP1
      INTEGER ITN,ITS,LOW,MP1,ENM1,IEND,IERR
      REAL HR(NM,*),HI(NM,*),WR(*),WI(*),ZR(NM,*),ZI(NM,*)
      REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,S1,S2
      INTEGER INT(*)
C
C***FIRST EXECUTABLE STATEMENT  COMLR2
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
C                FROM THE INFORMATION LEFT BY COMHES ..........
      IEND = IGH - LOW - 1
      IF (IEND .LE. 0) GO TO 180
C     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 160 II = 1, IEND
         I = IGH - II
         IP1 = I + 1
C
         DO 120 K = IP1, IGH
            ZR(K,I) = HR(K,I-1)
            ZI(K,I) = HI(K,I-1)
  120    CONTINUE
C
         J = INT(I)
         IF (I .EQ. J) GO TO 160
C
         DO 140 K = I, IGH
            ZR(I,K) = ZR(J,K)
            ZI(I,K) = ZI(J,K)
            ZR(J,K) = 0.0E0
            ZI(J,K) = 0.0E0
  140    CONTINUE
C
         ZR(J,I) = 1.0E0
  160 CONTINUE
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
     1             + ABS(HR(L,L)) + ABS(HI(L,L))
         S2 = S1 + ABS(HR(L,L-1)) + ABS(HI(L,L-1))
         IF (S2 .EQ. S1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1) - HI(ENM1,EN) * HI(EN,ENM1)
      XI = HR(ENM1,EN) * HI(EN,ENM1) + HI(ENM1,EN) * HR(EN,ENM1)
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
      SI = ABS(HI(EN,ENM1)) + ABS(HI(ENM1,EN-2))
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
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS ..........
      XR = ABS(HR(ENM1,ENM1)) + ABS(HI(ENM1,ENM1))
      YR = ABS(HR(EN,ENM1)) + ABS(HI(EN,ENM1))
      ZZR = ABS(HR(EN,EN)) + ABS(HI(EN,EN))
C     .......... FOR M=EN-1 STEP -1 UNTIL L DO -- ..........
      DO 380 MM = L, ENM1
         M = ENM1 + L - MM
         IF (M .EQ. L) GO TO 420
         YI = YR
         YR = ABS(HR(M,M-1)) + ABS(HI(M,M-1))
         XI = ZZR
         ZZR = XR
         XR = ABS(HR(M-1,M-1)) + ABS(HI(M-1,M-1))
         S1 = ZZR / YI * (ZZR + XR + XI)
         S2 = S1 + YR
         IF (S2 .EQ. S1) GO TO 420
  380 CONTINUE
C     .......... TRIANGULAR DECOMPOSITION H=L*R ..........
  420 MP1 = M + 1
C
      DO 520 I = MP1, EN
         IM1 = I - 1
         XR = HR(IM1,IM1)
         XI = HI(IM1,IM1)
         YR = HR(I,IM1)
         YI = HI(I,IM1)
         IF (ABS(XR) + ABS(XI) .GE. ABS(YR) + ABS(YI)) GO TO 460
C     .......... INTERCHANGE ROWS OF HR AND HI ..........
         DO 440 J = IM1, N
            ZZR = HR(IM1,J)
            HR(IM1,J) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(IM1,J)
            HI(IM1,J) = HI(I,J)
            HI(I,J) = ZZI
  440    CONTINUE
C
         CALL CDIV(XR,XI,YR,YI,ZZR,ZZI)
         WR(I) = 1.0E0
         GO TO 480
  460    CALL CDIV(YR,YI,XR,XI,ZZR,ZZI)
         WR(I) = -1.0E0
  480    HR(I,IM1) = ZZR
         HI(I,IM1) = ZZI
C
         DO 500 J = I, N
            HR(I,J) = HR(I,J) - ZZR * HR(IM1,J) + ZZI * HI(IM1,J)
            HI(I,J) = HI(I,J) - ZZR * HI(IM1,J) - ZZI * HR(IM1,J)
  500    CONTINUE
C
  520 CONTINUE
C     .......... COMPOSITION R*L=H ..........
      DO 640 J = MP1, EN
         XR = HR(J,J-1)
         XI = HI(J,J-1)
         HR(J,J-1) = 0.0E0
         HI(J,J-1) = 0.0E0
C     .......... INTERCHANGE COLUMNS OF HR, HI, ZR, AND ZI,
C                IF NECESSARY ..........
         IF (WR(J) .LE. 0.0E0) GO TO 580
C
         DO 540 I = 1, J
            ZZR = HR(I,J-1)
            HR(I,J-1) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(I,J-1)
            HI(I,J-1) = HI(I,J)
            HI(I,J) = ZZI
  540    CONTINUE
C
         DO 560 I = LOW, IGH
            ZZR = ZR(I,J-1)
            ZR(I,J-1) = ZR(I,J)
            ZR(I,J) = ZZR
            ZZI = ZI(I,J-1)
            ZI(I,J-1) = ZI(I,J)
            ZI(I,J) = ZZI
  560    CONTINUE
C
  580    DO 600 I = 1, J
            HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)
            HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)
  600    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
         DO 620 I = LOW, IGH
            ZR(I,J-1) = ZR(I,J-1) + XR * ZR(I,J) - XI * ZI(I,J)
            ZI(I,J-1) = ZI(I,J-1) + XR * ZI(I,J) + XI * ZR(I,J)
  620    CONTINUE
C
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
      HR(1,1) = NORM
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
      DO 840 I = 1, ENM1
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
