*DECK COMLR
      SUBROUTINE COMLR (NM, N, LOW, IGH, HR, HI, WR, WI, IERR)
C***BEGIN PROLOGUE  COMLR
C***PURPOSE  Compute the eigenvalues of a complex upper Hessenberg
C            matrix using the modified LR method.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C2B
C***TYPE      COMPLEX (COMLR-C)
C***KEYWORDS  EIGENVALUES, EISPACK, LR METHOD
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure COMLR,
C     NUM. MATH. 12, 369-376(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
C
C     This subroutine finds the eigenvalues of a COMPLEX
C     UPPER Hessenberg matrix by the modified LR method.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, HR and HI, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix H=(HR,HI).  N is an INTEGER
C          variable.  N must be less than or equal to NM.
C
C        LOW and IGH are two INTEGER variables determined by the
C          balancing subroutine  CBAL.  If  CBAL  has not been used,
C          set LOW=1 and IGH equal to the order of the matrix, N.
C
C        HR and HI contain the real and imaginary parts, respectively,
C          of the complex upper Hessenberg matrix.  Their lower
C          triangles below the subdiagonal contain the multipliers
C          which were used in the reduction by  COMHES, if performed.
C          HR and HI are two-dimensional REAL arrays, dimensioned
C          HR(NM,N) and HI(NM,N).
C
C     On OUTPUT
C
C        The upper Hessenberg portions of HR and HI have been
C          destroyed.  Therefore, they must be saved before calling
C          COMLR  if subsequent calculation of eigenvectors is to
C          be performed.
C
C        WR and WI contain the real and imaginary parts, respectively,
C          of the eigenvalues of the upper Hessenberg matrix.  If an
C          error exit is made, the eigenvalues should be correct for
C          indices IERR+1, IERR+2, ..., N.  WR and WI are one-
C          dimensional REAL arrays, dimensioned WR(N) and WI(N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after a total of 30*N iterations.
C                     The eigenvalues should be correct for indices
C                     IERR+1, IERR+2, ..., N.
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
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  COMLR
C
      INTEGER I,J,L,M,N,EN,LL,MM,NM,IGH,IM1,ITN,ITS,LOW,MP1,ENM1,IERR
      REAL HR(NM,*),HI(NM,*),WR(*),WI(*)
      REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,S1,S2
C
C***FIRST EXECUTABLE STATEMENT  COMLR
      IERR = 0
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
      DO 200 I = 1, N
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
  220 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW E0 -- ..........
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
         DO 440 J = IM1, EN
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
         DO 500 J = I, EN
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
C     .......... INTERCHANGE COLUMNS OF HR AND HI,
C                IF NECESSARY ..........
         IF (WR(J) .LE. 0.0E0) GO TO 580
C
         DO 540 I = L, J
            ZZR = HR(I,J-1)
            HR(I,J-1) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(I,J-1)
            HI(I,J-1) = HI(I,J)
            HI(I,J) = ZZI
  540    CONTINUE
C
  580    DO 600 I = L, J
            HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)
            HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)
  600    CONTINUE
C
  640 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 WR(EN) = HR(EN,EN) + TR
      WI(EN) = HI(EN,EN) + TI
      EN = ENM1
      GO TO 220
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
