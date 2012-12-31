*DECK DCOEF
      SUBROUTINE DCOEF (YH, YP, NCOMP, NROWB, NFC, NIC, B, BETA, COEF,
     +   INHOMO, RE, AE, BY, CVEC, WORK, IWORK, IFLAG, NFCC)
C***BEGIN PROLOGUE  DCOEF
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (SCOEF-S, DCOEF-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C INPUT to DCOEF
C **********************************************************************
C
C     YH = matrix of homogeneous solutions.
C     YP = vector containing particular solution.
C     NCOMP = number of components per solution vector.
C     NROWB = first dimension of B in calling program.
C     NFC = number of base solution vectors.
C     NFCC = 2*NFC for the special treatment of COMPLEX*16 valued
C            equations. Otherwise, NFCC=NFC.
C     NIC = number of specified initial conditions.
C     B = boundary condition matrix at X = XFINAL.
C     BETA = vector of nonhomogeneous boundary conditions at X = XFINAL.
C              1 - nonzero particular solution
C     INHOMO = 2 - zero particular solution
C              3 - eigenvalue problem
C     RE = relative error tolerance.
C     AE = absolute error tolerance.
C     BY = storage space for the matrix  B*YH
C     CVEC = storage space for the vector  BETA-B*YP
C     WORK = double precision array of internal storage. Dimension must
C     be GE
C            NFCC*(NFCC+4)
C     IWORK = integer array of internal storage. Dimension must be GE
C             3+NFCC
C
C **********************************************************************
C OUTPUT from DCOEF
C **********************************************************************
C
C     COEF = array containing superposition constants.
C     IFLAG = indicator of success from DSUDS in solving the
C             boundary equations.
C           = 0 boundary equations are solved.
C           = 1 boundary equations appear to have many solutions.
C           = 2 boundary equations appear to be inconsistent.
C           = 3 for this value of an eigenparameter, the boundary
C               equations have only the zero solution.
C
C **********************************************************************
C
C     Subroutine DCOEF solves for the superposition constants from the
C     linear equations defined by the boundary conditions at X = XFINAL.
C
C                          B*YP + B*YH*COEF = BETA
C
C **********************************************************************
C
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  DDOT, DSUDS, XGETF, XSETF
C***COMMON BLOCKS    DML5MC
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   890921  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DCOEF
C
      DOUBLE PRECISION DDOT
      INTEGER I, IFLAG, INHOMO, IWORK(*), J, K, KFLAG, KI, L, LPAR,
     1     MLSO, NCOMP, NCOMP2, NF, NFC, NFCC, NFCCM1, NIC,
     2     NROWB
      DOUBLE PRECISION AE, B(NROWB,*), BBN, BETA(*), BN, BRN,
     1     BY(NFCC,*), BYKL, BYS, COEF(*), CONS, CVEC(*), EPS,
     2     FOURU, GAM, RE, SQOVFL, SRU, TWOU, UN, URO, WORK(*),
     3     YH(NCOMP,*), YP(*), YPN
C
      COMMON /DML5MC/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
C***FIRST EXECUTABLE STATEMENT  DCOEF
C
C     SET UP MATRIX  B*YH  AND VECTOR  BETA - B*YP
C
      NCOMP2 = NCOMP/2
      DO 80 K = 1, NFCC
         DO 10 J = 1, NFC
            L = J
            IF (NFC .NE. NFCC) L = 2*J - 1
            BY(K,L) = DDOT(NCOMP,B(K,1),NROWB,YH(1,J),1)
   10    CONTINUE
         IF (NFC .EQ. NFCC) GO TO 30
            DO 20 J = 1, NFC
               L = 2*J
               BYKL = DDOT(NCOMP2,B(K,1),NROWB,YH(NCOMP2+1,J),1)
               BY(K,L) = DDOT(NCOMP2,B(K,NCOMP2+1),NROWB,YH(1,J),1)
     1                   - BYKL
   20       CONTINUE
   30    CONTINUE
         GO TO (40,50,60), INHOMO
C        CASE 1
   40    CONTINUE
            CVEC(K) = BETA(K) - DDOT(NCOMP,B(K,1),NROWB,YP,1)
         GO TO 70
C        CASE 2
   50    CONTINUE
            CVEC(K) = BETA(K)
         GO TO 70
C        CASE 3
   60    CONTINUE
            CVEC(K) = 0.0D0
   70    CONTINUE
   80 CONTINUE
      CONS = ABS(CVEC(1))
      BYS = ABS(BY(1,1))
C
C     ******************************************************************
C         SOLVE LINEAR SYSTEM
C
      IFLAG = 0
      MLSO = 0
      IF (INHOMO .EQ. 3) MLSO = 1
      KFLAG = 0.5D0 * LOG10(EPS)
      CALL XGETF(NF)
      CALL XSETF(0)
   90 CONTINUE
         CALL DSUDS(BY,COEF,CVEC,NFCC,NFCC,NFCC,KFLAG,MLSO,WORK,IWORK)
         IF (KFLAG .NE. 3) GO TO 100
         KFLAG = 1
         IFLAG = 1
      GO TO 90
  100 CONTINUE
      IF (KFLAG .EQ. 4) IFLAG = 2
      CALL XSETF(NF)
      IF (NFCC .EQ. 1) GO TO 180
         IF (INHOMO .NE. 3) GO TO 170
            IF (IWORK(1) .LT. NFCC) GO TO 140
               IFLAG = 3
               DO 110 K = 1, NFCC
                  COEF(K) = 0.0D0
  110          CONTINUE
               COEF(NFCC) = 1.0D0
               NFCCM1 = NFCC - 1
               DO 130 K = 1, NFCCM1
                  J = NFCC - K
                  L = NFCC - J + 1
                  GAM = DDOT(L,BY(J,J),NFCC,COEF(J),1)/(WORK(J)*BY(J,J))
                  DO 120 I = J, NFCC
                     COEF(I) = COEF(I) + GAM*BY(J,I)
  120             CONTINUE
  130          CONTINUE
            GO TO 160
  140       CONTINUE
               DO 150 K = 1, NFCC
                  KI = 4*NFCC + K
                  COEF(K) = WORK(KI)
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
      GO TO 220
  180 CONTINUE
C
C        ***************************************************************
C            TESTING FOR EXISTENCE AND UNIQUENESS OF BOUNDARY-VALUE
C            PROBLEM SOLUTION IN A SCALAR CASE
C
         BN = 0.0D0
         UN = 0.0D0
         YPN = 0.0D0
         DO 190 K = 1, NCOMP
            UN = MAX(UN,ABS(YH(K,1)))
            YPN = MAX(YPN,ABS(YP(K)))
            BN = MAX(BN,ABS(B(1,K)))
  190    CONTINUE
         BBN = MAX(BN,ABS(BETA(1)))
         IF (BYS .GT. 10.0D0*(RE*UN + AE)*BN) GO TO 200
            BRN = BBN/BN*BYS
            IF (CONS .GE. 0.1D0*BRN .AND. CONS .LE. 10.0D0*BRN)
     1         IFLAG = 1
            IF (CONS .GT. 10.0D0*BRN) IFLAG = 2
            IF (CONS .LE. RE*ABS(BETA(1)) + AE + (RE*YPN + AE)*BN)
     1         IFLAG = 1
            IF (INHOMO .EQ. 3) COEF(1) = 1.0D0
         GO TO 210
  200    CONTINUE
         IF (INHOMO .NE. 3) GO TO 210
            IFLAG = 3
            COEF(1) = 1.0D0
  210    CONTINUE
  220 CONTINUE
      RETURN
      END
