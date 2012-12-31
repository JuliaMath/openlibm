*DECK SCOEF
      SUBROUTINE SCOEF (YH, YP, NCOMP, NROWB, NFC, NIC, B, BETA, COEF,
     +   INHOMO, RE, AE, BY, CVEC, WORK, IWORK, IFLAG, NFCC)
C***BEGIN PROLOGUE  SCOEF
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SCOEF-S, DCOEF-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C INPUT TO SCOEF
C **********************************************************************
C
C     YH = Matrix of homogeneous solutions.
C     YP = Vector containing particular solution.
C     NCOMP = Number of components per solution vector.
C     NROWB = First dimension of B in calling program.
C     NFC = Number of base solution vectors.
C     NFCC = 2*NFC for the special treatment of complex valued
C            equations. Otherwise, NFCC=NFC.
C     NIC = Number of specified initial conditions.
C     B = Boundary condition matrix at X = Xfinal.
C     BETA = Vector of nonhomogeneous boundary conditions at X = Xfinal.
C              1 - Nonzero particular solution
C     INHOMO = 2 - Zero particular solution
C              3 - Eigenvalue problem
C     RE = Relative error tolerance
C     AE = Absolute error tolerance
C     BY = Storage space for the matrix  B*YH
C     CVEC = Storage space for the vector  BETA-B*YP
C     WORK = Real array of internal storage. Dimension must be .GE.
C            NFCC*(NFCC+4)
C     IWORK = Integer array of internal storage. Dimension must be .GE.
C             3+NFCC
C
C **********************************************************************
C OUTPUT FROM SCOEF
C **********************************************************************
C
C     COEF = Array containing superposition constants.
C     IFLAG = Indicator of success from SUDS in solving the
C             boundary equations
C           = 0 Boundary equations are solved
C           = 1 Boundary equations appear to have many solutions
C           = 2 Boundary equations appear to be inconsistent
C           = 3 For this value of an eigenparameter, the boundary
C               equations have only the zero solution.
C
C **********************************************************************
C
C     Subroutine SCOEF solves for the superposition constants from the
C     linear equations defined by the boundary conditions at X = Xfinal.
C
C                          B*YP + B*YH*COEF = BETA
C
C **********************************************************************
C
C***SEE ALSO  BVSUP
C***ROUTINES CALLED  SDOT, SUDS, XGETF, XSETF
C***COMMON BLOCKS    ML5MCO
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  SCOEF
C
      DIMENSION YH(NCOMP,*),YP(*),B(NROWB,*),BETA(*),
     1          COEF(*),BY(NFCC,*),CVEC(*),WORK(*),IWORK(*)
C
      COMMON /ML5MCO/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
C
C     SET UP MATRIX  B*YH  AND VECTOR  BETA - B*YP
C
C***FIRST EXECUTABLE STATEMENT  SCOEF
      NCOMP2=NCOMP/2
      DO 7 K = 1,NFCC
      DO 1 J = 1,NFC
      L=J
      IF (NFC .NE. NFCC) L=2*J-1
    1 BY(K,L) = SDOT(NCOMP,B(K,1),NROWB,YH(1,J),1)
      IF (NFC .EQ. NFCC) GO TO 3
      DO 2 J=1,NFC
      L=2*J
      BYKL=SDOT(NCOMP2,B(K,1),NROWB,YH(NCOMP2+1,J),1)
      BY(K,L)=SDOT(NCOMP2,B(K,NCOMP2+1),NROWB,YH(1,J),1) - BYKL
    2 CONTINUE
    3 GO TO (4,5,6), INHOMO
C     CASE 1
    4 CVEC(K) = BETA(K) - SDOT(NCOMP,B(K,1),NROWB,YP,1)
      GO TO 7
C     CASE 2
    5 CVEC(K) = BETA(K)
      GO TO 7
C     CASE 3
    6 CVEC(K) = 0.
    7 CONTINUE
      CONS=ABS(CVEC(1))
      BYS=ABS(BY(1,1))
C
C **********************************************************************
C     SOLVE LINEAR SYSTEM
C
      IFLAG=0
      MLSO=0
      IF (INHOMO .EQ. 3) MLSO=1
      KFLAG = 0.5 * LOG10(EPS)
      CALL XGETF(NF)
      CALL XSETF(0)
   10 CALL SUDS(BY,COEF,CVEC,NFCC,NFCC,NFCC,KFLAG,MLSO,WORK,IWORK)
      IF (KFLAG .NE. 3) GO TO 13
      KFLAG=1
      IFLAG=1
      GO TO 10
   13 IF (KFLAG .EQ. 4) IFLAG=2
      CALL XSETF(NF)
      IF (NFCC .EQ. 1) GO TO 25
      IF (INHOMO .NE. 3) RETURN
      IF (IWORK(1) .LT. NFCC) GO TO 17
      IFLAG=3
      DO 14 K=1,NFCC
   14 COEF(K)=0.
      COEF(NFCC)=1.
      NFCCM1=NFCC-1
      DO 15 K=1,NFCCM1
      J=NFCC-K
      L=NFCC-J+1
      GAM=SDOT(L,BY(J,J),NFCC,COEF(J),1)/(WORK(J)*BY(J,J))
      DO 15 I=J,NFCC
   15 COEF(I)=COEF(I)+GAM*BY(J,I)
      RETURN
   17 DO 20 K=1,NFCC
      KI=4*NFCC+K
   20 COEF(K)=WORK(KI)
      RETURN
C
C **********************************************************************
C     TESTING FOR EXISTENCE AND UNIQUENESS OF BOUNDARY-VALUE PROBLEM
C     SOLUTION IN A SCALAR CASE
C
   25 BN = 0.
      UN = 0.
      YPN=0.
      DO 30 K = 1,NCOMP
      UN = MAX(UN,ABS(YH(K,1)))
      YPN=MAX(YPN,ABS(YP(K)))
   30 BN = MAX(BN,ABS(B(1,K)))
      BBN = MAX(BN,ABS(BETA(1)))
      IF (BYS .GT. 10.*(RE*UN + AE)*BN)  GO TO 35
      BRN = BBN / BN * BYS
      IF (CONS .GE. 0.1*BRN  .AND.  CONS .LE. 10.*BRN) IFLAG=1
      IF (CONS .GT. 10.*BRN) IFLAG=2
      IF (CONS  .LE.  RE*ABS(BETA(1))+AE + (RE*YPN+AE)*BN) IFLAG=1
      IF (INHOMO .EQ. 3) COEF(1)=1.
      RETURN
   35 IF (INHOMO .NE. 3) RETURN
      IFLAG=3
      COEF(1)=1.
      RETURN
      END
