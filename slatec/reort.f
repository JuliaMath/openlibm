*DECK REORT
      SUBROUTINE REORT (NCOMP, Y, YP, YHP, NIV, W, S, P, IP, STOWA,
     +   IFLAG)
C***BEGIN PROLOGUE  REORT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (REORT-S, DREORT-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C   INPUT
C *********
C     Y, YP and YHP = homogeneous solution matrix and particular
C                     solution vector to be orthonormalized.
C     IFLAG = 1 --  store YHP into Y and YP, test for
C                   reorthonormalization, orthonormalize if needed,
C                   save restart data.
C             2 --  store YHP into Y and YP, reorthonormalization,
C                   no restarts.
C                   (preset orthonormalization mode)
C             3 --  store YHP into Y and YP, reorthonormalization
C                   (when INHOMO=3 and X=XEND).
C **********************************************************************
C   OUTPUT
C *********
C     Y, YP = orthonormalized solutions.
C     NIV = number of independent vectors returned from DMGSBV.
C     IFLAG = 0 --  reorthonormalization was performed.
C            10 --  solution process must be restarted at the last
C                   orthonormalization point.
C            30 --  solutions are linearly dependent, problem must
C                   be restarted from the beginning.
C     W, P, IP = orthonormalization information.
C **********************************************************************
C
C***SEE ALSO  BVSUP
C***ROUTINES CALLED  MGSBV, SDOT, STOR1, STWAY
C***COMMON BLOCKS    ML15TO, ML18JR, ML8SZ
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  REORT
C
      DIMENSION Y(NCOMP,*),YP(*),W(*),S(*),P(*),IP(*),
     1          STOWA(*),YHP(NCOMP,*)
C
C **********************************************************************
C
      COMMON /ML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFC
      COMMON /ML15TO/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP,
     1                KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
      COMMON /ML18JR/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ,
     1                INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC,
     2                ICOCO
C
C **********************************************************************
C***FIRST EXECUTABLE STATEMENT  REORT
      NFCP=NFC+1
C
C     CHECK TO SEE IF ORTHONORMALIZATION TEST IS TO BE PERFORMED
C
      IF (IFLAG .NE. 1) GO TO 5
      KNSWOT=KNSWOT+1
      IF (KNSWOT .GE. NSWOT) GO TO 5
      IF ((XEND-X)*(X-XOT) .LT. 0.) RETURN
    5 CALL STOR1(Y,YHP,YP,YHP(1,NFCP),1,0,0)
C
C     ****************************************
C
C     ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y
C     AND PARTICULAR SOLUTION YP.
C
      NIV=NFC
      CALL MGSBV(NCOMP,NFC,Y,NCOMP,NIV,MFLAG,S,P,IP,INHOMO,YP,W,WCND)
C
C     ****************************************
C
C  CHECK FOR LINEAR DEPENDENCE OF THE SOLUTIONS.
C
      IF (MFLAG .EQ. 0)  GO TO 25
      IF (IFLAG .EQ. 2) GO TO 15
      IF (NSWOT .GT. 1  .OR.  LOTJP .EQ. 0) GO TO 20
   15 IFLAG=30
      RETURN
C
C     RETRIEVE DATA FOR A RESTART AT LAST ORTHONORMALIZATION POINT
C
   20 CALL STWAY(Y,YP,YHP,1,STOWA)
      LOTJP=1
      NSWOT=1
      KNSWOT=0
      MNSWOT=MNSWOT/2
      TND=TND+1.
      IFLAG=10
      RETURN
C
C     ****************************************
C
   25 IF (IFLAG .NE. 1) GO TO 60
C
C     TEST FOR ORTHONORMALIZATION
C
      IF (WCND .LT. 50.*TOL) GO TO 60
      DO 30 IJK=1,NFCP
      IF (S(IJK) .GT. 1.0E+20) GO TO 60
   30 CONTINUE
C
C     USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE NORM
C     DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION CHECKPOINT.
C     OTHER CONTROLS ON THE NUMBER OF STEPS TO THE NEXT CHECKPOINT
C     ARE ADDED FOR SAFETY PURPOSES.
C
      NSWOT=KNSWOT
      KNSWOT=0
      LOTJP=0
      WCND=LOG10(WCND)
      IF (WCND .GT. TND+3.) NSWOT=2*NSWOT
      IF (WCND .GE. PWCND) GO TO 40
      DX=X-PX
      DND=PWCND-WCND
      IF (DND .GE. 4) NSWOT=NSWOT/2
      DNDT=WCND-TND
      IF (ABS(DX*DNDT) .GT. DND*ABS(XEND-X)) GO TO 40
      XOT=X+DX*DNDT/DND
      GO TO 50
   40 XOT=XEND
   50 NSWOT=MIN(MNSWOT,NSWOT)
      PWCND=WCND
      PX=X
      RETURN
C
C     ****************************************
C
C     ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE HOMOGENEOUS
C     SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
C
   60 NSWOT=1
      KNSWOT=0
      LOTJP=1
      KK = 1
      L=1
      DO 70 K = 1,NFCC
      SRP=SQRT(P(KK))
      IF (INHOMO .EQ. 1) W(K)=SRP*W(K)
      VNORM=1./SRP
      P(KK)=VNORM
      KK = KK + NFCC + 1 - K
      IF (NFC .EQ. NFCC) GO TO 63
      IF (L .NE. K/2) GO TO 70
   63 DO 65 J = 1,NCOMP
   65 Y(J,L) = Y(J,L)*VNORM
      L=L+1
   70 CONTINUE
C
      IF (INHOMO .NE. 1  .OR.  NPS .EQ. 1)  GO TO 100
C
C     NORMALIZE THE PARTICULAR SOLUTION
C
      YPNM=SDOT(NCOMP,YP,1,YP,1)
      IF (YPNM .EQ. 0.0)  YPNM = 1.0
      YPNM = SQRT(YPNM)
      S(NFCP) = YPNM
      DO 80 J = 1,NCOMP
   80 YP(J) = YP(J) / YPNM
      DO 90 J = 1,NFCC
   90 W(J) = C * W(J)
C
  100 IF (IFLAG .EQ. 1) CALL STWAY(Y,YP,YHP,0,STOWA)
      IFLAG=0
      RETURN
      END
