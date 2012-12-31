*DECK DREORT
      SUBROUTINE DREORT (NCOMP, Y, YP, YHP, NIV, W, S, P, IP, STOWA,
     +   IFLAG)
C***BEGIN PROLOGUE  DREORT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (REORT-S, DREORT-D)
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
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  DDOT, DMGSBV, DSTOR1, DSTWAY
C***COMMON BLOCKS    DML15T, DML18J, DML8SZ
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DREORT
C
      DOUBLE PRECISION DDOT
      INTEGER ICOCO, IFLAG, IGOFX, IJK, INDPVT, INFO, INHOMO, INTEG,
     1     IP(*), ISTKOP, IVP, J, K, KK, KNSWOT, KOP, L, LOTJP, MFLAG,
     2     MNSWOT, MXNON, NCOMP, NCOMPD, NDISK, NEQ, NEQIVP, NFC,
     3     NFCC, NFCP, NIC, NIV, NOPG, NPS, NSWOT, NTAPE, NTP, NUMORT,
     4     NXPTS
      DOUBLE PRECISION AE, C, DND, DNDT, DX, P(*), PWCND, PX, RE, S(*),
     1     SRP, STOWA(*), TND, TOL, VNORM, W(*), WCND, X, XBEG, XEND,
     2     XOP, XOT, XSAV, Y(NCOMP,*), YHP(NCOMP,*), YP(*), YPNM
C
C     ******************************************************************
C
      COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFC
      COMMON /DML15T/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP,
     1                KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
      COMMON /DML18J/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ,
     1                INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC,
     2                ICOCO
C
C **********************************************************************
C     BEGIN BLOCK PERMITTING ...EXITS TO 210
C        BEGIN BLOCK PERMITTING ...EXITS TO 10
C***FIRST EXECUTABLE STATEMENT  DREORT
            NFCP = NFC + 1
C
C           CHECK TO SEE IF ORTHONORMALIZATION TEST IS TO BE PERFORMED
C
C        ...EXIT
            IF (IFLAG .NE. 1) GO TO 10
            KNSWOT = KNSWOT + 1
C        ...EXIT
            IF (KNSWOT .GE. NSWOT) GO TO 10
C     ......EXIT
            IF ((XEND - X)*(X - XOT) .LT. 0.0D0) GO TO 210
   10    CONTINUE
         CALL DSTOR1(Y,YHP,YP,YHP(1,NFCP),1,0,0)
C
C        ***************************************************************
C
C        ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y
C        AND PARTICULAR SOLUTION YP.
C
         NIV = NFC
         CALL DMGSBV(NCOMP,NFC,Y,NCOMP,NIV,MFLAG,S,P,IP,INHOMO,YP,W,
     1               WCND)
C
C           ************************************************************
C
C        CHECK FOR LINEAR DEPENDENCE OF THE SOLUTIONS.
C
         IF (MFLAG .EQ. 0) GO TO 50
C           BEGIN BLOCK PERMITTING ...EXITS TO 40
               IF (IFLAG .EQ. 2) GO TO 30
                  IF (NSWOT .LE. 1 .AND. LOTJP .NE. 0) GO TO 20
C
C                    RETRIEVE DATA FOR A RESTART AT LAST
C                    ORTHONORMALIZATION POINT
C
                     CALL DSTWAY(Y,YP,YHP,1,STOWA)
                     LOTJP = 1
                     NSWOT = 1
                     KNSWOT = 0
                     MNSWOT = MNSWOT/2
                     TND = TND + 1.0D0
                     IFLAG = 10
C           .........EXIT
                     GO TO 40
   20             CONTINUE
   30          CONTINUE
               IFLAG = 30
   40       CONTINUE
         GO TO 200
   50    CONTINUE
C           BEGIN BLOCK PERMITTING ...EXITS TO 190
C              BEGIN BLOCK PERMITTING ...EXITS TO 110
C
C                 ******************************************************
C
C              ...EXIT
                  IF (IFLAG .NE. 1) GO TO 110
C
C                 TEST FOR ORTHONORMALIZATION
C
C              ...EXIT
                  IF (WCND .LT. 50.0D0*TOL) GO TO 110
                  DO 60 IJK = 1, NFCP
C              ......EXIT
                     IF (S(IJK) .GT. 1.0D20) GO TO 110
   60             CONTINUE
C
C                 USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE
C                 NORM DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION
C                 CHECKPOINT.  OTHER CONTROLS ON THE NUMBER OF STEPS TO
C                 THE NEXT CHECKPOINT ARE ADDED FOR SAFETY PURPOSES.
C
                  NSWOT = KNSWOT
                  KNSWOT = 0
                  LOTJP = 0
                  WCND = LOG10(WCND)
                  IF (WCND .GT. TND + 3.0D0) NSWOT = 2*NSWOT
                  IF (WCND .LT. PWCND) GO TO 70
                     XOT = XEND
                     NSWOT = MIN(MNSWOT,NSWOT)
                     PWCND = WCND
                     PX = X
                  GO TO 100
   70             CONTINUE
                     DX = X - PX
                     DND = PWCND - WCND
                     IF (DND .GE. 4) NSWOT = NSWOT/2
                     DNDT = WCND - TND
                     IF (ABS(DX*DNDT) .LE. DND*ABS(XEND-X)) GO TO 80
                        XOT = XEND
                        NSWOT = MIN(MNSWOT,NSWOT)
                        PWCND = WCND
                        PX = X
                     GO TO 90
   80                CONTINUE
                        XOT = X + DX*DNDT/DND
                        NSWOT = MIN(MNSWOT,NSWOT)
                        PWCND = WCND
                        PX = X
   90                CONTINUE
  100             CONTINUE
C           ......EXIT
                  GO TO 190
  110          CONTINUE
C
C              *********************************************************
C
C              ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE
C              HOMOGENEOUS SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
C
               NSWOT = 1
               KNSWOT = 0
               LOTJP = 1
               KK = 1
               L = 1
               DO 150 K = 1, NFCC
C                 BEGIN BLOCK PERMITTING ...EXITS TO 140
                     SRP = SQRT(P(KK))
                     IF (INHOMO .EQ. 1) W(K) = SRP*W(K)
                     VNORM = 1.0D0/SRP
                     P(KK) = VNORM
                     KK = KK + NFCC + 1 - K
                     IF (NFC .EQ. NFCC) GO TO 120
C                 ......EXIT
                        IF (L .NE. K/2) GO TO 140
  120                CONTINUE
                     DO 130 J = 1, NCOMP
                        Y(J,L) = Y(J,L)*VNORM
  130                CONTINUE
                     L = L + 1
  140             CONTINUE
  150          CONTINUE
C
               IF (INHOMO .NE. 1 .OR. NPS .EQ. 1) GO TO 180
C
C                 NORMALIZE THE PARTICULAR SOLUTION
C
                  YPNM = DDOT(NCOMP,YP,1,YP,1)
                  IF (YPNM .EQ. 0.0D0) YPNM = 1.0D0
                  YPNM = SQRT(YPNM)
                  S(NFCP) = YPNM
                  DO 160 J = 1, NCOMP
                     YP(J) = YP(J)/YPNM
  160             CONTINUE
                  DO 170 J = 1, NFCC
                     W(J) = C*W(J)
  170             CONTINUE
  180          CONTINUE
C
               IF (IFLAG .EQ. 1) CALL DSTWAY(Y,YP,YHP,0,STOWA)
               IFLAG = 0
  190       CONTINUE
  200    CONTINUE
  210 CONTINUE
      RETURN
      END
