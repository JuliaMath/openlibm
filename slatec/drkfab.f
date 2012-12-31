*DECK DRKFAB
      SUBROUTINE DRKFAB (NCOMP, XPTS, NXPTS, NFC, IFLAG, Z, MXNON, P,
     +   NTP, IP, YHP, NIV, U, V, W, S, STOWA, G, WORK, IWORK, NFCC)
C***BEGIN PROLOGUE  DRKFAB
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (RKFAB-S, DRKFAB-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C
C     Subroutine DRKFAB integrates the initial value equations using
C     the variable-step Runge-Kutta-Fehlberg integration scheme or
C     the variable-order Adams method and orthonormalization
C     determined by a linear dependence test.
C
C **********************************************************************
C
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  DBVDER, DDEABM, DDERKF, DREORT, DSTOR1
C***COMMON BLOCKS    DML15T, DML17B, DML18J, DML8SZ
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DRKFAB
C
      INTEGER ICOCO, IDID, IFLAG, IGOFX, INDPVT, INFO, INHOMO, INTEG,
     1     IPAR, ISTKOP, IVP, J, JFLAG, JON,
     2     K1, K10, K11, K2, K3, K4, K5, K6, K7, K8, K9, KKKINT,
     3     KKKZPW, KNSWOT, KOD, KOP, KOPP, L1, L2, LLLINT, LOTJP,
     4     MNSWOT, MXNON, MXNOND, NCOMP, NCOMPD, NDISK, NEEDIW, NEEDW,
     5     NEQ, NEQIVP, NFC, NFCC, NFCCD, NFCD, NFCP1, NIC, NIV, NON,
     6     NOPG, NPS, NSWOT, NTAPE, NTP, NTPD, NUMORT, NXPTS, NXPTSD,
     7     IP(NFCC,*), IWORK(*)
      DOUBLE PRECISION AE, C, G(*), P(NTP,*), PWCND, PX, RE,
     1     S(*), STOWA(*), TND, TOL, U(NCOMP,NFC,*),
     2     V(NCOMP,*), W(NFCC,*), WORK(*), X, XBEG, XEND, XOP,
     3     XOT, XPTS(*), XSAV, XXOP, YHP(NCOMP,*), Z(*)
C
C     ******************************************************************
C
      COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFCD
      COMMON /DML15T/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP,
     1                KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
      COMMON /DML18J/ AE,RE,TOL,NXPTSD,NIC,NOPG,MXNOND,NDISK,NTAPE,NEQ,
     1                INDPVT,INTEG,NPS,NTPD,NEQIVP,NUMORT,NFCCD,
     2                ICOCO
      COMMON /DML17B/ KKKZPW,NEEDW,NEEDIW,K1,K2,K3,K4,K5,K6,K7,K8,K9,
     1                K10,K11,L1,L2,KKKINT,LLLINT
C
      EXTERNAL DBVDER
C
C      *****************************************************************
C       INITIALIZATION OF COUNTERS AND VARIABLES.
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 220
C        BEGIN BLOCK PERMITTING ...EXITS TO 10
C***FIRST EXECUTABLE STATEMENT  DRKFAB
            KOD = 1
            NON = 1
            X = XBEG
            JON = 1
            INFO(1) = 0
            INFO(2) = 0
            INFO(3) = 1
            INFO(4) = 1
            WORK(1) = XEND
C        ...EXIT
            IF (NOPG .EQ. 0) GO TO 10
            INFO(3) = 0
            IF (X .EQ. Z(1)) JON = 2
   10    CONTINUE
         NFCP1 = NFC + 1
C
C        ***************************************************************
C        *****BEGINNING OF INTEGRATION LOOP AT OUTPUT
C        POINTS.******************
C        ***************************************************************
C
         DO 210 KOPP = 2, NXPTS
            KOP = KOPP
            XOP = XPTS(KOP)
            IF (NDISK .EQ. 0) KOD = KOP
C
   20       CONTINUE
C
C              STEP BY STEP INTEGRATION LOOP BETWEEN OUTPUT POINTS.
C
C              BEGIN BLOCK PERMITTING ...EXITS TO 190
C                 BEGIN BLOCK PERMITTING ...EXITS TO 30
                     XXOP = XOP
C                 ...EXIT
                     IF (NOPG .EQ. 0) GO TO 30
                     IF (XEND .GT. XBEG .AND. XOP .GT. Z(JON))
     1                  XXOP = Z(JON)
                     IF (XEND .LT. XBEG .AND. XOP .LT. Z(JON))
     1                  XXOP = Z(JON)
   30             CONTINUE
C
C                 ******************************************************
   40             CONTINUE
C                    BEGIN BLOCK PERMITTING ...EXITS TO 170
                        GO TO (50,60), INTEG
C                       DDERKF INTEGRATOR
C
   50                   CONTINUE
                           CALL DDERKF(DBVDER,NEQ,X,YHP,XXOP,INFO,RE,AE,
     1                                 IDID,WORK,KKKINT,IWORK,LLLINT,G,
     2                                 IPAR)
                        GO TO 70
C                       DDEABM INTEGRATOR
C
   60                   CONTINUE
                           CALL DDEABM(DBVDER,NEQ,X,YHP,XXOP,INFO,RE,AE,
     1                                 IDID,WORK,KKKINT,IWORK,LLLINT,G,
     2                                 IPAR)
   70                   CONTINUE
                        IF (IDID .GE. 1) GO TO 80
                           INFO(1) = 1
C                    ......EXIT
                           IF (IDID .EQ. -1) GO TO 170
                           IFLAG = 20 - IDID
C     .....................EXIT
                           GO TO 220
   80                   CONTINUE
C
C                       ************************************************
C                           GRAM-SCHMIDT ORTHOGONALIZATION TEST FOR
C                           ORTHONORMALIZATION (TEMPORARILY USING U AND
C                           V IN THE TEST)
C
                        IF (NOPG .EQ. 0) GO TO 100
                           IF (XXOP .EQ. Z(JON)) GO TO 90
C
C                             ******************************************
C                                 CONTINUE INTEGRATION IF WE ARE NOT AT
C                                 AN OUTPUT POINT.
C
C           ..................EXIT
                              IF (IDID .NE. 1) GO TO 200
C                    .........EXIT
                              GO TO 170
   90                      CONTINUE
                           JFLAG = 2
                        GO TO 110
  100                   CONTINUE
                           JFLAG = 1
                           IF (INHOMO .EQ. 3 .AND. X .EQ. XEND)
     1                        JFLAG = 3
  110                   CONTINUE
C
                        IF (NDISK .EQ. 0) NON = NUMORT + 1
                        CALL DREORT(NCOMP,U(1,1,KOD),V(1,KOD),YHP,NIV,
     1                              W(1,NON),S,P(1,NON),IP(1,NON),STOWA,
     2                              JFLAG)
C
                        IF (JFLAG .NE. 30) GO TO 120
                           IFLAG = 30
C     .....................EXIT
                           GO TO 220
  120                   CONTINUE
C
                        IF (JFLAG .NE. 10) GO TO 130
                           XOP = XPTS(KOP)
                           IF (NDISK .EQ. 0) KOD = KOP
C              ............EXIT
                           GO TO 190
  130                   CONTINUE
C
                        IF (JFLAG .EQ. 0) GO TO 140
C
C                          *********************************************
C                              CONTINUE INTEGRATION IF WE ARE NOT AT AN
C                              OUTPUT POINT.
C
C           ...............EXIT
                           IF (IDID .NE. 1) GO TO 200
C                    ......EXIT
                           GO TO 170
  140                   CONTINUE
C
C                       ************************************************
C                           STORE ORTHONORMALIZED VECTORS INTO SOLUTION
C                           VECTORS.
C
                        IF (NUMORT .LT. MXNON) GO TO 150
                        IF (X .EQ. XEND) GO TO 150
                           IFLAG = 13
C     .....................EXIT
                           GO TO 220
  150                   CONTINUE
C
                        NUMORT = NUMORT + 1
                        CALL DSTOR1(YHP,U(1,1,KOD),YHP(1,NFCP1),
     1                              V(1,KOD),1,NDISK,NTAPE)
C
C                       ************************************************
C                           STORE ORTHONORMALIZATION INFORMATION,
C                           INITIALIZE INTEGRATION FLAG, AND CONTINUE
C                           INTEGRATION TO THE NEXT ORTHONORMALIZATION
C                           POINT OR OUTPUT POINT.
C
                        Z(NUMORT) = X
                        IF (INHOMO .EQ. 1 .AND. NPS .EQ. 0)
     1                     C = S(NFCP1)*C
                        IF (NDISK .EQ. 0) GO TO 160
                           IF (INHOMO .EQ. 1)
     1                        WRITE (NTAPE) (W(J,1), J = 1, NFCC)
                           WRITE (NTAPE)
     1                           (IP(J,1), J = 1, NFCC),
     2                           (P(J,1), J = 1, NTP)
  160                   CONTINUE
                        INFO(1) = 0
                        JON = JON + 1
C                 ......EXIT
                        IF (NOPG .EQ. 1 .AND. X .NE. XOP) GO TO 180
C
C                       ************************************************
C                           CONTINUE INTEGRATION IF WE ARE NOT AT AN
C                           OUTPUT POINT.
C
C           ............EXIT
                        IF (IDID .NE. 1) GO TO 200
  170                CONTINUE
                  GO TO 40
  180             CONTINUE
  190          CONTINUE
            GO TO 20
  200       CONTINUE
C
C           STORAGE OF HOMOGENEOUS SOLUTIONS IN U AND THE PARTICULAR
C           SOLUTION IN V AT THE OUTPUT POINTS.
C
            CALL DSTOR1(U(1,1,KOD),YHP,V(1,KOD),YHP(1,NFCP1),0,NDISK,
     1                  NTAPE)
  210    CONTINUE
C        ***************************************************************
C        ***************************************************************
C
         IFLAG = 0
  220 CONTINUE
      RETURN
      END
