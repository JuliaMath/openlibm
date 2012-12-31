*DECK DSTWAY
      SUBROUTINE DSTWAY (U, V, YHP, INOUT, STOWA)
C***BEGIN PROLOGUE  DSTWAY
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (STWAY-S, DSTWAY-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C  This subroutine stores (recalls) integration data in the event
C  that a restart is needed (the homogeneous solution vectors become
C  too dependent to continue).
C
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  DSTOR1
C***COMMON BLOCKS    DML15T, DML18J, DML8SZ
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DSTWAY
C
      INTEGER ICOCO, IGOFX, INDPVT, INFO, INHOMO, INOUT, INTEG, ISTKOP,
     1     IVP, J, K, KNSWOT, KO, KOP, KS, KSJ, LOTJP, MNSWOT, MXNON,
     2     NCOMP, NDISK, NEQ, NEQIVP, NFC, NFCC, NIC, NOPG, NPS, NSWOT,
     3     NTAPE, NTP, NUMORT, NXPTS
      DOUBLE PRECISION AE, C, PWCND, PX, RE, STOWA(*), TND, TOL, U(*),
     1     V(*), X, XBEG, XEND, XOP, XOT, XSAV, YHP(*)
C
      COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
      COMMON /DML15T/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP,
     1                KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
      COMMON /DML18J/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ,
     1                INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC,
     2                ICOCO
C
C***FIRST EXECUTABLE STATEMENT  DSTWAY
      IF (INOUT .EQ. 1) GO TO 30
C
C        SAVE IN STOWA ARRAY AND ISTKOP
C
         KS = NFC*NCOMP
         CALL DSTOR1(STOWA,U,STOWA(KS+1),V,1,0,0)
         KS = KS + NCOMP
         IF (NEQIVP .LT. 1) GO TO 20
         DO 10 J = 1, NEQIVP
            KSJ = KS + J
            STOWA(KSJ) = YHP(KSJ)
   10    CONTINUE
   20    CONTINUE
         KS = KS + NEQIVP
         STOWA(KS+1) = X
         ISTKOP = KOP
         IF (XOP .EQ. X) ISTKOP = KOP + 1
      GO TO 80
   30 CONTINUE
C
C        RECALL FROM STOWA ARRAY AND ISTKOP
C
         KS = NFC*NCOMP
         CALL DSTOR1(YHP,STOWA,YHP(KS+1),STOWA(KS+1),1,0,0)
         KS = KS + NCOMP
         IF (NEQIVP .LT. 1) GO TO 50
         DO 40 J = 1, NEQIVP
            KSJ = KS + J
            YHP(KSJ) = STOWA(KSJ)
   40    CONTINUE
   50    CONTINUE
         KS = KS + NEQIVP
         X = STOWA(KS+1)
         INFO(1) = 0
         KO = KOP - ISTKOP
         KOP = ISTKOP
         IF (NDISK .EQ. 0 .OR. KO .EQ. 0) GO TO 70
            DO 60 K = 1, KO
               BACKSPACE NTAPE
   60       CONTINUE
   70    CONTINUE
   80 CONTINUE
      RETURN
      END
