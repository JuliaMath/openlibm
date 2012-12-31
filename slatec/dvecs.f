*DECK DVECS
      SUBROUTINE DVECS (NCOMP, LNFC, YHP, WORK, IWORK, INHOMO, IFLAG)
C***BEGIN PROLOGUE  DVECS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (SVECS-S, DVECS-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C  This subroutine is used for the special structure of COMPLEX*16
C  valued problems. DMGSBV is called upon to obtain LNFC vectors from an
C  original set of 2*LNFC independent vectors so that the resulting
C  LNFC vectors together with their imaginary product or mate vectors
C  form an independent set.
C
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  DMGSBV
C***COMMON BLOCKS    DML18J
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891009  Removed unreferenced statement label.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DVECS
C
      INTEGER ICOCO, IDP, IFLAG, INDPVT, INHOMO, INTEG, IWORK(*), K,
     1     KP, LNFC, LNFCC, MXNON, NCOMP, NDISK, NEQ, NEQIVP, NIC, NIV,
     2     NOPG, NPS, NTAPE, NTP, NUMORT, NXPTS
      DOUBLE PRECISION AE, DUM, RE, TOL, WORK(*), YHP(NCOMP,*)
      COMMON /DML18J/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ,
     1                INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,LNFCC,
     2                ICOCO
C***FIRST EXECUTABLE STATEMENT  DVECS
         IF (LNFC .NE. 1) GO TO 20
            DO 10 K = 1, NCOMP
               YHP(K,LNFC+1) = YHP(K,LNFCC+1)
   10       CONTINUE
            IFLAG = 1
         GO TO 60
   20    CONTINUE
            NIV = LNFC
            LNFC = 2*LNFC
            LNFCC = 2*LNFCC
            KP = LNFC + 2 + LNFCC
            IDP = INDPVT
            INDPVT = 0
            CALL DMGSBV(NCOMP,LNFC,YHP,NCOMP,NIV,IFLAG,WORK(1),WORK(KP),
     1                  IWORK(1),INHOMO,YHP(1,LNFC+1),WORK(LNFC+2),DUM)
            LNFC = LNFC/2
            LNFCC = LNFCC/2
            INDPVT = IDP
            IF (IFLAG .NE. 0 .OR. NIV .NE. LNFC) GO TO 40
               DO 30 K = 1, NCOMP
                  YHP(K,LNFC+1) = YHP(K,LNFCC+1)
   30          CONTINUE
               IFLAG = 1
            GO TO 50
   40       CONTINUE
               IFLAG = 99
   50       CONTINUE
   60    CONTINUE
      CONTINUE
      RETURN
      END
