*DECK SVECS
      SUBROUTINE SVECS (NCOMP, LNFC, YHP, WORK, IWORK, INHOMO, IFLAG)
C***BEGIN PROLOGUE  SVECS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SVECS-S, DVECS-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C  This subroutine is used for the special structure of complex valued
C  problems. MGSBV is called upon to obtain LNFC vectors from an
C  original set of 2*LNFC independent vectors so that the resulting
C  LNFC vectors together with their imaginary product or mate vectors
C  form an independent set.
C
C***SEE ALSO  BVSUP
C***ROUTINES CALLED  MGSBV
C***COMMON BLOCKS    ML18JR
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  SVECS
C
      DIMENSION YHP(NCOMP,*),WORK(*),IWORK(*)
      COMMON /ML18JR/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ,
     1                INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,LNFCC,
     2                ICOCO
C***FIRST EXECUTABLE STATEMENT  SVECS
      IF (LNFC .EQ. 1) GO TO 5
      NIV=LNFC
      LNFC=2*LNFC
      LNFCC=2*LNFCC
      KP=LNFC+2+LNFCC
      IDP=INDPVT
      INDPVT=0
      CALL MGSBV(NCOMP,LNFC,YHP,NCOMP,NIV,IFLAG,WORK(1),WORK(KP),
     1         IWORK(1),INHOMO,YHP(1,LNFC+1),WORK(LNFC+2),DUM)
      LNFC=LNFC/2
      LNFCC=LNFCC/2
      INDPVT=IDP
      IF (IFLAG .EQ. 0  .AND.  NIV .EQ. LNFC) GO TO 5
      IFLAG=99
      RETURN
    5 DO 6 K=1,NCOMP
    6 YHP(K,LNFC+1)=YHP(K,LNFCC+1)
      IFLAG=1
      RETURN
      END
