*DECK XPMUP
      SUBROUTINE XPMUP (NU1, NU2, MU1, MU2, PQA, IPQA, IERROR)
C***BEGIN PROLOGUE  XPMUP
C***SUBSIDIARY
C***PURPOSE  To compute the values of Legendre functions for XLEGF.
C            This subroutine transforms an array of Legendre functions
C            of the first kind of negative order stored in array PQA
C            into Legendre functions of the first kind of positive
C            order stored in array PQA. The original array is destroyed.
C***LIBRARY   SLATEC
C***CATEGORY  C3A2, C9
C***TYPE      SINGLE PRECISION (XPMUP-S, DXPMUP-D)
C***KEYWORDS  LEGENDRE FUNCTIONS
C***AUTHOR  Smith, John M., (NBS and George Mason University)
C***ROUTINES CALLED  XADJ
C***REVISION HISTORY  (YYMMDD)
C   820728  DATE WRITTEN
C   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  XPMUP
      REAL DMU,NU,NU1,NU2,PQA,PROD
      DIMENSION PQA(*),IPQA(*)
C***FIRST EXECUTABLE STATEMENT  XPMUP
      IERROR=0
      NU=NU1
      MU=MU1
      DMU=MU
      N=INT(NU2-NU1+.1)+(MU2-MU1)+1
      J=1
      IF(MOD(NU,1.).NE.0.) GO TO 210
  200 IF(DMU.LT.NU+1.) GO TO 210
      PQA(J)=0.
      IPQA(J)=0
      J=J+1
      IF(J.GT.N) RETURN
C        INCREMENT EITHER MU OR NU AS APPROPRIATE.
      IF(NU2-NU1.GT..5) NU=NU+1.
      IF(MU2.GT.MU1) MU=MU+1
      GO TO 200
C
C        TRANSFORM P(-MU,NU,X) TO P(MU,NU,X) USING
C        P(MU,NU,X)=(NU-MU+1)*(NU-MU+2)*...*(NU+MU)*P(-MU,NU,X)*(-1)**MU
C
  210 PROD=1.
      IPROD=0
      K=2*MU
      IF(K.EQ.0) GO TO 222
      DO 220 L=1,K
      PROD=PROD*(DMU-NU-L)
  220 CALL XADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
  222 CONTINUE
      DO 240 I=J,N
      IF(MU.EQ.0) GO TO 225
      PQA(I)=PQA(I)*PROD*(-1)**MU
      IPQA(I)=IPQA(I)+IPROD
      CALL XADJ(PQA(I),IPQA(I),IERROR)
      IF (IERROR.NE.0) RETURN
  225 IF(NU2-NU1.GT..5) GO TO 230
      PROD=(DMU-NU)*PROD*(-DMU-NU-1.)
      CALL XADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
      MU=MU+1
      DMU=DMU+1.
      GO TO 240
  230 PROD=PROD*(-DMU-NU-1.)/(DMU-NU-1.)
      CALL XADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
      NU=NU+1.
  240 CONTINUE
      RETURN
      END
