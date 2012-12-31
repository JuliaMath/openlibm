*DECK DXPNRM
      SUBROUTINE DXPNRM (NU1, NU2, MU1, MU2, PQA, IPQA, IERROR)
C***BEGIN PROLOGUE  DXPNRM
C***SUBSIDIARY
C***PURPOSE  To compute the values of Legendre functions for DXLEGF.
C            This subroutine transforms an array of Legendre functions
C            of the first kind of negative order stored in array PQA
C            into normalized Legendre polynomials stored in array PQA.
C            The original array is destroyed.
C***LIBRARY   SLATEC
C***CATEGORY  C3A2, C9
C***TYPE      DOUBLE PRECISION (XPNRM-S, DXPNRM-D)
C***KEYWORDS  LEGENDRE FUNCTIONS
C***AUTHOR  Smith, John M., (NBS and George Mason University)
C***ROUTINES CALLED  DXADJ
C***REVISION HISTORY  (YYMMDD)
C   820728  DATE WRITTEN
C   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  DXPNRM
      DOUBLE PRECISION C1,DMU,NU,NU1,NU2,PQA,PROD
      DIMENSION PQA(*),IPQA(*)
C***FIRST EXECUTABLE STATEMENT  DXPNRM
      IERROR=0
      L=(MU2-MU1)+(NU2-NU1+1.5D0)
      MU=MU1
      DMU=MU1
      NU=NU1
C
C         IF MU .GT.NU, NORM P =0.
C
      J=1
  500 IF(DMU.LE.NU) GO TO 505
      PQA(J)=0.D0
      IPQA(J)=0
      J=J+1
      IF(J.GT.L) RETURN
C
C        INCREMENT EITHER MU OR NU AS APPROPRIATE.
C
      IF(MU2.GT.MU1) DMU=DMU+1.D0
      IF(NU2-NU1.GT..5D0) NU=NU+1.D0
      GO TO 500
C
C         TRANSFORM P(-MU,NU,X) INTO NORMALIZED P(MU,NU,X) USING
C              NORM P(MU,NU,X)=
C                 SQRT((NU+.5)*FACTORIAL(NU+MU)/FACTORIAL(NU-MU))
C                              *P(-MU,NU,X)
C
  505 PROD=1.D0
      IPROD=0
      K=2*MU
      IF(K.LE.0) GO TO 520
      DO 510 I=1,K
      PROD=PROD*SQRT(NU+DMU+1.D0-I)
  510 CALL DXADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
  520 DO 540 I=J,L
      C1=PROD*SQRT(NU+.5D0)
      PQA(I)=PQA(I)*C1
      IPQA(I)=IPQA(I)+IPROD
      CALL DXADJ(PQA(I),IPQA(I),IERROR)
      IF (IERROR.NE.0) RETURN
      IF(NU2-NU1.GT..5D0) GO TO 530
      IF(DMU.GE.NU) GO TO 525
      PROD=SQRT(NU+DMU+1.D0)*PROD
      IF(NU.GT.DMU) PROD=PROD*SQRT(NU-DMU)
      CALL DXADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
      MU=MU+1
      DMU=DMU+1.D0
      GO TO 540
  525 PROD=0.D0
      IPROD=0
      MU=MU+1
      DMU=DMU+1.D0
      GO TO 540
  530 PROD=SQRT(NU+DMU+1.D0)*PROD
      IF(NU.NE.DMU-1.D0) PROD=PROD/SQRT(NU-DMU+1.D0)
      CALL DXADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
      NU=NU+1.D0
  540 CONTINUE
      RETURN
      END
