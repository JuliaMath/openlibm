*DECK XQMU
      SUBROUTINE XQMU (NU1, NU2, MU1, MU2, THETA, X, SX, ID, PQA, IPQA,
     1   IERROR)
C***BEGIN PROLOGUE  XQMU
C***SUBSIDIARY
C***PURPOSE  To compute the values of Legendre functions for XLEGF.
C            Method: forward mu-wise recurrence for Q(MU,NU,X) for fixed
C            nu to obtain Q(MU1,NU,X), Q(MU1+1,NU,X), ..., Q(MU2,NU,X).
C***LIBRARY   SLATEC
C***CATEGORY  C3A2, C9
C***TYPE      SINGLE PRECISION (XQMU-S, DXQMU-D)
C***KEYWORDS  LEGENDRE FUNCTIONS
C***AUTHOR  Smith, John M., (NBS and George Mason University)
C***ROUTINES CALLED  XADD, XADJ, XPQNU
C***REVISION HISTORY  (YYMMDD)
C   820728  DATE WRITTEN
C   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  XQMU
      DIMENSION PQA(*),IPQA(*)
      REAL DMU,NU,NU1,NU2,PQ,PQA,PQ1,PQ2,SX,X,X1,X2
      REAL THETA
C***FIRST EXECUTABLE STATEMENT  XQMU
      IERROR=0
      MU=0
C
C        CALL XPQNU TO OBTAIN Q(0.,NU1,X)
C
      CALL XPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      PQ2=PQA(1)
      IPQ2=IPQA(1)
      MU=1
C
C        CALL XPQNU TO OBTAIN Q(1.,NU1,X)
C
      CALL XPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      NU=NU1
      K=0
      MU=1
      DMU=1.
      PQ1=PQA(1)
      IPQ1=IPQA(1)
      IF(MU1.GT.0) GO TO 310
      K=K+1
      PQA(K)=PQ2
      IPQA(K)=IPQ2
      IF(MU2.LT.1) GO TO 330
  310 IF(MU1.GT.1) GO TO 320
      K=K+1
      PQA(K)=PQ1
      IPQA(K)=IPQ1
      IF(MU2.LE.1) GO TO 330
  320 CONTINUE
C
C        FORWARD RECURRENCE IN MU TO OBTAIN
C                  Q(MU1,NU,X),Q(MU1+1,NU,X),....,Q(MU2,NU,X) USING
C             Q(MU+1,NU,X)=-2.*MU*X*SQRT(1./(1.-X**2))*Q(MU,NU,X)
C                               -(NU+MU)*(NU-MU+1.)*Q(MU-1,NU,X)
C
      X1=-2.*DMU*X*SX*PQ1
      X2=(NU+DMU)*(NU-DMU+1.)*PQ2
      CALL XADD(X1,IPQ1,-X2,IPQ2,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      CALL XADJ(PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      PQ2=PQ1
      IPQ2=IPQ1
      PQ1=PQ
      IPQ1=IPQ
      MU=MU+1
      DMU=DMU+1.
      IF(MU.LT.MU1) GO TO 320
      K=K+1
      PQA(K)=PQ
      IPQA(K)=IPQ
      IF(MU2.GT.MU) GO TO 320
  330 RETURN
      END
