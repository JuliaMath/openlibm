*DECK DXQNU
      SUBROUTINE DXQNU (NU1, NU2, MU1, THETA, X, SX, ID, PQA, IPQA,
     1   IERROR)
C***BEGIN PROLOGUE  DXQNU
C***SUBSIDIARY
C***PURPOSE  To compute the values of Legendre functions for DXLEGF.
C            Method: backward nu-wise recurrence for Q(MU,NU,X) for
C            fixed mu to obtain Q(MU1,NU1,X), Q(MU1,NU1+1,X), ...,
C            Q(MU1,NU2,X).
C***LIBRARY   SLATEC
C***CATEGORY  C3A2, C9
C***TYPE      DOUBLE PRECISION (XQNU-S, DXQNU-D)
C***KEYWORDS  LEGENDRE FUNCTIONS
C***AUTHOR  Smith, John M., (NBS and George Mason University)
C***ROUTINES CALLED  DXADD, DXADJ, DXPQNU
C***REVISION HISTORY  (YYMMDD)
C   820728  DATE WRITTEN
C   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  DXQNU
      DIMENSION PQA(*),IPQA(*)
      DOUBLE PRECISION DMU,NU,NU1,NU2,PQ,PQA,PQ1,PQ2,SX,X,X1,X2
      DOUBLE PRECISION THETA,PQL1,PQL2
C***FIRST EXECUTABLE STATEMENT  DXQNU
      IERROR=0
      K=0
      PQ2=0.0D0
      IPQ2=0
      PQL2=0.0D0
      IPQL2=0
      IF(MU1.EQ.1) GO TO 290
      MU=0
C
C        CALL DXPQNU TO OBTAIN Q(0.,NU2,X) AND Q(0.,NU2-1,X)
C
      CALL DXPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      IF(MU1.EQ.0) RETURN
      K=(NU2-NU1+1.5D0)
      PQ2=PQA(K)
      IPQ2=IPQA(K)
      PQL2=PQA(K-1)
      IPQL2=IPQA(K-1)
  290 MU=1
C
C        CALL DXPQNU TO OBTAIN Q(1.,NU2,X) AND Q(1.,NU2-1,X)
C
      CALL DXPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      IF(MU1.EQ.1) RETURN
      NU=NU2
      PQ1=PQA(K)
      IPQ1=IPQA(K)
      PQL1=PQA(K-1)
      IPQL1=IPQA(K-1)
  300 MU=1
      DMU=1.D0
  320 CONTINUE
C
C        FORWARD RECURRENCE IN MU TO OBTAIN Q(MU1,NU2,X) AND
C              Q(MU1,NU2-1,X) USING
C              Q(MU+1,NU,X)=-2.*MU*X*SQRT(1./(1.-X**2))*Q(MU,NU,X)
C                   -(NU+MU)*(NU-MU+1.)*Q(MU-1,NU,X)
C
C              FIRST FOR NU=NU2
C
      X1=-2.D0*DMU*X*SX*PQ1
      X2=(NU+DMU)*(NU-DMU+1.D0)*PQ2
      CALL DXADD(X1,IPQ1,-X2,IPQ2,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      CALL DXADJ(PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      PQ2=PQ1
      IPQ2=IPQ1
      PQ1=PQ
      IPQ1=IPQ
      MU=MU+1
      DMU=DMU+1.D0
      IF(MU.LT.MU1) GO TO 320
      PQA(K)=PQ
      IPQA(K)=IPQ
      IF(K.EQ.1) RETURN
      IF(NU.LT.NU2) GO TO 340
C
C              THEN FOR NU=NU2-1
C
      NU=NU-1.D0
      PQ2=PQL2
      IPQ2=IPQL2
      PQ1=PQL1
      IPQ1=IPQL1
      K=K-1
      GO TO 300
C
C         BACKWARD RECURRENCE IN NU TO OBTAIN
C              Q(MU1,NU1,X),Q(MU1,NU1+1,X),....,Q(MU1,NU2,X)
C              USING
C              (NU-MU+1.)*Q(MU,NU+1,X)=
C                       (2.*NU+1.)*X*Q(MU,NU,X)-(NU+MU)*Q(MU,NU-1,X)
C
  340 PQ1=PQA(K)
      IPQ1=IPQA(K)
      PQ2=PQA(K+1)
      IPQ2=IPQA(K+1)
  350 IF(NU.LE.NU1) RETURN
      K=K-1
      X1=(2.D0*NU+1.D0)*X*PQ1/(NU+DMU)
      X2=-(NU-DMU+1.D0)*PQ2/(NU+DMU)
      CALL DXADD(X1,IPQ1,X2,IPQ2,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      CALL DXADJ(PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      PQ2=PQ1
      IPQ2=IPQ1
      PQ1=PQ
      IPQ1=IPQ
      PQA(K)=PQ
      IPQA(K)=IPQ
      NU=NU-1.D0
      GO TO 350
      END
