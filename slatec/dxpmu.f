*DECK DXPMU
      SUBROUTINE DXPMU (NU1, NU2, MU1, MU2, THETA, X, SX, ID, PQA, IPQA,
     1   IERROR)
C***BEGIN PROLOGUE  DXPMU
C***SUBSIDIARY
C***PURPOSE  To compute the values of Legendre functions for DXLEGF.
C            Method: backward mu-wise recurrence for P(-MU,NU,X) for
C            fixed nu to obtain P(-MU2,NU1,X), P(-(MU2-1),NU1,X), ...,
C            P(-MU1,NU1,X) and store in ascending mu order.
C***LIBRARY   SLATEC
C***CATEGORY  C3A2, C9
C***TYPE      DOUBLE PRECISION (XPMU-S, DXPMU-D)
C***KEYWORDS  LEGENDRE FUNCTIONS
C***AUTHOR  Smith, John M., (NBS and George Mason University)
C***ROUTINES CALLED  DXADD, DXADJ, DXPQNU
C***REVISION HISTORY  (YYMMDD)
C   820728  DATE WRITTEN
C   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  DXPMU
      DOUBLE PRECISION PQA,NU1,NU2,P0,X,SX,THETA,X1,X2
      DIMENSION PQA(*),IPQA(*)
C
C        CALL DXPQNU TO OBTAIN P(-MU2,NU,X)
C
C***FIRST EXECUTABLE STATEMENT  DXPMU
      IERROR=0
      CALL DXPQNU(NU1,NU2,MU2,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      P0=PQA(1)
      IP0=IPQA(1)
      MU=MU2-1
C
C        CALL DXPQNU TO OBTAIN P(-MU2-1,NU,X)
C
      CALL DXPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      N=MU2-MU1+1
      PQA(N)=P0
      IPQA(N)=IP0
      IF(N.EQ.1) GO TO 300
      PQA(N-1)=PQA(1)
      IPQA(N-1)=IPQA(1)
      IF(N.EQ.2) GO TO 300
      J=N-2
  290 CONTINUE
C
C        BACKWARD RECURRENCE IN MU TO OBTAIN
C              P(-MU2,NU1,X),P(-(MU2-1),NU1,X),....P(-MU1,NU1,X)
C              USING
C              (NU-MU)*(NU+MU+1.)*P(-(MU+1),NU,X)=
C                2.*MU*X*SQRT((1./(1.-X**2))*P(-MU,NU,X)-P(-(MU-1),NU,X)
C
      X1=2.D0*MU*X*SX*PQA(J+1)
      X2=-(NU1-MU)*(NU1+MU+1.D0)*PQA(J+2)
      CALL DXADD(X1,IPQA(J+1),X2,IPQA(J+2),PQA(J),IPQA(J),IERROR)
      IF (IERROR.NE.0) RETURN
      CALL DXADJ(PQA(J),IPQA(J),IERROR)
      IF (IERROR.NE.0) RETURN
      IF(J.EQ.1) GO TO 300
      J=J-1
      MU=MU-1
      GO TO 290
  300 RETURN
      END
