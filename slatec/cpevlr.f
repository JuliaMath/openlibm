*DECK CPEVLR
      SUBROUTINE CPEVLR (N, M, A, X, C)
C***BEGIN PROLOGUE  CPEVLR
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CPZERO
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (CPEVLR-S)
C***AUTHOR  (UNKNOWN)
C***SEE ALSO  CPZERO
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810223  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CPEVLR
      REAL A(*),C(*)
C***FIRST EXECUTABLE STATEMENT  CPEVLR
      NP1=N+1
      DO 1 J=1,NP1
            CI=0.0
            CIM1=A(J)
            MINI=MIN(M+1,N+2-J)
            DO 1 I=1,MINI
               IF(J .NE. 1) CI=C(I)
               IF(I .NE. 1) CIM1=C(I-1)
               C(I)=CIM1+X*CI
    1 CONTINUE
      RETURN
      END
