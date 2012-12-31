*DECK CDIV
      SUBROUTINE CDIV (AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  CDIV
C***SUBSIDIARY
C***PURPOSE  Compute the complex quotient of two complex numbers.
C***LIBRARY   SLATEC
C***TYPE      COMPLEX (CDIV-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Complex division, (CR,CI) = (AR,AI)/(BR,BI)
C
C***SEE ALSO  EISDOC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   811101  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CDIV
      REAL AR,AI,BR,BI,CR,CI
C
      REAL S,ARS,AIS,BRS,BIS
C***FIRST EXECUTABLE STATEMENT  CDIV
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END
