*DECK BCRH
      FUNCTION BCRH (XLL, XRR, IZ, C, A, BH, F, SGN)
C***BEGIN PROLOGUE  BCRH
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CBLKTR
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BCRH-S, BSRH-S)
C***AUTHOR  (UNKNOWN)
C***SEE ALSO  CBLKTR
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    CCBLK
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  BCRH
      DIMENSION       A(*)       ,C(*)       ,BH(*)
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
C***FIRST EXECUTABLE STATEMENT  BCRH
      XL = XLL
      XR = XRR
      DX = .5*ABS(XR-XL)
  101 X = .5*(XL+XR)
      IF (SGN*F(X,IZ,C,A,BH)) 103,105,102
  102 XR = X
      GO TO 104
  103 XL = X
  104 DX = .5*DX
      IF (DX-CNV) 105,105,101
  105 BCRH = .5*(XL+XR)
      RETURN
      END
