*DECK BSRH
      FUNCTION BSRH (XLL, XRR, IZ, C, A, BH, F, SGN)
C***BEGIN PROLOGUE  BSRH
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BLKTRI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BCRH-S, BSRH-S)
C***AUTHOR  (UNKNOWN)
C***SEE ALSO  BLKTRI
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    CBLKT
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  BSRH
      DIMENSION       A(*)       ,C(*)       ,BH(*)
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
C***FIRST EXECUTABLE STATEMENT  BSRH
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
  105 BSRH = .5*(XL+XR)
      RETURN
      END
