*DECK INDXB
      SUBROUTINE INDXB (I, IR, IDX, IDP)
C***BEGIN PROLOGUE  INDXB
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BLKTRI
C***LIBRARY   SLATEC
C***TYPE      INTEGER (INDXB-I)
C***AUTHOR  (UNKNOWN)
C***SEE ALSO  BLKTRI
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    CBLKT
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   920422  Added statement so IDX would always be defined.  (WRB)
C***END PROLOGUE  INDXB
C
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
C***FIRST EXECUTABLE STATEMENT  INDXB
      IDX = I
      IDP = 0
      IF (IR) 107,101,103
  101 IF (I-NM) 102,102,107
  102 IDX = I
      IDP = 1
      RETURN
  103 IZH = 2**IR
      ID = I-IZH-IZH
      IDX = ID+ID+(IR-1)*IK+IR+(IK-I)/IZH+4
      IPL = IZH-1
      IDP = IZH+IZH-1
      IF (I-IPL-NM) 105,105,104
  104 IDP = 0
      RETURN
  105 IF (I+IPL-NM) 107,107,106
  106 IDP = NM+IPL-I+1
  107 RETURN
      END
