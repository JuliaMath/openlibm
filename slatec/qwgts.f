*DECK QWGTS
      REAL FUNCTION QWGTS (X, A, B, ALFA, BETA, INTEGR)
C***BEGIN PROLOGUE  QWGTS
C***SUBSIDIARY
C***PURPOSE  This function subprogram is used together with the
C            routine QAWS and defines the WEIGHT function.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (QWGTS-S, DQWGTS-D)
C***KEYWORDS  ALGEBRAICO-LOGARITHMIC, END POINT SINGULARITIES,
C             WEIGHT FUNCTION
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***SEE ALSO  QK15W
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  QWGTS
C
      REAL A,ALFA,B,BETA,BMX,X,XMA
      INTEGER INTEGR
C***FIRST EXECUTABLE STATEMENT  QWGTS
      XMA = X-A
      BMX = B-X
      QWGTS = XMA**ALFA*BMX**BETA
      GO TO (40,10,20,30),INTEGR
   10 QWGTS = QWGTS*LOG(XMA)
      GO TO 40
   20 QWGTS = QWGTS*LOG(BMX)
      GO TO 40
   30 QWGTS = QWGTS*LOG(XMA)*LOG(BMX)
   40 RETURN
      END
