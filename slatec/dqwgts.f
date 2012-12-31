*DECK DQWGTS
      DOUBLE PRECISION FUNCTION DQWGTS (X, A, B, ALFA, BETA, INTEGR)
C***BEGIN PROLOGUE  DQWGTS
C***SUBSIDIARY
C***PURPOSE  This function subprogram is used together with the
C            routine DQAWS and defines the WEIGHT function.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (QWGTS-S, DQWGTS-D)
C***KEYWORDS  ALGEBRAICO-LOGARITHMIC, END POINT SINGULARITIES,
C             WEIGHT FUNCTION
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***SEE ALSO  DQK15W
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DQWGTS
C
      DOUBLE PRECISION A,ALFA,B,BETA,BMX,X,XMA
      INTEGER INTEGR
C***FIRST EXECUTABLE STATEMENT  DQWGTS
      XMA = X-A
      BMX = B-X
      DQWGTS = XMA**ALFA*BMX**BETA
      GO TO (40,10,20,30),INTEGR
   10 DQWGTS = DQWGTS*LOG(XMA)
      GO TO 40
   20 DQWGTS = DQWGTS*LOG(BMX)
      GO TO 40
   30 DQWGTS = DQWGTS*LOG(XMA)*LOG(BMX)
   40 RETURN
      END
