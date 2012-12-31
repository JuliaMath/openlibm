*DECK VNWRMS
      REAL FUNCTION VNWRMS (N, V, W)
C***BEGIN PROLOGUE  VNWRMS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DEBDF
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (VNWRMS-S, DVNRMS-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   VNWRMS computes a weighted root-mean-square vector norm for the
C   integrator package DEBDF.
C
C***SEE ALSO  DEBDF
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  VNWRMS
C
C
CLLL. OPTIMIZE
C-----------------------------------------------------------------------
C THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED ROOT-MEAN-SQUARE NORM
C OF THE VECTOR OF LENGTH N CONTAINED IN THE ARRAY V, WITH WEIGHTS
C CONTAINED IN THE ARRAY W OF LENGTH N..
C   VNWRMS = SQRT( (1/N) * SUM( V(I)/W(I) )**2 )
C-----------------------------------------------------------------------
      INTEGER N, I
      REAL V, W, SUM
      DIMENSION V(*), W(*)
C***FIRST EXECUTABLE STATEMENT  VNWRMS
      SUM = 0.0E0
      DO 10 I = 1,N
 10     SUM = SUM + (V(I)/W(I))**2
      VNWRMS = SQRT(SUM/N)
      RETURN
C----------------------- END OF FUNCTION VNWRMS ------------------------
      END
