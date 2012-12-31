*DECK DVNRMS
      DOUBLE PRECISION FUNCTION DVNRMS (N, V, W)
C***BEGIN PROLOGUE  DVNRMS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (VNWRMS-S, DVNRMS-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   DVNRMS computes a weighted root-mean-square vector norm for the
C   integrator package DDEBDF.
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DVNRMS
      INTEGER I, N
      DOUBLE PRECISION SUM, V, W
      DIMENSION V(*),W(*)
C***FIRST EXECUTABLE STATEMENT  DVNRMS
      SUM = 0.0D0
      DO 10 I = 1, N
         SUM = SUM + (V(I)/W(I))**2
   10 CONTINUE
      DVNRMS = SQRT(SUM/N)
      RETURN
C     ----------------------- END OF FUNCTION DVNRMS
C     ------------------------
      END
