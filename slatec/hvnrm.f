*DECK HVNRM
      FUNCTION HVNRM (V, NCOMP)
C***BEGIN PROLOGUE  HVNRM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DEABM, DEBDF and DERKF
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (HVNRM-S, DHVNRM-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     Compute the maximum norm of the vector V(*) of length NCOMP and
C     return the result as HVNRM.
C
C***SEE ALSO  DEABM, DEBDF, DERKF
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800501  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891024  Changed routine name from VNORM to HVNRM.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  HVNRM
      DIMENSION V(*)
C***FIRST EXECUTABLE STATEMENT  HVNRM
      HVNRM=0.
      DO 10 K=1,NCOMP
   10   HVNRM=MAX(HVNRM,ABS(V(K)))
      RETURN
      END
