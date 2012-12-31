*DECK DHVNRM
      DOUBLE PRECISION FUNCTION DHVNRM (V, NCOMP)
C***BEGIN PROLOGUE  DHVNRM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEABM, DDEBDF and DDERKF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (HVNRM-S, DHVNRM-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     Compute the maximum norm of the vector V(*) of length NCOMP and
C     return the result as DHVNRM
C
C***SEE ALSO  DDEABM, DDEBDF, DDERKF
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891024  Changed references from DVNORM to DHVNRM.  (WRB)
C   891024  Changed routine name from DVNORM to DHVNRM.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DHVNRM
C
      INTEGER K, NCOMP
      DOUBLE PRECISION V
      DIMENSION V(*)
C***FIRST EXECUTABLE STATEMENT  DHVNRM
      DHVNRM = 0.0D0
      DO 10 K = 1, NCOMP
         DHVNRM = MAX(DHVNRM,ABS(V(K)))
   10 CONTINUE
      RETURN
      END
