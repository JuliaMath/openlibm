*DECK DWNLT3
      SUBROUTINE DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
C***BEGIN PROLOGUE  DWNLT3
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIT
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLT3-S, DWNLT3-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     Perform column interchange.
C     Exchange elements of permuted index vector and perform column
C     interchanges.
C
C***SEE ALSO  DWNLIT
C***ROUTINES CALLED  DSWAP
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
C   900604  DP version created from SP version.  (RWC)
C***END PROLOGUE  DWNLT3
      INTEGER I, IMAX, IPIVOT(*), M, MDW
      DOUBLE PRECISION H(*), W(MDW,*)
C
      EXTERNAL DSWAP
C
      DOUBLE PRECISION T
      INTEGER ITEMP
C
C***FIRST EXECUTABLE STATEMENT  DWNLT3
      IF (IMAX.NE.I) THEN
         ITEMP        = IPIVOT(I)
         IPIVOT(I)    = IPIVOT(IMAX)
         IPIVOT(IMAX) = ITEMP
C
         CALL DSWAP(M, W(1,IMAX), 1, W(1,I), 1)
C
         T       = H(IMAX)
         H(IMAX) = H(I)
         H(I)    = T
      ENDIF
      RETURN
      END
