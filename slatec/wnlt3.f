*DECK WNLT3
      SUBROUTINE WNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
C***BEGIN PROLOGUE  WNLT3
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIT
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (WNLT3-S, DWNLT3-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     Perform column interchange.
C     Exchange elements of permuted index vector and perform column
C     interchanges.
C
C***SEE ALSO  WNLIT
C***ROUTINES CALLED  SSWAP
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLT and made a subroutine.  (RWC))
C***END PROLOGUE  WNLT3
      INTEGER I, IMAX, IPIVOT(*), M, MDW
      REAL             H(*), W(MDW,*)
C
      EXTERNAL SSWAP
C
      REAL             T
      INTEGER ITEMP
C
C***FIRST EXECUTABLE STATEMENT  WNLT3
      IF (IMAX.NE.I) THEN
         ITEMP        = IPIVOT(I)
         IPIVOT(I)    = IPIVOT(IMAX)
         IPIVOT(IMAX) = ITEMP
C
         CALL SSWAP(M, W(1,IMAX), 1, W(1,I), 1)
C
         T       = H(IMAX)
         H(IMAX) = H(I)
         H(I)    = T
      ENDIF
      RETURN
      END
