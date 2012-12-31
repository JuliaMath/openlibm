*DECK IDLOC
      INTEGER FUNCTION IDLOC (LOC, SX, IX)
C***BEGIN PROLOGUE  IDLOC
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DSPLP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (IPLOC-S, IDLOC-D)
C***KEYWORDS  RELATIVE ADDRESS DETERMINATION FUNCTION, SLATEC
C***AUTHOR  Boland, W. Robert, (LANL)
C           Nicol, Tom, (University of British Columbia)
C***DESCRIPTION
C
C   Given a "virtual" location,  IDLOC returns the relative working
C   address of the vector component stored in SX, IX.  Any necessary
C   page swaps are performed automatically for the user in this
C   function subprogram.
C
C   LOC       is the "virtual" address of the data to be retrieved.
C   SX ,IX    represent the matrix where the data is stored.
C
C***SEE ALSO  DSPLP
C***ROUTINES CALLED  DPRWPG, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   890606  DATE WRITTEN
C   890606  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   910731  Added code to set IDLOC to 0 if LOC is non-positive.  (WRB)
C***END PROLOGUE  IDLOC
      DOUBLE PRECISION SX(*)
      INTEGER IX(*)
C***FIRST EXECUTABLE STATEMENT  IDLOC
      IF (LOC.LE.0) THEN
         CALL XERMSG ('SLATEC', 'IDLOC',
     +     'A value of LOC, the first argument, .LE. 0 was encountered',
     +     55, 1)
         IDLOC = 0
         RETURN
      ENDIF
C
C     Two cases exist:  (1.LE.LOC.LE.K) .OR. (LOC.GT.K).
C
      K = IX(3) + 4
      LMX = IX(1)
      LMXM1 = LMX - 1
      IF (LOC.LE.K) THEN
         IDLOC = LOC
         RETURN
      ENDIF
C
C     Compute length of the page, starting address of the page, page
C     number and relative working address.
C
      LPG = LMX-K
      ITEMP = LOC - K - 1
      IPAGE = ITEMP/LPG + 1
      IDLOC = MOD(ITEMP,LPG) + K + 1
      NP = ABS(IX(LMXM1))
C
C     Determine if a page fault has occurred.  If so, write page NP
C     and read page IPAGE.  Write the page only if it has been
C     modified.
C
      IF (IPAGE.NE.NP) THEN
         IF (SX(LMX).EQ.1.0) THEN
            SX(LMX) = 0.0
            KEY = 2
            CALL DPRWPG (KEY, NP, LPG, SX, IX)
         ENDIF
         KEY = 1
         CALL DPRWPG (KEY, IPAGE, LPG, SX, IX)
      ENDIF
      RETURN
      END
