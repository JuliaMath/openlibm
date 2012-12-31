*DECK SVCO
      SUBROUTINE SVCO (RSAV, ISAV)
C***BEGIN PROLOGUE  SVCO
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DEBDF
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SVCO-S, DSVCO-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   SVCO transfers data from a common block to arrays within the
C   integrator package DEBDF.
C
C***SEE ALSO  DEBDF
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DEBDF1
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  SVCO
C
C
C-----------------------------------------------------------------------
C THIS ROUTINE STORES IN RSAV AND ISAV THE CONTENTS OF COMMON BLOCK
C DEBDF1  , WHICH IS USED INTERNALLY IN THE DEBDF PACKAGE.
C
C RSAV = REAL ARRAY OF LENGTH 218 OR MORE.
C ISAV = INTEGER ARRAY OF LENGTH 33 OR MORE.
C-----------------------------------------------------------------------
      INTEGER ISAV, I,      ILS, LENILS, LENRLS
      REAL RSAV, RLS
      DIMENSION RSAV(*), ISAV(*)
      COMMON /DEBDF1/ RLS(218), ILS(33)
      SAVE LENRLS, LENILS
      DATA LENRLS/218/, LENILS/33/
C
C***FIRST EXECUTABLE STATEMENT  SVCO
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      RETURN
C----------------------- END OF SUBROUTINE SVCO -----------------------
      END
