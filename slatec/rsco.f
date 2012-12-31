*DECK RSCO
      SUBROUTINE RSCO (RSAV, ISAV)
C***BEGIN PROLOGUE  RSCO
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DEBDF
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (RSCO-S, DRSCO-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   RSCO transfers data from arrays to a common block within the
C   integrator package DEBDF.
C
C***SEE ALSO  DEBDF
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DEBDF1
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  RSCO
C
C
C-----------------------------------------------------------------------
C THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON
C BLOCK DEBDF1  , WHICH IS USED INTERNALLY IN THE DEBDF
C PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS
C OF SUBROUTINE SVCO OR THE EQUIVALENT.
C-----------------------------------------------------------------------
      INTEGER ISAV, I,      ILS, LENILS, LENRLS
      REAL RSAV, RLS
      DIMENSION RSAV(*), ISAV(*)
      COMMON /DEBDF1/ RLS(218), ILS(33)
      SAVE LENRLS, LENILS
      DATA LENRLS/218/, LENILS/33/
C
C***FIRST EXECUTABLE STATEMENT  RSCO
      DO 10 I = 1,LENRLS
 10     RLS(I) = RSAV(I)
      DO 20 I = 1,LENILS
 20     ILS(I) = ISAV(I)
      RETURN
C----------------------- END OF SUBROUTINE RSCO -----------------------
      END
