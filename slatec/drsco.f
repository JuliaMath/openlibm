*DECK DRSCO
      SUBROUTINE DRSCO (RSAV, ISAV)
C***BEGIN PROLOGUE  DRSCO
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (RSCO-S, DRSCO-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   DRSCO transfers data from arrays to a common block within the
C   integrator package DDEBDF.
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DRSCO
C-----------------------------------------------------------------------
C THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON
C BLOCK DDEBD1  , WHICH IS USED INTERNALLY IN THE DDEBDF
C PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS
C OF SUBROUTINE DSVCO OR THE EQUIVALENT.
C-----------------------------------------------------------------------
C
      INTEGER I, ILS, ISAV, LENILS, LENRLS
      DOUBLE PRECISION RLS, RSAV
      DIMENSION RSAV(*),ISAV(*)
      SAVE LENRLS, LENILS
      COMMON /DDEBD1/ RLS(218),ILS(33)
      DATA LENRLS /218/, LENILS /33/
C
C***FIRST EXECUTABLE STATEMENT  DRSCO
      DO 10 I = 1, LENRLS
         RLS(I) = RSAV(I)
   10 CONTINUE
      DO 20 I = 1, LENILS
         ILS(I) = ISAV(I)
   20 CONTINUE
      RETURN
C     ----------------------- END OF SUBROUTINE DRSCO
C     -----------------------
      END
