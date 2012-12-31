*DECK MACON
      SUBROUTINE MACON
C***BEGIN PROLOGUE  MACON
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (MACON-S, DMACON-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C    Sets up machine constants using R1MACH
C
C***SEE ALSO  BVSUP
C***ROUTINES CALLED  R1MACH
C***COMMON BLOCKS    ML5MCO
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  MACON
      COMMON /ML5MCO/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
C***FIRST EXECUTABLE STATEMENT  MACON
      URO=R1MACH(4)
      SRU=SQRT(URO)
      DD=-LOG10(URO)
      LPAR=0.5*DD
      KE=0.5+0.75*DD
      EPS=10.**(-2*KE)
      SQOVFL=SQRT(R1MACH(2))
      TWOU=2.0*URO
      FOURU=4.0*URO
      RETURN
      END
