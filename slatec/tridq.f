*DECK TRIDQ
      SUBROUTINE TRIDQ (MR, A, B, C, Y, D)
C***BEGIN PROLOGUE  TRIDQ
C***SUBSIDIARY
C***PURPOSE  Subsidiary to POIS3D
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (TRIDQ-S)
C***AUTHOR  (UNKNOWN)
C***SEE ALSO  POIS3D
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900308  Renamed routine from TRID to TRIDQ.  (WRB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  TRIDQ
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,Y(*)       ,
     1                D(*)
C***FIRST EXECUTABLE STATEMENT  TRIDQ
      M = MR
      MM1 = M-1
      Z = 1./B(1)
      D(1) = C(1)*Z
      Y(1) = Y(1)*Z
      DO 101 I=2,MM1
         Z = 1./(B(I)-A(I)*D(I-1))
         D(I) = C(I)*Z
         Y(I) = (Y(I)-A(I)*Y(I-1))*Z
  101 CONTINUE
      Z = B(M)-A(M)*D(MM1)
      IF (Z .NE. 0.) GO TO 102
      Y(M) = 0.
      GO TO 103
  102 Y(M) = (Y(M)-A(M)*Y(MM1))/Z
  103 CONTINUE
      DO 104 IP=1,MM1
         I = M-IP
         Y(I) = Y(I)-D(I)*Y(I+1)
  104 CONTINUE
      RETURN
      END
