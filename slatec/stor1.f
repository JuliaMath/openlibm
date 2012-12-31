*DECK STOR1
      SUBROUTINE STOR1 (U, YH, V, YP, NTEMP, NDISK, NTAPE)
C***BEGIN PROLOGUE  STOR1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (STOR1-S, DSTOR1-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C             0 -- Storage at output points.
C     NTEMP =
C             1 -- Temporary storage
C **********************************************************************
C
C***SEE ALSO  BVSUP
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    ML8SZ
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  STOR1
      DIMENSION U(*),YH(*),V(*),YP(*)
C
C **********************************************************************
C
      COMMON /ML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
C
C **********************************************************************
C
C***FIRST EXECUTABLE STATEMENT  STOR1
      NCTNF = NCOMP * NFC
      DO 10 J = 1,NCTNF
   10 U(J) = YH(J)
      IF (INHOMO .EQ. 1)  GO TO 30
C
C   ZERO PARTICULAR SOLUTION
C
      IF (NTEMP .EQ. 1)  RETURN
      DO 20 J = 1,NCOMP
   20 V(J) = 0.
      GO TO 70
C
C   NONZERO PARTICULAR SOLUTION
C
   30 IF (NTEMP .EQ. 0)  GO TO 50
C
      DO 40 J = 1,NCOMP
   40 V(J) = YP(J)
      RETURN
C
   50 DO 60 J = 1,NCOMP
   60 V(J) = C * YP(J)
C
C  IS OUTPUT INFORMATION TO BE WRITTEN TO DISK
C
   70 IF (NDISK .EQ. 1)  WRITE (NTAPE) (V(J),J=1,NCOMP),(U(J),J=1,NCTNF)
C
      RETURN
      END
