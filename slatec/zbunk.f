*DECK ZBUNK
      SUBROUTINE ZBUNK (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZBUNK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CBUNI-A, ZBUNI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
C     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
C
C***SEE ALSO  ZBESH, ZBESK
C***ROUTINES CALLED  ZUNK1, ZUNK2
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZBUNK
C     COMPLEX Y,Z
      DOUBLE PRECISION ALIM, AX, AY, ELIM, FNU, TOL, YI, YR, ZI, ZR
      INTEGER KODE, MR, N, NZ
      DIMENSION YR(N), YI(N)
C***FIRST EXECUTABLE STATEMENT  ZBUNK
      NZ = 0
      AX = ABS(ZR)*1.7321D0
      AY = ABS(ZI)
      IF (AY.GT.AX) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
      RETURN
      END
