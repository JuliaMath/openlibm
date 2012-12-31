*DECK CBUNK
      SUBROUTINE CBUNK (Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  CBUNK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CBESH and CBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CBUNK-A, ZBUNK-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     CBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
C     IN CUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN CUNK2
C
C***SEE ALSO  CBESH, CBESK
C***ROUTINES CALLED  CUNK1, CUNK2
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CBUNK
      COMPLEX Y, Z
      REAL ALIM, AX, AY, ELIM, FNU, TOL, XX, YY
      INTEGER KODE, MR, N, NZ
      DIMENSION Y(N)
C***FIRST EXECUTABLE STATEMENT  CBUNK
      NZ = 0
      XX = REAL(Z)
      YY = AIMAG(Z)
      AX = ABS(XX)*1.7321E0
      AY = ABS(YY)
      IF (AY.GT.AX) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL CUNK1(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL CUNK2(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
      RETURN
      END
