*DECK CBINU
      SUBROUTINE CBINU (Z, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  CBINU
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CAIRY, CBESH, CBESI, CBESJ, CBESK and CBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (CBINU-A, ZBINU-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     CBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
C
C***SEE ALSO  CAIRY, CBESH, CBESI, CBESJ, CBESK, CBIRY
C***ROUTINES CALLED  CASYI, CBUNI, CMLRI, CSERI, CUOIK, CWRSK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CBINU
      COMPLEX CW, CY, CZERO, Z
      REAL ALIM, AZ, DFNU, ELIM, FNU, FNUL, RL, TOL
      INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
      DIMENSION CY(N), CW(2)
      DATA CZERO / (0.0E0,0.0E0) /
C***FIRST EXECUTABLE STATEMENT  CBINU
      NZ = 0
      AZ = ABS(Z)
      NN = N
      DFNU = FNU + (N-1)
      IF (AZ.LE.2.0E0) GO TO 10
      IF (AZ*AZ*0.25E0.GT.DFNU+1.0E0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES
C-----------------------------------------------------------------------
      CALL CSERI(Z, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
      INW = ABS(NW)
      NZ = NZ + INW
      NN = NN - INW
      IF (NN.EQ.0) RETURN
      IF (NW.GE.0) GO TO 120
      DFNU = FNU + (NN-1)
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 40
      IF (DFNU.LE.1.0E0) GO TO 30
      IF (AZ+AZ.LT.DFNU*DFNU) GO TO 50
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z
C-----------------------------------------------------------------------
   30 CONTINUE
      CALL CASYI(Z, FNU, KODE, NN, CY, NW, RL, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
   40 CONTINUE
      IF (DFNU.LE.1.0E0) GO TO 70
   50 CONTINUE
C-----------------------------------------------------------------------
C     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      CALL CUOIK(Z, FNU, KODE, 1, NN, CY, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      NN = NN - NW
      IF (NN.EQ.0) RETURN
      DFNU = FNU+(NN-1)
      IF (DFNU.GT.FNUL) GO TO 110
      IF (AZ.GT.FNUL) GO TO 110
   60 CONTINUE
      IF (AZ.GT.RL) GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES
C-----------------------------------------------------------------------
      CALL CMLRI(Z, FNU, KODE, NN, CY, NW, TOL)
      IF(NW.LT.0) GO TO 130
      GO TO 120
   80 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
C-----------------------------------------------------------------------
      CALL CUOIK(Z, FNU, KODE, 2, 2, CW, NW, TOL, ELIM, ALIM)
      IF (NW.GE.0) GO TO 100
      NZ = NN
      DO 90 I=1,NN
        CY(I) = CZERO
   90 CONTINUE
      RETURN
  100 CONTINUE
      IF (NW.GT.0) GO TO 130
      CALL CWRSK(Z, FNU, KODE, NN, CY, NW, CW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
  110 CONTINUE
C-----------------------------------------------------------------------
C     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
C-----------------------------------------------------------------------
      NUI = FNUL-DFNU + 1
      NUI = MAX(NUI,0)
      CALL CBUNI(Z, FNU, KODE, NN, CY, NW, NUI, NLAST, FNUL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      IF (NLAST.EQ.0) GO TO 120
      NN = NLAST
      GO TO 60
  120 CONTINUE
      RETURN
  130 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
