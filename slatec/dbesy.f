*DECK DBESY
      SUBROUTINE DBESY (X, FNU, N, Y)
C***BEGIN PROLOGUE  DBESY
C***PURPOSE  Implement forward recursion on the three term recursion
C            relation for a sequence of non-negative order Bessel
C            functions Y/SUB(FNU+I-1)/(X), I=1,...,N for real, positive
C            X and non-negative orders FNU.
C***LIBRARY   SLATEC
C***CATEGORY  C10A3
C***TYPE      DOUBLE PRECISION (BESY-S, DBESY-D)
C***KEYWORDS  SPECIAL FUNCTIONS, Y BESSEL FUNCTION
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract  **** a double precision routine ****
C         DBESY implements forward recursion on the three term
C         recursion relation for a sequence of non-negative order Bessel
C         functions Y/sub(FNU+I-1)/(X), I=1,N for real X .GT. 0.0D0 and
C         non-negative orders FNU.  If FNU .LT. NULIM, orders FNU and
C         FNU+1 are obtained from DBSYNU which computes by a power
C         series for X .LE. 2, the K Bessel function of an imaginary
C         argument for 2 .LT. X .LE. 20 and the asymptotic expansion for
C         X .GT. 20.
C
C         If FNU .GE. NULIM, the uniform asymptotic expansion is coded
C         in DASYJY for orders FNU and FNU+1 to start the recursion.
C         NULIM is 70 or 100 depending on whether N=1 or N .GE. 2.  An
C         overflow test is made on the leading term of the asymptotic
C         expansion before any extensive computation is done.
C
C         The maximum number of significant digits obtainable
C         is the smaller of 14 and the number of digits carried in
C         double precision arithmetic.
C
C     Description of Arguments
C
C         Input
C           X      - X .GT. 0.0D0
C           FNU    - order of the initial Y function, FNU .GE. 0.0D0
C           N      - number of members in the sequence, N .GE. 1
C
C         Output
C           Y      - a vector whose first N components contain values
C                    for the sequence Y(I)=Y/sub(FNU+I-1)/(X), I=1,N.
C
C     Error Conditions
C         Improper input arguments - a fatal error
C         Overflow - a fatal error
C
C***REFERENCES  F. W. J. Olver, Tables of Bessel Functions of Moderate
C                 or Large Orders, NPL Mathematical Tables 6, Her
C                 Majesty's Stationery Office, London, 1962.
C               N. M. Temme, On the numerical evaluation of the modified
C                 Bessel function of the third kind, Journal of
C                 Computational Physics 19, (1975), pp. 324-337.
C               N. M. Temme, On the numerical evaluation of the ordinary
C                 Bessel function of the second kind, Journal of
C                 Computational Physics 21, (1976), pp. 343-350.
C***ROUTINES CALLED  D1MACH, DASYJY, DBESY0, DBESY1, DBSYNU, DYAIRY,
C                    I1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800501  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBESY
C
      EXTERNAL DYAIRY
      INTEGER I, IFLW, J, N, NB, ND, NN, NUD, NULIM
      INTEGER I1MACH
      DOUBLE PRECISION AZN,CN,DNU,ELIM,FLGJY,FN,FNU,RAN,S,S1,S2,TM,TRX,
     1           W,WK,W2N,X,XLIM,XXN,Y
      DOUBLE PRECISION DBESY0, DBESY1, D1MACH
      DIMENSION W(2), NULIM(2), Y(*), WK(7)
      SAVE NULIM
      DATA NULIM(1),NULIM(2) / 70 , 100 /
C***FIRST EXECUTABLE STATEMENT  DBESY
      NN = -I1MACH(15)
      ELIM = 2.303D0*(NN*D1MACH(5)-3.0D0)
      XLIM = D1MACH(1)*1.0D+3
      IF (FNU.LT.0.0D0) GO TO 140
      IF (X.LE.0.0D0) GO TO 150
      IF (X.LT.XLIM) GO TO 170
      IF (N.LT.1) GO TO 160
C
C     ND IS A DUMMY VARIABLE FOR N
C
      ND = N
      NUD = INT(FNU)
      DNU = FNU - NUD
      NN = MIN(2,ND)
      FN = FNU + N - 1
      IF (FN.LT.2.0D0) GO TO 100
C
C     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
C     FOR THE LAST ORDER, FNU+N-1.GE.NULIM
C
      XXN = X/FN
      W2N = 1.0D0-XXN*XXN
      IF(W2N.LE.0.0D0) GO TO 10
      RAN = SQRT(W2N)
      AZN = LOG((1.0D0+RAN)/XXN) - RAN
      CN = FN*AZN
      IF(CN.GT.ELIM) GO TO 170
   10 CONTINUE
      IF (NUD.LT.NULIM(NN)) GO TO 20
C
C     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1.GE.NULIM
C
      FLGJY = -1.0D0
      CALL DASYJY(DYAIRY,X,FNU,FLGJY,NN,Y,WK,IFLW)
      IF(IFLW.NE.0) GO TO 170
      IF (NN.EQ.1) RETURN
      TRX = 2.0D0/X
      TM = (FNU+FNU+2.0D0)/X
      GO TO 80
C
   20 CONTINUE
      IF (DNU.NE.0.0D0) GO TO 30
      S1 = DBESY0(X)
      IF (NUD.EQ.0 .AND. ND.EQ.1) GO TO 70
      S2 = DBESY1(X)
      GO TO 40
   30 CONTINUE
      NB = 2
      IF (NUD.EQ.0 .AND. ND.EQ.1) NB = 1
      CALL DBSYNU(X, DNU, NB, W)
      S1 = W(1)
      IF (NB.EQ.1) GO TO 70
      S2 = W(2)
   40 CONTINUE
      TRX = 2.0D0/X
      TM = (DNU+DNU+2.0D0)/X
C     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2)
      IF (ND.EQ.1) NUD = NUD - 1
      IF (NUD.GT.0) GO TO 50
      IF (ND.GT.1) GO TO 70
      S1 = S2
      GO TO 70
   50 CONTINUE
      DO 60 I=1,NUD
        S = S2
        S2 = TM*S2 - S1
        S1 = S
        TM = TM + TRX
   60 CONTINUE
      IF (ND.EQ.1) S1 = S2
   70 CONTINUE
      Y(1) = S1
      IF (ND.EQ.1) RETURN
      Y(2) = S2
   80 CONTINUE
      IF (ND.EQ.2) RETURN
C     FORWARD RECUR FROM FNU+2 TO FNU+N-1
      DO 90 I=3,ND
        Y(I) = TM*Y(I-1) - Y(I-2)
        TM = TM + TRX
   90 CONTINUE
      RETURN
C
  100 CONTINUE
C     OVERFLOW TEST
      IF (FN.LE.1.0D0) GO TO 110
      IF (-FN*(LOG(X)-0.693D0).GT.ELIM) GO TO 170
  110 CONTINUE
      IF (DNU.EQ.0.0D0) GO TO 120
      CALL DBSYNU(X, FNU, ND, Y)
      RETURN
  120 CONTINUE
      J = NUD
      IF (J.EQ.1) GO TO 130
      J = J + 1
      Y(J) = DBESY0(X)
      IF (ND.EQ.1) RETURN
      J = J + 1
  130 CONTINUE
      Y(J) = DBESY1(X)
      IF (ND.EQ.1) RETURN
      TRX = 2.0D0/X
      TM = TRX
      GO TO 80
C
C
C
  140 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESY', 'ORDER, FNU, LESS THAN ZERO', 2,
     +   1)
      RETURN
  150 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESY', 'X LESS THAN OR EQUAL TO ZERO',
     +   2, 1)
      RETURN
  160 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESY', 'N LESS THAN ONE', 2, 1)
      RETURN
  170 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESY',
     +   'OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL', 6, 1)
      RETURN
      END
