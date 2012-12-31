*DECK DASYIK
      SUBROUTINE DASYIK (X, FNU, KODE, FLGIK, RA, ARG, IN, Y)
C***BEGIN PROLOGUE  DASYIK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBESI and DBESK
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (ASYIK-S, DASYIK-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C                    DASYIK computes Bessel functions I and K
C                  for arguments X.GT.0.0 and orders FNU.GE.35
C                  on FLGIK = 1 and FLGIK = -1 respectively.
C
C                                    INPUT
C
C      X    - Argument, X.GT.0.0D0
C      FNU  - Order of first Bessel function
C      KODE - A parameter to indicate the scaling option
C             KODE=1 returns Y(I)=        I/SUB(FNU+I-1)/(X), I=1,IN
C                    or      Y(I)=        K/SUB(FNU+I-1)/(X), I=1,IN
C                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
C             KODE=2 returns Y(I)=EXP(-X)*I/SUB(FNU+I-1)/(X), I=1,IN
C                    or      Y(I)=EXP( X)*K/SUB(FNU+I-1)/(X), I=1,IN
C                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
C     FLGIK - Selection parameter for I or K FUNCTION
C             FLGIK =  1.0D0 gives the I function
C             FLGIK = -1.0D0 gives the K function
C        RA - SQRT(1.+Z*Z), Z=X/FNU
C       ARG - Argument of the leading exponential
C        IN - Number of functions desired, IN=1 or 2
C
C                                    OUTPUT
C
C         Y - A vector whose first IN components contain the sequence
C
C     Abstract  **** A double precision routine ****
C         DASYIK implements the uniform asymptotic expansion of
C         the I and K Bessel functions for FNU.GE.35 and real
C         X.GT.0.0D0. The forms are identical except for a change
C         in sign of some of the terms. This change in sign is
C         accomplished by means of the FLAG FLGIK = 1 or -1.
C
C***SEE ALSO  DBESI, DBESK
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR section.  (WRB)
C***END PROLOGUE  DASYIK
C
      INTEGER IN, J, JN, K, KK, KODE, L
      DOUBLE PRECISION AK,AP,ARG,C,COEF,CON,ETX,FLGIK,FN,FNU,GLN,RA,
     1 S1, S2, T, TOL, T2, X, Y, Z
      DOUBLE PRECISION D1MACH
      DIMENSION Y(*), C(65), CON(2)
      SAVE CON, C
      DATA CON(1), CON(2)  /
     1        3.98942280401432678D-01,    1.25331413731550025D+00/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3       -2.08333333333333D-01,        1.25000000000000D-01,
     4        3.34201388888889D-01,       -4.01041666666667D-01,
     5        7.03125000000000D-02,       -1.02581259645062D+00,
     6        1.84646267361111D+00,       -8.91210937500000D-01,
     7        7.32421875000000D-02,        4.66958442342625D+00,
     8       -1.12070026162230D+01,        8.78912353515625D+00,
     9       -2.36408691406250D+00,        1.12152099609375D-01,
     1       -2.82120725582002D+01,        8.46362176746007D+01,
     2       -9.18182415432400D+01,        4.25349987453885D+01,
     3       -7.36879435947963D+00,        2.27108001708984D-01,
     4        2.12570130039217D+02,       -7.65252468141182D+02,
     5        1.05999045252800D+03,       -6.99579627376133D+02/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3        2.18190511744212D+02,       -2.64914304869516D+01,
     4        5.72501420974731D-01,       -1.91945766231841D+03,
     5        8.06172218173731D+03,       -1.35865500064341D+04,
     6        1.16553933368645D+04,       -5.30564697861340D+03,
     7        1.20090291321635D+03,       -1.08090919788395D+02,
     8        1.72772750258446D+00,        2.02042913309661D+04,
     9       -9.69805983886375D+04,        1.92547001232532D+05,
     1       -2.03400177280416D+05,        1.22200464983017D+05,
     2       -4.11926549688976D+04,        7.10951430248936D+03,
     3       -4.93915304773088D+02,        6.07404200127348D+00,
     4       -2.42919187900551D+05,        1.31176361466298D+06,
     5       -2.99801591853811D+06,        3.76327129765640D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65)/
     3       -2.81356322658653D+06,        1.26836527332162D+06,
     4       -3.31645172484564D+05,        4.52187689813627D+04,
     5       -2.49983048181121D+03,        2.43805296995561D+01,
     6        3.28446985307204D+06,       -1.97068191184322D+07,
     7        5.09526024926646D+07,       -7.41051482115327D+07,
     8        6.63445122747290D+07,       -3.75671766607634D+07,
     9        1.32887671664218D+07,       -2.78561812808645D+06,
     1        3.08186404612662D+05,       -1.38860897537170D+04,
     2        1.10017140269247D+02/
C***FIRST EXECUTABLE STATEMENT  DASYIK
      TOL = D1MACH(3)
      TOL = MAX(TOL,1.0D-15)
      FN = FNU
      Z  = (3.0D0-FLGIK)/2.0D0
      KK = INT(Z)
      DO 50 JN=1,IN
        IF (JN.EQ.1) GO TO 10
        FN = FN - FLGIK
        Z = X/FN
        RA = SQRT(1.0D0+Z*Z)
        GLN = LOG((1.0D0+RA)/Z)
        ETX = KODE - 1
        T = RA*(1.0D0-ETX) + ETX/(Z+RA)
        ARG = FN*(T-GLN)*FLGIK
   10   COEF = EXP(ARG)
        T = 1.0D0/RA
        T2 = T*T
        T = T/FN
        T = SIGN(T,FLGIK)
        S2 = 1.0D0
        AP = 1.0D0
        L = 0
        DO 30 K=2,11
          L = L + 1
          S1 = C(L)
          DO 20 J=2,K
            L = L + 1
            S1 = S1*T2 + C(L)
   20     CONTINUE
          AP = AP*T
          AK = AP*S1
          S2 = S2 + AK
          IF (MAX(ABS(AK),ABS(AP)) .LT.TOL) GO TO 40
   30   CONTINUE
   40   CONTINUE
      T = ABS(T)
      Y(JN) = S2*COEF*SQRT(T)*CON(KK)
   50 CONTINUE
      RETURN
      END
