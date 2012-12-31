*DECK ASYIK
      SUBROUTINE ASYIK (X, FNU, KODE, FLGIK, RA, ARG, IN, Y)
C***BEGIN PROLOGUE  ASYIK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BESI and BESK
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (ASYIK-S, DASYIK-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C                    ASYIK computes Bessel functions I and K
C                  for arguments X.GT.0.0 and orders FNU.GE.35
C                  on FLGIK = 1 and FLGIK = -1 respectively.
C
C                                    INPUT
C
C      X    - argument, X.GT.0.0E0
C      FNU  - order of first Bessel function
C      KODE - a parameter to indicate the scaling option
C             KODE=1 returns Y(I)=        I/SUB(FNU+I-1)/(X), I=1,IN
C                    or      Y(I)=        K/SUB(FNU+I-1)/(X), I=1,IN
C                    on FLGIK = 1.0E0 or FLGIK = -1.0E0
C             KODE=2 returns Y(I)=EXP(-X)*I/SUB(FNU+I-1)/(X), I=1,IN
C                    or      Y(I)=EXP( X)*K/SUB(FNU+I-1)/(X), I=1,IN
C                    on FLGIK = 1.0E0 or FLGIK = -1.0E0
C     FLGIK - selection parameter for I or K function
C             FLGIK =  1.0E0 gives the I function
C             FLGIK = -1.0E0 gives the K function
C        RA - SQRT(1.+Z*Z), Z=X/FNU
C       ARG - argument of the leading exponential
C        IN - number of functions desired, IN=1 or 2
C
C                                    OUTPUT
C
C         Y - a vector whose first in components contain the sequence
C
C     Abstract
C         ASYIK implements the uniform asymptotic expansion of
C         the I and K Bessel functions for FNU.GE.35 and real
C         X.GT.0.0E0. The forms are identical except for a change
C         in sign of some of the terms. This change in sign is
C         accomplished by means of the flag FLGIK = 1 or -1.
C
C***SEE ALSO  BESI, BESK
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR section.  (WRB)
C***END PROLOGUE  ASYIK
C
      INTEGER IN, J, JN, K, KK, KODE, L
      REAL AK,AP,ARG,C, COEF,CON,ETX,FLGIK,FN, FNU,GLN,RA,S1,S2,
     1 T, TOL, T2, X, Y, Z
      REAL R1MACH
      DIMENSION Y(*), C(65), CON(2)
      SAVE CON, C
      DATA CON(1), CON(2)  /
     1        3.98942280401432678E-01,    1.25331413731550025E+00/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3       -2.08333333333333E-01,        1.25000000000000E-01,
     4        3.34201388888889E-01,       -4.01041666666667E-01,
     5        7.03125000000000E-02,       -1.02581259645062E+00,
     6        1.84646267361111E+00,       -8.91210937500000E-01,
     7        7.32421875000000E-02,        4.66958442342625E+00,
     8       -1.12070026162230E+01,        8.78912353515625E+00,
     9       -2.36408691406250E+00,        1.12152099609375E-01,
     1       -2.82120725582002E+01,        8.46362176746007E+01,
     2       -9.18182415432400E+01,        4.25349987453885E+01,
     3       -7.36879435947963E+00,        2.27108001708984E-01,
     4        2.12570130039217E+02,       -7.65252468141182E+02,
     5        1.05999045252800E+03,       -6.99579627376133E+02/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3        2.18190511744212E+02,       -2.64914304869516E+01,
     4        5.72501420974731E-01,       -1.91945766231841E+03,
     5        8.06172218173731E+03,       -1.35865500064341E+04,
     6        1.16553933368645E+04,       -5.30564697861340E+03,
     7        1.20090291321635E+03,       -1.08090919788395E+02,
     8        1.72772750258446E+00,        2.02042913309661E+04,
     9       -9.69805983886375E+04,        1.92547001232532E+05,
     1       -2.03400177280416E+05,        1.22200464983017E+05,
     2       -4.11926549688976E+04,        7.10951430248936E+03,
     3       -4.93915304773088E+02,        6.07404200127348E+00,
     4       -2.42919187900551E+05,        1.31176361466298E+06,
     5       -2.99801591853811E+06,        3.76327129765640E+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65)/
     3       -2.81356322658653E+06,        1.26836527332162E+06,
     4       -3.31645172484564E+05,        4.52187689813627E+04,
     5       -2.49983048181121E+03,        2.43805296995561E+01,
     6        3.28446985307204E+06,       -1.97068191184322E+07,
     7        5.09526024926646E+07,       -7.41051482115327E+07,
     8        6.63445122747290E+07,       -3.75671766607634E+07,
     9        1.32887671664218E+07,       -2.78561812808645E+06,
     1        3.08186404612662E+05,       -1.38860897537170E+04,
     2        1.10017140269247E+02/
C***FIRST EXECUTABLE STATEMENT  ASYIK
      TOL = R1MACH(3)
      TOL = MAX(TOL,1.0E-15)
      FN = FNU
      Z  = (3.0E0-FLGIK)/2.0E0
      KK = INT(Z)
      DO 50 JN=1,IN
        IF (JN.EQ.1) GO TO 10
        FN = FN - FLGIK
        Z = X/FN
        RA = SQRT(1.0E0+Z*Z)
        GLN = LOG((1.0E0+RA)/Z)
        ETX = KODE - 1
        T = RA*(1.0E0-ETX) + ETX/(Z+RA)
        ARG = FN*(T-GLN)*FLGIK
   10   COEF = EXP(ARG)
        T = 1.0E0/RA
        T2 = T*T
        T = T/FN
        T = SIGN(T,FLGIK)
        S2 = 1.0E0
        AP = 1.0E0
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
          IF (MAX(ABS(AK),ABS(AP)) .LT. TOL) GO TO 40
   30   CONTINUE
   40   CONTINUE
      T = ABS(T)
      Y(JN) = S2*COEF*SQRT(T)*CON(KK)
   50 CONTINUE
      RETURN
      END
