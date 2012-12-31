*DECK SPENC
      FUNCTION SPENC (X)
C***BEGIN PROLOGUE  SPENC
C***PURPOSE  Compute a form of Spence's integral due to K. Mitchell.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C5
C***TYPE      SINGLE PRECISION (SPENC-S, DSPENC-D)
C***KEYWORDS  FNLIB, SPECIAL FUNCTIONS, SPENCE'S INTEGRAL
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate a form of Spence's function defined by
C        integral from 0 to X of  -LOG(1-Y)/Y  DY.
C For ABS(X) .LE. 1, the uniformly convergent expansion
C        SPENC = sum K=1,infinity  X**K / K**2     is valid.
C
C Spence's function can be used to evaluate much more general integral
C forms.  For example,
C        integral from 0 to Z of  LOG(A*X+B)/(C*X+D)  DX  =
C             LOG(ABS(B-A*D/C))*LOG(ABS(A*(C*X+D)/(A*D-B*C)))/C
C             - SPENC (A*(C*Z+D)/(A*D-B*C)) / C.
C
C Ref -- K. Mitchell, Philosophical Magazine, 40, p. 351 (1949).
C        Stegun and Abromowitz, AMS 55, p. 1004.
C
C
C Series for SPEN       on the interval  0.          to  5.00000D-01
C                                        with weighted error   6.82E-17
C                                         log weighted error  16.17
C                               significant figures required  15.22
C                                    decimal places required  16.81
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   780201  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  SPENC
      DIMENSION SPENCS(19)
      LOGICAL FIRST
      SAVE SPENCS, PI26, NSPENC, XBIG, FIRST
      DATA SPENCS( 1) /    .1527365598 892406E0 /
      DATA SPENCS( 2) /    .0816965805 8051014E0 /
      DATA SPENCS( 3) /    .0058141571 4077873E0 /
      DATA SPENCS( 4) /    .0005371619 8145415E0 /
      DATA SPENCS( 5) /    .0000572470 4675185E0 /
      DATA SPENCS( 6) /    .0000066745 4612164E0 /
      DATA SPENCS( 7) /    .0000008276 4673397E0 /
      DATA SPENCS( 8) /    .0000001073 3156730E0 /
      DATA SPENCS( 9) /    .0000000144 0077294E0 /
      DATA SPENCS(10) /    .0000000019 8444202E0 /
      DATA SPENCS(11) /    .0000000002 7940058E0 /
      DATA SPENCS(12) /    .0000000000 4003991E0 /
      DATA SPENCS(13) /    .0000000000 0582346E0 /
      DATA SPENCS(14) /    .0000000000 0085767E0 /
      DATA SPENCS(15) /    .0000000000 0012768E0 /
      DATA SPENCS(16) /    .0000000000 0001918E0 /
      DATA SPENCS(17) /    .0000000000 0000290E0 /
      DATA SPENCS(18) /    .0000000000 0000044E0 /
      DATA SPENCS(19) /    .0000000000 0000006E0 /
      DATA PI26 / 1.644934066 848226E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  SPENC
      IF (FIRST) THEN
         NSPENC = INITS (SPENCS, 19, 0.1*R1MACH(3))
         XBIG = 1.0/R1MACH(3)
      ENDIF
      FIRST = .FALSE.
C
      IF (X.GT.2.0) GO TO 60
      IF (X.GT.1.0) GO TO 50
      IF (X.GT.0.5) GO TO 40
      IF (X.GE.0.0) GO TO 30
      IF (X.GT.(-1.)) GO TO 20
C
C HERE IF X .LE. -1.0
C
      ALN = LOG(1.0-X)
      SPENC = -PI26 - 0.5*ALN*(2.0*LOG(-X)-ALN)
      IF (X.GT.(-XBIG)) SPENC = SPENC
     1  + (1.0 + CSEVL (4.0/(1.0-X)-1.0, SPENCS, NSPENC)) / (1.0-X)
      RETURN
C
C -1.0 .LT. X .LT. 0.0
C
 20   SPENC = -0.5*LOG(1.0-X)**2
     1  - X*(1.0 + CSEVL (4.0*X/(X-1.0)-1.0, SPENCS, NSPENC)) / (X-1.0)
      RETURN
C
C 0.0 .LE. X .LE. 0.5
C
 30   SPENC = X*(1.0 + CSEVL (4.0*X-1.0, SPENCS, NSPENC))
      RETURN
C
C 0.5 .LT. X .LE. 1.0
C
 40   SPENC = PI26
      IF (X.NE.1.0) SPENC = PI26 - LOG(X)*LOG(1.0-X)
     1  - (1.0-X)*(1.0 + CSEVL (4.0*(1.0-X)-1.0, SPENCS, NSPENC))
      RETURN
C
C 1.0 .LT. X .LE. 2.0
C
 50   SPENC = PI26 - 0.5*LOG(X)*LOG((X-1.0)**2/X)
     1  + (X-1.)*(1.0 + CSEVL (4.0*(X-1.)/X-1.0, SPENCS, NSPENC))/X
      RETURN
C
C X .GT. 2.0
C
 60   SPENC = 2.0*PI26 - 0.5*LOG(X)**2
      IF (X.LT.XBIG) SPENC = SPENC
     1  - (1.0 + CSEVL (4.0/X-1.0, SPENCS, NSPENC))/X
      RETURN
C
      END
