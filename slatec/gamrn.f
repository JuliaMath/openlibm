*DECK GAMRN
      REAL FUNCTION GAMRN (X)
C***BEGIN PROLOGUE  GAMRN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BSKIN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (GAMRN-S, DGAMRN-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C         GAMRN computes the GAMMA function ratio GAMMA(X)/GAMMA(X+0.5)
C         for real X.gt.0. If X.ge.XMIN, an asymptotic expansion is
C         evaluated. If X.lt.XMIN, an integer is added to X to form a
C         new value of X.ge.XMIN and the asymptotic expansion is eval-
C         uated for this new value of X. Successive application of the
C         recurrence relation
C
C                      W(X)=W(X+1)*(1+0.5/X)
C
C         reduces the argument to its original value. XMIN and comp-
C         utational tolerances are computed as a function of the number
C         of digits carried in a word by calls to I1MACH and R1MACH.
C         However, the computational accuracy is limited to the max-
C         imum of unit roundoff (=R1MACH(4)) and 1.0E-18 since critical
C         constants are given to only 18 digits.
C
C         Input
C           X      - Argument, X.gt.0.0
C
C         OUTPUT
C           GAMRN  - Ratio  GAMMA(X)/GAMMA(X+0.5)
C
C***SEE ALSO  BSKIN
C***REFERENCES  Y. L. Luke, The Special Functions and Their
C                 Approximations, Vol. 1, Math In Sci. And
C                 Eng. Series 53, Academic Press, New York, 1969,
C                 pp. 34-35.
C***ROUTINES CALLED  I1MACH, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C   920520  Added REFERENCES section.  (WRB)
C***END PROLOGUE  GAMRN
      INTEGER I, I1M11, K, MX, NX
      INTEGER I1MACH
      REAL FLN, GR, RLN, S, TOL, TRM, X, XDMY, XINC, XM, XMIN, XP, XSQ
      REAL R1MACH
      DIMENSION GR(12)
      SAVE GR
C
      DATA GR(1), GR(2), GR(3), GR(4), GR(5), GR(6), GR(7), GR(8),
     * GR(9), GR(10), GR(11), GR(12) /1.00000000000000000E+00,
     * -1.56250000000000000E-02,2.56347656250000000E-03,
     * -1.27983093261718750E-03,1.34351104497909546E-03,
     * -2.43289663922041655E-03,6.75423753364157164E-03,
     * -2.66369606131178216E-02,1.41527455519564332E-01,
     * -9.74384543032201613E-01,8.43686251229783675E+00,
     * -8.97258321640552515E+01/
C
C***FIRST EXECUTABLE STATEMENT  GAMRN
      NX = INT(X)
      TOL = MAX(R1MACH(4),1.0E-18)
      I1M11 = I1MACH(11)
      RLN = R1MACH(5)*I1M11
      FLN = MIN(RLN,20.0E0)
      FLN = MAX(FLN,3.0E0)
      FLN = FLN - 3.0E0
      XM = 2.0E0 + FLN*(0.2366E0+0.01723E0*FLN)
      MX = INT(XM) + 1
      XMIN = MX
      XDMY = X - 0.25E0
      XINC = 0.0E0
      IF (X.GE.XMIN) GO TO 10
      XINC = XMIN - NX
      XDMY = XDMY + XINC
   10 CONTINUE
      S = 1.0E0
      IF (XDMY*TOL.GT.1.0E0) GO TO 30
      XSQ = 1.0E0/(XDMY*XDMY)
      XP = XSQ
      DO 20 K=2,12
        TRM = GR(K)*XP
        IF (ABS(TRM).LT.TOL) GO TO 30
        S = S + TRM
        XP = XP*XSQ
   20 CONTINUE
   30 CONTINUE
      S = S/SQRT(XDMY)
      IF (XINC.NE.0.0E0) GO TO 40
      GAMRN = S
      RETURN
   40 CONTINUE
      NX = INT(XINC)
      XP = 0.0E0
      DO 50 I=1,NX
        S = S*(1.0E0+0.5E0/(X+XP))
        XP = XP + 1.0E0
   50 CONTINUE
      GAMRN = S
      RETURN
      END
