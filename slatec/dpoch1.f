*DECK DPOCH1
      DOUBLE PRECISION FUNCTION DPOCH1 (A, X)
C***BEGIN PROLOGUE  DPOCH1
C***PURPOSE  Calculate a generalization of Pochhammer's symbol starting
C            from first order.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C1, C7A
C***TYPE      DOUBLE PRECISION (POCH1-S, DPOCH1-D)
C***KEYWORDS  FIRST ORDER, FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate a double precision generalization of Pochhammer's symbol
C for double precision A and X for special situations that require
C especially accurate values when X is small in
C        POCH1(A,X) = (POCH(A,X)-1)/X
C                   = (GAMMA(A+X)/GAMMA(A) - 1.0)/X .
C This specification is particularly suited for stably computing
C expressions such as
C        (GAMMA(A+X)/GAMMA(A) - GAMMA(B+X)/GAMMA(B))/X
C             = POCH1(A,X) - POCH1(B,X)
C Note that POCH1(A,0.0) = PSI(A)
C
C When ABS(X) is so small that substantial cancellation will occur if
C the straightforward formula is used, we use an expansion due
C to Fields and discussed by Y. L. Luke, The Special Functions and Their
C Approximations, Vol. 1, Academic Press, 1969, page 34.
C
C The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
C        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) .
C In order to maintain significance in POCH1, we write for positive a
C        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q)
C                       = 1.0 + Q*EXPREL(Q) .
C Likewise the polynomial is written
C        POLY = 1.0 + X*POLY1(A,X) .
C Thus,
C        POCH1(A,X) = (POCH(A,X) - 1) / X
C                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCOT, DEXPRL, DPOCH, DPSI, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  DPOCH1
      DOUBLE PRECISION A, X, ABSA, ABSX, ALNEPS, ALNVAR, B, BERN(20),
     1  BINV, BP, GBERN(21), GBK, PI, POLY1, Q, RHO, SINPXX, SINPX2,
     2  SQTBIG, TERM, TRIG, VAR, VAR2, D1MACH, DPSI, DEXPRL, DCOT, DPOCH
      LOGICAL FIRST
      EXTERNAL DCOT
      SAVE BERN, PI, SQTBIG, ALNEPS, FIRST
      DATA BERN  (  1) / +.8333333333 3333333333 3333333333 333 D-1    /
      DATA BERN  (  2) / -.1388888888 8888888888 8888888888 888 D-2    /
      DATA BERN  (  3) / +.3306878306 8783068783 0687830687 830 D-4    /
      DATA BERN  (  4) / -.8267195767 1957671957 6719576719 576 D-6    /
      DATA BERN  (  5) / +.2087675698 7868098979 2100903212 014 D-7    /
      DATA BERN  (  6) / -.5284190138 6874931848 4768220217 955 D-9    /
      DATA BERN  (  7) / +.1338253653 0684678832 8269809751 291 D-10   /
      DATA BERN  (  8) / -.3389680296 3225828668 3019539124 944 D-12   /
      DATA BERN  (  9) / +.8586062056 2778445641 3590545042 562 D-14   /
      DATA BERN  ( 10) / -.2174868698 5580618730 4151642386 591 D-15   /
      DATA BERN  ( 11) / +.5509002828 3602295152 0265260890 225 D-17   /
      DATA BERN  ( 12) / -.1395446468 5812523340 7076862640 635 D-18   /
      DATA BERN  ( 13) / +.3534707039 6294674716 9322997780 379 D-20   /
      DATA BERN  ( 14) / -.8953517427 0375468504 0261131811 274 D-22   /
      DATA BERN  ( 15) / +.2267952452 3376830603 1095073886 816 D-23   /
      DATA BERN  ( 16) / -.5744724395 2026452383 4847971943 400 D-24   /
      DATA BERN  ( 17) / +.1455172475 6148649018 6626486727 132 D-26   /
      DATA BERN  ( 18) / -.3685994940 6653101781 8178247990 866 D-28   /
      DATA BERN  ( 19) / +.9336734257 0950446720 3255515278 562 D-30   /
      DATA BERN  ( 20) / -.2365022415 7006299345 5963519636 983 D-31   /
      DATA PI / 3.1415926535 8979323846 2643383279 503 D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DPOCH1
      IF (FIRST) THEN
         SQTBIG = 1.0D0/SQRT(24.0D0*D1MACH(1))
         ALNEPS = LOG(D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (X.EQ.0.0D0) DPOCH1 = DPSI(A)
      IF (X.EQ.0.0D0) RETURN
C
      ABSX = ABS(X)
      ABSA = ABS(A)
      IF (ABSX.GT.0.1D0*ABSA) GO TO 70
      IF (ABSX*LOG(MAX(ABSA,2.0D0)).GT.0.1D0) GO TO 70
C
      BP = A
      IF (A.LT.(-0.5D0)) BP = 1.0D0 - A - X
      INCR = 0
      IF (BP.LT.10.0D0) INCR = 11.0D0 - BP
      B = BP + INCR
C
      VAR = B + 0.5D0*(X-1.0D0)
      ALNVAR = LOG(VAR)
      Q = X*ALNVAR
C
      POLY1 = 0.0D0
      IF (VAR.GE.SQTBIG) GO TO 40
      VAR2 = (1.0D0/VAR)**2
C
      RHO = 0.5D0*(X+1.0D0)
      GBERN(1) = 1.0D0
      GBERN(2) = -RHO/12.0D0
      TERM = VAR2
      POLY1 = GBERN(2)*TERM
C
      NTERMS = -0.5D0*ALNEPS/ALNVAR + 1.0D0
      IF (NTERMS .GT. 20) CALL XERMSG ('SLATEC', 'DPOCH1',
     +   'NTERMS IS TOO BIG, MAYBE D1MACH(3) IS BAD', 1, 2)
      IF (NTERMS.LT.2) GO TO 40
C
      DO 30 K=2,NTERMS
        GBK = 0.0D0
        DO 20 J=1,K
          NDX = K - J + 1
          GBK = GBK + BERN(NDX)*GBERN(J)
 20     CONTINUE
        GBERN(K+1) = -RHO*GBK/K
C
        TERM = TERM * (2*K-2-X)*(2*K-1-X)*VAR2
        POLY1 = POLY1 + GBERN(K+1)*TERM
 30   CONTINUE
C
 40   POLY1 = (X-1.0D0)*POLY1
      DPOCH1 = DEXPRL(Q)*(ALNVAR+Q*POLY1) + POLY1
C
      IF (INCR.EQ.0) GO TO 60
C
C WE HAVE DPOCH1(B,X), BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
C TO OBTAIN DPOCH1(BP,X).
C
      DO 50 II=1,INCR
        I = INCR - II
        BINV = 1.0D0/(BP+I)
        DPOCH1 = (DPOCH1 - BINV) / (1.0D0 + X*BINV)
 50   CONTINUE
C
 60   IF (BP.EQ.A) RETURN
C
C WE HAVE DPOCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A REFLECTION
C FORMULA TO OBTAIN DPOCH1(A,X).
C
      SINPXX = SIN(PI*X)/X
      SINPX2 = SIN(0.5D0*PI*X)
      TRIG = SINPXX*DCOT(PI*B) - 2.0D0*SINPX2*(SINPX2/X)
C
      DPOCH1 = TRIG + (1.0D0 + X*TRIG)*DPOCH1
      RETURN
C
 70   DPOCH1 = (DPOCH(A,X) - 1.0D0) / X
      RETURN
C
      END
