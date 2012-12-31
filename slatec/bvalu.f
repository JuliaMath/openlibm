*DECK BVALU
      FUNCTION BVALU (T, A, N, K, IDERIV, X, INBV, WORK)
C***BEGIN PROLOGUE  BVALU
C***PURPOSE  Evaluate the B-representation of a B-spline at X for the
C            function value or any of its derivatives.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (BVALU-S, DBVALU-D)
C***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract
C         BVALU is the BVALUE function of the reference.
C
C         BVALU evaluates the B-representation (T,A,N,K) of a B-spline
C         at X for the function value on IDERIV = 0 or any of its
C         derivatives on IDERIV = 1,2,...,K-1.  Right limiting values
C         (right derivatives) are returned except at the right end
C         point X=T(N+1) where left limiting values are computed.  The
C         spline is defined on T(K) .LE. X .LE. T(N+1).  BVALU returns
C         a fatal error message when X is outside of this interval.
C
C         To compute left derivatives or left limiting values at a
C         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
C
C         BVALU calls INTRV
C
C     Description of Arguments
C         Input
C          T       - knot vector of length N+K
C          A       - B-spline coefficient vector of length N
C          N       - number of B-spline coefficients
C                    N = sum of knot multiplicities-K
C          K       - order of the B-spline, K .GE. 1
C          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
C                    IDERIV=0 returns the B-spline value
C          X       - argument, T(K) .LE. X .LE. T(N+1)
C          INBV    - an initialization parameter which must be set
C                    to 1 the first time BVALU is called.
C
C         Output
C          INBV    - INBV contains information for efficient process-
C                    ing after the initial call and INBV must not
C                    be changed by the user.  Distinct splines require
C                    distinct INBV parameters.
C          WORK    - work vector of length 3*K.
C          BVALU   - value of the IDERIV-th derivative at X
C
C     Error Conditions
C         An improper input is a fatal error
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  INTRV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BVALU
C
      INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ,
     1 IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
      REAL A, FKMJ, T, WORK, X
C     DIMENSION T(N+K), WORK(3*K)
      DIMENSION T(*), A(*), WORK(*)
C***FIRST EXECUTABLE STATEMENT  BVALU
      BVALU = 0.0E0
      IF(K.LT.1) GO TO 102
      IF(N.LT.K) GO TO 101
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
      KMIDER = K - IDERIV
C
C *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
C     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K - 1
      CALL INTRV(T, N+1, X, INBV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 120
      IF (MFLAG.EQ.0) GO TO 20
      IF (X.GT.T(I)) GO TO 130
   10 IF (I.EQ.K) GO TO 140
      I = I - 1
      IF (X.EQ.T(I)) GO TO 10
C
C *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
C     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
C
   20 IMK = I - K
      DO 30 J=1,K
        IMKPJ = IMK + J
        WORK(J) = A(IMKPJ)
   30 CONTINUE
      IF (IDERIV.EQ.0) GO TO 60
      DO 50 J=1,IDERIV
        KMJ = K - J
        FKMJ = KMJ
        DO 40 JJ=1,KMJ
          IHI = I + JJ
          IHMKMJ = IHI - KMJ
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
C
C *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
C     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60 IF (IDERIV.EQ.KM1) GO TO 100
      IP1 = I + 1
      KPK = K + K
      J1 = K + 1
      J2 = KPK + 1
      DO 70 J=1,KMIDER
        IPJ = I + J
        WORK(J1) = T(IPJ) - X
        IP1MJ = IP1 - J
        WORK(J2) = X - T(IP1MJ)
        J1 = J1 + 1
        J2 = J2 + 1
   70 CONTINUE
      IDERP1 = IDERIV + 1
      DO 90 J=IDERP1,KM1
        KMJ = K - J
        ILO = KMJ
        DO 80 JJ=1,KMJ
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ)
     1              *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
          ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
  100 BVALU = WORK(1)
      RETURN
C
C
  101 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU', 'N DOES NOT SATISFY N.GE.K', 2,
     +   1)
      RETURN
  102 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU',
     +   'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', 2, 1)
      RETURN
  120 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU',
     +   'X IS N0T GREATER THAN OR EQUAL TO T(K)', 2, 1)
      RETURN
  130 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU',
     +   'X IS NOT LESS THAN OR EQUAL TO T(N+1)', 2, 1)
      RETURN
  140 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU',
     +   'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)
      RETURN
      END
