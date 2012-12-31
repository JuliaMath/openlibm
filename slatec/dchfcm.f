*DECK DCHFCM
      INTEGER FUNCTION DCHFCM (D1, D2, DELTA)
C***BEGIN PROLOGUE  DCHFCM
C***SUBSIDIARY
C***PURPOSE  Check a single cubic for monotonicity.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      DOUBLE PRECISION (CHFCM-S, DCHFCM-D)
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C *Usage:
C
C        DOUBLE PRECISION  D1, D2, DELTA
C        INTEGER  ISMON, DCHFCM
C
C        ISMON = DCHFCM (D1, D2, DELTA)
C
C *Arguments:
C
C     D1,D2:IN  are the derivative values at the ends of an interval.
C
C     DELTA:IN  is the data slope over that interval.
C
C *Function Return Values:
C     ISMON : indicates the monotonicity of the cubic segment:
C             ISMON = -3  if function is probably decreasing;
C             ISMON = -1  if function is strictly decreasing;
C             ISMON =  0  if function is constant;
C             ISMON =  1  if function is strictly increasing;
C             ISMON =  2  if function is non-monotonic;
C             ISMON =  3  if function is probably increasing.
C           If ABS(ISMON)=3, the derivative values are too close to the
C           boundary of the monotonicity region to declare monotonicity
C           in the presence of roundoff error.
C
C *Description:
C
C          DCHFCM:  Cubic Hermite Function -- Check Monotonicity.
C
C    Called by  DPCHCM  to determine the monotonicity properties of the
C    cubic with boundary derivative values D1,D2 and chord slope DELTA.
C
C *Cautions:
C     This is essentially the same as old DCHFMC, except that a
C     new output value, -3, was added February 1989.  (Formerly, -3
C     and +3 were lumped together in the single value 3.)  Codes that
C     flag nonmonotonicity by "IF (ISMON.EQ.2)" need not be changed.
C     Codes that check via "IF (ISMON.GE.3)" should change the test to
C     "IF (IABS(ISMON).GE.3)".  Codes that declare monotonicity via
C     "IF (ISMON.LE.1)" should change to "IF (IABS(ISMON).LE.1)".
C
C   REFER TO  DPCHCM
C
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   820518  DATE WRITTEN
C   820805  Converted to SLATEC library version.
C   831201  Changed from  ISIGN  to SIGN  to correct bug that
C           produced wrong sign when -1 .LT. DELTA .LT. 0 .
C   890206  Added SAVE statements.
C   890209  Added sign to returned value ISMON=3 and corrected
C           argument description accordingly.
C   890306  Added caution about changed output.
C   890407  Changed name from DCHFMC to DCHFCM, as requested at the
C           March 1989 SLATEC CML meeting, and made a few other
C           minor modifications necessitated by this change.
C   890407  Converted to new SLATEC format.
C   890407  Modified DESCRIPTION to LDOC format.
C   891214  Moved SAVE statements.  (WRB)
C***END PROLOGUE  DCHFCM
C
C  Fortran intrinsics used:  DSIGN.
C  Other routines used:  D1MACH.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     TEN is actually a tuning parameter, which determines the width of
C     the fuzz around the elliptical boundary.
C
C     To produce a single precision version, simply:
C        a. Change DCHFCM to CHFCM wherever it occurs,
C        b. Change the double precision declarations to real, and
C        c. Change the constants ZERO, ONE, ... to single precision.
C
C  DECLARE ARGUMENTS.
C
      DOUBLE PRECISION  D1, D2, DELTA, D1MACH
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER ISMON, ITRUE
      DOUBLE PRECISION  A, B, EPS, FOUR, ONE, PHI, TEN, THREE, TWO,
     * ZERO
      SAVE ZERO, ONE, TWO, THREE, FOUR
      SAVE TEN
C
C  INITIALIZE.
C
      DATA ZERO /0.D0/, ONE/1.D0/, TWO/2.D0/, THREE/3.D0/, FOUR/4.D0/,
     1      TEN /10.D0/
C
C        MACHINE-DEPENDENT PARAMETER -- SHOULD BE ABOUT 10*UROUND.
C***FIRST EXECUTABLE STATEMENT  DCHFCM
      EPS = TEN*D1MACH(4)
C
C  MAKE THE CHECK.
C
      IF (DELTA .EQ. ZERO)  THEN
C        CASE OF CONSTANT DATA.
         IF ((D1.EQ.ZERO) .AND. (D2.EQ.ZERO))  THEN
            ISMON = 0
         ELSE
            ISMON = 2
         ENDIF
      ELSE
C        DATA IS NOT CONSTANT -- PICK UP SIGN.
         ITRUE = DSIGN (ONE, DELTA)
         A = D1/DELTA
         B = D2/DELTA
         IF ((A.LT.ZERO) .OR. (B.LT.ZERO))  THEN
            ISMON = 2
         ELSE IF ((A.LE.THREE-EPS) .AND. (B.LE.THREE-EPS))  THEN
C           INSIDE SQUARE (0,3)X(0,3)  IMPLIES   OK.
            ISMON = ITRUE
         ELSE IF ((A.GT.FOUR+EPS) .AND. (B.GT.FOUR+EPS))  THEN
C           OUTSIDE SQUARE (0,4)X(0,4)  IMPLIES   NONMONOTONIC.
            ISMON = 2
         ELSE
C           MUST CHECK AGAINST BOUNDARY OF ELLIPSE.
            A = A - TWO
            B = B - TWO
            PHI = ((A*A + B*B) + A*B) - THREE
            IF (PHI .LT. -EPS)  THEN
               ISMON = ITRUE
            ELSE IF (PHI .GT. EPS)  THEN
               ISMON = 2
            ELSE
C              TO CLOSE TO BOUNDARY TO TELL,
C                  IN THE PRESENCE OF ROUND-OFF ERRORS.
               ISMON = 3*ITRUE
            ENDIF
         ENDIF
      ENDIF
C
C  RETURN VALUE.
C
      DCHFCM = ISMON
      RETURN
C------------- LAST LINE OF DCHFCM FOLLOWS -----------------------------
      END
