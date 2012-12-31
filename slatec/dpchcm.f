*DECK DPCHCM
      SUBROUTINE DPCHCM (N, X, F, D, INCFD, SKIP, ISMON, IERR)
C***BEGIN PROLOGUE  DPCHCM
C***PURPOSE  Check a cubic Hermite function for monotonicity.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E3
C***TYPE      DOUBLE PRECISION (PCHCM-S, DPCHCM-D)
C***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
C             PCHIP, PIECEWISE CUBIC INTERPOLATION, UTILITY ROUTINE
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Computing & Mathematics Research Division
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C *Usage:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, ISMON(N), IERR
C        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
C        LOGICAL  SKIP
C
C        CALL  DPCHCM (N, X, F, D, INCFD, SKIP, ISMON, IERR)
C
C *Arguments:
C
C     N:IN  is the number of data points.  (Error return if N.LT.2 .)
C
C     X:IN  is a real*8 array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F:IN  is a real*8 array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D:IN  is a real*8 array of derivative values.  D(1+(I-1)*INCFD) is
C           is the value corresponding to X(I).
C
C     INCFD:IN  is the increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP:INOUT  is a logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed.
C           SKIP will be set to .TRUE. on normal return.
C
C     ISMON:OUT  is an integer array indicating on which intervals the
C           PCH function defined by  N, X, F, D  is monotonic.
C           For data interval [X(I),X(I+1)],
C             ISMON(I) = -3  if function is probably decreasing;
C             ISMON(I) = -1  if function is strictly decreasing;
C             ISMON(I) =  0  if function is constant;
C             ISMON(I) =  1  if function is strictly increasing;
C             ISMON(I) =  2  if function is non-monotonic;
C             ISMON(I) =  3  if function is probably increasing.
C                If ABS(ISMON)=3, this means that the D-values are near
C                the boundary of the monotonicity region.  A small
C                increase produces non-monotonicity; decrease, strict
C                monotonicity.
C           The above applies to I=1(1)N-1.  ISMON(N) indicates whether
C              the entire function is monotonic on [X(1),X(N)].
C
C     IERR:OUT  is an error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C          (The ISMON-array has not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C *Description:
C
C          DPCHCM:  Piecewise Cubic Hermite -- Check Monotonicity.
C
C     Checks the piecewise cubic Hermite function defined by  N,X,F,D
C     for monotonicity.
C
C     To provide compatibility with DPCHIM and DPCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C *Cautions:
C     This provides the same capability as old DPCHMC, except that a
C     new output value, -3, was added February 1989.  (Formerly, -3
C     and +3 were lumped together in the single value 3.)  Codes that
C     flag nonmonotonicity by "IF (ISMON.EQ.2)" need not be changed.
C     Codes that check via "IF (ISMON.GE.3)" should change the test to
C     "IF (IABS(ISMON).GE.3)".  Codes that declare monotonicity via
C     "IF (ISMON.LE.1)" should change to "IF (IABS(ISMON).LE.1)".
C
C***REFERENCES  F. N. Fritsch and R. E. Carlson, Monotone piecewise
C                 cubic interpolation, SIAM Journal on Numerical Ana-
C                 lysis 17, 2 (April 1980), pp. 238-246.
C***ROUTINES CALLED  DCHFCM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820518  DATE WRITTEN
C   820804  Converted to SLATEC library version.
C   831201  Reversed order of subscripts of F and D, so that the
C           routine will work properly when INCFD.GT.1 .  (Bug!)
C   870707  Corrected XERROR calls for d.p. name(s).
C   890206  Corrected XERROR calls.
C   890209  Added possible ISMON value of -3 and modified code so
C           that 1,3,-1 produces ISMON(N)=2, rather than 3.
C   890306  Added caution about changed output.
C   890407  Changed name from DPCHMC to DPCHCM, as requested at the
C           March 1989 SLATEC CML meeting, and made a few other
C           minor modifications necessitated by this change.
C   890407  Converted to new SLATEC format.
C   890407  Modified DESCRIPTION to LDOC format.
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920429  Revised format and order of references.  (WRB,FNF)
C***END PROLOGUE  DPCHCM
C
C  Fortran intrinsics used:  ISIGN.
C  Other routines used:  CHFCM, XERMSG.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     An alternate organization would have separate loops for computing
C     ISMON(i), i=1,...,NSEG, and for the computation of ISMON(N).  The
C     first loop can be readily parallelized, since the NSEG calls to
C     CHFCM are independent.  The second loop can be cut short if
C     ISMON(N) is ever equal to 2, for it cannot be changed further.
C
C     To produce a single precision version, simply:
C        a. Change DPCHCM to PCHCM wherever it occurs,
C        b. Change DCHFCM to CHFCM wherever it occurs, and
C        c. Change the double precision declarations to real.
C
C  DECLARE ARGUMENTS.
C
      INTEGER N, INCFD, ISMON(N), IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER I, NSEG
      DOUBLE PRECISION  DELTA
      INTEGER DCHFCM
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DPCHCM
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
      SKIP = .TRUE.
C
C  FUNCTION DEFINITION IS OK -- GO ON.
C
    5 CONTINUE
      NSEG = N - 1
      DO 90  I = 1, NSEG
         DELTA = (F(1,I+1)-F(1,I))/(X(I+1)-X(I))
C                   -------------------------------
         ISMON(I) = DCHFCM (D(1,I), D(1,I+1), DELTA)
C                   -------------------------------
         IF (I .EQ. 1)  THEN
            ISMON(N) = ISMON(1)
         ELSE
C           Need to figure out cumulative monotonicity from following
C           "multiplication table":
C
C                    +        I S M O N (I)
C                     +  -3  -1   0   1   3   2
C                      +------------------------+
C               I   -3 I -3  -3  -3   2   2   2 I
C               S   -1 I -3  -1  -1   2   2   2 I
C               M    0 I -3  -1   0   1   3   2 I
C               O    1 I  2   2   1   1   3   2 I
C               N    3 I  2   2   3   3   3   2 I
C              (N)   2 I  2   2   2   2   2   2 I
C                      +------------------------+
C           Note that the 2 row and column are out of order so as not
C           to obscure the symmetry in the rest of the table.
C
C           No change needed if equal or constant on this interval or
C           already declared nonmonotonic.
            IF ( (ISMON(I).NE.ISMON(N)) .AND. (ISMON(I).NE.0)
     .                                  .AND. (ISMON(N).NE.2) )  THEN
               IF ( (ISMON(I).EQ.2) .OR. (ISMON(N).EQ.0) )  THEN
                  ISMON(N) =  ISMON(I)
               ELSE IF (ISMON(I)*ISMON(N) .LT. 0)  THEN
C                 This interval has opposite sense from curve so far.
                  ISMON(N) = 2
               ELSE
C                 At this point, both are nonzero with same sign, and
C                 we have already eliminated case both +-1.
                  ISMON(N) = ISIGN (3, ISMON(N))
               ENDIF
            ENDIF
         ENDIF
   90 CONTINUE
C
C  NORMAL RETURN.
C
      IERR = 0
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'DPCHCM',
     +   'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'DPCHCM', 'INCREMENT LESS THAN ONE', IERR,
     +   1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERMSG ('SLATEC', 'DPCHCM',
     +   'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
      RETURN
C------------- LAST LINE OF DPCHCM FOLLOWS -----------------------------
      END
