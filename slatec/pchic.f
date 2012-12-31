*DECK PCHIC
      SUBROUTINE PCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK,
     +   IERR)
C***BEGIN PROLOGUE  PCHIC
C***PURPOSE  Set derivatives needed to determine a piecewise monotone
C            piecewise cubic Hermite interpolant to given data.
C            User control is available over boundary conditions and/or
C            treatment of points where monotonicity switches direction.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E1A
C***TYPE      SINGLE PRECISION (PCHIC-S, DPCHIC-D)
C***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
C             PCHIP, PIECEWISE CUBIC INTERPOLATION,
C             SHAPE-PRESERVING INTERPOLATION
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C         PCHIC:  Piecewise Cubic Hermite Interpolation Coefficients.
C
C     Sets derivatives needed to determine a piecewise monotone piece-
C     wise cubic interpolant to the data given in X and F satisfying the
C     boundary conditions specified by IC and VC.
C
C     The treatment of points where monotonicity switches direction is
C     controlled by argument SWITCH.
C
C     To facilitate two-dimensional applications, includes an increment
C     between successive values of the F- and D-arrays.
C
C     The resulting piecewise cubic Hermite function may be evaluated
C     by PCHFE or PCHFD.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  IC(2), N, NWK, IERR
C        REAL  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
C
C        CALL  PCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR)
C
C   Parameters:
C
C     IC -- (input) integer array of length 2 specifying desired
C           boundary conditions:
C           IC(1) = IBEG, desired condition at beginning of data.
C           IC(2) = IEND, desired condition at end of data.
C
C           IBEG = 0  for the default boundary condition (the same as
C                     used by PCHIM).
C           If IBEG.NE.0, then its sign indicates whether the boundary
C                     derivative is to be adjusted, if necessary, to be
C                     compatible with monotonicity:
C              IBEG.GT.0  if no adjustment is to be performed.
C              IBEG.LT.0  if the derivative is to be adjusted for
C                     monotonicity.
C
C           Allowable values for the magnitude of IBEG are:
C           IBEG = 1  if first derivative at X(1) is given in VC(1).
C           IBEG = 2  if second derivative at X(1) is given in VC(1).
C           IBEG = 3  to use the 3-point difference formula for D(1).
C                     (Reverts to the default b.c. if N.LT.3 .)
C           IBEG = 4  to use the 4-point difference formula for D(1).
C                     (Reverts to the default b.c. if N.LT.4 .)
C           IBEG = 5  to set D(1) so that the second derivative is con-
C              tinuous at X(2). (Reverts to the default b.c. if N.LT.4.)
C              This option is somewhat analogous to the "not a knot"
C              boundary condition provided by PCHSP.
C
C          NOTES (IBEG):
C           1. An error return is taken if ABS(IBEG).GT.5 .
C           2. Only in case  IBEG.LE.0  is it guaranteed that the
C              interpolant will be monotonic in the first interval.
C              If the returned value of D(1) lies between zero and
C              3*SLOPE(1), the interpolant will be monotonic.  This
C              is **NOT** checked if IBEG.GT.0 .
C           3. If IBEG.LT.0 and D(1) had to be changed to achieve mono-
C              tonicity, a warning error is returned.
C
C           IEND may take on the same values as IBEG, but applied to
C           derivative at X(N).  In case IEND = 1 or 2, the value is
C           given in VC(2).
C
C          NOTES (IEND):
C           1. An error return is taken if ABS(IEND).GT.5 .
C           2. Only in case  IEND.LE.0  is it guaranteed that the
C              interpolant will be monotonic in the last interval.
C              If the returned value of D(1+(N-1)*INCFD) lies between
C              zero and 3*SLOPE(N-1), the interpolant will be monotonic.
C              This is **NOT** checked if IEND.GT.0 .
C           3. If IEND.LT.0 and D(1+(N-1)*INCFD) had to be changed to
C              achieve monotonicity, a warning error is returned.
C
C     VC -- (input) real array of length 2 specifying desired boundary
C           values, as indicated above.
C           VC(1) need be set only if IC(1) = 1 or 2 .
C           VC(2) need be set only if IC(2) = 1 or 2 .
C
C     SWITCH -- (input) indicates desired treatment of points where
C           direction of monotonicity switches:
C           Set SWITCH to zero if interpolant is required to be mono-
C           tonic in each interval, regardless of monotonicity of data.
C             NOTES:
C              1. This will cause D to be set to zero at all switch
C                 points, thus forcing extrema there.
C              2. The result of using this option with the default boun-
C                 dary conditions will be identical to using PCHIM, but
C                 will generally cost more compute time.
C                 This option is provided only to facilitate comparison
C                 of different switch and/or boundary conditions.
C           Set SWITCH nonzero to use a formula based on the 3-point
C              difference formula in the vicinity of switch points.
C           If SWITCH is positive, the interpolant on each interval
C              containing an extremum is controlled to not deviate from
C              the data by more than SWITCH*DFLOC, where DFLOC is the
C              maximum of the change of F on this interval and its two
C              immediate neighbors.
C           If SWITCH is negative, no such control is to be imposed.
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of dependent variable values to be inter-
C           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
C
C     D -- (output) real array of derivative values at the data points.
C           These values will determine a monotone cubic Hermite func-
C           tion on each subinterval on which the data are monotonic,
C           except possibly adjacent to switches in monotonicity.
C           The value corresponding to X(I) is stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C           No other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in F and D.
C           This argument is provided primarily for 2-D applications.
C           (Error return if  INCFD.LT.1 .)
C
C     WK -- (scratch) real array of working storage.  The user may wish
C           to know that the returned values are:
C              WK(I)     = H(I)     = X(I+1) - X(I) ;
C              WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)
C           for  I = 1(1)N-1.
C
C     NWK -- (input) length of work array.
C           (Error return if  NWK.LT.2*(N-1) .)
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning errors:
C              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
C                        monotonicity.
C              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
C                        adjusted for monotonicity.
C              IERR = 3  if both of the above are true.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if ABS(IBEG).GT.5 .
C              IERR = -5  if ABS(IEND).GT.5 .
C              IERR = -6  if both of the above are true.
C              IERR = -7  if NWK.LT.2*(N-1) .
C             (The D-array has not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  1. F. N. Fritsch, Piecewise Cubic Hermite Interpolation
C                 Package, Report UCRL-87285, Lawrence Livermore Nation-
C                 al Laboratory, July 1982.  [Poster presented at the
C                 SIAM 30th Anniversary Meeting, 19-23 July 1982.]
C               2. F. N. Fritsch and J. Butland, A method for construc-
C                 ting local monotone piecewise cubic interpolants, SIAM
C                 Journal on Scientific and Statistical Computing 5, 2
C                 (June 1984), pp. 300-304.
C               3. F. N. Fritsch and R. E. Carlson, Monotone piecewise
C                 cubic interpolation, SIAM Journal on Numerical Ana-
C                 lysis 17, 2 (April 1980), pp. 238-246.
C***ROUTINES CALLED  PCHCE, PCHCI, PCHCS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820218  DATE WRITTEN
C   820804  Converted to SLATEC library version.
C   870813  Updated Reference 2.
C   890411  Added SAVE statements (Vers. 3.2).
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890703  Corrected category record.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920429  Revised format and order of references.  (WRB,FNF)
C***END PROLOGUE  PCHIC
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change PCHIC to DPCHIC wherever it occurs,
C        b. Change PCHCE to DPCHCE wherever it occurs,
C        c. Change PCHCI to DPCHCI wherever it occurs,
C        d. Change PCHCS to DPCHCS wherever it occurs,
C        e. Change the real declarations to double precision, and
C        f. Change the constant  ZERO  to double precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  IC(2), N, INCFD, NWK, IERR
      REAL  VC(2), SWITCH, X(*), F(INCFD,*), D(INCFD,*), WK(NWK)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IBEG, IEND, NLESS1
      REAL  ZERO
      SAVE ZERO
      DATA  ZERO /0./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHIC
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
      IF (ABS(IBEG) .GT. 5)  IERR = IERR - 1
      IF (ABS(IEND) .GT. 5)  IERR = IERR - 2
      IF (IERR .LT. 0)  GO TO 5004
C
C  FUNCTION DEFINITION IS OK -- GO ON.
C
      NLESS1 = N - 1
      IF ( NWK .LT. 2*NLESS1 )  GO TO 5007
C
C  SET UP H AND SLOPE ARRAYS.
C
      DO 20  I = 1, NLESS1
         WK(I) = X(I+1) - X(I)
         WK(NLESS1+I) = (F(1,I+1) - F(1,I)) / WK(I)
   20 CONTINUE
C
C  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
C
      IF (NLESS1 .GT. 1)  GO TO 1000
      D(1,1) = WK(2)
      D(1,N) = WK(2)
      GO TO 3000
C
C  NORMAL CASE  (N .GE. 3) .
C
 1000 CONTINUE
C
C  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.
C
C     --------------------------------------
      CALL PCHCI (N, WK(1), WK(N), D, INCFD)
C     --------------------------------------
C
C  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.
C
      IF (SWITCH .EQ. ZERO)  GO TO 3000
C     ----------------------------------------------------
      CALL PCHCS (SWITCH, N, WK(1), WK(N), D, INCFD, IERR)
C     ----------------------------------------------------
      IF (IERR .NE. 0)  GO TO 5008
C
C  SET END CONDITIONS.
C
 3000 CONTINUE
      IF ( (IBEG.EQ.0) .AND. (IEND.EQ.0) )  GO TO 5000
C     -------------------------------------------------------
      CALL PCHCE (IC, VC, N, X, WK(1), WK(N), D, INCFD, IERR)
C     -------------------------------------------------------
      IF (IERR .LT. 0)  GO TO 5009
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'PCHIC',
     +   'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'PCHIC', 'INCREMENT LESS THAN ONE', IERR,
     +   1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERMSG ('SLATEC', 'PCHIC', 'X-ARRAY NOT STRICTLY INCREASING'
     +   , IERR, 1)
      RETURN
C
 5004 CONTINUE
C     IC OUT OF RANGE RETURN.
      IERR = IERR - 3
      CALL XERMSG ('SLATEC', 'PCHIC', 'IC OUT OF RANGE', IERR, 1)
      RETURN
C
 5007 CONTINUE
C     NWK .LT. 2*(N-1)  RETURN.
      IERR = -7
      CALL XERMSG ('SLATEC', 'PCHIC', 'WORK ARRAY TOO SMALL', IERR, 1)
      RETURN
C
 5008 CONTINUE
C     ERROR RETURN FROM PCHCS.
      IERR = -8
      CALL XERMSG ('SLATEC', 'PCHIC', 'ERROR RETURN FROM PCHCS', IERR,
     +   1)
      RETURN
C
 5009 CONTINUE
C     ERROR RETURN FROM PCHCE.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
      CALL XERMSG ('SLATEC', 'PCHIC', 'ERROR RETURN FROM PCHCE', IERR,
     +   1)
      RETURN
C------------- LAST LINE OF PCHIC FOLLOWS ------------------------------
      END
