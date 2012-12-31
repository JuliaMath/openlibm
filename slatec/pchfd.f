*DECK PCHFD
      SUBROUTINE PCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
C***BEGIN PROLOGUE  PCHFD
C***PURPOSE  Evaluate a piecewise cubic Hermite function and its first
C            derivative at an array of points.  May be used by itself
C            for Hermite interpolation, or as an evaluator for PCHIM
C            or PCHIC.  If only function values are required, use
C            PCHFE instead.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E3, H1
C***TYPE      SINGLE PRECISION (PCHFD-S, DPCHFD-D)
C***KEYWORDS  CUBIC HERMITE DIFFERENTIATION, CUBIC HERMITE EVALUATION,
C             HERMITE INTERPOLATION, PCHIP, PIECEWISE CUBIC EVALUATION
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C          PCHFD:  Piecewise Cubic Hermite Function and Derivative
C                  evaluator
C
C     Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
C     gether with its first derivative, at the points  XE(J), J=1(1)NE.
C
C     If only function values are required, use PCHFE, instead.
C
C     To provide compatibility with PCHIM and PCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, NE, IERR
C        REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE), DE(NE)
C        LOGICAL  SKIP
C
C        CALL  PCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in PCHIM or PCHIC).
C           SKIP will be set to .TRUE. on normal return.
C
C     NE -- (input) number of evaluation points.  (Error return if
C           NE.LT.1 .)
C
C     XE -- (input) real array of points at which the functions are to
C           be evaluated.
C
C
C          NOTES:
C           1. The evaluation will be most efficient if the elements
C              of XE are increasing relative to X;
C              that is,   XE(J) .GE. X(I)
C              implies    XE(K) .GE. X(I),  all K.GE.J .
C           2. If any of the XE are outside the interval [X(1),X(N)],
C              values are extrapolated from the nearest extreme cubic,
C              and a warning error is returned.
C
C     FE -- (output) real array of values of the cubic Hermite function
C           defined by  N, X, F, D  at the points  XE.
C
C     DE -- (output) real array of values of the first derivative of
C           the same function at the points  XE.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  means that extrapolation was performed at
C                 IERR points.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if NE.LT.1 .
C           (Output arrays have not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C              IERR = -5  if an error has occurred in the lower-level
C                         routine CHFDV.  NB: this should never happen.
C                         Notify the author **IMMEDIATELY** if it does.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CHFDV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   811020  DATE WRITTEN
C   820803  Minor cosmetic changes for release 1.
C   870707  Minor cosmetic changes to prologue.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  PCHFD
C  Programming notes:
C
C     1. To produce a double precision version, simply:
C        a. Change PCHFD to DPCHFD, and CHFDV to DCHFDV, wherever they
C           occur,
C        b. Change the real declaration to double precision,
C
C     2. Most of the coding between the call to CHFDV and the end of
C        the IR-loop could be eliminated if it were permissible to
C        assume that XE is ordered relative to X.
C
C     3. CHFDV does not assume that X1 is less than X2.  thus, it would
C        be possible to write a version of PCHFD that assumes a strict-
C        ly decreasing X-array by simply running the IR-loop backwards
C        (and reversing the order of appropriate tests).
C
C     4. The present code has a minor bug, which I have decided is not
C        worth the effort that would be required to fix it.
C        If XE contains points in [X(N-1),X(N)], followed by points .LT.
C        X(N-1), followed by points .GT.X(N), the extrapolation points
C        will be counted (at least) twice in the total returned in IERR.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, NE, IERR
      REAL  X(*), F(INCFD,*), D(INCFD,*), XE(*), FE(*), DE(*)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IERC, IR, J, JFIRST, NEXT(2), NJ
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHFD
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      IF ( NE.LT.1 )  GO TO 5004
      IERR = 0
      SKIP = .TRUE.
C
C  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
C                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
      JFIRST = 1
      IR = 2
   10 CONTINUE
C
C     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
C
         IF (JFIRST .GT. NE)  GO TO 5000
C
C     LOCATE ALL POINTS IN INTERVAL.
C
         DO 20  J = JFIRST, NE
            IF (XE(J) .GE. X(IR))  GO TO 30
   20    CONTINUE
         J = NE + 1
         GO TO 40
C
C     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
C
   30    CONTINUE
         IF (IR .EQ. N)  J = NE + 1
C
   40    CONTINUE
         NJ = J - JFIRST
C
C     SKIP EVALUATION IF NO POINTS IN INTERVAL.
C
         IF (NJ .EQ. 0)  GO TO 50
C
C     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
C
C       ----------------------------------------------------------------
        CALL CHFDV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR),
     *              NJ, XE(JFIRST), FE(JFIRST), DE(JFIRST), NEXT, IERC)
C       ----------------------------------------------------------------
         IF (IERC .LT. 0)  GO TO 5005
C
         IF (NEXT(2) .EQ. 0)  GO TO 42
C        IF (NEXT(2) .GT. 0)  THEN
C           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
C           RIGHT OF X(IR).
C
            IF (IR .LT. N)  GO TO 41
C           IF (IR .EQ. N)  THEN
C              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(2)
               GO TO 42
   41       CONTINUE
C           ELSE
C              WE SHOULD NEVER HAVE GOTTEN HERE.
               GO TO 5005
C           ENDIF
C        ENDIF
   42    CONTINUE
C
         IF (NEXT(1) .EQ. 0)  GO TO 49
C        IF (NEXT(1) .GT. 0)  THEN
C           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
C           LEFT OF X(IR-1).
C
            IF (IR .GT. 2)  GO TO 43
C           IF (IR .EQ. 2)  THEN
C              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(1)
               GO TO 49
   43       CONTINUE
C           ELSE
C              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
C              EVALUATION INTERVAL.
C
C              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
               DO 44  I = JFIRST, J-1
                  IF (XE(I) .LT. X(IR-1))  GO TO 45
   44          CONTINUE
C              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
C                     IN CHFDV.
               GO TO 5005
C
   45          CONTINUE
C              RESET J.  (THIS WILL BE THE NEW JFIRST.)
               J = I
C
C              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
               DO 46  I = 1, IR-1
                  IF (XE(J) .LT. X(I)) GO TO 47
   46          CONTINUE
C              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
C
   47          CONTINUE
C              AT THIS POINT, EITHER  XE(J) .LT. X(1)
C                 OR      X(I-1) .LE. XE(J) .LT. X(I) .
C              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
C              CYCLING.
               IR = MAX(1, I-1)
C           ENDIF
C        ENDIF
   49    CONTINUE
C
         JFIRST = J
C
C     END OF IR-LOOP.
C
   50 CONTINUE
      IR = IR + 1
      IF (IR .LE. N)  GO TO 10
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
      CALL XERMSG ('SLATEC', 'PCHFD',
     +   'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'PCHFD', 'INCREMENT LESS THAN ONE', IERR,
     +   1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERMSG ('SLATEC', 'PCHFD', 'X-ARRAY NOT STRICTLY INCREASING'
     +   , IERR, 1)
      RETURN
C
 5004 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -4
      CALL XERMSG ('SLATEC', 'PCHFD',
     +   'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
      RETURN
C
 5005 CONTINUE
C     ERROR RETURN FROM CHFDV.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -5
      CALL XERMSG ('SLATEC', 'PCHFD',
     +   'ERROR RETURN FROM CHFDV -- FATAL', IERR, 2)
      RETURN
C------------- LAST LINE OF PCHFD FOLLOWS ------------------------------
      END
