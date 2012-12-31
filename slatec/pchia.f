*DECK PCHIA
      REAL FUNCTION PCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
C***BEGIN PROLOGUE  PCHIA
C***PURPOSE  Evaluate the definite integral of a piecewise cubic
C            Hermite function over an arbitrary interval.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E3, H2A1B2
C***TYPE      SINGLE PRECISION (PCHIA-S, DPCHIA-D)
C***KEYWORDS  CUBIC HERMITE INTERPOLATION, NUMERICAL INTEGRATION, PCHIP,
C             QUADRATURE
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C          PCHIA:  Piecewise Cubic Hermite Integrator, Arbitrary limits
C
C     Evaluates the definite integral of the cubic Hermite function
C     defined by  N, X, F, D  over the interval [A, B].
C
C     To provide compatibility with PCHIM and PCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IERR
C        REAL  X(N), F(INCFD,N), D(INCFD,N), A, B
C        REAL  VALUE, PCHIA
C        LOGICAL  SKIP
C
C        VALUE = PCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
C
C   Parameters:
C
C     VALUE -- (output) value of the requested integral.
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
C           SKIP will be set to .TRUE. on return with IERR.GE.0 .
C
C     A,B -- (input) the limits of integration.
C           NOTE:  There is no requirement that [A,B] be contained in
C                  [X(1),X(N)].  However, the resulting integral value
C                  will be highly suspect, if not.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning errors:
C              IERR = 1  if  A  is outside the interval [X(1),X(N)].
C              IERR = 2  if  B  is outside the interval [X(1),X(N)].
C              IERR = 3  if both of the above are true.  (Note that this
C                        means that either [A,B] contains data interval
C                        or the intervals do not intersect at all.)
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C                (VALUE will be zero in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C              IERR = -4  in case of an error return from PCHID (which
C                         should never occur).
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CHFIE, PCHID, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820730  DATE WRITTEN
C   820804  Converted to SLATEC library version.
C   870707  Corrected double precision conversion instructions.
C   870813  Minor cosmetic changes.
C   890411  Added SAVE statements (Vers. 3.2).
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890703  Corrected category record.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   930503  Corrected to set VALUE=0 when IERR.lt.0.  (FNF)
C   930504  Changed CHFIV to CHFIE.  (FNF)
C***END PROLOGUE  PCHIA
C
C  Programming notes:
C  1. The error flag from PCHID is tested, because a logic flaw
C     could conceivably result in IERD=-4, which should be reported.
C**End
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, IERR
      REAL  X(*), F(INCFD,*), D(INCFD,*), A, B
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IA, IB, IERD, IL, IR
      REAL  VALUE, XA, XB, ZERO
      SAVE ZERO
      REAL  CHFIE, PCHID
C
C  INITIALIZE.
C
      DATA  ZERO /0./
C***FIRST EXECUTABLE STATEMENT  PCHIA
      VALUE = ZERO
C
C  VALIDITY-CHECK ARGUMENTS.
C
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
      SKIP = .TRUE.
      IERR = 0
      IF ( (A.LT.X(1)) .OR. (A.GT.X(N)) )  IERR = IERR + 1
      IF ( (B.LT.X(1)) .OR. (B.GT.X(N)) )  IERR = IERR + 2
C
C  COMPUTE INTEGRAL VALUE.
C
      IF (A .NE. B)  THEN
         XA = MIN (A, B)
         XB = MAX (A, B)
         IF (XB .LE. X(2))  THEN
C           INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
C                   --------------------------------------
            VALUE = CHFIE (X(1),X(2), F(1,1),F(1,2),
     +                                D(1,1),D(1,2), A, B)
C                   --------------------------------------
         ELSE IF (XA .GE. X(N-1))  THEN
C           INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
C                   -----------------------------------------
            VALUE = CHFIE(X(N-1),X(N), F(1,N-1),F(1,N),
     +                                 D(1,N-1),D(1,N), A, B)
C                   -----------------------------------------
         ELSE
C           'NORMAL' CASE -- XA.LT.XB, XA.LT.X(N-1), XB.GT.X(2).
C      ......LOCATE IA AND IB SUCH THAT
C               X(IA-1).LT.XA.LE.X(IA).LE.X(IB).LE.XB.LE.X(IB+1)
            IA = 1
            DO 10  I = 1, N-1
               IF (XA .GT. X(I))  IA = I + 1
   10       CONTINUE
C             IA = 1 IMPLIES XA.LT.X(1) .  OTHERWISE,
C             IA IS LARGEST INDEX SUCH THAT X(IA-1).LT.XA,.
C
            IB = N
            DO 20  I = N, IA, -1
               IF (XB .LT. X(I))  IB = I - 1
   20       CONTINUE
C             IB = N IMPLIES XB.GT.X(N) .  OTHERWISE,
C             IB IS SMALLEST INDEX SUCH THAT XB.LT.X(IB+1) .
C
C     ......COMPUTE THE INTEGRAL.
            IF (IB .LT. IA)  THEN
C              THIS MEANS IB = IA-1 AND
C                 (A,B) IS A SUBSET OF (X(IB),X(IA)).
C                      ------------------------------------------
               VALUE = CHFIE (X(IB),X(IA), F(1,IB),F(1,IA),
     +                                     D(1,IB),D(1,IA), A, B)
C                      ------------------------------------------
            ELSE
C
C              FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
C                (Case (IB .EQ. IA) is taken care of by initialization
C                 of VALUE to ZERO.)
               IF (IB .GT. IA)  THEN
C                         ---------------------------------------------
                  VALUE = PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERD)
C                         ---------------------------------------------
                  IF (IERD .LT. 0)  GO TO 5004
               ENDIF
C
C              THEN ADD ON INTEGRAL OVER (XA,X(IA)).
               IF (XA .LT. X(IA))  THEN
                  IL = MAX(1, IA-1)
                  IR = IL + 1
C                                 -------------------------------------
                  VALUE = VALUE + CHFIE (X(IL),X(IR), F(1,IL),F(1,IR),
     +                                      D(1,IL),D(1,IR), XA, X(IA))
C                                 -------------------------------------
               ENDIF
C
C              THEN ADD ON INTEGRAL OVER (X(IB),XB).
               IF (XB .GT. X(IB))  THEN
                  IR = MIN (IB+1, N)
                  IL = IR - 1
C                                 -------------------------------------
                  VALUE = VALUE + CHFIE (X(IL),X(IR), F(1,IL),F(1,IR),
     +                                      D(1,IL),D(1,IR), X(IB), XB)
C                                 -------------------------------------
               ENDIF
C
C              FINALLY, ADJUST SIGN IF NECESSARY.
               IF (A .GT. B)  VALUE = -VALUE
            ENDIF
         ENDIF
      ENDIF
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      PCHIA = VALUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'PCHIA',
     +   'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
      GO TO 5000
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'PCHIA', 'INCREMENT LESS THAN ONE', IERR,
     +   1)
      GO TO 5000
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERMSG ('SLATEC', 'PCHIA',
     +   'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
      GO TO 5000
C
 5004 CONTINUE
C     TROUBLE IN PCHID.  (SHOULD NEVER OCCUR.)
      IERR = -4
      CALL XERMSG ('SLATEC', 'PCHIA', 'TROUBLE IN PCHID', IERR, 1)
      GO TO 5000
C------------- LAST LINE OF PCHIA FOLLOWS ------------------------------
      END
