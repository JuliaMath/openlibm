*DECK DPCHID
      DOUBLE PRECISION FUNCTION DPCHID (N, X, F, D, INCFD, SKIP, IA, IB,
     +   IERR)
C***BEGIN PROLOGUE  DPCHID
C***PURPOSE  Evaluate the definite integral of a piecewise cubic
C            Hermite function over an interval whose endpoints are data
C            points.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E3, H2A1B2
C***TYPE      DOUBLE PRECISION (PCHID-S, DPCHID-D)
C***KEYWORDS  CUBIC HERMITE INTERPOLATION, NUMERICAL INTEGRATION, PCHIP,
C             QUADRATURE
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C          DPCHID:  Piecewise Cubic Hermite Integrator, Data limits
C
C     Evaluates the definite integral of the cubic Hermite function
C     defined by  N, X, F, D  over the interval [X(IA), X(IB)].
C
C     To provide compatibility with DPCHIM and DPCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IA, IB, IERR
C        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
C        LOGICAL  SKIP
C
C        VALUE = DPCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
C
C   Parameters:
C
C     VALUE -- (output) value of the requested integral.
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real*8 array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
C           is the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in DPCHIM or DPCHIC).
C           SKIP will be set to .TRUE. on return with IERR = 0 or -4.
C
C     IA,IB -- (input) indices in X-array for the limits of integration.
C           both must be in the range [1,N].  (Error return if not.)
C           No restrictions on their relative values.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if IA or IB is out of range.
C                (VALUE will be zero in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820723  DATE WRITTEN
C   820804  Converted to SLATEC library version.
C   870707  Corrected XERROR calls for d.p. name(s).
C   870813  Minor cosmetic changes.
C   890206  Corrected XERROR calls.
C   890411  Added SAVE statements (Vers. 3.2).
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890703  Corrected category record.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   930504  Corrected to set VALUE=0 when IERR.ne.0.  (FNF)
C***END PROLOGUE  DPCHID
C
C  Programming notes:
C  1. This routine uses a special formula that is valid only for
C     integrals whose limits coincide with data values.  This is
C     mathematically equivalent to, but much more efficient than,
C     calls to DCHFIE.
C**End
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, IA, IB, IERR
      DOUBLE PRECISION  X(*), F(INCFD,*), D(INCFD,*)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IUP, LOW
      DOUBLE PRECISION  H, HALF, SIX, SUM, VALUE, ZERO
      SAVE ZERO, HALF, SIX
C
C  INITIALIZE.
C
      DATA  ZERO /0.D0/,  HALF/.5D0/, SIX/6.D0/
C***FIRST EXECUTABLE STATEMENT  DPCHID
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
      IF ((IA.LT.1) .OR. (IA.GT.N))  GO TO 5004
      IF ((IB.LT.1) .OR. (IB.GT.N))  GO TO 5004
      IERR = 0
C
C  COMPUTE INTEGRAL VALUE.
C
      IF (IA .NE. IB)  THEN
         LOW = MIN(IA, IB)
         IUP = MAX(IA, IB) - 1
         SUM = ZERO
         DO 10  I = LOW, IUP
            H = X(I+1) - X(I)
            SUM = SUM + H*( (F(1,I) + F(1,I+1)) +
     *                      (D(1,I) - D(1,I+1))*(H/SIX) )
   10    CONTINUE
         VALUE = HALF * SUM
         IF (IA .GT. IB)  VALUE = -VALUE
      ENDIF
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      DPCHID = VALUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'DPCHID',
     +   'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
      GO TO 5000
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'DPCHID', 'INCREMENT LESS THAN ONE', IERR,
     +   1)
      GO TO 5000
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERMSG ('SLATEC', 'DPCHID',
     +   'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
      GO TO 5000
C
 5004 CONTINUE
C     IA OR IB OUT OF RANGE RETURN.
      IERR = -4
      CALL XERMSG ('SLATEC', 'DPCHID', 'IA OR IB OUT OF RANGE', IERR,
     +   1)
      GO TO 5000
C------------- LAST LINE OF DPCHID FOLLOWS -----------------------------
      END
