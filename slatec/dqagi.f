*DECK DQAGI
      SUBROUTINE DQAGI (F, BOUND, INF, EPSABS, EPSREL, RESULT, ABSERR,
     +   NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
C***BEGIN PROLOGUE  DQAGI
C***PURPOSE  The routine calculates an approximation result to a given
C            INTEGRAL   I = Integral of F over (BOUND,+INFINITY)
C            OR I = Integral of F over (-INFINITY,BOUND)
C            OR I = Integral of F over (-INFINITY,+INFINITY)
C            Hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A3A1, H2A4A1
C***TYPE      DOUBLE PRECISION (QAGI-S, DQAGI-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE,
C             GLOBALLY ADAPTIVE, INFINITE INTERVALS, QUADPACK,
C             QUADRATURE, TRANSFORMATION
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Integration over infinite intervals
C        Standard fortran subroutine
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            BOUND  - Double precision
C                     Finite bound of integration range
C                     (has no meaning if interval is doubly-infinite)
C
C            INF    - Integer
C                     indicating the kind of integration range involved
C                     INF = 1 corresponds to  (BOUND,+INFINITY),
C                     INF = -1            to  (-INFINITY,BOUND),
C                     INF = 2             to (-INFINITY,+INFINITY).
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     the routine will end with IER = 6.
C
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate of the modulus of the absolute error,
C                     which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                   - IER.GT.0 abnormal termination of the routine. The
C                             estimates for result and error are less
C                             reliable. It is assumed that the requested
C                             accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account). However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties. If
C                             the position of a local difficulty can be
C                             determined (e.g. SINGULARITY,
C                             DISCONTINUITY within the interval) one
C                             will probably gain from splitting up the
C                             interval at this point and calling the
C                             integrator on the subranges. If possible,
C                             an appropriate special-purpose integrator
C                             should be used, which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                             The error may be under-estimated.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 4 The algorithm does not converge.
C                             Roundoff error is detected in the
C                             extrapolation table.
C                             It is assumed that the requested tolerance
C                             cannot be achieved, and that the returned
C                             RESULT is the best which can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                              or LIMIT.LT.1 or LENIW.LT.LIMIT*4.
C                             RESULT, ABSERR, NEVAL, LAST are set to
C                             zero.  Except when LIMIT or LENIW is
C                             invalid, IWORK(1), WORK(LIMIT*2+1) and
C                             WORK(LIMIT*3+1) are set to ZERO, WORK(1)
C                             is set to A and WORK(LIMIT+1) to B.
C
C         DIMENSIONING PARAMETERS
C            LIMIT - Integer
C                    Dimensioning parameter for IWORK
C                    LIMIT determines the maximum number of subintervals
C                    in the partition of the given integration interval
C                    (A,B), LIMIT.GE.1.
C                    If LIMIT.LT.1, the routine will end with IER = 6.
C
C            LENW  - Integer
C                    Dimensioning parameter for WORK
C                    LENW must be at least LIMIT*4.
C                    If LENW.LT.LIMIT*4, the routine will end
C                    with IER = 6.
C
C            LAST  - Integer
C                    On return, LAST equals the number of subintervals
C                    produced in the subdivision process, which
C                    determines the number of significant elements
C                    actually in the WORK ARRAYS.
C
C         WORK ARRAYS
C            IWORK - Integer
C                    Vector of dimension at least LIMIT, the first
C                    K elements of which contain pointers
C                    to the error estimates over the subintervals,
C                    such that WORK(LIMIT*3+IWORK(1)),... ,
C                    WORK(LIMIT*3+IWORK(K)) form a decreasing
C                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2), and
C                    K = LIMIT+1-LAST otherwise
C
C            WORK  - Double precision
C                    Vector of dimension at least LENW
C                    on return
C                    WORK(1), ..., WORK(LAST) contain the left
C                     end points of the subintervals in the
C                     partition of (A,B),
C                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) Contain
C                     the right end points,
C                    WORK(LIMIT*2+1), ...,WORK(LIMIT*2+LAST) contain the
C                     integral approximations over the subintervals,
C                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3)
C                     contain the error estimates.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAGIE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DQAGI
C
      DOUBLE PRECISION ABSERR,BOUND,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,INF,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(*),WORK(*)
C
      EXTERNAL F
C
C         CHECK VALIDITY OF LIMIT AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAGI
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LIMIT.LT.1.OR.LENW.LT.LIMIT*4) GO TO 10
C
C         PREPARE CALL FOR DQAGIE.
C
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
C
      CALL DQAGIE(F,BOUND,INF,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,
     1  NEVAL,IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
C         CALL ERROR HANDLER IF NECESSARY.
C
       LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAGI',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
