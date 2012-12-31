*DECK DQAGP
      SUBROUTINE DQAGP (F, A, B, NPTS2, POINTS, EPSABS, EPSREL, RESULT,
     +   ABSERR, NEVAL, IER, LENIW, LENW, LAST, IWORK, WORK)
C***BEGIN PROLOGUE  DQAGP
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral I = Integral of F over (A,B),
C            hopefully satisfying following claim for accuracy
C            break points of the integration interval, where local
C            difficulties of the integrand may occur (e.g.
C            SINGULARITIES, DISCONTINUITIES), are provided by the user.
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A1
C***TYPE      DOUBLE PRECISION (QAGP-S, DQAGP-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE,
C             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE,
C             SINGULARITIES AT USER SPECIFIED POINTS
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Computation of a definite integral
C        Standard fortran subroutine
C        Double precision version
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     Function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Lower limit of integration
C
C            B      - Double precision
C                     Upper limit of integration
C
C            NPTS2  - Integer
C                     Number equal to two more than the number of
C                     user-supplied break points within the integration
C                     range, NPTS.GE.2.
C                     If NPTS2.LT.2, The routine will end with IER = 6.
C
C            POINTS - Double precision
C                     Vector of dimension NPTS2, the first (NPTS2-2)
C                     elements of which are the user provided break
C                     points. If these points do not constitute an
C                     ascending sequence there will be an automatic
C                     sorting.
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     The routine will end with IER = 6.
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
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine.
C                             The estimates for integral and error are
C                             less reliable. it is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. one can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account). However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties. If
C                             the position of a local difficulty can be
C                             determined (i.e. SINGULARITY,
C                             DISCONTINUITY within the interval), it
C                             should be supplied to the routine as an
C                             element of the vector points. If necessary
C                             an appropriate special-purpose integrator
C                             must be used, which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                             The error may be under-estimated.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 4 The algorithm does not converge.
C                             roundoff error is detected in the
C                             extrapolation table.
C                             It is presumed that the requested
C                             tolerance cannot be achieved, and that
C                             the returned RESULT is the best which
C                             can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. it must be noted that
C                             divergence can occur with any other value
C                             of IER.GT.0.
C                         = 6 The input is invalid because
C                             NPTS2.LT.2 or
C                             break points are specified outside
C                             the integration range or
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             RESULT, ABSERR, NEVAL, LAST are set to
C                             zero.  Except when LENIW or LENW or NPTS2
C                             is invalid, IWORK(1), IWORK(LIMIT+1),
C                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1)
C                             are set to zero.
C                             WORK(1) is set to A and WORK(LIMIT+1)
C                             to B (where LIMIT = (LENIW-NPTS2)/2).
C
C         DIMENSIONING PARAMETERS
C            LENIW - Integer
C                    Dimensioning parameter for IWORK
C                    LENIW determines LIMIT = (LENIW-NPTS2)/2,
C                    which is the maximum number of subintervals in the
C                    partition of the given integration interval (A,B),
C                    LENIW.GE.(3*NPTS2-2).
C                    If LENIW.LT.(3*NPTS2-2), the routine will end with
C                    IER = 6.
C
C            LENW  - Integer
C                    Dimensioning parameter for WORK
C                    LENW must be at least LENIW*2-NPTS2.
C                    If LENW.LT.LENIW*2-NPTS2, the routine will end
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
C                    Vector of dimension at least LENIW. on return,
C                    the first K elements of which contain
C                    pointers to the error estimates over the
C                    subintervals, such that WORK(LIMIT*3+IWORK(1)),...,
C                    WORK(LIMIT*3+IWORK(K)) form a decreasing
C                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2), and
C                    K = LIMIT+1-LAST otherwise
C                    IWORK(LIMIT+1), ...,IWORK(LIMIT+LAST) Contain the
C                     subdivision levels of the subintervals, i.e.
C                     if (AA,BB) is a subinterval of (P1,P2)
C                     where P1 as well as P2 is a user-provided
C                     break point or integration LIMIT, then (AA,BB) has
C                     level L if ABS(BB-AA) = ABS(P2-P1)*2**(-L),
C                    IWORK(LIMIT*2+1), ..., IWORK(LIMIT*2+NPTS2) have
C                     no significance for the user,
C                    note that LIMIT = (LENIW-NPTS2)/2.
C
C            WORK  - Double precision
C                    Vector of dimension at least LENW
C                    on return
C                    WORK(1), ..., WORK(LAST) contain the left
C                     end points of the subintervals in the
C                     partition of (A,B),
C                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
C                     the right end points,
C                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
C                     the integral approximations over the subintervals,
C                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
C                     contain the corresponding error estimates,
C                    WORK(LIMIT*4+1), ..., WORK(LIMIT*4+NPTS2)
C                     contain the integration limits and the
C                     break points sorted in an ascending sequence.
C                    note that LIMIT = (LENIW-NPTS2)/2.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAGPE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DQAGP
C
      DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,POINTS,RESULT,WORK
      INTEGER IER,IWORK,LAST,LENIW,LENW,LIMIT,LVL,L1,L2,L3,L4,NEVAL,
     1  NPTS2
C
      DIMENSION IWORK(*),POINTS(*),WORK(*)
C
      EXTERNAL F
C
C         CHECK VALIDITY OF LIMIT AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAGP
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LENIW.LT.(3*NPTS2-2).OR.LENW.LT.(LENIW*2-NPTS2).OR.NPTS2.LT.2)
     1  GO TO 10
C
C         PREPARE CALL FOR DQAGPE.
C
      LIMIT = (LENIW-NPTS2)/2
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
      L4 = LIMIT+L3
C
      CALL DQAGPE(F,A,B,NPTS2,POINTS,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,
     1  NEVAL,IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),WORK(L4),
     2  IWORK(1),IWORK(L1),IWORK(L2),LAST)
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAGP',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
