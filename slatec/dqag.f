*DECK DQAG
      SUBROUTINE DQAG (F, A, B, EPSABS, EPSREL, KEY, RESULT, ABSERR,
     +   NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
C***BEGIN PROLOGUE  DQAG
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral I = integral of F over (A,B),
C            hopefully satisfying following claim for accuracy
C            ABS(I-RESULT)LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A1
C***TYPE      DOUBLE PRECISION (QAG-S, DQAG-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD RULES,
C             GENERAL-PURPOSE, GLOBALLY ADAPTIVE, INTEGRAND EXAMINATOR,
C             QUADPACK, QUADRATURE
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
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     Function F(X). The actual name for F needs to be
C                     Declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Lower limit of integration
C
C            B      - Double precision
C                     Upper limit of integration
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     The routine will end with IER = 6.
C
C            KEY    - Integer
C                     Key for choice of local integration rule
C                     A GAUSS-KRONROD PAIR is used with
C                       7 - 15 POINTS If KEY.LT.2,
C                      10 - 21 POINTS If KEY = 2,
C                      15 - 31 POINTS If KEY = 3,
C                      20 - 41 POINTS If KEY = 4,
C                      25 - 51 POINTS If KEY = 5,
C                      30 - 61 POINTS If KEY.GT.5.
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate of the modulus of the absolute error,
C                     Which should EQUAL or EXCEED ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine
C                             The estimates for RESULT and ERROR are
C                             Less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C                      ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account). HOWEVER, If
C                             this yield no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties.
C                             If the position of a local difficulty can
C                             be determined (I.E. SINGULARITY,
C                             DISCONTINUITY WITHIN THE INTERVAL) One
C                             will probably gain from splitting up the
C                             interval at this point and calling the
C                             INTEGRATOR on the SUBRANGES. If possible,
C                             AN APPROPRIATE SPECIAL-PURPOSE INTEGRATOR
C                             should be used which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 AND
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
C                             RESULT, ABSERR, NEVAL, LAST are set
C                             to zero.
C                             EXCEPT when LENW is invalid, IWORK(1),
C                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1) are
C                             set to zero, WORK(1) is set to A and
C                             WORK(LIMIT+1) to B.
C
C         DIMENSIONING PARAMETERS
C            LIMIT - Integer
C                    Dimensioning parameter for IWORK
C                    Limit determines the maximum number of subintervals
C                    in the partition of the given integration interval
C                    (A,B), LIMIT.GE.1.
C                    If LIMIT.LT.1, the routine will end with IER = 6.
C
C            LENW  - Integer
C                    Dimensioning parameter for work
C                    LENW must be at least LIMIT*4.
C                    IF LENW.LT.LIMIT*4, the routine will end with
C                    IER = 6.
C
C            LAST  - Integer
C                    On return, LAST equals the number of subintervals
C                    produced in the subdivision process, which
C                    determines the number of significant elements
C                    actually in the WORK ARRAYS.
C
C         WORK ARRAYS
C            IWORK - Integer
C                    Vector of dimension at least limit, the first K
C                    elements of which contain pointers to the error
C                    estimates over the subintervals, such that
C                    WORK(LIMIT*3+IWORK(1)),... , WORK(LIMIT*3+IWORK(K))
C                    form a decreasing sequence with K = LAST If
C                    LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST otherwise
C
C            WORK  - Double precision
C                    Vector of dimension at least LENW
C                    on return
C                    WORK(1), ..., WORK(LAST) contain the left end
C                    points of the subintervals in the partition of
C                     (A,B),
C                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain the
C                     right end points,
C                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
C                     the integral approximations over the subintervals,
C                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) contain
C                     the error estimates.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAGE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DQAG
      DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(*),WORK(*)
C
      EXTERNAL F
C***FIRST EXECUTABLE STATEMENT  DQAG
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF (LIMIT.GE.1 .AND. LENW.GE.LIMIT*4) THEN
C
C        PREPARE CALL FOR DQAGE.
C
         L1 = LIMIT+1
         L2 = LIMIT+L1
         L3 = LIMIT+L2
C
         CALL DQAGE(F,A,B,EPSABS,EPSREL,KEY,LIMIT,RESULT,ABSERR,NEVAL,
     1     IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
C        CALL ERROR HANDLER IF NECESSARY.
C
         LVL = 0
      ENDIF
C
      IF (IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAG',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
