*DECK DQAGS
      SUBROUTINE DQAGS (F, A, B, EPSABS, EPSREL, RESULT, ABSERR, NEVAL,
     +   IER, LIMIT, LENW, LAST, IWORK, WORK)
C***BEGIN PROLOGUE  DQAGS
C***PURPOSE  The routine calculates an approximation result to a given
C            Definite integral  I = Integral of F over (A,B),
C            Hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A1
C***TYPE      DOUBLE PRECISION (QAGS-S, DQAGS-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, END POINT SINGULARITIES,
C             EXTRAPOLATION, GENERAL-PURPOSE, GLOBALLY ADAPTIVE,
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
C
C        PARAMETERS
C         ON ENTRY
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
C                     IER.GT.0 Abnormal termination of the routine
C                             The estimates for integral and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more sub-
C                             divisions by increasing the value of LIMIT
C                             (and taking the according dimension
C                             adjustments into account. However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties. If
C                             the position of a local difficulty can be
C                             determined (E.G. SINGULARITY,
C                             DISCONTINUITY WITHIN THE INTERVAL) one
C                             will probably gain from splitting up the
C                             interval at this point and calling the
C                             integrator on the subranges. If possible,
C                             an appropriate special-purpose integrator
C                             should be used, which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is detec-
C                             ted, which prevents the requested
C                             tolerance from being achieved.
C                             The error may be under-estimated.
C                         = 3 Extremely bad integrand behaviour
C                             occurs at some points of the integration
C                             interval.
C                         = 4 The algorithm does not converge.
C                             Roundoff error is detected in the
C                             Extrapolation table. It is presumed that
C                             the requested tolerance cannot be
C                             achieved, and that the returned result is
C                             the best which can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 AND
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)
C                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
C                             RESULT, ABSERR, NEVAL, LAST are set to
C                             zero.  Except when LIMIT or LENW is
C                             invalid, IWORK(1), WORK(LIMIT*2+1) and
C                             WORK(LIMIT*3+1) are set to zero, WORK(1)
C                             is set to A and WORK(LIMIT+1) TO B.
C
C         DIMENSIONING PARAMETERS
C            LIMIT - Integer
C                    DIMENSIONING PARAMETER FOR IWORK
C                    LIMIT determines the maximum number of subintervals
C                    in the partition of the given integration interval
C                    (A,B), LIMIT.GE.1.
C                    IF LIMIT.LT.1, the routine will end with IER = 6.
C
C            LENW  - Integer
C                    DIMENSIONING PARAMETER FOR WORK
C                    LENW must be at least LIMIT*4.
C                    If LENW.LT.LIMIT*4, the routine will end
C                    with IER = 6.
C
C            LAST  - Integer
C                    On return, LAST equals the number of subintervals
C                    produced in the subdivision process, determines the
C                    number of significant elements actually in the WORK
C                    Arrays.
C
C         WORK ARRAYS
C            IWORK - Integer
C                    Vector of dimension at least LIMIT, the first K
C                    elements of which contain pointers
C                    to the error estimates over the subintervals
C                    such that WORK(LIMIT*3+IWORK(1)),... ,
C                    WORK(LIMIT*3+IWORK(K)) form a decreasing
C                    sequence, with K = LAST IF LAST.LE.(LIMIT/2+2),
C                    and K = LIMIT+1-LAST otherwise
C
C            WORK  - Double precision
C                    Vector of dimension at least LENW
C                    on return
C                    WORK(1), ..., WORK(LAST) contain the left
C                     end-points of the subintervals in the
C                     partition of (A,B),
C                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
C                     the right end-points,
C                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
C                     the integral approximations over the subintervals,
C                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
C                     contain the error estimates.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAGSE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DQAGS
C
C
      DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(*),WORK(*)
C
      EXTERNAL F
C
C         CHECK VALIDITY OF LIMIT AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAGS
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LIMIT.LT.1.OR.LENW.LT.LIMIT*4) GO TO 10
C
C         PREPARE CALL FOR DQAGSE.
C
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
C
      CALL DQAGSE(F,A,B,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,NEVAL,
     1  IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAGS',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
