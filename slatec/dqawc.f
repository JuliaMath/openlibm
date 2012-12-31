*DECK DQAWC
      SUBROUTINE DQAWC (F, A, B, C, EPSABS, EPSREL, RESULT, ABSERR,
     +   NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
C***BEGIN PROLOGUE  DQAWC
C***PURPOSE  The routine calculates an approximation result to a
C            Cauchy principal value I = INTEGRAL of F*W over (A,B)
C            (W(X) = 1/((X-C), C.NE.A, C.NE.B), hopefully satisfying
C            following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABE,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A1, J4
C***TYPE      DOUBLE PRECISION (QAWC-S, DQAWC-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, CAUCHY PRINCIPAL VALUE,
C             CLENSHAW-CURTIS METHOD, GLOBALLY ADAPTIVE, QUADPACK,
C             QUADRATURE, SPECIAL-PURPOSE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Computation of a Cauchy principal value
C        Standard fortran subroutine
C        Double precision version
C
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     Function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Under limit of integration
C
C            B      - Double precision
C                     Upper limit of integration
C
C            C      - Parameter in the weight function, C.NE.A, C.NE.B.
C                     If C = A or C = B, the routine will end with
C                     IER = 6 .
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     the routine will end with IER = 6.
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate or the modulus of the absolute error,
C                     Which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine
C                             the estimates for integral and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more sub-
C                             divisions by increasing the value of LIMIT
C                             (and taking the according dimension
C                             adjustments into account). However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties.
C                             If the position of a local difficulty
C                             can be determined (e.g. SINGULARITY,
C                             DISCONTINUITY within the interval) one
C                             will probably gain from splitting up the
C                             interval at this point and calling
C                             appropriate integrators on the subranges.
C                         = 2 The occurrence of roundoff error is detec-
C                             ted, which prevents the requested
C                             tolerance from being achieved.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 6 The input is invalid, because
C                             C = A or C = B or
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             or LIMIT.LT.1 or LENW.LT.LIMIT*4.
C                             RESULT, ABSERR, NEVAL, LAST are set to
C                             zero.  Except when LENW or LIMIT is
C                             invalid, IWORK(1), WORK(LIMIT*2+1) and
C                             WORK(LIMIT*3+1) are set to zero, WORK(1)
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
C           LENW   - Integer
C                    Dimensioning parameter for WORK
C                    LENW must be at least LIMIT*4.
C                    If LENW.LT.LIMIT*4, the routine will end with
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
C                    Vector of dimension at least LIMIT, the first K
C                    elements of which contain pointers
C                    to the error estimates over the subintervals,
C                    such that WORK(LIMIT*3+IWORK(1)), ... ,
C                    WORK(LIMIT*3+IWORK(K)) form a decreasing
C                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2),
C                    and K = LIMIT+1-LAST otherwise
C
C            WORK  - Double precision
C                    Vector of dimension at least LENW
C                    On return
C                    WORK(1), ..., WORK(LAST) contain the left
C                     end points of the subintervals in the
C                     partition of (A,B),
C                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
C                     the right end points,
C                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
C                     the integral approximations over the subintervals,
C                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
C                     contain the error estimates.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAWCE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DQAWC
C
      DOUBLE PRECISION A,ABSERR,B,C,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(*),WORK(*)
C
      EXTERNAL F
C
C         CHECK VALIDITY OF LIMIT AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAWC
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LIMIT.LT.1.OR.LENW.LT.LIMIT*4) GO TO 10
C
C         PREPARE CALL FOR DQAWCE.
C
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
      CALL DQAWCE(F,A,B,C,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,NEVAL,IER,
     1  WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAWC',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
