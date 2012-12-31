*DECK DQAWS
      SUBROUTINE DQAWS (F, A, B, ALFA, BETA, INTEGR, EPSABS, EPSREL,
     +   RESULT, ABSERR, NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
C***BEGIN PROLOGUE  DQAWS
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral I = Integral of F*W over (A,B),
C            (where W shows a singular behaviour at the end points
C            see parameter INTEGR).
C            Hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A1
C***TYPE      DOUBLE PRECISION (QAWS-S, DQAWS-D)
C***KEYWORDS  ALGEBRAIC-LOGARITHMIC END POINT SINGULARITIES,
C             AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD,
C             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE, SPECIAL-PURPOSE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Integration of functions having algebraico-logarithmic
C        end point singularities
C        Standard fortran subroutine
C        Double precision version
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Lower limit of integration
C
C            B      - Double precision
C                     Upper limit of integration, B.GT.A
C                     If B.LE.A, the routine will end with IER = 6.
C
C            ALFA   - Double precision
C                     Parameter in the integrand function, ALFA.GT.(-1)
C                     If ALFA.LE.(-1), the routine will end with
C                     IER = 6.
C
C            BETA   - Double precision
C                     Parameter in the integrand function, BETA.GT.(-1)
C                     If BETA.LE.(-1), the routine will end with
C                     IER = 6.
C
C            INTEGR - Integer
C                     Indicates which WEIGHT function is to be used
C                     = 1  (X-A)**ALFA*(B-X)**BETA
C                     = 2  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
C                     = 3  (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
C                     = 4  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*LOG(B-X)
C                     If INTEGR.LT.1 or INTEGR.GT.4, the routine
C                     will end with IER = 6.
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
C                     Estimate of the modulus of the absolute error,
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
C                             The estimates for the integral and error
C                             are less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account). However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand, in order to
C                             determine the integration difficulties
C                             which prevent the requested tolerance from
C                             being achieved. In case of a jump
C                             discontinuity or a local singularity
C                             of algebraico-logarithmic type at one or
C                             more interior points of the integration
C                             range, one should proceed by splitting up
C                             the interval at these points and calling
C                             the integrator on the subranges.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 6 The input is invalid, because
C                             B.LE.A or ALFA.LE.(-1) or BETA.LE.(-1) or
C                             or INTEGR.LT.1 or INTEGR.GT.4 or
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             or LIMIT.LT.2 or LENW.LT.LIMIT*4.
C                             RESULT, ABSERR, NEVAL, LAST are set to
C                             zero. Except when LENW or LIMIT is invalid
C                             IWORK(1), WORK(LIMIT*2+1) and
C                             WORK(LIMIT*3+1) are set to zero, WORK(1)
C                             is set to A and WORK(LIMIT+1) to B.
C
C         DIMENSIONING PARAMETERS
C            LIMIT  - Integer
C                     Dimensioning parameter for IWORK
C                     LIMIT determines the maximum number of
C                     subintervals in the partition of the given
C                     integration interval (A,B), LIMIT.GE.2.
C                     If LIMIT.LT.2, the routine will end with IER = 6.
C
C            LENW   - Integer
C                     Dimensioning parameter for WORK
C                     LENW must be at least LIMIT*4.
C                     If LENW.LT.LIMIT*4, the routine will end
C                     with IER = 6.
C
C            LAST   - Integer
C                     On return, LAST equals the number of
C                     subintervals produced in the subdivision process,
C                     which determines the significant number of
C                     elements actually in the WORK ARRAYS.
C
C         WORK ARRAYS
C            IWORK  - Integer
C                     Vector of dimension LIMIT, the first K
C                     elements of which contain pointers
C                     to the error estimates over the subintervals,
C                     such that WORK(LIMIT*3+IWORK(1)), ...,
C                     WORK(LIMIT*3+IWORK(K)) form a decreasing
C                     sequence with K = LAST if LAST.LE.(LIMIT/2+2),
C                     and K = LIMIT+1-LAST otherwise
C
C            WORK   - Double precision
C                     Vector of dimension LENW
C                     On return
C                     WORK(1), ..., WORK(LAST) contain the left
C                      end points of the subintervals in the
C                      partition of (A,B),
C                     WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
C                      the right end points,
C                     WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST)
C                      contain the integral approximations over
C                      the subintervals,
C                     WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
C                      contain the error estimates.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAWSE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DQAWS
C
      DOUBLE PRECISION A,ABSERR,ALFA,B,BETA,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,INTEGR,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(*),WORK(*)
C
      EXTERNAL F
C
C         CHECK VALIDITY OF LIMIT AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAWS
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LIMIT.LT.2.OR.LENW.LT.LIMIT*4) GO TO 10
C
C         PREPARE CALL FOR DQAWSE.
C
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
C
      CALL DQAWSE(F,A,B,ALFA,BETA,INTEGR,EPSABS,EPSREL,LIMIT,RESULT,
     1  ABSERR,NEVAL,IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAWS',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
