*DECK DQAWO
      SUBROUTINE DQAWO (F, A, B, OMEGA, INTEGR, EPSABS, EPSREL, RESULT,
     +   ABSERR, NEVAL, IER, LENIW, MAXP1, LENW, LAST, IWORK, WORK)
C***BEGIN PROLOGUE  DQAWO
C***PURPOSE  Calculate an approximation to a given definite integral
C            I= Integral of F(X)*W(X) over (A,B), where
C                   W(X) = COS(OMEGA*X)
C               or  W(X) = SIN(OMEGA*X),
C            hopefully satisfying the following claim for accuracy
C                ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A1
C***TYPE      DOUBLE PRECISION (QAWO-S, DQAWO-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD,
C             EXTRAPOLATION, GLOBALLY ADAPTIVE,
C             INTEGRAND WITH OSCILLATORY COS OR SIN FACTOR, QUADPACK,
C             QUADRATURE, SPECIAL-PURPOSE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Computation of oscillatory integrals
C        Standard fortran subroutine
C        Double precision version
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the function
C                     F(X).  The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Lower limit of integration
C
C            B      - Double precision
C                     Upper limit of integration
C
C            OMEGA  - Double precision
C                     Parameter in the integrand weight function
C
C            INTEGR - Integer
C                     Indicates which of the weight functions is used
C                     INTEGR = 1      W(X) = COS(OMEGA*X)
C                     INTEGR = 2      W(X) = SIN(OMEGA*X)
C                     If INTEGR.NE.1.AND.INTEGR.NE.2, the routine will
C                     end with IER = 6.
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If EPSABS.LE.0 and
C                     EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     the routine will end with IER = 6.
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
C                   - IER.GT.0 Abnormal termination of the routine.
C                             The estimates for integral and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved (= LENIW/2). One can
C                             allow more subdivisions by increasing the
C                             value of LENIW (and taking the according
C                             dimension adjustments into account).
C                             However, if this yields no improvement it
C                             is advised to analyze the integrand in
C                             order to determine the integration
C                             difficulties. If the position of a local
C                             difficulty can be determined (e.g.
C                             SINGULARITY, DISCONTINUITY within the
C                             interval) one will probably gain from
C                             splitting up the interval at this point
C                             and calling the integrator on the
C                             subranges. If possible, an appropriate
C                             special-purpose integrator should be used
C                             which is designed for handling the type of
C                             difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                             The error may be under-estimated.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some interior points of the
C                             integration interval.
C                         = 4 The algorithm does not converge.
C                             Roundoff error is detected in the
C                             extrapolation table. It is presumed that
C                             the requested tolerance cannot be achieved
C                             due to roundoff in the extrapolation
C                             table, and that the returned result is
C                             the best which can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             or (INTEGR.NE.1 AND INTEGR.NE.2),
C                             or LENIW.LT.2 OR MAXP1.LT.1 or
C                             LENW.LT.LENIW*2+MAXP1*25.
C                             RESULT, ABSERR, NEVAL, LAST are set to
C                             zero. Except when LENIW, MAXP1 or LENW are
C                             invalid, WORK(LIMIT*2+1), WORK(LIMIT*3+1),
C                             IWORK(1), IWORK(LIMIT+1) are set to zero,
C                             WORK(1) is set to A and WORK(LIMIT+1) to
C                             B.
C
C         DIMENSIONING PARAMETERS
C            LENIW  - Integer
C                     Dimensioning parameter for IWORK.
C                     LENIW/2 equals the maximum number of subintervals
C                     allowed in the partition of the given integration
C                     interval (A,B), LENIW.GE.2.
C                     If LENIW.LT.2, the routine will end with IER = 6.
C
C            MAXP1  - Integer
C                     Gives an upper bound on the number of Chebyshev
C                     moments which can be stored, i.e. for the
C                     intervals of lengths ABS(B-A)*2**(-L),
C                     L=0,1, ..., MAXP1-2, MAXP1.GE.1
C                     If MAXP1.LT.1, the routine will end with IER = 6.
C
C            LENW   - Integer
C                     Dimensioning parameter for WORK
C                     LENW must be at least LENIW*2+MAXP1*25.
C                     If LENW.LT.(LENIW*2+MAXP1*25), the routine will
C                     end with IER = 6.
C
C            LAST   - Integer
C                     On return, LAST equals the number of subintervals
C                     produced in the subdivision process, which
C                     determines the number of significant elements
C                     actually in the WORK ARRAYS.
C
C         WORK ARRAYS
C            IWORK  - Integer
C                     Vector of dimension at least LENIW
C                     on return, the first K elements of which contain
C                     pointers to the error estimates over the
C                     subintervals, such that WORK(LIMIT*3+IWORK(1)), ..
C                     WORK(LIMIT*3+IWORK(K)) form a decreasing
C                     sequence, with LIMIT = LENW/2 , and K = LAST
C                     if LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
C                     otherwise.
C                     Furthermore, IWORK(LIMIT+1), ..., IWORK(LIMIT+
C                     LAST) indicate the subdivision levels of the
C                     subintervals, such that IWORK(LIMIT+I) = L means
C                     that the subinterval numbered I is of length
C                     ABS(B-A)*2**(1-L).
C
C            WORK   - Double precision
C                     Vector of dimension at least LENW
C                     On return
C                     WORK(1), ..., WORK(LAST) contain the left
C                      end points of the subintervals in the
C                      partition of (A,B),
C                     WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
C                      the right end points,
C                     WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
C                      the integral approximations over the
C                      subintervals,
C                     WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
C                      contain the error estimates.
C                     WORK(LIMIT*4+1), ..., WORK(LIMIT*4+MAXP1*25)
C                      Provide space for storing the Chebyshev moments.
C                     Note that LIMIT = LENW/2.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAWOE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DQAWO
C
       DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,OMEGA,RESULT,WORK
       INTEGER IER,INTEGR,IWORK,LAST,LIMIT,LENW,LENIW,LVL,L1,L2,L3,L4,
     1  MAXP1,MOMCOM,NEVAL
C
       DIMENSION IWORK(*),WORK(*)
C
       EXTERNAL F
C
C         CHECK VALIDITY OF LENIW, MAXP1 AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAWO
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LENIW.LT.2.OR.MAXP1.LT.1.OR.LENW.LT.(LENIW*2+MAXP1*25))
     1  GO TO 10
C
C         PREPARE CALL FOR DQAWOE
C
      LIMIT = LENIW/2
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
      L4 = LIMIT+L3
      CALL DQAWOE(F,A,B,OMEGA,INTEGR,EPSABS,EPSREL,LIMIT,1,MAXP1,RESULT,
     1   ABSERR,NEVAL,IER,LAST,WORK(1),WORK(L1),WORK(L2),WORK(L3),
     2   IWORK(1),IWORK(L1),MOMCOM,WORK(L4))
C
C         CALL ERROR HANDLER IF NECESSARY
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 0
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAWO',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
