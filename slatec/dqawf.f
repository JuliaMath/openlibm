*DECK DQAWF
      SUBROUTINE DQAWF (F, A, OMEGA, INTEGR, EPSABS, RESULT, ABSERR,
     +   NEVAL, IER, LIMLST, LST, LENIW, MAXP1, LENW, IWORK, WORK)
C***BEGIN PROLOGUE  DQAWF
C***PURPOSE  The routine calculates an approximation result to a given
C            Fourier integral I=Integral of F(X)*W(X) over (A,INFINITY)
C            where W(X) = COS(OMEGA*X) or W(X) = SIN(OMEGA*X).
C            Hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.EPSABS.
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A3A1
C***TYPE      DOUBLE PRECISION (QAWF-S, DQAWF-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, CONVERGENCE ACCELERATION,
C             FOURIER INTEGRALS, INTEGRATION BETWEEN ZEROS, QUADPACK,
C             QUADRATURE, SPECIAL-PURPOSE INTEGRAL
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Computation of Fourier integrals
C        Standard fortran subroutine
C        Double precision version
C
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
C            OMEGA  - Double precision
C                     Parameter in the integrand WEIGHT function
C
C            INTEGR - Integer
C                     Indicates which of the WEIGHT functions is used
C                     INTEGR = 1      W(X) = COS(OMEGA*X)
C                     INTEGR = 2      W(X) = SIN(OMEGA*X)
C                     IF INTEGR.NE.1.AND.INTEGR.NE.2, the routine
C                     will end with IER = 6.
C
C            EPSABS - Double precision
C                     Absolute accuracy requested, EPSABS.GT.0.
C                     If EPSABS.LE.0, the routine will end with IER = 6.
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
C                     IER.GT.0 Abnormal termination of the routine.
C                             The estimates for integral and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                    If OMEGA.NE.0
C                     IER = 1 Maximum number of cycles allowed
C                             has been achieved, i.e. of subintervals
C                             (A+(K-1)C,A+KC) where
C                             C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
C                             FOR K = 1, 2, ..., LST.
C                             One can allow more cycles by increasing
C                             the value of LIMLST (and taking the
C                             according dimension adjustments into
C                             account). Examine the array IWORK which
C                             contains the error flags on the cycles, in
C                             order to look for eventual local
C                             integration difficulties.
C                             If the position of a local difficulty
C                             can be determined (e.g. singularity,
C                             discontinuity within the interval) one
C                             will probably gain from splitting up the
C                             interval at this point and calling
C                             appropriate integrators on the subranges.
C                         = 4 The extrapolation table constructed for
C                             convergence acceleration of the series
C                             formed by the integral contributions over
C                             the cycles, does not converge to within
C                             the requested accuracy.
C                             As in the case of IER = 1, it is advised
C                             to examine the array IWORK which contains
C                             the error flags on the cycles.
C                         = 6 The input is invalid because
C                             (INTEGR.NE.1 AND INTEGR.NE.2) or
C                              EPSABS.LE.0 or LIMLST.LT.1 or
C                              LENIW.LT.(LIMLST+2) or MAXP1.LT.1 or
C                              LENW.LT.(LENIW*2+MAXP1*25).
C                              RESULT, ABSERR, NEVAL, LST are set to
C                              zero.
C                         = 7 Bad integrand behaviour occurs within
C                             one or more of the cycles. Location and
C                             type of the difficulty involved can be
C                             determined from the first LST elements of
C                             vector IWORK.  Here LST is the number of
C                             cycles actually needed (see below).
C                             IWORK(K) = 1 The maximum number of
C                                          subdivisions (=(LENIW-LIMLST)
C                                          /2) has been achieved on the
C                                          K th cycle.
C                                      = 2 Occurrence of roundoff error
C                                          is detected and prevents the
C                                          tolerance imposed on the K th
C                                          cycle, from being achieved
C                                          on this cycle.
C                                      = 3 Extremely bad integrand
C                                          behaviour occurs at some
C                                          points of the K th cycle.
C                                      = 4 The integration procedure
C                                          over the K th cycle does
C                                          not converge (to within the
C                                          required accuracy) due to
C                                          roundoff in the extrapolation
C                                          procedure invoked on this
C                                          cycle. It is assumed that the
C                                          result on this interval is
C                                          the best which can be
C                                          obtained.
C                                      = 5 The integral over the K th
C                                          cycle is probably divergent
C                                          or slowly convergent. It must
C                                          be noted that divergence can
C                                          occur with any other value of
C                                          IWORK(K).
C                    If OMEGA = 0 and INTEGR = 1,
C                    The integral is calculated by means of DQAGIE,
C                    and IER = IWORK(1) (with meaning as described
C                    for IWORK(K),K = 1).
C
C         DIMENSIONING PARAMETERS
C            LIMLST - Integer
C                     LIMLST gives an upper bound on the number of
C                     cycles, LIMLST.GE.3.
C                     If LIMLST.LT.3, the routine will end with IER = 6.
C
C            LST    - Integer
C                     On return, LST indicates the number of cycles
C                     actually needed for the integration.
C                     If OMEGA = 0, then LST is set to 1.
C
C            LENIW  - Integer
C                     Dimensioning parameter for IWORK. On entry,
C                     (LENIW-LIMLST)/2 equals the maximum number of
C                     subintervals allowed in the partition of each
C                     cycle, LENIW.GE.(LIMLST+2).
C                     If LENIW.LT.(LIMLST+2), the routine will end with
C                     IER = 6.
C
C            MAXP1  - Integer
C                     MAXP1 gives an upper bound on the number of
C                     Chebyshev moments which can be stored, i.e. for
C                     the intervals of lengths ABS(B-A)*2**(-L),
C                     L = 0,1, ..., MAXP1-2, MAXP1.GE.1.
C                     If MAXP1.LT.1, the routine will end with IER = 6.
C            LENW   - Integer
C                     Dimensioning parameter for WORK
C                     LENW must be at least LENIW*2+MAXP1*25.
C                     If LENW.LT.(LENIW*2+MAXP1*25), the routine will
C                     end with IER = 6.
C
C         WORK ARRAYS
C            IWORK  - Integer
C                     Vector of dimension at least LENIW
C                     On return, IWORK(K) FOR K = 1, 2, ..., LST
C                     contain the error flags on the cycles.
C
C            WORK   - Double precision
C                     Vector of dimension at least
C                     On return,
C                     WORK(1), ..., WORK(LST) contain the integral
C                      approximations over the cycles,
C                     WORK(LIMLST+1), ..., WORK(LIMLST+LST) contain
C                      the error estimates over the cycles.
C                     further elements of WORK have no specific
C                     meaning for the user.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAWFE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced variable.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DQAWF
C
       DOUBLE PRECISION A,ABSERR,EPSABS,F,OMEGA,RESULT,WORK
       INTEGER IER,INTEGR,IWORK,LENIW,LENW,LIMIT,LIMLST,LL2,LVL,
     1  LST,L1,L2,L3,L4,L5,L6,MAXP1,NEVAL
C
       DIMENSION IWORK(*),WORK(*)
C
       EXTERNAL F
C
C         CHECK VALIDITY OF LIMLST, LENIW, MAXP1 AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAWF
      IER = 6
      NEVAL = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LIMLST.LT.3.OR.LENIW.LT.(LIMLST+2).OR.MAXP1.LT.1.OR.LENW.LT.
     1   (LENIW*2+MAXP1*25)) GO TO 10
C
C         PREPARE CALL FOR DQAWFE
C
      LIMIT = (LENIW-LIMLST)/2
      L1 = LIMLST+1
      L2 = LIMLST+L1
      L3 = LIMIT+L2
      L4 = LIMIT+L3
      L5 = LIMIT+L4
      L6 = LIMIT+L5
      LL2 = LIMIT+L1
      CALL DQAWFE(F,A,OMEGA,INTEGR,EPSABS,LIMLST,LIMIT,MAXP1,RESULT,
     1  ABSERR,NEVAL,IER,WORK(1),WORK(L1),IWORK(1),LST,WORK(L2),
     2  WORK(L3),WORK(L4),WORK(L5),IWORK(L1),IWORK(LL2),WORK(L6))
C
C         CALL ERROR HANDLER IF NECESSARY
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAWF',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
