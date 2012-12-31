*DECK DQAGIE
      SUBROUTINE DQAGIE (F, BOUND, INF, EPSABS, EPSREL, LIMIT, RESULT,
     +   ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
C***BEGIN PROLOGUE  DQAGIE
C***PURPOSE  The routine calculates an approximation result to a given
C            integral   I = Integral of F over (BOUND,+INFINITY)
C            or I = Integral of F over (-INFINITY,BOUND)
C            or I = Integral of F over (-INFINITY,+INFINITY),
C            hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I))
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A3A1, H2A4A1
C***TYPE      DOUBLE PRECISION (QAGIE-S, DQAGIE-D)
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
C Integration over infinite intervals
C Standard fortran subroutine
C
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            BOUND  - Double precision
C                     Finite bound of integration range
C                     (has no meaning if interval is doubly-infinite)
C
C            INF    - Double precision
C                     Indicating the kind of integration range involved
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
C            LIMIT  - Integer
C                     Gives an upper bound on the number of subintervals
C                     in the partition of (A,B), LIMIT.GE.1
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
C                   - IER.GT.0 Abnormal termination of the routine. The
C                             estimates for result and error are less
C                             reliable. It is assumed that the requested
C                             accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account).  However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties.
C                             If the position of a local difficulty can
C                             be determined (e.g. SINGULARITY,
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
C                             result is the best which can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
C                             ELIST(1) and IORD(1) are set to zero.
C                             ALIST(1) and BLIST(1) are set to 0
C                             and 1 respectively.
C
C            ALIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the left
C                     end points of the subintervals in the partition
C                     of the transformed integration range (0,1).
C
C            BLIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the right
C                     end points of the subintervals in the partition
C                     of the transformed integration range (0,1).
C
C            RLIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the integral
C                     approximations on the subintervals
C
C            ELIST  - Double precision
C                     Vector of dimension at least LIMIT,  the first
C                     LAST elements of which are the moduli of the
C                     absolute error estimates on the subintervals
C
C            IORD   - Integer
C                     Vector of dimension LIMIT, the first K
C                     elements of which are pointers to the
C                     error estimates over the subintervals,
C                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
C                     form a decreasing sequence, with K = LAST
C                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
C                     otherwise
C
C            LAST   - Integer
C                     Number of subintervals actually produced
C                     in the subdivision process
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DQELG, DQK15I, DQPSRT
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQAGIE
      DOUBLE PRECISION ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,
     1  A2,BLIST,BOUN,BOUND,B1,B2,CORREC,DEFABS,DEFAB1,DEFAB2,
     2  DRES,D1MACH,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,
     3  ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,ERTEST,F,OFLOW,RESABS,
     4  RESEPS,RESULT,RES3LA,RLIST,RLIST2,SMALL,UFLOW
      INTEGER ID,IER,IERRO,INF,IORD,IROFF1,IROFF2,IROFF3,JUPBND,K,KSGN,
     1  KTMIN,LAST,LIMIT,MAXERR,NEVAL,NRES,NRMAX,NUMRL2
      LOGICAL EXTRAP,NOEXT
C
      DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*),
     1  RES3LA(3),RLIST(*),RLIST2(52)
C
      EXTERNAL F
C
C            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
C            LIMEXP IN SUBROUTINE DQELG.
C
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                       (ALIST(I),BLIST(I))
C           RLIST2    - ARRAY OF DIMENSION AT LEAST (LIMEXP+2),
C                       CONTAINING THE PART OF THE EPSILON TABLE
C                       WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
C                       ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
C                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
C           NUMRL2    - NUMBER OF ELEMENTS CURRENTLY IN RLIST2. IF AN
C                       APPROPRIATE APPROXIMATION TO THE COMPOUNDED
C                       INTEGRAL HAS BEEN OBTAINED, IT IS PUT IN
C                       RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
C                       BY ONE.
C           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED UP
C                       TO NOW, MULTIPLIED BY 1.5
C           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
C                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
C           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
C                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
C                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
C                       TRY TO DECREASE THE VALUE OF ERLARG.
C           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
C                       IS NO LONGER ALLOWED (TRUE-VALUE)
C
C            MACHINE DEPENDENT CONSTANTS
C            ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQAGIE
      EPMACH = D1MACH(4)
C
C           TEST ON VALIDITY OF PARAMETERS
C           -----------------------------
C
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = 0.0D+00
      BLIST(1) = 0.1D+01
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      IF(EPSABS.LE.0.0D+00.AND.EPSREL.LT.MAX(0.5D+02*EPMACH,0.5D-28))
     1  IER = 6
       IF(IER.EQ.6) GO TO 999
C
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
C           DETERMINE THE INTERVAL TO BE MAPPED ONTO (0,1).
C           IF INF = 2 THE INTEGRAL IS COMPUTED AS I = I1+I2, WHERE
C           I1 = INTEGRAL OF F OVER (-INFINITY,0),
C           I2 = INTEGRAL OF F OVER (0,+INFINITY).
C
      BOUN = BOUND
      IF(INF.EQ.2) BOUN = 0.0D+00
      CALL DQK15I(F,BOUN,INF,0.0D+00,0.1D+01,RESULT,ABSERR,
     1  DEFABS,RESABS)
C
C           TEST ON ACCURACY
C
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
      DRES = ABS(RESULT)
      ERRBND = MAX(EPSABS,EPSREL*DRES)
      IF(ABSERR.LE.1.0D+02*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS).OR.
     1  ABSERR.EQ.0.0D+00) GO TO 130
C
C           INITIALIZATION
C           --------------
C
      UFLOW = D1MACH(1)
      OFLOW = D1MACH(2)
      RLIST2(1) = RESULT
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      ABSERR = OFLOW
      NRMAX = 1
      NRES = 0
      KTMIN = 0
      NUMRL2 = 2
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      IERRO = 0
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      KSGN = -1
      IF(DRES.GE.(0.1D+01-0.5D+02*EPMACH)*DEFABS) KSGN = 1
C
C           MAIN DO-LOOP
C           ------------
C
      DO 90 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE.
C
        A1 = ALIST(MAXERR)
        B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL DQK15I(F,BOUN,INF,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        CALL DQK15I(F,BOUN,INF,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2)GO TO 15
        IF(ABS(RLIST(MAXERR)-AREA12).GT.0.1D-04*ABS(AREA12)
     1  .OR.ERRO12.LT.0.99D+00*ERRMAX) GO TO 10
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   10   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   15   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
        IF(IROFF1+IROFF2.GE.10.OR.IROFF3.GE.20) IER = 2
        IF(IROFF2.GE.5) IERRO = 3
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
C           SUBINTERVALS EQUALS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT SOME POINTS OF THE INTEGRATION RANGE.
C
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1D+01+0.1D+03*EPMACH)*
     1  (ABS(A2)+0.1D+04*UFLOW)) IER = 4
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
        IF(ERROR2.GT.ERROR1) GO TO 20
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 30
   20   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
   30   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
        IF(ERRSUM.LE.ERRBND) GO TO 115
        IF(IER.NE.0) GO TO 100
        IF(LAST.EQ.2) GO TO 80
        IF(NOEXT) GO TO 90
        ERLARG = ERLARG-ERLAST
        IF(ABS(B1-A1).GT.SMALL) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 40
C
C           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
C           SMALLEST INTERVAL.
C
        IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
        EXTRAP = .TRUE.
        NRMAX = 2
   40   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 60
C
C           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
C           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER THE
C           LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
C
        ID = NRMAX
        JUPBND = LAST
        IF(LAST.GT.(2+LIMIT/2)) JUPBND = LIMIT+3-LAST
        DO 50 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
          IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
          NRMAX = NRMAX+1
   50   CONTINUE
C
C           PERFORM EXTRAPOLATION.
C
   60   NUMRL2 = NUMRL2+1
        RLIST2(NUMRL2) = AREA
        CALL DQELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
        KTMIN = KTMIN+1
        IF(KTMIN.GT.5.AND.ABSERR.LT.0.1D-02*ERRSUM) IER = 5
        IF(ABSEPS.GE.ABSERR) GO TO 70
        KTMIN = 0
        ABSERR = ABSEPS
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = MAX(EPSABS,EPSREL*ABS(RESEPS))
        IF(ABSERR.LE.ERTEST) GO TO 100
C
C            PREPARE BISECTION OF THE SMALLEST INTERVAL.
C
   70   IF(NUMRL2.EQ.1) NOEXT = .TRUE.
        IF(IER.EQ.5) GO TO 100
        MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        SMALL = SMALL*0.5D+00
        ERLARG = ERRSUM
        GO TO 90
   80   SMALL = 0.375D+00
        ERLARG = ERRSUM
        ERTEST = ERRBND
        RLIST2(2) = AREA
   90 CONTINUE
C
C           SET FINAL RESULT AND ERROR ESTIMATE.
C           ------------------------------------
C
  100 IF(ABSERR.EQ.OFLOW) GO TO 115
      IF((IER+IERRO).EQ.0) GO TO 110
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0D+00.AND.AREA.NE.0.0D+00)GO TO 105
      IF(ABSERR.GT.ERRSUM)GO TO 115
      IF(AREA.EQ.0.0D+00) GO TO 130
      GO TO 110
  105 IF(ABSERR/ABS(RESULT).GT.ERRSUM/ABS(AREA))GO TO 115
C
C           TEST ON DIVERGENCE
C
  110 IF(KSGN.EQ.(-1).AND.MAX(ABS(RESULT),ABS(AREA)).LE.
     1 DEFABS*0.1D-01) GO TO 130
      IF (0.1D-01.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1D+03
     1    .OR.ERRSUM.GT.ABS(AREA)) IER = 6
      GO TO 130
C
C           COMPUTE GLOBAL INTEGRAL SUM.
C
  115 RESULT = 0.0D+00
      DO 120 K = 1,LAST
        RESULT = RESULT+RLIST(K)
  120 CONTINUE
      ABSERR = ERRSUM
  130 NEVAL = 30*LAST-15
      IF(INF.EQ.2) NEVAL = 2*NEVAL
      IF(IER.GT.2) IER=IER-1
  999 RETURN
      END
