*DECK QAGPE
      SUBROUTINE QAGPE (F, A, B, NPTS2, POINTS, EPSABS, EPSREL, LIMIT,
     +   RESULT, ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, PTS,
     +   IORD, LEVEL, NDIN, LAST)
C***BEGIN PROLOGUE  QAGPE
C***PURPOSE  Approximate a given definite integral I = Integral of F
C            over (A,B), hopefully satisfying the accuracy claim:
C                  ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C            Break points of the integration interval, where local
C            difficulties of the integrand may occur (e.g. singularities
C            or discontinuities) are provided by the user.
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A1
C***TYPE      SINGLE PRECISION (QAGPE-S, DQAGPE-D)
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
C        Real version
C
C        PARAMETERS
C         ON ENTRY
C            F      - Real
C                     Function subprogram defining the integrand
C                     function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            A      - Real
C                     Lower limit of integration
C
C            B      - Real
C                     Upper limit of integration
C
C            NPTS2  - Integer
C                     Number equal to two more than the number of
C                     user-supplied break points within the integration
C                     range, NPTS2.GE.2.
C                     If NPTS2.LT.2, the routine will end with IER = 6.
C
C            POINTS - Real
C                     Vector of dimension NPTS2, the first (NPTS2-2)
C                     elements of which are the user provided break
C                     POINTS. If these POINTS do not constitute an
C                     ascending sequence there will be an automatic
C                     sorting.
C
C            EPSABS - Real
C                     Absolute accuracy requested
C            EPSREL - Real
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     the routine will end with IER = 6.
C
C            LIMIT  - Integer
C                     Gives an upper bound on the number of subintervals
C                     in the partition of (A,B), LIMIT.GE.NPTS2
C                     If LIMIT.LT.NPTS2, the routine will end with
C                     IER = 6.
C
C         ON RETURN
C            RESULT - Real
C                     Approximation to the integral
C
C            ABSERR - Real
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
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
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
C                             At some points of the integration
C                             interval.
C                         = 4 The algorithm does not converge.
C                             Roundoff error is detected in the
C                             extrapolation table. It is presumed that
C                             the requested tolerance cannot be
C                             achieved, and that the returned result is
C                             the best which can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.GT.0.
C                         = 6 The input is invalid because
C                             NPTS2.LT.2 or
C                             Break points are specified outside
C                             the integration range or
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             or LIMIT.LT.NPTS2.
C                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
C                             and ELIST(1) are set to zero. ALIST(1) and
C                             BLIST(1) are set to A and B respectively.
C
C            ALIST  - Real
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the left end points
C                     of the subintervals in the partition of the given
C                     integration range (A,B)
C
C            BLIST  - Real
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the right end points
C                     of the subintervals in the partition of the given
C                     integration range (A,B)
C
C            RLIST  - Real
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the integral
C                     approximations on the subintervals
C
C            ELIST  - Real
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the moduli of the
C                     absolute error estimates on the subintervals
C
C            PTS    - Real
C                     Vector of dimension at least NPTS2, containing the
C                     integration limits and the break points of the
C                     interval in ascending sequence.
C
C            LEVEL  - Integer
C                     Vector of dimension at least LIMIT, containing the
C                     subdivision levels of the subinterval, i.e. if
C                     (AA,BB) is a subinterval of (P1,P2) where P1 as
C                     well as P2 is a user-provided break point or
C                     integration limit, then (AA,BB) has level L if
C                     ABS(BB-AA) = ABS(P2-P1)*2**(-L).
C
C            NDIN   - Integer
C                     Vector of dimension at least NPTS2, after first
C                     integration over the intervals (PTS(I)),PTS(I+1),
C                     I = 0,1, ..., NPTS2-2, the error estimates over
C                     some of the intervals may have been increased
C                     artificially, in order to put their subdivision
C                     forward. If this happens for the subinterval
C                     numbered K, NDIN(K) is put to 1, otherwise
C                     NDIN(K) = 0.
C
C            IORD   - Integer
C                     Vector of dimension at least LIMIT, the first K
C                     elements of which are pointers to the
C                     error estimates over the subintervals,
C                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
C                     form a decreasing sequence, with K = LAST
C                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
C                     otherwise
C
C            LAST   - Integer
C                     Number of subintervals actually produced in the
C                     subdivisions process
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  QELG, QK21, QPSRT, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QAGPE
      REAL A,ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,
     1  A2,B,BLIST,B1,B2,CORREC,DEFABS,DEFAB1,DEFAB2,
     2  DRES,R1MACH,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND,
     3  ERRMAX,ERROR1,ERRO12,ERROR2,ERRSUM,ERTEST,F,OFLOW,POINTS,PTS,
     4  RESA,RESABS,RESEPS,RESULT,RES3LA,RLIST,RLIST2,SIGN,TEMP,
     5  UFLOW
      INTEGER I,ID,IER,IERRO,IND1,IND2,IORD,IP1,IROFF1,IROFF2,
     1  IROFF3,J,JLOW,JUPBND,K,KSGN,KTMIN,LAST,LEVCUR,LEVEL,LEVMAX,
     2  LIMIT,MAXERR,NDIN,NEVAL,NINT,NINTP1,NPTS,NPTS2,NRES,
     3  NRMAX,NUMRL2
      LOGICAL EXTRAP,NOEXT
C
C
      DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*),
     1  LEVEL(*),NDIN(*),POINTS(*),PTS(*),RES3LA(3),
     2  RLIST(*),RLIST2(52)
C
      EXTERNAL F
C
C            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
C            LIMEXP IN SUBROUTINE EPSALG (RLIST2 SHOULD BE OF DIMENSION
C            (LIMEXP+2) AT LEAST).
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
C           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2
C                       CONTAINING THE PART OF THE EPSILON TABLE WHICH
C                       IS STILL NEEDED FOR FURTHER COMPUTATIONS
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
C           NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN
C                       APPROPRIATE APPROXIMATION TO THE COMPOUNDED
C                       INTEGRAL HAS BEEN OBTAINED, IT IS PUT IN
C                       RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
C                       BY ONE.
C           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
C                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
C           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
C                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
C                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
C                       TRY TO DECREASE THE VALUE OF ERLARG.
C           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION IS
C                       NO LONGER ALLOWED (TRUE-VALUE)
C
C            MACHINE DEPENDENT CONSTANTS
C            ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QAGPE
      EPMACH = R1MACH(4)
C
C            TEST ON VALIDITY OF PARAMETERS
C            -----------------------------
C
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0E+00
      ABSERR = 0.0E+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0E+00
      ELIST(1) = 0.0E+00
      IORD(1) = 0
      LEVEL(1) = 0
      NPTS = NPTS2-2
      IF(NPTS2.LT.2.OR.LIMIT.LE.NPTS.OR.(EPSABS.LE.0.0E+00.AND.
     1  EPSREL.LT.MAX(0.5E+02*EPMACH,0.5E-14))) IER = 6
      IF(IER.EQ.6) GO TO 999
C
C            IF ANY BREAK POINTS ARE PROVIDED, SORT THEM INTO AN
C            ASCENDING SEQUENCE.
C
      SIGN = 1.0E+00
      IF(A.GT.B) SIGN = -1.0E+00
      PTS(1) = MIN(A,B)
      IF(NPTS.EQ.0) GO TO 15
      DO 10 I = 1,NPTS
        PTS(I+1) = POINTS(I)
   10 CONTINUE
   15 PTS(NPTS+2) = MAX(A,B)
      NINT = NPTS+1
      A1 = PTS(1)
      IF(NPTS.EQ.0) GO TO 40
      NINTP1 = NINT+1
      DO 20 I = 1,NINT
        IP1 = I+1
        DO 20 J = IP1,NINTP1
          IF(PTS(I).LE.PTS(J)) GO TO 20
          TEMP = PTS(I)
          PTS(I) = PTS(J)
          PTS(J) = TEMP
   20 CONTINUE
      IF(PTS(1).NE.MIN(A,B).OR.PTS(NINTP1).NE.
     1  MAX(A,B)) IER = 6
      IF(IER.EQ.6) GO TO 999
C
C            COMPUTE FIRST INTEGRAL AND ERROR APPROXIMATIONS.
C            ------------------------------------------------
C
   40 RESABS = 0.0E+00
      DO 50 I = 1,NINT
        B1 = PTS(I+1)
        CALL QK21(F,A1,B1,AREA1,ERROR1,DEFABS,RESA)
        ABSERR = ABSERR+ERROR1
        RESULT = RESULT+AREA1
        NDIN(I) = 0
        IF(ERROR1.EQ.RESA.AND.ERROR1.NE.0.0E+00) NDIN(I) = 1
        RESABS = RESABS+DEFABS
        LEVEL(I) = 0
        ELIST(I) = ERROR1
        ALIST(I) = A1
        BLIST(I) = B1
        RLIST(I) = AREA1
        IORD(I) = I
        A1 = B1
   50 CONTINUE
      ERRSUM = 0.0E+00
      DO 55 I = 1,NINT
        IF(NDIN(I).EQ.1) ELIST(I) = ABSERR
        ERRSUM = ERRSUM+ELIST(I)
   55 CONTINUE
C
C           TEST ON ACCURACY.
C
      LAST = NINT
      NEVAL = 21*NINT
      DRES = ABS(RESULT)
      ERRBND = MAX(EPSABS,EPSREL*DRES)
      IF(ABSERR.LE.0.1E+03*EPMACH*RESABS.AND.ABSERR.GT.
     1  ERRBND) IER = 2
      IF(NINT.EQ.1) GO TO 80
      DO 70 I = 1,NPTS
        JLOW = I+1
        IND1 = IORD(I)
        DO 60 J = JLOW,NINT
          IND2 = IORD(J)
          IF(ELIST(IND1).GT.ELIST(IND2)) GO TO 60
          IND1 = IND2
          K = J
   60   CONTINUE
        IF(IND1.EQ.IORD(I)) GO TO 70
        IORD(K) = IORD(I)
        IORD(I) = IND1
   70 CONTINUE
      IF(LIMIT.LT.NPTS2) IER = 1
   80 IF(IER.NE.0.OR.ABSERR.LE.ERRBND) GO TO 999
C
C           INITIALIZATION
C           --------------
C
      RLIST2(1) = RESULT
      MAXERR = IORD(1)
      ERRMAX = ELIST(MAXERR)
      AREA = RESULT
      NRMAX = 1
      NRES = 0
      NUMRL2 = 1
      KTMIN = 0
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      ERLARG = ERRSUM
      ERTEST = ERRBND
      LEVMAX = 1
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      IERRO = 0
      UFLOW = R1MACH(1)
      OFLOW = R1MACH(2)
      ABSERR = OFLOW
      KSGN = -1
      IF(DRES.GE.(0.1E+01-0.5E+02*EPMACH)*RESABS) KSGN = 1
C
C           MAIN DO-LOOP
C           ------------
C
      DO 160 LAST = NPTS2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST
C           ERROR ESTIMATE.
C
        LEVCUR = LEVEL(MAXERR)+1
        A1 = ALIST(MAXERR)
        B1 = 0.5E+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL QK21(F,A1,B1,AREA1,ERROR1,RESA,DEFAB1)
        CALL QK21(F,A2,B2,AREA2,ERROR2,RESA,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        NEVAL = NEVAL+42
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 95
        IF(ABS(RLIST(MAXERR)-AREA12).GT.0.1E-04*ABS(AREA12)
     1  .OR.ERRO12.LT.0.99E+00*ERRMAX) GO TO 90
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   90   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   95   LEVEL(MAXERR) = LEVCUR
        LEVEL(LAST) = LEVCUR
        RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY
C           SET ERROR FLAG.
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
C           AT A POINT OF THE INTEGRATION RANGE
C
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1E+01+0.1E+03*EPMACH)*
     1  (ABS(A2)+0.1E+04*UFLOW)) IER = 4
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
        IF(ERROR2.GT.ERROR1) GO TO 100
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 110
  100   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE QPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE
C           SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE
C           BISECTED NEXT).
C
  110   CALL QPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C ***JUMP OUT OF DO-LOOP
        IF(ERRSUM.LE.ERRBND) GO TO 190
C ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0) GO TO 170
        IF(NOEXT) GO TO 160
        ERLARG = ERLARG-ERLAST
        IF(LEVCUR+1.LE.LEVMAX) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 120
C
C           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
C           SMALLEST INTERVAL.
C
        IF(LEVEL(MAXERR)+1.LE.LEVMAX) GO TO 160
        EXTRAP = .TRUE.
        NRMAX = 2
  120   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 140
C
C           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
C           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS
C           OVER THE LARGER INTERVALS (ERLARG) AND PERFORM
C           EXTRAPOLATION.
C
        ID = NRMAX
        JUPBND = LAST
        IF(LAST.GT.(2+LIMIT/2)) JUPBND = LIMIT+3-LAST
        DO 130 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
C ***JUMP OUT OF DO-LOOP
          IF(LEVEL(MAXERR)+1.LE.LEVMAX) GO TO 160
          NRMAX = NRMAX+1
  130   CONTINUE
C
C           PERFORM EXTRAPOLATION.
C
  140   NUMRL2 = NUMRL2+1
        RLIST2(NUMRL2) = AREA
        IF(NUMRL2.LE.2) GO TO 155
        CALL QELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
        KTMIN = KTMIN+1
        IF(KTMIN.GT.5.AND.ABSERR.LT.0.1E-02*ERRSUM) IER = 5
        IF(ABSEPS.GE.ABSERR) GO TO 150
        KTMIN = 0
        ABSERR = ABSEPS
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = MAX(EPSABS,EPSREL*ABS(RESEPS))
C ***JUMP OUT OF DO-LOOP
        IF(ABSERR.LT.ERTEST) GO TO 170
C
C           PREPARE BISECTION OF THE SMALLEST INTERVAL.
C
  150   IF(NUMRL2.EQ.1) NOEXT = .TRUE.
        IF(IER.GE.5) GO TO 170
  155   MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        LEVMAX = LEVMAX+1
        ERLARG = ERRSUM
  160 CONTINUE
C
C           SET THE FINAL RESULT.
C           ---------------------
C
C
  170 IF(ABSERR.EQ.OFLOW) GO TO 190
      IF((IER+IERRO).EQ.0) GO TO 180
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0E+00.AND.AREA.NE.0.0E+00)GO TO 175
      IF(ABSERR.GT.ERRSUM)GO TO 190
      IF(AREA.EQ.0.0E+00) GO TO 210
      GO TO 180
  175 IF(ABSERR/ABS(RESULT).GT.ERRSUM/ABS(AREA))GO TO 190
C
C           TEST ON DIVERGENCE.
C
  180 IF(KSGN.EQ.(-1).AND.MAX(ABS(RESULT),ABS(AREA)).LE.
     1  DEFABS*0.1E-01) GO TO 210
      IF(0.1E-01.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1E+03.OR.
     1  ERRSUM.GT.ABS(AREA)) IER = 6
      GO TO 210
C
C           COMPUTE GLOBAL INTEGRAL SUM.
C
  190 RESULT = 0.0E+00
      DO 200 K = 1,LAST
        RESULT = RESULT+RLIST(K)
  200 CONTINUE
      ABSERR = ERRSUM
  210 IF(IER.GT.2) IER = IER - 1
      RESULT = RESULT*SIGN
 999  RETURN
      END
