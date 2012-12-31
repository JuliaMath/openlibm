*DECK DQAWOE
      SUBROUTINE DQAWOE (F, A, B, OMEGA, INTEGR, EPSABS, EPSREL, LIMIT,
     +   ICALL, MAXP1, RESULT, ABSERR, NEVAL, IER, LAST, ALIST, BLIST,
     +   RLIST, ELIST, IORD, NNLOG, MOMCOM, CHEBMO)
C***BEGIN PROLOGUE  DQAWOE
C***PURPOSE  Calculate an approximation to a given definite integral
C            I = Integral of F(X)*W(X) over (A,B), where
C                     W(X) = COS(OMEGA*X)
C                 or  W(X)=SIN(OMEGA*X),
C            hopefully satisfying the following claim for accuracy
C                 ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A1
C***TYPE      DOUBLE PRECISION (QAWOE-S, DQAWOE-D)
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
C        Computation of Oscillatory integrals
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
C                     Upper limit of integration
C
C            OMEGA  - Double precision
C                     Parameter in the integrand weight function
C
C            INTEGR - Integer
C                     Indicates which of the WEIGHT functions is to be
C                     used
C                     INTEGR = 1      W(X) = COS(OMEGA*X)
C                     INTEGR = 2      W(X) = SIN(OMEGA*X)
C                     If INTEGR.NE.1 and INTEGR.NE.2, the routine
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
C            LIMIT  - Integer
C                     Gives an upper bound on the number of subdivisions
C                     in the partition of (A,B), LIMIT.GE.1.
C
C            ICALL  - Integer
C                     If DQAWOE is to be used only once, ICALL must
C                     be set to 1.  Assume that during this call, the
C                     Chebyshev moments (for CLENSHAW-CURTIS integration
C                     of degree 24) have been computed for intervals of
C                     lengths (ABS(B-A))*2**(-L), L=0,1,2,...MOMCOM-1.
C                     If ICALL.GT.1 this means that DQAWOE has been
C                     called twice or more on intervals of the same
C                     length ABS(B-A). The Chebyshev moments already
C                     computed are then re-used in subsequent calls.
C                     If ICALL.LT.1, the routine will end with IER = 6.
C
C            MAXP1  - Integer
C                     Gives an upper bound on the number of Chebyshev
C                     moments which can be stored, i.e. for the
C                     intervals of lengths ABS(B-A)*2**(-L),
C                     L=0,1, ..., MAXP1-2, MAXP1.GE.1.
C                     If MAXP1.LT.1, the routine will end with IER = 6.
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
C                             routine. It is assumed that the
C                             requested accuracy has been achieved.
C                   - IER.GT.0 Abnormal termination of the routine.
C                             The estimates for integral and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking according dimension
C                             adjustments into account). However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand, in order to
C                             determine the integration difficulties.
C                             If the position of a local difficulty can
C                             be determined (e.g. SINGULARITY,
C                             DISCONTINUITY within the interval) one
C                             will probably gain from splitting up the
C                             interval at this point and calling the
C                             integrator on the subranges. If possible,
C                             an appropriate special-purpose integrator
C                             should be used which is designed for
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
C                             It is presumed that the requested
C                             tolerance cannot be achieved due to
C                             roundoff in the extrapolation table,
C                             and that the returned result is the
C                             best which can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.GT.0.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             or (INTEGR.NE.1 and INTEGR.NE.2) or
C                             ICALL.LT.1 or MAXP1.LT.1.
C                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
C                             ELIST(1), IORD(1) and NNLOG(1) are set
C                             to ZERO. ALIST(1) and BLIST(1) are set
C                             to A and B respectively.
C
C            LAST  -  Integer
C                     On return, LAST equals the number of
C                     subintervals produces in the subdivision
C                     process, which determines the number of
C                     significant elements actually in the
C                     WORK ARRAYS.
C            ALIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the left
C                     end points of the subintervals in the partition
C                     of the given integration range (A,B)
C
C            BLIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the right
C                     end points of the subintervals in the partition
C                     of the given integration range (A,B)
C
C            RLIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the integral
C                     approximations on the subintervals
C
C            ELIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the moduli of the
C                     absolute error estimates on the subintervals
C
C            IORD   - Integer
C                     Vector of dimension at least LIMIT, the first K
C                     elements of which are pointers to the error
C                     estimates over the subintervals,
C                     such that ELIST(IORD(1)), ...,
C                     ELIST(IORD(K)) form a decreasing sequence, with
C                     K = LAST if LAST.LE.(LIMIT/2+2), and
C                     K = LIMIT+1-LAST otherwise.
C
C            NNLOG  - Integer
C                     Vector of dimension at least LIMIT, containing the
C                     subdivision levels of the subintervals, i.e.
C                     IWORK(I) = L means that the subinterval
C                     numbered I is of length ABS(B-A)*2**(1-L)
C
C         ON ENTRY AND RETURN
C            MOMCOM - Integer
C                     Indicating that the Chebyshev moments
C                     have been computed for intervals of lengths
C                     (ABS(B-A))*2**(-L), L=0,1,2, ..., MOMCOM-1,
C                     MOMCOM.LT.MAXP1
C
C            CHEBMO - Double precision
C                     Array of dimension (MAXP1,25) containing the
C                     Chebyshev moments
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DQC25F, DQELG, DQPSRT
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQAWOE
C
      DOUBLE PRECISION A,ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,
     1  A2,B,BLIST,B1,B2,CHEBMO,CORREC,DEFAB1,DEFAB2,DEFABS,
     2  DOMEGA,D1MACH,DRES,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,
     3  ERRBND,ERRMAX,ERROR1,ERRO12,ERROR2,ERRSUM,ERTEST,F,OFLOW,
     4  OMEGA,RESABS,RESEPS,RESULT,RES3LA,RLIST,RLIST2,SMALL,UFLOW,WIDTH
      INTEGER ICALL,ID,IER,IERRO,INTEGR,IORD,IROFF1,IROFF2,IROFF3,
     1  JUPBND,K,KSGN,KTMIN,LAST,LIMIT,MAXERR,MAXP1,MOMCOM,NEV,NEVAL,
     2  NNLOG,NRES,NRMAX,NRMOM,NUMRL2
      LOGICAL EXTRAP,NOEXT,EXTALL
C
      DIMENSION ALIST(*),BLIST(*),RLIST(*),ELIST(*),
     1  IORD(*),RLIST2(52),RES3LA(3),CHEBMO(MAXP1,25),NNLOG(*)
C
      EXTERNAL F
C
C            THE DIMENSION OF RLIST2 IS DETERMINED BY  THE VALUE OF
C            LIMEXP IN SUBROUTINE DQELG (RLIST2 SHOULD BE OF
C            DIMENSION (LIMEXP+2) AT LEAST).
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
C                       CONTAINING THE PART OF THE EPSILON TABLE
C                       WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
C                       ERROR ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
C           NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN APPROPRIATE
C                       APPROXIMATION TO THE COMPOUNDED INTEGRAL HAS
C                       BEEN OBTAINED IT IS PUT IN RLIST2(NUMRL2) AFTER
C                       NUMRL2 HAS BEEN INCREASED BY ONE
C           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED
C                       UP TO NOW, MULTIPLIED BY 1.5
C           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
C                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
C           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE IS
C                       ATTEMPTING TO PERFORM EXTRAPOLATION, I.E. BEFORE
C                       SUBDIVIDING THE SMALLEST INTERVAL WE TRY TO
C                       DECREASE THE VALUE OF ERLARG
C           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
C                       IS NO LONGER ALLOWED (TRUE  VALUE)
C
C            MACHINE DEPENDENT CONSTANTS
C            ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQAWOE
      EPMACH = D1MACH(4)
C
C         TEST ON VALIDITY OF PARAMETERS
C         ------------------------------
C
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      NNLOG(1) = 0
      IF((INTEGR.NE.1.AND.INTEGR.NE.2).OR.(EPSABS.LE.0.0D+00.AND.
     1  EPSREL.LT.MAX(0.5D+02*EPMACH,0.5D-28)).OR.ICALL.LT.1.OR.
     2  MAXP1.LT.1) IER = 6
      IF(IER.EQ.6) GO TO 999
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      DOMEGA = ABS(OMEGA)
      NRMOM = 0
      IF (ICALL.GT.1) GO TO 5
      MOMCOM = 0
    5 CALL DQC25F(F,A,B,DOMEGA,INTEGR,NRMOM,MAXP1,0,RESULT,ABSERR,
     1  NEVAL,DEFABS,RESABS,MOMCOM,CHEBMO)
C
C           TEST ON ACCURACY.
C
      DRES = ABS(RESULT)
      ERRBND = MAX(EPSABS,EPSREL*DRES)
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
      IF(ABSERR.LE.0.1D+03*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.ABSERR.LE.ERRBND) GO TO 200
C
C           INITIALIZATIONS
C           ---------------
C
      UFLOW = D1MACH(1)
      OFLOW = D1MACH(2)
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      ABSERR = OFLOW
      NRMAX = 1
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      IERRO = 0
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      KTMIN = 0
      SMALL = ABS(B-A)*0.75D+00
      NRES = 0
      NUMRL2 = 0
      EXTALL = .FALSE.
      IF(0.5D+00*ABS(B-A)*DOMEGA.GT.0.2D+01) GO TO 10
      NUMRL2 = 1
      EXTALL = .TRUE.
      RLIST2(1) = RESULT
   10 IF(0.25D+00*ABS(B-A)*DOMEGA.LE.0.2D+01) EXTALL = .TRUE.
      KSGN = -1
      IF(DRES.GE.(0.1D+01-0.5D+02*EPMACH)*DEFABS) KSGN = 1
C
C           MAIN DO-LOOP
C           ------------
C
      DO 140 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST
C           ERROR ESTIMATE.
C
        NRMOM = NNLOG(MAXERR)+1
        A1 = ALIST(MAXERR)
        B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL DQC25F(F,A1,B1,DOMEGA,INTEGR,NRMOM,MAXP1,0,
     1  AREA1,ERROR1,NEV,RESABS,DEFAB1,MOMCOM,CHEBMO)
        NEVAL = NEVAL+NEV
        CALL DQC25F(F,A2,B2,DOMEGA,INTEGR,NRMOM,MAXP1,1,
     1  AREA2,ERROR2,NEV,RESABS,DEFAB2,MOMCOM,CHEBMO)
        NEVAL = NEVAL+NEV
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 25
        IF(ABS(RLIST(MAXERR)-AREA12).GT.0.1D-04*ABS(AREA12)
     1  .OR.ERRO12.LT.0.99D+00*ERRMAX) GO TO 20
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   20   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   25   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        NNLOG(MAXERR) = NRMOM
        NNLOG(LAST) = NRMOM
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
C           AT A POINT OF THE INTEGRATION RANGE.
C
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1D+01+0.1D+03*EPMACH)
     1  *(ABS(A2)+0.1D+04*UFLOW)) IER = 4
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
        IF(ERROR2.GT.ERROR1) GO TO 30
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 40
   30   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BISECTED NEXT).
C
   40   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C ***JUMP OUT OF DO-LOOP
      IF(ERRSUM.LE.ERRBND) GO TO 170
      IF(IER.NE.0) GO TO 150
        IF(LAST.EQ.2.AND.EXTALL) GO TO 120
        IF(NOEXT) GO TO 140
        IF(.NOT.EXTALL) GO TO 50
        ERLARG = ERLARG-ERLAST
        IF(ABS(B1-A1).GT.SMALL) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 70
C
C           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
C           SMALLEST INTERVAL.
C
   50   WIDTH = ABS(BLIST(MAXERR)-ALIST(MAXERR))
        IF(WIDTH.GT.SMALL) GO TO 140
        IF(EXTALL) GO TO 60
C
C           TEST WHETHER WE CAN START WITH THE EXTRAPOLATION PROCEDURE
C           (WE DO THIS IF WE INTEGRATE OVER THE NEXT INTERVAL WITH
C           USE OF A GAUSS-KRONROD RULE - SEE SUBROUTINE DQC25F).
C
        SMALL = SMALL*0.5D+00
        IF(0.25D+00*WIDTH*DOMEGA.GT.0.2D+01) GO TO 140
        EXTALL = .TRUE.
        GO TO 130
   60   EXTRAP = .TRUE.
        NRMAX = 2
   70   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 90
C
C           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
C           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER
C           THE LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
C
        JUPBND = LAST
        IF (LAST.GT.(LIMIT/2+2)) JUPBND = LIMIT+3-LAST
        ID = NRMAX
        DO 80 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
          IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 140
          NRMAX = NRMAX+1
   80   CONTINUE
C
C           PERFORM EXTRAPOLATION.
C
   90   NUMRL2 = NUMRL2+1
        RLIST2(NUMRL2) = AREA
        IF(NUMRL2.LT.3) GO TO 110
        CALL DQELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
        KTMIN = KTMIN+1
        IF(KTMIN.GT.5.AND.ABSERR.LT.0.1D-02*ERRSUM) IER = 5
        IF(ABSEPS.GE.ABSERR) GO TO 100
        KTMIN = 0
        ABSERR = ABSEPS
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = MAX(EPSABS,EPSREL*ABS(RESEPS))
C ***JUMP OUT OF DO-LOOP
        IF(ABSERR.LE.ERTEST) GO TO 150
C
C           PREPARE BISECTION OF THE SMALLEST INTERVAL.
C
  100   IF(NUMRL2.EQ.1) NOEXT = .TRUE.
        IF(IER.EQ.5) GO TO 150
  110   MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        SMALL = SMALL*0.5D+00
        ERLARG = ERRSUM
        GO TO 140
  120   SMALL = SMALL*0.5D+00
        NUMRL2 = NUMRL2+1
        RLIST2(NUMRL2) = AREA
  130   ERTEST = ERRBND
        ERLARG = ERRSUM
  140 CONTINUE
C
C           SET THE FINAL RESULT.
C           ---------------------
C
  150 IF(ABSERR.EQ.OFLOW.OR.NRES.EQ.0) GO TO 170
      IF(IER+IERRO.EQ.0) GO TO 165
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0D+00.AND.AREA.NE.0.0D+00) GO TO 160
      IF(ABSERR.GT.ERRSUM) GO TO 170
      IF(AREA.EQ.0.0D+00) GO TO 190
      GO TO 165
  160 IF(ABSERR/ABS(RESULT).GT.ERRSUM/ABS(AREA)) GO TO 170
C
C           TEST ON DIVERGENCE.
C
  165 IF(KSGN.EQ.(-1).AND.MAX(ABS(RESULT),ABS(AREA)).LE.
     1 DEFABS*0.1D-01) GO TO 190
      IF(0.1D-01.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1D+03
     1 .OR.ERRSUM.GE.ABS(AREA)) IER = 6
      GO TO 190
C
C           COMPUTE GLOBAL INTEGRAL SUM.
C
  170 RESULT = 0.0D+00
      DO 180 K=1,LAST
        RESULT = RESULT+RLIST(K)
  180 CONTINUE
      ABSERR = ERRSUM
  190 IF (IER.GT.2) IER=IER-1
  200 IF (INTEGR.EQ.2.AND.OMEGA.LT.0.0D+00) RESULT=-RESULT
  999 RETURN
      END
