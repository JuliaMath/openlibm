*DECK DQAWSE
      SUBROUTINE DQAWSE (F, A, B, ALFA, BETA, INTEGR, EPSABS, EPSREL,
     +   LIMIT, RESULT, ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST,
     +   IORD, LAST)
C***BEGIN PROLOGUE  DQAWSE
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral I = Integral of F*W over (A,B),
C            (where W shows a singular behaviour at the end points,
C            see parameter INTEGR).
C            Hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A1
C***TYPE      DOUBLE PRECISION (QAWSE-S, DQAWSE-D)
C***KEYWORDS  ALGEBRAIC-LOGARITHMIC END POINT SINGULARITIES,
C             AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD, QUADPACK,
C             QUADRATURE, SPECIAL-PURPOSE
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
C                     Parameter in the WEIGHT function, ALFA.GT.(-1)
C                     If ALFA.LE.(-1), the routine will end with
C                     IER = 6.
C
C            BETA   - Double precision
C                     Parameter in the WEIGHT function, BETA.GT.(-1)
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
C            LIMIT  - Integer
C                     Gives an upper bound on the number of subintervals
C                     in the partition of (A,B), LIMIT.GE.2
C                     If LIMIT.LT.2, the routine will end with IER = 6.
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
C                             the estimates for the integral and error
C                             are less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                         = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT. However, if this yields no
C                             improvement, it is advised to analyze the
C                             integrand in order to determine the
C                             integration difficulties which prevent the
C                             requested tolerance from being achieved.
C                             In case of a jump DISCONTINUITY or a local
C                             SINGULARITY of algebraico-logarithmic type
C                             at one or more interior points of the
C                             integration range, one should proceed by
C                             splitting up the interval at these
C                             points and calling the integrator on the
C                             subranges.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 6 The input is invalid, because
C                             B.LE.A or ALFA.LE.(-1) or BETA.LE.(-1), or
C                             INTEGR.LT.1 or INTEGR.GT.4, or
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                             or LIMIT.LT.2.
C                             RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1),
C                             IORD(1) and LAST are set to zero. ALIST(1)
C                             and BLIST(1) are set to A and B
C                             respectively.
C
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
C                     of which are pointers to the error
C                     estimates over the subintervals, so that
C                     ELIST(IORD(1)), ..., ELIST(IORD(K)) with K = LAST
C                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
C                     otherwise form a decreasing sequence
C
C            LAST   - Integer
C                     Number of subintervals actually produced in
C                     the subdivision process
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DQC25S, DQMOMO, DQPSRT
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQAWSE
C
      DOUBLE PRECISION A,ABSERR,ALFA,ALIST,AREA,AREA1,AREA12,AREA2,A1,
     1  A2,B,BETA,BLIST,B1,B2,CENTRE,D1MACH,ELIST,EPMACH,
     2  EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERRO12,ERROR2,ERRSUM,F,
     3  RESAS1,RESAS2,RESULT,RG,RH,RI,RJ,RLIST,UFLOW
      INTEGER IER,INTEGR,IORD,IROFF1,IROFF2,K,LAST,LIMIT,MAXERR,NEV,
     1  NEVAL,NRMAX
C
      EXTERNAL F
C
      DIMENSION ALIST(*),BLIST(*),RLIST(*),ELIST(*),
     1  IORD(*),RI(25),RJ(25),RH(25),RG(25)
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
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
C                       ERROR ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C
C
C            MACHINE DEPENDENT CONSTANTS
C            ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQAWSE
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
      IER = 6
      NEVAL = 0
      LAST = 0
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF (B.LE.A.OR.(EPSABS.EQ.0.0D+00.AND.
     1    EPSREL.LT.MAX(0.5D+02*EPMACH,0.5D-28)).OR.ALFA.LE.(-0.1D+01)
     2    .OR.BETA.LE.(-0.1D+01).OR.INTEGR.LT.1.OR.INTEGR.GT.4.OR.
     3    LIMIT.LT.2) GO TO 999
      IER = 0
C
C           COMPUTE THE MODIFIED CHEBYSHEV MOMENTS.
C
      CALL DQMOMO(ALFA,BETA,RI,RJ,RG,RH,INTEGR)
C
C           INTEGRATE OVER THE INTERVALS (A,(A+B)/2) AND ((A+B)/2,B).
C
      CENTRE = 0.5D+00*(B+A)
      CALL DQC25S(F,A,B,A,CENTRE,ALFA,BETA,RI,RJ,RG,RH,AREA1,
     1  ERROR1,RESAS1,INTEGR,NEV)
      NEVAL = NEV
      CALL DQC25S(F,A,B,CENTRE,B,ALFA,BETA,RI,RJ,RG,RH,AREA2,
     1  ERROR2,RESAS2,INTEGR,NEV)
      LAST = 2
      NEVAL = NEVAL+NEV
      RESULT = AREA1+AREA2
      ABSERR = ERROR1+ERROR2
C
C           TEST ON ACCURACY.
C
      ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
C
C           INITIALIZATION
C           --------------
C
      IF(ERROR2.GT.ERROR1) GO TO 10
      ALIST(1) = A
      ALIST(2) = CENTRE
      BLIST(1) = CENTRE
      BLIST(2) = B
      RLIST(1) = AREA1
      RLIST(2) = AREA2
      ELIST(1) = ERROR1
      ELIST(2) = ERROR2
      GO TO 20
   10 ALIST(1) = CENTRE
      ALIST(2) = A
      BLIST(1) = B
      BLIST(2) = CENTRE
      RLIST(1) = AREA2
      RLIST(2) = AREA1
      ELIST(1) = ERROR2
      ELIST(2) = ERROR1
   20 IORD(1) = 1
      IORD(2) = 2
      IF(LIMIT.EQ.2) IER = 1
      IF(ABSERR.LE.ERRBND.OR.IER.EQ.1) GO TO 999
      ERRMAX = ELIST(1)
      MAXERR = 1
      NRMAX = 1
      AREA = RESULT
      ERRSUM = ABSERR
      IROFF1 = 0
      IROFF2 = 0
C
C            MAIN DO-LOOP
C            ------------
C
      DO 60 LAST = 3,LIMIT
C
C           BISECT THE SUBINTERVAL WITH LARGEST ERROR ESTIMATE.
C
        A1 = ALIST(MAXERR)
        B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
C
        CALL DQC25S(F,A,B,A1,B1,ALFA,BETA,RI,RJ,RG,RH,AREA1,
     1  ERROR1,RESAS1,INTEGR,NEV)
        NEVAL = NEVAL+NEV
        CALL DQC25S(F,A,B,A2,B2,ALFA,BETA,RI,RJ,RG,RH,AREA2,
     1  ERROR2,RESAS2,INTEGR,NEV)
        NEVAL = NEVAL+NEV
C
C           IMPROVE PREVIOUS APPROXIMATIONS INTEGRAL AND ERROR
C           AND TEST FOR ACCURACY.
C
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(A.EQ.A1.OR.B.EQ.B2) GO TO 30
        IF(RESAS1.EQ.ERROR1.OR.RESAS2.EQ.ERROR2) GO TO 30
C
C           TEST FOR ROUNDOFF ERROR.
C
        IF(ABS(RLIST(MAXERR)-AREA12).LT.0.1D-04*ABS(AREA12)
     1  .AND.ERRO12.GE.0.99D+00*ERRMAX) IROFF1 = IROFF1+1
        IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF2 = IROFF2+1
   30   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
C
C           TEST ON ACCURACY.
C
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
        IF(ERRSUM.LE.ERRBND) GO TO 35
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF INTERVAL
C           BISECTIONS EXCEEDS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C
C           SET ERROR FLAG IN THE CASE OF ROUNDOFF ERROR.
C
        IF(IROFF1.GE.6.OR.IROFF2.GE.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT INTERIOR POINTS OF INTEGRATION RANGE.
C
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1D+01+0.1D+03*EPMACH)*
     1  (ABS(A2)+0.1D+04*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
   35   IF(ERROR2.GT.ERROR1) GO TO 40
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 50
   40   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
   50   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C ***JUMP OUT OF DO-LOOP
        IF (IER.NE.0.OR.ERRSUM.LE.ERRBND) GO TO 70
   60 CONTINUE
C
C           COMPUTE FINAL RESULT.
C           ---------------------
C
   70 RESULT = 0.0D+00
      DO 80 K=1,LAST
        RESULT = RESULT+RLIST(K)
   80 CONTINUE
      ABSERR = ERRSUM
  999 RETURN
      END
