*DECK DQAWCE
      SUBROUTINE DQAWCE (F, A, B, C, EPSABS, EPSREL, LIMIT, RESULT,
     +   ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
C***BEGIN PROLOGUE  DQAWCE
C***PURPOSE  The routine calculates an approximation result to a
C            CAUCHY PRINCIPAL VALUE I = Integral of F*W over (A,B)
C            (W(X) = 1/(X-C), (C.NE.A, C.NE.B), hopefully satisfying
C            following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I))
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A1, J4
C***TYPE      DOUBLE PRECISION (QAWCE-S, DQAWCE-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, CAUCHY PRINCIPAL VALUE,
C             CLENSHAW-CURTIS METHOD, QUADPACK, QUADRATURE,
C             SPECIAL-PURPOSE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Computation of a CAUCHY PRINCIPAL VALUE
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
C            C      - Double precision
C                     Parameter in the WEIGHT function, C.NE.A, C.NE.B
C                     If C = A OR C = B, the routine will end with
C                     IER = 6.
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
C                     IER.GT.0 Abnormal termination of the routine
C                             the estimates for integral and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more sub-
C                             divisions by increasing the value of
C                             LIMIT. However, if this yields no
C                             improvement it is advised to analyze the
C                             the integrand, in order to determine the
C                             the integration difficulties. If the
C                             position of a local difficulty can be
C                             determined (e.g. SINGULARITY,
C                             DISCONTINUITY within the interval) one
C                             will probably gain from splitting up the
C                             interval at this point and calling
C                             appropriate integrators on the subranges.
C                         = 2 The occurrence of roundoff error is detec-
C                             ted, which prevents the requested
C                             tolerance from being achieved.
C                         = 3 Extremely bad integrand behaviour
C                             occurs at some interior points of
C                             the integration interval.
C                         = 6 The input is invalid, because
C                             C = A or C = B or
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             or LIMIT.LT.1.
C                             RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1),
C                             IORD(1) and LAST are set to zero. ALIST(1)
C                             and BLIST(1) are set to A and B
C                             respectively.
C
C            ALIST   - Double precision
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the left
C                      end points of the subintervals in the partition
C                      of the given integration range (A,B)
C
C            BLIST   - Double precision
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the right
C                      end points of the subintervals in the partition
C                      of the given integration range (A,B)
C
C            RLIST   - Double precision
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the integral
C                      approximations on the subintervals
C
C            ELIST   - Double precision
C                      Vector of dimension LIMIT, the first  LAST
C                      elements of which are the moduli of the absolute
C                      error estimates on the subintervals
C
C            IORD    - Integer
C                      Vector of dimension at least LIMIT, the first K
C                      elements of which are pointers to the error
C                      estimates over the subintervals, so that
C                      ELIST(IORD(1)), ..., ELIST(IORD(K)) with K = LAST
C                      If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
C                      otherwise, form a decreasing sequence
C
C            LAST    - Integer
C                      Number of subintervals actually produced in
C                      the subdivision process
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DQC25C, DQPSRT
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQAWCE
C
      DOUBLE PRECISION A,AA,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2,
     1  B,BB,BLIST,B1,B2,C,D1MACH,ELIST,EPMACH,EPSABS,EPSREL,
     2  ERRBND,ERRMAX,ERROR1,ERRO12,ERROR2,ERRSUM,F,RESULT,RLIST,UFLOW
      INTEGER IER,IORD,IROFF1,IROFF2,K,KRULE,LAST,LIMIT,MAXERR,NEV,
     1  NEVAL,NRMAX
C
      DIMENSION ALIST(*),BLIST(*),RLIST(*),ELIST(*),
     1  IORD(*)
C
      EXTERNAL F
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
C***FIRST EXECUTABLE STATEMENT  DQAWCE
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
      IER = 6
      NEVAL = 0
      LAST = 0
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF (C.EQ.A.OR.C.EQ.B.OR.(EPSABS.LE.0.0D+00.AND.
     1    EPSREL.LT.MAX(0.5D+02*EPMACH,0.5D-28))) GO TO 999
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      AA=A
      BB=B
      IF (A.LE.B) GO TO 10
      AA=B
      BB=A
10    IER=0
      KRULE = 1
      CALL DQC25C(F,AA,BB,C,RESULT,ABSERR,KRULE,NEVAL)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
      ALIST(1) = A
      BLIST(1) = B
C
C           TEST ON ACCURACY
C
      ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
      IF(LIMIT.EQ.1) IER = 1
      IF(ABSERR.LT.MIN(0.1D-01*ABS(RESULT),ERRBND)
     1  .OR.IER.EQ.1) GO TO 70
C
C           INITIALIZATION
C           --------------
C
      ALIST(1) = AA
      BLIST(1) = BB
      RLIST(1) = RESULT
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      NRMAX = 1
      IROFF1 = 0
      IROFF2 = 0
C
C           MAIN DO-LOOP
C           ------------
C
      DO 40 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST
C           ERROR ESTIMATE.
C
        A1 = ALIST(MAXERR)
        B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
        B2 = BLIST(MAXERR)
        IF(C.LE.B1.AND.C.GT.A1) B1 = 0.5D+00*(C+B2)
        IF(C.GT.B1.AND.C.LT.B2) B1 = 0.5D+00*(A1+C)
        A2 = B1
        KRULE = 2
        CALL DQC25C(F,A1,B1,C,AREA1,ERROR1,KRULE,NEV)
        NEVAL = NEVAL+NEV
        CALL DQC25C(F,A2,B2,C,AREA2,ERROR2,KRULE,NEV)
        NEVAL = NEVAL+NEV
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(ABS(RLIST(MAXERR)-AREA12).LT.0.1D-04*ABS(AREA12)
     1    .AND.ERRO12.GE.0.99D+00*ERRMAX.AND.KRULE.EQ.0)
     2    IROFF1 = IROFF1+1
        IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX.AND.KRULE.EQ.0)
     1    IROFF2 = IROFF2+1
        RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
        IF(ERRSUM.LE.ERRBND) GO TO 15
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
        IF(IROFF1.GE.6.AND.IROFF2.GT.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE THAT NUMBER OF INTERVAL
C           BISECTIONS EXCEEDS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE.
C
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1D+01+0.1D+03*EPMACH)
     1    *(ABS(A2)+0.1D+04*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
   15   IF(ERROR2.GT.ERROR1) GO TO 20
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
   30    CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0.OR.ERRSUM.LE.ERRBND) GO TO 50
   40 CONTINUE
C
C           COMPUTE FINAL RESULT.
C           ---------------------
C
   50 RESULT = 0.0D+00
      DO 60 K=1,LAST
        RESULT = RESULT+RLIST(K)
   60 CONTINUE
      ABSERR = ERRSUM
   70 IF (AA.EQ.B) RESULT=-RESULT
  999 RETURN
      END
