*DECK DDSTP
      SUBROUTINE DDSTP (EPS, F, FA, HMAX, IMPL, IERROR, JACOBN, MATDIM,
     8   MAXORD, MINT, MITER, ML, MU, N, NDE, YWT, UROUND, USERS, AVGH,
     8   AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD, NFE, NJE, NQUSED,
     8   NSTEP, T, Y, YH, A, CONVRG, DFDY, EL, FAC, HOLD, IPVT, JSTATE,
     8   JSTEPL, NQ, NWAIT, RC, RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG,
     8   MTRSV, MXRDSV)
C***BEGIN PROLOGUE  DDSTP
C***SUBSIDIARY
C***PURPOSE  DDSTP performs one step of the integration of an initial
C            value problem for a system of ordinary differential
C            equations.
C***LIBRARY   SLATEC (SDRIVE)
C***TYPE      DOUBLE PRECISION (SDSTP-S, DDSTP-D, CDSTP-C)
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             Mail Stop D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  Communication with DDSTP is done with the following variables:
C
C    YH      An N by MAXORD+1 array containing the dependent variables
C              and their scaled derivatives.  MAXORD, the maximum order
C              used, is currently 12 for the Adams methods and 5 for the
C              Gear methods.  YH(I,J+1) contains the J-th derivative of
C              Y(I), scaled by H**J/factorial(J).  Only Y(I),
C              1 .LE. I .LE. N, need be set by the calling program on
C              the first entry.  The YH array should not be altered by
C              the calling program.  When referencing YH as a
C              2-dimensional array, use a column length of N, as this is
C              the value used in DDSTP.
C    DFDY    A block of locations used for partial derivatives if MITER
C              is not 0.  If MITER is 1 or 2 its length must be at least
C              N*N.  If MITER is 4 or 5 its length must be at least
C              (2*ML+MU+1)*N.
C    YWT     An array of N locations used in convergence and error tests
C    SAVE1
C    SAVE2   Arrays of length N used for temporary storage.
C    IPVT    An integer array of length N used by the linear system
C              solvers for the storage of row interchange information.
C    A       A block of locations used to store the matrix A, when using
C              the implicit method.  If IMPL is 1, A is a MATDIM by N
C              array.  If MITER is 1 or 2 MATDIM is N, and if MITER is 4
C              or 5 MATDIM is 2*ML+MU+1.  If IMPL is 2 its length is N.
C              If IMPL is 3, A is a MATDIM by NDE array.
C    JTASK   An integer used on input.
C              It has the following values and meanings:
C                 .EQ. 0  Perform the first step.  This value enables
C                         the subroutine to initialize itself.
C                .GT. 0  Take a new step continuing from the last.
C                         Assumes the last step was successful and
C                         user has not changed any parameters.
C                 .LT. 0  Take a new step with a new value of H and/or
C                         MINT and/or MITER.
C    JSTATE  A completion code with the following meanings:
C                1  The step was successful.
C                2  A solution could not be obtained with H .NE. 0.
C                3  A solution was not obtained in MXTRY attempts.
C                4  For IMPL .NE. 0, the matrix A is singular.
C              On a return with JSTATE .GT. 1, the values of T and
C              the YH array are as of the beginning of the last
C              step, and H is the last step size attempted.
C
C***ROUTINES CALLED  DDCOR, DDCST, DDNTL, DDPSC, DDPST, DDSCL, DNRM2
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDSTP
      EXTERNAL F, JACOBN, FA, USERS
      INTEGER I, IERROR, IMPL, IPVT(*), ISWFLG, ITER, J, JSTATE, JSTEPL,
     8        JTASK, MATDIM, MAXORD, MINT, MITER, ML, MNTOLD, MTROLD,
     8        MTRSV, MU, MXFAIL, MXITER, MXRDSV, MXTRY, N, NDE, NDJSTP,
     8        NFAIL, NFE, NJE, NQ, NQUSED, NSTEP, NSV, NTRY, NWAIT
      DOUBLE PRECISION A(MATDIM,*), AVGH, AVGORD, BIAS1, BIAS2, BIAS3,
     8     BND, CTEST, D, DENOM, DFDY(MATDIM,*), D1, EL(13,12), EPS,
     8     ERDN, ERUP, ETEST, FAC(*), H, HMAX, HN, HOLD, HS, HUSED,
     8     NUMER, RC, RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL, RMNORM,
     8     SAVE1(*), SAVE2(*), DNRM2, T, TOLD, TQ(3,12), TREND, TRSHLD,
     8     UROUND, Y(*), YH(N,*), YWT(*), Y0NRM
      LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH
      PARAMETER(BIAS1 = 1.3D0, BIAS2 = 1.2D0, BIAS3 = 1.4D0, MXFAIL = 3,
     8          MXITER = 3, MXTRY = 50, RCTEST = .3D0, RMFAIL = 2.D0,
     8          RMNORM = 10.D0, TRSHLD = 1.D0)
      PARAMETER (NDJSTP = 10)
      DATA IER /.FALSE./
C***FIRST EXECUTABLE STATEMENT  DDSTP
      NSV = N
      BND = 0.D0
      SWITCH = .FALSE.
      NTRY = 0
      TOLD = T
      NFAIL = 0
      IF (JTASK .LE. 0) THEN
        CALL DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC,
     8              YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH,
     8              RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
        IF (N .EQ. 0) GO TO 440
        IF (H .EQ. 0.D0) GO TO 400
        IF (IER) GO TO 420
      END IF
 100  NTRY = NTRY + 1
      IF (NTRY .GT. MXTRY) GO TO 410
      T = T + H
      CALL DDPSC (1, N, NQ,  YH)
      EVALJC = (((ABS(RC - 1.D0) .GT. RCTEST) .OR.
     8  (NSTEP .GE. JSTEPL + NDJSTP)) .AND. (MITER .NE. 0))
      EVALFA = .NOT. EVALJC
C
 110  ITER = 0
      DO 115 I = 1,N
 115    Y(I) = YH(I,1)
      CALL F (N, T, Y, SAVE2)
      IF (N .EQ. 0) THEN
        JSTATE = 6
        GO TO 430
      END IF
      NFE = NFE + 1
      IF (EVALJC .OR. IER) THEN
        CALL DDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,
     8              MU, N, NDE, NQ, SAVE2, T, USERS, Y, YH, YWT, UROUND,
     8              NFE, NJE,  A, DFDY, FAC, IER, IPVT, SAVE1, ISWFLG,
     8              BND, JSTATE)
        IF (N .EQ. 0) GO TO 430
        IF (IER) GO TO 160
        CONVRG = .FALSE.
        RC = 1.D0
        JSTEPL = NSTEP
      END IF
      DO 125 I = 1,N
 125    SAVE1(I) = 0.D0
C                      Up to MXITER corrector iterations are taken.
C                      Convergence is tested by requiring the r.m.s.
C                      norm of changes to be less than EPS.  The sum of
C                      the corrections is accumulated in the vector
C                      SAVE1(I).  It is approximately equal to the L-th
C                      derivative of Y multiplied by
C                      H**L/(factorial(L-1)*EL(L,NQ)), and is thus
C                      proportional to the actual errors to the lowest
C                      power of H present (H**L).  The YH array is not
C                      altered in the correction loop.  The norm of the
C                      iterate difference is stored in D.  If
C                      ITER .GT. 0, an estimate of the convergence rate
C                      constant is stored in TREND, and this is used in
C                      the convergence test.
C
 130  CALL DDCOR (DFDY, EL, FA, H, IERROR, IMPL, IPVT, MATDIM, MITER,
     8            ML, MU, N, NDE, NQ, T, USERS, Y, YH, YWT,  EVALFA,
     8            SAVE1, SAVE2,  A, D, JSTATE)
        IF (N .EQ. 0) GO TO 430
      IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
        IF (ITER .EQ. 0) THEN
          NUMER = DNRM2(N, SAVE1, 1)
          DO 132 I = 1,N
 132        DFDY(1,I) = SAVE1(I)
          Y0NRM = DNRM2(N, YH, 1)
        ELSE
          DENOM = NUMER
          DO 134 I = 1,N
 134        DFDY(1,I) = SAVE1(I) - DFDY(1,I)
          NUMER = DNRM2(N, DFDY, MATDIM)
          IF (EL(1,NQ)*NUMER .LE. 100.D0*UROUND*Y0NRM) THEN
            IF (RMAX .EQ. RMFAIL) THEN
              SWITCH = .TRUE.
              GO TO 170
            END IF
          END IF
          DO 136 I = 1,N
 136        DFDY(1,I) = SAVE1(I)
          IF (DENOM .NE. 0.D0)
     8    BND = MAX(BND, NUMER/(DENOM*ABS(H)*EL(1,NQ)))
        END IF
      END IF
      IF (ITER .GT. 0) TREND = MAX(.9D0*TREND, D/D1)
      D1 = D
      CTEST = MIN(2.D0*TREND, 1.D0)*D
      IF (CTEST .LE. EPS) GO TO 170
      ITER = ITER + 1
      IF (ITER .LT. MXITER) THEN
        DO 140 I = 1,N
 140      Y(I) = YH(I,1) + EL(1,NQ)*SAVE1(I)
        CALL F (N, T, Y, SAVE2)
        IF (N .EQ. 0) THEN
          JSTATE = 6
          GO TO 430
        END IF
        NFE = NFE + 1
        GO TO 130
      END IF
C                     The corrector iteration failed to converge in
C                     MXITER tries.  If partials are involved but are
C                     not up to date, they are reevaluated for the next
C                     try.  Otherwise the YH array is retracted to its
C                     values before prediction, and H is reduced, if
C                     possible.  If not, a no-convergence exit is taken.
      IF (CONVRG) THEN
        EVALJC = .TRUE.
        EVALFA = .FALSE.
        GO TO 110
      END IF
 160  T = TOLD
      CALL DDPSC (-1, N, NQ,  YH)
      NWAIT = NQ + 2
      IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
      IF (ITER .EQ. 0) THEN
        RH = .3D0
      ELSE
        RH = .9D0*(EPS/CTEST)**(.2D0)
      END IF
      IF (RH*H .EQ. 0.D0) GO TO 400
      CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
      GO TO 100
C                          The corrector has converged.  CONVRG is set
C                          to .TRUE. if partial derivatives were used,
C                          to indicate that they may need updating on
C                          subsequent steps.  The error test is made.
 170  CONVRG = (MITER .NE. 0)
      IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
        DO 180 I = 1,NDE
 180      SAVE2(I) = SAVE1(I)/YWT(I)
      ELSE
        DO 185 I = 1,NDE
 185      SAVE2(I) = SAVE1(I)/MAX(ABS(Y(I)), YWT(I))
      END IF
      ETEST = DNRM2(NDE, SAVE2, 1)/(TQ(2,NQ)*SQRT(DBLE(NDE)))
C
C                           The error test failed.  NFAIL keeps track of
C                           multiple failures.  Restore T and the YH
C                           array to their previous values, and prepare
C                           to try the step again.  Compute the optimum
C                           step size for this or one lower order.
      IF (ETEST .GT. EPS) THEN
        T = TOLD
        CALL DDPSC (-1, N, NQ,  YH)
        NFAIL = NFAIL + 1
        IF (NFAIL .LT. MXFAIL .OR. NQ .EQ. 1) THEN
          IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
          RH2 = 1.D0/(BIAS2*(ETEST/EPS)**(1.D0/(NQ+1)))
          IF (NQ .GT. 1) THEN
            IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
              DO 190 I = 1,NDE
 190            SAVE2(I) = YH(I,NQ+1)/YWT(I)
            ELSE
              DO 195 I = 1,NDE
 195            SAVE2(I) = YH(I,NQ+1)/MAX(ABS(Y(I)), YWT(I))
            END IF
            ERDN = DNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(DBLE(NDE)))
            RH1 = 1.D0/MAX(1.D0, BIAS1*(ERDN/EPS)**(1.D0/NQ))
            IF (RH2 .LT. RH1) THEN
              NQ = NQ - 1
              RC = RC*EL(1,NQ)/EL(1,NQ+1)
              RH = RH1
            ELSE
              RH = RH2
            END IF
          ELSE
            RH = RH2
          END IF
          NWAIT = NQ + 2
          IF (RH*H .EQ. 0.D0) GO TO 400
          CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
          GO TO 100
        END IF
C                Control reaches this section if the error test has
C                failed MXFAIL or more times.  It is assumed that the
C                derivatives that have accumulated in the YH array have
C                errors of the wrong order.  Hence the first derivative
C                is recomputed, the order is set to 1, and the step is
C                retried.
        NFAIL = 0
        JTASK = 2
        DO 215 I = 1,N
 215      Y(I) = YH(I,1)
        CALL DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC,
     8              YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH,
     8              RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
        RMAX = RMNORM
        IF (N .EQ. 0) GO TO 440
        IF (H .EQ. 0.D0) GO TO 400
        IF (IER) GO TO 420
        GO TO 100
      END IF
C                          After a successful step, update the YH array.
      NSTEP = NSTEP + 1
      HUSED = H
      NQUSED = NQ
      AVGH = ((NSTEP-1)*AVGH + H)/NSTEP
      AVGORD = ((NSTEP-1)*AVGORD + NQ)/NSTEP
      DO 230 J = 1,NQ+1
        DO 230 I = 1,N
 230      YH(I,J) = YH(I,J) + EL(J,NQ)*SAVE1(I)
      DO 235 I = 1,N
 235    Y(I) = YH(I,1)
C                                          If ISWFLG is 3, consider
C                                          changing integration methods.
      IF (ISWFLG .EQ. 3) THEN
        IF (BND .NE. 0.D0) THEN
          IF (MINT .EQ. 1 .AND. NQ .LE. 5) THEN
            HN = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.D0/(NQ+1)))
            HN = MIN(HN, 1.D0/(2.D0*EL(1,NQ)*BND))
            HS = ABS(H)/MAX(UROUND,
     8      (ETEST/(EPS*EL(NQ+1,1)))**(1.D0/(NQ+1)))
            IF (HS .GT. 1.2D0*HN) THEN
              MINT = 2
              MNTOLD = MINT
              MITER = MTRSV
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 5)
              RC = 0.D0
              RMAX = RMNORM
              TREND = 1.D0
              CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          ELSE IF (MINT .EQ. 2) THEN
            HS = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.D0/(NQ+1)))
            HN = ABS(H)/MAX(UROUND,
     8      (ETEST*EL(NQ+1,1)/EPS)**(1.D0/(NQ+1)))
            HN = MIN(HN, 1.D0/(2.D0*EL(1,NQ)*BND))
            IF (HN .GE. HS) THEN
              MINT = 1
              MNTOLD = MINT
              MITER = 0
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 12)
              RMAX = RMNORM
              TREND = 1.D0
              CONVRG = .FALSE.
              CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          END IF
        END IF
      END IF
      IF (SWITCH) THEN
        MINT = 2
        MNTOLD = MINT
        MITER = MTRSV
        MTROLD = MITER
        MAXORD = MIN(MXRDSV, 5)
        NQ = MIN(NQ, MAXORD)
        RC = 0.D0
        RMAX = RMNORM
        TREND = 1.D0
        CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
        NWAIT = NQ + 2
      END IF
C                           Consider changing H if NWAIT = 1.  Otherwise
C                           decrease NWAIT by 1.  If NWAIT is then 1 and
C                           NQ.LT.MAXORD, then SAVE1 is saved for use in
C                           a possible order increase on the next step.
C
      IF (JTASK .EQ. 0 .OR. JTASK .EQ. 2) THEN
        RH = 1.D0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.D0/(NQ+1)))
        IF (RH.GT.TRSHLD) CALL DDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
      ELSE IF (NWAIT .GT. 1) THEN
        NWAIT = NWAIT - 1
        IF (NWAIT .EQ. 1 .AND. NQ .LT. MAXORD) THEN
          DO 250 I = 1,NDE
 250        YH(I,MAXORD+1) = SAVE1(I)
        END IF
C             If a change in H is considered, an increase or decrease in
C             order by one is considered also.  A change in H is made
C             only if it is by a factor of at least TRSHLD.  Factors
C             RH1, RH2, and RH3 are computed, by which H could be
C             multiplied at order NQ - 1, order NQ, or order NQ + 1,
C             respectively.  The largest of these is determined and the
C             new order chosen accordingly.  If the order is to be
C             increased, we compute one additional scaled derivative.
C             If there is a change of order, reset NQ and the
C             coefficients.  In any case H is reset according to RH and
C             the YH array is rescaled.
      ELSE
        IF (NQ .EQ. 1) THEN
          RH1 = 0.D0
        ELSE
          IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
            DO 270 I = 1,NDE
 270          SAVE2(I) = YH(I,NQ+1)/YWT(I)
          ELSE
            DO 275 I = 1,NDE
 275          SAVE2(I) = YH(I,NQ+1)/MAX(ABS(Y(I)), YWT(I))
          END IF
          ERDN = DNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(DBLE(NDE)))
          RH1 = 1.D0/MAX(UROUND, BIAS1*(ERDN/EPS)**(1.D0/NQ))
        END IF
        RH2 = 1.D0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.D0/(NQ+1)))
        IF (NQ .EQ. MAXORD) THEN
          RH3 = 0.D0
        ELSE
          IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
            DO 290 I = 1,NDE
 290          SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/YWT(I)
          ELSE
            DO 295 I = 1,NDE
              SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/
     8        MAX(ABS(Y(I)), YWT(I))
 295          CONTINUE
          END IF
          ERUP = DNRM2(NDE, SAVE2, 1)/(TQ(3,NQ)*SQRT(DBLE(NDE)))
          RH3 = 1.D0/MAX(UROUND, BIAS3*(ERUP/EPS)**(1.D0/(NQ+2)))
        END IF
        IF (RH1 .GT. RH2 .AND. RH1 .GE. RH3) THEN
          RH = RH1
          IF (RH .LE. TRSHLD) GO TO 380
          NQ = NQ - 1
          RC = RC*EL(1,NQ)/EL(1,NQ+1)
        ELSE IF (RH2 .GE. RH1 .AND. RH2 .GE. RH3) THEN
          RH = RH2
          IF (RH .LE. TRSHLD) GO TO 380
        ELSE
          RH = RH3
          IF (RH .LE. TRSHLD) GO TO 380
          DO 360 I = 1,N
 360        YH(I,NQ+2) = SAVE1(I)*EL(NQ+1,NQ)/(NQ+1)
          NQ = NQ + 1
          RC = RC*EL(1,NQ)/EL(1,NQ-1)
        END IF
        IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
          IF (BND.NE.0.D0) RH = MIN(RH, 1.D0/(2.D0*EL(1,NQ)*BND*ABS(H)))
        END IF
        CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
        RMAX = RMNORM
 380    NWAIT = NQ + 2
      END IF
C               All returns are made through this section.  H is saved
C               in HOLD to allow the caller to change H on the next step
      JSTATE = 1
      HOLD = H
      RETURN
C
 400  JSTATE = 2
      HOLD = H
      DO 405 I = 1,N
 405    Y(I) = YH(I,1)
      RETURN
C
 410  JSTATE = 3
      HOLD = H
      RETURN
C
 420  JSTATE = 4
      HOLD = H
      RETURN
C
 430  T = TOLD
      CALL DDPSC (-1, NSV, NQ,  YH)
      DO 435 I = 1,NSV
 435    Y(I) = YH(I,1)
 440  HOLD = H
      RETURN
      END
