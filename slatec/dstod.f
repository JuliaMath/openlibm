*DECK DSTOD
      SUBROUTINE DSTOD (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, WM, IWM,
     +   DF, DJAC, RPAR, IPAR)
C***BEGIN PROLOGUE  DSTOD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (STOD-S, DSTOD-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   DSTOD integrates a system of first order odes over one step in the
C   integrator package DDEBDF.
C ----------------------------------------------------------------------
C DSTOD  performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C Note.. DSTOD  is independent of the value of the iteration method
C indicator MITER, when this is .NE. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTOD  is done with the following variables..
C
C Y      = An array of length .GE. N used as the Y argument in
C          all calls to DF and DJAC.
C NEQ    = Integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to DF and DJAC.
C YH     = An NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate
C          J-th derivative of Y(I), scaled by H**J/FACTORIAL(J)
C          (J = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = A constant integer .GE. N, the first dimension of YH.
C YH1    = A one-dimensional array occupying the same space as YH.
C EWT    = An array of N elements with which the estimated local
C          errors in YH are compared.
C SAVF   = An array of working storage, of length N.
C ACOR   = A work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(I) contains
C          the estimated one-step local error in Y(I).
C WM,IWM = DOUBLE PRECISION and INTEGER work arrays associated with
C          matrix operations in chord iteration (MITER .NE. 0).
C DPJAC   = Name of routine to evaluate and preprocess Jacobian matrix
C          if a chord method is being used.
C DSLVS   = Name of routine to solve linear system in chord iteration.
C H      = The step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = The minimum absolute value of the step size H to be used.
C HMXI   = Inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = The independent variable. TN is updated on each step taken.
C JSTART = An integer used for input only, with the following
C          values and meanings..
C               0  Perform the first step.
C           .GT.0  Take a new step continuing from the last.
C              -1  Take the next step with a new value of H, MAXORD,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  Take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings..
C               0  The step was successful.
C              -1  The requested error could not be achieved.
C              -2  Corrector convergence could not be achieved.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = The maximum order of integration method to be allowed.
C METH/MITER = The method flags.  See description in driver.
C N      = The number of first-order differential equations.
C ----------------------------------------------------------------------
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  DCFOD, DPJAC, DSLVS, DVNRMS
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C   920422  Changed DIMENSION statement.  (WRB)
C***END PROLOGUE  DSTOD
C
      INTEGER I, I1, IALTH, IER, IOD, IOWND, IPAR, IPUP, IREDO, IRET,
     1      IWM, J, JB, JSTART, KFLAG, KSTEPS, L, LMAX, M, MAXORD,
     2      MEO, METH, MITER, N, NCF, NEQ, NEWQ, NFE, NJE, NQ, NQNYH,
     3      NQU, NST, NSTEPJ, NYH
      DOUBLE PRECISION ACOR, CONIT, CRATE, DCON, DDN,
     1      DEL, DELP, DSM, DUP, DVNRMS, EL, EL0, ELCO,
     2      EWT, EXDN, EXSM, EXUP, H, HMIN, HMXI, HOLD, HU, R, RC,
     3      RH, RHDN, RHSM, RHUP, RMAX, ROWND, RPAR, SAVF, TESCO,
     4      TN, TOLD, UROUND, WM, Y, YH, YH1
      EXTERNAL DF, DJAC
C
      DIMENSION Y(*),YH(NYH,*),YH1(*),EWT(*),SAVF(*),ACOR(*),WM(*),
     1          IWM(*),RPAR(*),IPAR(*)
      COMMON /DDEBD1/ ROWND,CONIT,CRATE,EL(13),ELCO(13,12),HOLD,RC,RMAX,
     1                TESCO(3,12),EL0,H,HMIN,HMXI,HU,TN,UROUND,IOWND(7),
     2                KSTEPS,IOD(6),IALTH,IPUP,LMAX,MEO,NQNYH,NSTEPJ,
     3                IER,JSTART,KFLAG,L,METH,MITER,MAXORD,N,NQ,NST,NFE,
     4                NJE,NQU
C
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 690
C        BEGIN BLOCK PERMITTING ...EXITS TO 60
C***FIRST EXECUTABLE STATEMENT  DSTOD
            KFLAG = 0
            TOLD = TN
            NCF = 0
            IF (JSTART .GT. 0) GO TO 160
            IF (JSTART .EQ. -1) GO TO 10
               IF (JSTART .EQ. -2) GO TO 90
C              ---------------------------------------------------------
C               ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER
C               VARIABLES ARE INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY
C               WHICH H CAN BE INCREASED IN A SINGLE STEP.  IT IS
C               INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL INITIAL H,
C               BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE OCCURS
C               (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT
C               2 FOR THE NEXT INCREASE.
C              ---------------------------------------------------------
               LMAX = MAXORD + 1
               NQ = 1
               L = 2
               IALTH = 2
               RMAX = 10000.0D0
               RC = 0.0D0
               EL0 = 1.0D0
               CRATE = 0.7D0
               DELP = 0.0D0
               HOLD = H
               MEO = METH
               NSTEPJ = 0
               IRET = 3
            GO TO 50
   10       CONTINUE
C              BEGIN BLOCK PERMITTING ...EXITS TO 30
C                 ------------------------------------------------------
C                  THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN
C                  JSTART = -1.  IPUP IS SET TO MITER TO FORCE A MATRIX
C                  UPDATE.  IF AN ORDER INCREASE IS ABOUT TO BE
C                  CONSIDERED (IALTH = 1), IALTH IS RESET TO 2 TO
C                  POSTPONE CONSIDERATION ONE MORE STEP.  IF THE CALLER
C                  HAS CHANGED METH, DCFOD  IS CALLED TO RESET THE
C                  COEFFICIENTS OF THE METHOD.  IF THE CALLER HAS
C                  CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT
C                  ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN
C                  ACCORDINGLY.  IF H IS TO BE CHANGED, YH MUST BE
C                  RESCALED.  IF H OR METH IS BEING CHANGED, IALTH IS
C                  RESET TO L = NQ + 1 TO PREVENT FURTHER CHANGES IN H
C                  FOR THAT MANY STEPS.
C                 ------------------------------------------------------
                  IPUP = MITER
                  LMAX = MAXORD + 1
                  IF (IALTH .EQ. 1) IALTH = 2
                  IF (METH .EQ. MEO) GO TO 20
                     CALL DCFOD(METH,ELCO,TESCO)
                     MEO = METH
C              ......EXIT
                     IF (NQ .GT. MAXORD) GO TO 30
                     IALTH = L
                     IRET = 1
C        ............EXIT
                     GO TO 60
   20             CONTINUE
                  IF (NQ .LE. MAXORD) GO TO 90
   30          CONTINUE
               NQ = MAXORD
               L = LMAX
               DO 40 I = 1, L
                  EL(I) = ELCO(I,NQ)
   40          CONTINUE
               NQNYH = NQ*NYH
               RC = RC*EL(1)/EL0
               EL0 = EL(1)
               CONIT = 0.5D0/(NQ+2)
               DDN = DVNRMS(N,SAVF,EWT)/TESCO(1,L)
               EXDN = 1.0D0/L
               RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
               RH = MIN(RHDN,1.0D0)
               IREDO = 3
               IF (H .EQ. HOLD) GO TO 660
               RH = MIN(RH,ABS(H/HOLD))
               H = HOLD
               GO TO 100
   50       CONTINUE
C           ------------------------------------------------------------
C            DCFOD  IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS
C            FOR THE CURRENT METH.  THEN THE EL VECTOR AND RELATED
C            CONSTANTS ARE RESET WHENEVER THE ORDER NQ IS CHANGED, OR AT
C            THE START OF THE PROBLEM.
C           ------------------------------------------------------------
            CALL DCFOD(METH,ELCO,TESCO)
   60    CONTINUE
   70    CONTINUE
C           BEGIN BLOCK PERMITTING ...EXITS TO 680
               DO 80 I = 1, L
                  EL(I) = ELCO(I,NQ)
   80          CONTINUE
               NQNYH = NQ*NYH
               RC = RC*EL(1)/EL0
               EL0 = EL(1)
               CONIT = 0.5D0/(NQ+2)
               GO TO (90,660,160), IRET
C              ---------------------------------------------------------
C               IF H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST
C               RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH
C               IS SET TO L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT
C               MANY STEPS, UNLESS FORCED BY A CONVERGENCE OR ERROR TEST
C               FAILURE.
C              ---------------------------------------------------------
   90          CONTINUE
               IF (H .EQ. HOLD) GO TO 160
               RH = H/HOLD
               H = HOLD
               IREDO = 3
  100          CONTINUE
  110          CONTINUE
                  RH = MIN(RH,RMAX)
                  RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
                  R = 1.0D0
                  DO 130 J = 2, L
                     R = R*RH
                     DO 120 I = 1, N
                        YH(I,J) = YH(I,J)*R
  120                CONTINUE
  130             CONTINUE
                  H = H*RH
                  RC = RC*RH
                  IALTH = L
                  IF (IREDO .NE. 0) GO TO 150
                     RMAX = 10.0D0
                     R = 1.0D0/TESCO(2,NQU)
                     DO 140 I = 1, N
                        ACOR(I) = ACOR(I)*R
  140                CONTINUE
C     ...............EXIT
                     GO TO 690
  150             CONTINUE
C                 ------------------------------------------------------
C                  THIS SECTION COMPUTES THE PREDICTED VALUES BY
C                  EFFECTIVELY MULTIPLYING THE YH ARRAY BY THE PASCAL
C                  TRIANGLE MATRIX.  RC IS THE RATIO OF NEW TO OLD
C                  VALUES OF THE COEFFICIENT  H*EL(1).  WHEN RC DIFFERS
C                  FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER
C                  TO FORCE DPJAC TO BE CALLED, IF A JACOBIAN IS
C                  INVOLVED.  IN ANY CASE, DPJAC IS CALLED AT LEAST
C                  EVERY 20-TH STEP.
C                 ------------------------------------------------------
  160             CONTINUE
  170             CONTINUE
C                    BEGIN BLOCK PERMITTING ...EXITS TO 610
C                       BEGIN BLOCK PERMITTING ...EXITS TO 490
                           IF (ABS(RC-1.0D0) .GT. 0.3D0) IPUP = MITER
                           IF (NST .GE. NSTEPJ + 20) IPUP = MITER
                           TN = TN + H
                           I1 = NQNYH + 1
                           DO 190 JB = 1, NQ
                              I1 = I1 - NYH
                              DO 180 I = I1, NQNYH
                                 YH1(I) = YH1(I) + YH1(I+NYH)
  180                         CONTINUE
  190                      CONTINUE
                           KSTEPS = KSTEPS + 1
C                          ---------------------------------------------
C                           UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A
C                           CONVERGENCE TEST IS MADE ON THE R.M.S. NORM
C                           OF EACH CORRECTION, WEIGHTED BY THE ERROR
C                           WEIGHT VECTOR EWT.  THE SUM OF THE
C                           CORRECTIONS IS ACCUMULATED IN THE VECTOR
C                           ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE
C                           CORRECTOR LOOP.
C                          ---------------------------------------------
  200                      CONTINUE
                              M = 0
                              DO 210 I = 1, N
                                 Y(I) = YH(I,1)
  210                         CONTINUE
                              CALL DF(TN,Y,SAVF,RPAR,IPAR)
                              NFE = NFE + 1
                              IF (IPUP .LE. 0) GO TO 220
C                                ---------------------------------------
C                                 IF INDICATED, THE MATRIX P = I -
C                                 H*EL(1)*J IS REEVALUATED AND
C                                 PREPROCESSED BEFORE STARTING THE
C                                 CORRECTOR ITERATION.  IPUP IS SET TO 0
C                                 AS AN INDICATOR THAT THIS HAS BEEN
C                                 DONE.
C                                ---------------------------------------
                                 IPUP = 0
                                 RC = 1.0D0
                                 NSTEPJ = NST
                                 CRATE = 0.7D0
                                 CALL DPJAC(NEQ,Y,YH,NYH,EWT,ACOR,SAVF,
     1                                      WM,IWM,DF,DJAC,RPAR,IPAR)
C                          ......EXIT
                                 IF (IER .NE. 0) GO TO 440
  220                         CONTINUE
                              DO 230 I = 1, N
                                 ACOR(I) = 0.0D0
  230                         CONTINUE
  240                         CONTINUE
                                 IF (MITER .NE. 0) GO TO 270
C                                   ------------------------------------
C                                    IN THE CASE OF FUNCTIONAL
C                                    ITERATION, UPDATE Y DIRECTLY FROM
C                                    THE RESULT OF THE LAST FUNCTION
C                                    EVALUATION.
C                                   ------------------------------------
                                    DO 250 I = 1, N
                                       SAVF(I) = H*SAVF(I) - YH(I,2)
                                       Y(I) = SAVF(I) - ACOR(I)
  250                               CONTINUE
                                    DEL = DVNRMS(N,Y,EWT)
                                    DO 260 I = 1, N
                                       Y(I) = YH(I,1) + EL(1)*SAVF(I)
                                       ACOR(I) = SAVF(I)
  260                               CONTINUE
                                 GO TO 300
  270                            CONTINUE
C                                   ------------------------------------
C                                    IN THE CASE OF THE CHORD METHOD,
C                                    COMPUTE THE CORRECTOR ERROR, AND
C                                    SOLVE THE LINEAR SYSTEM WITH THAT
C                                    AS RIGHT-HAND SIDE AND P AS
C                                    COEFFICIENT MATRIX.
C                                   ------------------------------------
                                    DO 280 I = 1, N
                                       Y(I) = H*SAVF(I)
     1                                        - (YH(I,2) + ACOR(I))
  280                               CONTINUE
                                    CALL DSLVS(WM,IWM,Y,SAVF)
C                             ......EXIT
                                    IF (IER .NE. 0) GO TO 430
                                    DEL = DVNRMS(N,Y,EWT)
                                    DO 290 I = 1, N
                                       ACOR(I) = ACOR(I) + Y(I)
                                       Y(I) = YH(I,1) + EL(1)*ACOR(I)
  290                               CONTINUE
  300                            CONTINUE
C                                ---------------------------------------
C                                 TEST FOR CONVERGENCE.  IF M.GT.0, AN
C                                 ESTIMATE OF THE CONVERGENCE RATE
C                                 CONSTANT IS STORED IN CRATE, AND THIS
C                                 IS USED IN THE TEST.
C                                ---------------------------------------
                                 IF (M .NE. 0)
     1                              CRATE = MAX(0.2D0*CRATE,DEL/DELP)
                                 DCON = DEL*MIN(1.0D0,1.5D0*CRATE)
     1                                  /(TESCO(2,NQ)*CONIT)
                                 IF (DCON .GT. 1.0D0) GO TO 420
C                                   ------------------------------------
C                                    THE CORRECTOR HAS CONVERGED.  IPUP
C                                    IS SET TO -1 IF MITER .NE. 0, TO
C                                    SIGNAL THAT THE JACOBIAN INVOLVED
C                                    MAY NEED UPDATING LATER.  THE LOCAL
C                                    ERROR TEST IS MADE AND CONTROL
C                                    PASSES TO STATEMENT 500 IF IT
C                                    FAILS.
C                                   ------------------------------------
                                    IF (MITER .NE. 0) IPUP = -1
                                    IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
                                    IF (M .GT. 0)
     1                                 DSM = DVNRMS(N,ACOR,EWT)
     2                                       /TESCO(2,NQ)
                                    IF (DSM .GT. 1.0D0) GO TO 380
C                                      BEGIN BLOCK
C                                      PERMITTING ...EXITS TO 360
C                                         ------------------------------
C                                          AFTER A SUCCESSFUL STEP,
C                                          UPDATE THE YH ARRAY.
C                                          CONSIDER CHANGING H IF IALTH
C                                          = 1.  OTHERWISE DECREASE
C                                          IALTH BY 1.  IF IALTH IS THEN
C                                          1 AND NQ .LT. MAXORD, THEN
C                                          ACOR IS SAVED FOR USE IN A
C                                          POSSIBLE ORDER INCREASE ON
C                                          THE NEXT STEP.  IF A CHANGE
C                                          IN H IS CONSIDERED, AN
C                                          INCREASE OR DECREASE IN ORDER
C                                          BY ONE IS CONSIDERED ALSO.  A
C                                          CHANGE IN H IS MADE ONLY IF
C                                          IT IS BY A FACTOR OF AT LEAST
C                                          1.1.  IF NOT, IALTH IS SET TO
C                                          3 TO PREVENT TESTING FOR THAT
C                                          MANY STEPS.
C                                         ------------------------------
                                          KFLAG = 0
                                          IREDO = 0
                                          NST = NST + 1
                                          HU = H
                                          NQU = NQ
                                          DO 320 J = 1, L
                                             DO 310 I = 1, N
                                                YH(I,J) = YH(I,J)
     1                                                    + EL(J)
     2                                                      *ACOR(I)
  310                                        CONTINUE
  320                                     CONTINUE
                                          IALTH = IALTH - 1
                                          IF (IALTH .NE. 0) GO TO 340
C                                            ---------------------------
C                                             REGARDLESS OF THE SUCCESS
C                                             OR FAILURE OF THE STEP,
C                                             FACTORS RHDN, RHSM, AND
C                                             RHUP ARE COMPUTED, BY
C                                             WHICH H COULD BE
C                                             MULTIPLIED AT ORDER NQ -
C                                             1, ORDER NQ, OR ORDER NQ +
C                                             1, RESPECTIVELY.  IN THE
C                                             CASE OF FAILURE, RHUP =
C                                             0.0 TO AVOID AN ORDER
C                                             INCREASE.  THE LARGEST OF
C                                             THESE IS DETERMINED AND
C                                             THE NEW ORDER CHOSEN
C                                             ACCORDINGLY.  IF THE ORDER
C                                             IS TO BE INCREASED, WE
C                                             COMPUTE ONE ADDITIONAL
C                                             SCALED DERIVATIVE.
C                                            ---------------------------
                                             RHUP = 0.0D0
C                       .....................EXIT
                                             IF (L .EQ. LMAX) GO TO 490
                                             DO 330 I = 1, N
                                                SAVF(I) = ACOR(I)
     1                                                    - YH(I,LMAX)
  330                                        CONTINUE
                                             DUP = DVNRMS(N,SAVF,EWT)
     1                                             /TESCO(3,NQ)
                                             EXUP = 1.0D0/(L+1)
                                             RHUP = 1.0D0
     1                                              /(1.4D0*DUP**EXUP
     2                                                + 0.0000014D0)
C                       .....................EXIT
                                             GO TO 490
  340                                     CONTINUE
C                                      ...EXIT
                                          IF (IALTH .GT. 1) GO TO 360
C                                      ...EXIT
                                          IF (L .EQ. LMAX) GO TO 360
                                          DO 350 I = 1, N
                                             YH(I,LMAX) = ACOR(I)
  350                                     CONTINUE
  360                                  CONTINUE
                                       R = 1.0D0/TESCO(2,NQU)
                                       DO 370 I = 1, N
                                          ACOR(I) = ACOR(I)*R
  370                                  CONTINUE
C     .................................EXIT
                                       GO TO 690
  380                               CONTINUE
C                                   ------------------------------------
C                                    THE ERROR TEST FAILED.  KFLAG KEEPS
C                                    TRACK OF MULTIPLE FAILURES.
C                                    RESTORE TN AND THE YH ARRAY TO
C                                    THEIR PREVIOUS VALUES, AND PREPARE
C                                    TO TRY THE STEP AGAIN.  COMPUTE THE
C                                    OPTIMUM STEP SIZE FOR THIS OR ONE
C                                    LOWER ORDER.  AFTER 2 OR MORE
C                                    FAILURES, H IS FORCED TO DECREASE
C                                    BY A FACTOR OF 0.2 OR LESS.
C                                   ------------------------------------
                                    KFLAG = KFLAG - 1
                                    TN = TOLD
                                    I1 = NQNYH + 1
                                    DO 400 JB = 1, NQ
                                       I1 = I1 - NYH
                                       DO 390 I = I1, NQNYH
                                          YH1(I) = YH1(I) - YH1(I+NYH)
  390                                  CONTINUE
  400                               CONTINUE
                                    RMAX = 2.0D0
                                    IF (ABS(H) .GT. HMIN*1.00001D0)
     1                                 GO TO 410
C                                      ---------------------------------
C                                       ALL RETURNS ARE MADE THROUGH
C                                       THIS SECTION.  H IS SAVED IN
C                                       HOLD TO ALLOW THE CALLER TO
C                                       CHANGE H ON THE NEXT STEP.
C                                      ---------------------------------
                                       KFLAG = -1
C     .................................EXIT
                                       GO TO 690
  410                               CONTINUE
C                    ...............EXIT
                                    IF (KFLAG .LE. -3) GO TO 610
                                    IREDO = 2
                                    RHUP = 0.0D0
C                       ............EXIT
                                    GO TO 490
  420                            CONTINUE
                                 M = M + 1
C                             ...EXIT
                                 IF (M .EQ. 3) GO TO 430
C                             ...EXIT
                                 IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP)
     1                              GO TO 430
                                 DELP = DEL
                                 CALL DF(TN,Y,SAVF,RPAR,IPAR)
                                 NFE = NFE + 1
                              GO TO 240
  430                         CONTINUE
C                             ------------------------------------------
C                              THE CORRECTOR ITERATION FAILED TO
C                              CONVERGE IN 3 TRIES.  IF MITER .NE. 0 AND
C                              THE JACOBIAN IS OUT OF DATE, DPJAC IS
C                              CALLED FOR THE NEXT TRY.  OTHERWISE THE
C                              YH ARRAY IS RETRACTED TO ITS VALUES
C                              BEFORE PREDICTION, AND H IS REDUCED, IF
C                              POSSIBLE.  IF H CANNOT BE REDUCED OR 10
C                              FAILURES HAVE OCCURRED, EXIT WITH KFLAG =
C                              -2.
C                             ------------------------------------------
C                          ...EXIT
                              IF (IPUP .EQ. 0) GO TO 440
                              IPUP = MITER
                           GO TO 200
  440                      CONTINUE
                           TN = TOLD
                           NCF = NCF + 1
                           RMAX = 2.0D0
                           I1 = NQNYH + 1
                           DO 460 JB = 1, NQ
                              I1 = I1 - NYH
                              DO 450 I = I1, NQNYH
                                 YH1(I) = YH1(I) - YH1(I+NYH)
  450                         CONTINUE
  460                      CONTINUE
                           IF (ABS(H) .GT. HMIN*1.00001D0) GO TO 470
                              KFLAG = -2
C     ........................EXIT
                              GO TO 690
  470                      CONTINUE
                           IF (NCF .NE. 10) GO TO 480
                              KFLAG = -2
C     ........................EXIT
                              GO TO 690
  480                      CONTINUE
                           RH = 0.25D0
                           IPUP = MITER
                           IREDO = 1
C                 .........EXIT
                           GO TO 650
  490                   CONTINUE
                        EXSM = 1.0D0/L
                        RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
                        RHDN = 0.0D0
                        IF (NQ .EQ. 1) GO TO 500
                           DDN = DVNRMS(N,YH(1,L),EWT)/TESCO(1,NQ)
                           EXDN = 1.0D0/NQ
                           RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
  500                   CONTINUE
                        IF (RHSM .GE. RHUP) GO TO 550
                           IF (RHUP .LE. RHDN) GO TO 540
                              NEWQ = L
                              RH = RHUP
                              IF (RH .GE. 1.1D0) GO TO 520
                                 IALTH = 3
                                 R = 1.0D0/TESCO(2,NQU)
                                 DO 510 I = 1, N
                                    ACOR(I) = ACOR(I)*R
  510                            CONTINUE
C     ...........................EXIT
                                 GO TO 690
  520                         CONTINUE
                              R = EL(L)/L
                              DO 530 I = 1, N
                                 YH(I,NEWQ+1) = ACOR(I)*R
  530                         CONTINUE
                              NQ = NEWQ
                              L = NQ + 1
                              IRET = 2
C           ..................EXIT
                              GO TO 680
  540                      CONTINUE
                        GO TO 580
  550                   CONTINUE
                        IF (RHSM .LT. RHDN) GO TO 580
                           NEWQ = NQ
                           RH = RHSM
                           IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0)
     1                        GO TO 560
                              IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C                             ------------------------------------------
C                              IF THERE IS A CHANGE OF ORDER, RESET NQ,
C                              L, AND THE COEFFICIENTS.  IN ANY CASE H
C                              IS RESET ACCORDING TO RH AND THE YH ARRAY
C                              IS RESCALED.  THEN EXIT FROM 680 IF THE
C                              STEP WAS OK, OR REDO THE STEP OTHERWISE.
C                             ------------------------------------------
C                 ............EXIT
                              IF (NEWQ .EQ. NQ) GO TO 650
                              NQ = NEWQ
                              L = NQ + 1
                              IRET = 2
C           ..................EXIT
                              GO TO 680
  560                      CONTINUE
                           IALTH = 3
                           R = 1.0D0/TESCO(2,NQU)
                           DO 570 I = 1, N
                              ACOR(I) = ACOR(I)*R
  570                      CONTINUE
C     .....................EXIT
                           GO TO 690
  580                   CONTINUE
                        NEWQ = NQ - 1
                        RH = RHDN
                        IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
                        IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0) GO TO 590
                           IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C                          ---------------------------------------------
C                           IF THERE IS A CHANGE OF ORDER, RESET NQ, L,
C                           AND THE COEFFICIENTS.  IN ANY CASE H IS
C                           RESET ACCORDING TO RH AND THE YH ARRAY IS
C                           RESCALED.  THEN EXIT FROM 680 IF THE STEP
C                           WAS OK, OR REDO THE STEP OTHERWISE.
C                          ---------------------------------------------
C                 .........EXIT
                           IF (NEWQ .EQ. NQ) GO TO 650
                           NQ = NEWQ
                           L = NQ + 1
                           IRET = 2
C           ...............EXIT
                           GO TO 680
  590                   CONTINUE
                        IALTH = 3
                        R = 1.0D0/TESCO(2,NQU)
                        DO 600 I = 1, N
                           ACOR(I) = ACOR(I)*R
  600                   CONTINUE
C     ..................EXIT
                        GO TO 690
  610                CONTINUE
C                    ---------------------------------------------------
C                     CONTROL REACHES THIS SECTION IF 3 OR MORE FAILURES
C                     HAVE OCCURRED.  IF 10 FAILURES HAVE OCCURRED, EXIT
C                     WITH KFLAG = -1.  IT IS ASSUMED THAT THE
C                     DERIVATIVES THAT HAVE ACCUMULATED IN THE YH ARRAY
C                     HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST
C                     DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO
C                     1.  THEN H IS REDUCED BY A FACTOR OF 10, AND THE
C                     STEP IS RETRIED, UNTIL IT SUCCEEDS OR H REACHES
C                     HMIN.
C                    ---------------------------------------------------
                     IF (KFLAG .NE. -10) GO TO 620
C                       ------------------------------------------------
C                        ALL RETURNS ARE MADE THROUGH THIS SECTION.  H
C                        IS SAVED IN HOLD TO ALLOW THE CALLER TO CHANGE
C                        H ON THE NEXT STEP.
C                       ------------------------------------------------
                        KFLAG = -1
C     ..................EXIT
                        GO TO 690
  620                CONTINUE
                     RH = 0.1D0
                     RH = MAX(HMIN/ABS(H),RH)
                     H = H*RH
                     DO 630 I = 1, N
                        Y(I) = YH(I,1)
  630                CONTINUE
                     CALL DF(TN,Y,SAVF,RPAR,IPAR)
                     NFE = NFE + 1
                     DO 640 I = 1, N
                        YH(I,2) = H*SAVF(I)
  640                CONTINUE
                     IPUP = MITER
                     IALTH = 5
C              ......EXIT
                     IF (NQ .NE. 1) GO TO 670
                  GO TO 170
  650             CONTINUE
  660             CONTINUE
                  RH = MAX(RH,HMIN/ABS(H))
               GO TO 110
  670          CONTINUE
               NQ = 1
               L = 2
               IRET = 3
  680       CONTINUE
         GO TO 70
  690 CONTINUE
      HOLD = H
      JSTART = 1
      RETURN
C     ----------------------- END OF SUBROUTINE DSTOD
C     -----------------------
      END
