*DECK STOD
      SUBROUTINE STOD (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, WM, IWM,
     +   F, JAC, RPAR, IPAR)
C***BEGIN PROLOGUE  STOD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DEBDF
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (STOD-S, DSTOD-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   STOD integrates a system of first order odes over one step in the
C   integrator package DEBDF.
C ----------------------------------------------------------------------
C STOD  performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C Note.. STOD  is independent of the value of the iteration method
C indicator MITER, when this is .NE. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with STOD  is done with the following variables..
C
C Y      = An array of length .GE. n used as the Y argument in
C          all calls to F and JAC.
C NEQ    = Integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C YH     = An NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate
C          J-th derivative of Y(I), scaled by H**J/Factorial(j)
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
C WM,IWM = Real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .NE. 0).
C PJAC   = Name of routine to evaluate and preprocess Jacobian matrix
C          if a chord method is being used.
C SLVS   = Name of routine to solve linear system in chord iteration.
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
C***SEE ALSO  DEBDF
C***ROUTINES CALLED  CFOD, PJAC, SLVS, VNWRMS
C***COMMON BLOCKS    DEBDF1
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C   920422  Changed DIMENSION statement.  (WRB)
C***END PROLOGUE  STOD
      EXTERNAL F, JAC
C
CLLL. OPTIMIZE
      INTEGER NEQ, NYH, IWM, I, I1, IALTH, IER, IOWND, IREDO, IRET,
     1   IPUP, J, JB, JSTART, KFLAG, L, LMAX, M, MAXORD, MEO, METH,
     2   MITER, N, NCF, NEWQ, NFE, NJE, NQ, NQNYH, NQU, NST, NSTEPJ
      REAL Y, YH, YH1, EWT, SAVF, ACOR, WM,
     1   ROWND, CONIT, CRATE, EL, ELCO, HOLD, RC, RMAX, TESCO,
     2   EL0, H, HMIN, HMXI, HU, TN, UROUND,
     3   DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,
     4   R, RH, RHDN, RHSM, RHUP, TOLD, VNWRMS
      DIMENSION         Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   ACOR(*), WM(*), IWM(*), RPAR(*), IPAR(*)
      COMMON /DEBDF1/ ROWND, CONIT, CRATE, EL(13), ELCO(13,12),
     1   HOLD, RC, RMAX, TESCO(3,12),
     2   EL0, H, HMIN, HMXI, HU, TN, UROUND, IOWND(7), KSTEPS, IOD(6),
     3   IALTH, IPUP, LMAX, MEO, NQNYH, NSTEPJ,
     4   IER, JSTART, KFLAG, L, METH, MITER, MAXORD, N, NQ, NST, NFE,
     5   NJE, NQU
C
C
C***FIRST EXECUTABLE STATEMENT  STOD
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
C ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER VARIABLES ARE
C INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE INCREASED
C IN A SINGLE STEP.  IT IS INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL
C INITIAL H, BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE
C OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT 2
C FOR THE NEXT INCREASE.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0E0
      RC = 0.0E0
      EL0 = 1.0E0
      CRATE = 0.7E0
      DELP = 0.0E0
      HOLD = H
      MEO = METH
      NSTEPJ = 0
      IRET = 3
      GO TO 140
C-----------------------------------------------------------------------
C THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN JSTART = -1.
C IPUP IS SET TO MITER TO FORCE A MATRIX UPDATE.
C IF AN ORDER INCREASE IS ABOUT TO BE CONSIDERED (IALTH = 1),
C IALTH IS RESET TO 2 TO POSTPONE CONSIDERATION ONE MORE STEP.
C IF THE CALLER HAS CHANGED METH, CFOD  IS CALLED TO RESET
C THE COEFFICIENTS OF THE METHOD.
C IF THE CALLER HAS CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT
C ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN ACCORDINGLY.
C IF H IS TO BE CHANGED, YH MUST BE RESCALED.
C IF H OR METH IS BEING CHANGED, IALTH IS RESET TO L = NQ + 1
C TO PREVENT FURTHER CHANGES IN H FOR THAT MANY STEPS.
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MEO) GO TO 110
      CALL CFOD  (METH, ELCO, TESCO)
      MEO = METH
      IF (NQ .GT. MAXORD) GO TO 120
      IALTH = L
      IRET = 1
      GO TO 150
 110  IF (NQ .LE. MAXORD) GO TO 160
 120  NQ = MAXORD
      L = LMAX
      DO 125 I = 1,L
 125    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5E0/(NQ+2)
      DDN = VNWRMS (N, SAVF, EWT)/TESCO(1,L)
      EXDN = 1.0E0/L
      RHDN = 1.0E0/(1.3E0*DDN**EXDN + 0.0000013E0)
      RH = MIN(RHDN,1.0E0)
      IREDO = 3
      IF (H .EQ. HOLD) GO TO 170
      RH = MIN(RH,ABS(H/HOLD))
      H = HOLD
      GO TO 175
C-----------------------------------------------------------------------
C CFOD  IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS FOR THE
C CURRENT METH.  THEN THE EL VECTOR AND RELATED CONSTANTS ARE RESET
C WHENEVER THE ORDER NQ IS CHANGED, OR AT THE START OF THE PROBLEM.
C-----------------------------------------------------------------------
 140  CALL CFOD  (METH, ELCO, TESCO)
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5E0/(NQ+2)
      GO TO (160, 170, 200), IRET
C-----------------------------------------------------------------------
C IF H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST
C RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH IS SET TO
C L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT MANY STEPS, UNLESS
C FORCED BY A CONVERGENCE OR ERROR TEST FAILURE.
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0E0,ABS(H)*HMXI*RH)
      R = 1.0E0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 680
C-----------------------------------------------------------------------
C THIS SECTION COMPUTES THE PREDICTED VALUES BY EFFECTIVELY
C MULTIPLYING THE YH ARRAY BY THE PASCAL TRIANGLE MATRIX.
C RC IS THE RATIO OF NEW TO OLD VALUES OF THE COEFFICIENT  H*EL(1).
C WHEN RC DIFFERS FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER
C TO FORCE PJAC TO BE CALLED, IF A JACOBIAN IS INVOLVED.
C IN ANY CASE, PJAC IS CALLED AT LEAST EVERY 20-TH STEP.
C-----------------------------------------------------------------------
 200  IF (ABS(RC-1.0E0) .GT. 0.3E0) IPUP = MITER
      IF (NST .GE. NSTEPJ+20) IPUP = MITER
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
      KSTEPS = KSTEPS + 1
C-----------------------------------------------------------------------
C UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS
C MADE ON THE R.M.S. NORM OF EACH CORRECTION, WEIGHTED BY THE ERROR
C WEIGHT VECTOR EWT.  THE SUM OF THE CORRECTIONS IS ACCUMULATED IN THE
C VECTOR ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE CORRECTOR LOOP.
C-----------------------------------------------------------------------
 220  M = 0
      DO 230 I = 1,N
 230    Y(I) = YH(I,1)
      CALL F (TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C IF INDICATED, THE MATRIX P = I - H*EL(1)*J IS REEVALUATED AND
C PREPROCESSED BEFORE STARTING THE CORRECTOR ITERATION.  IPUP IS SET
C TO 0 AS AN INDICATOR THAT THIS HAS BEEN DONE.
C-----------------------------------------------------------------------
      IPUP = 0
      RC = 1.0E0
      NSTEPJ = NST
      CRATE = 0.7E0
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC,
     1          RPAR, IPAR)
      IF (IER .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0E0
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C IN THE CASE OF FUNCTIONAL ITERATION, UPDATE Y DIRECTLY FROM
C THE RESULT OF THE LAST FUNCTION EVALUATION.
C-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = VNWRMS (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
      GO TO 400
C-----------------------------------------------------------------------
C IN THE CASE OF THE CHORD METHOD, COMPUTE THE CORRECTOR ERROR,
C AND SOLVE THE LINEAR SYSTEM WITH THAT AS RIGHT-HAND SIDE AND
C P AS COEFFICIENT MATRIX.
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      CALL SLVS (WM, IWM, Y, SAVF)
      IF (IER .NE. 0) GO TO 410
      DEL = VNWRMS (N, Y, EWT)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + Y(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
C-----------------------------------------------------------------------
C TEST FOR CONVERGENCE.  IF M.GT.0, AN ESTIMATE OF THE CONVERGENCE
C RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.
C-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = MAX(0.2E0*CRATE,DEL/DELP)
      DCON = DEL*MIN(1.0E0,1.5E0*CRATE)/(TESCO(2,NQ)*CONIT)
      IF (DCON .LE. 1.0E0) GO TO 450
      M = M + 1
      IF (M .EQ. 3) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0E0*DELP) GO TO 410
      DELP = DEL
      CALL F (TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      GO TO 270
C-----------------------------------------------------------------------
C THE CORRECTOR ITERATION FAILED TO CONVERGE IN 3 TRIES.
C IF MITER .NE. 0 AND THE JACOBIAN IS OUT OF DATE, PJAC IS CALLED FOR
C THE NEXT TRY.  OTHERWISE THE YH ARRAY IS RETRACTED TO ITS VALUES
C BEFORE PREDICTION, AND H IS REDUCED, IF POSSIBLE.  IF H CANNOT BE
C REDUCED OR 10 FAILURES HAVE OCCURRED, EXIT WITH KFLAG = -2.
C-----------------------------------------------------------------------
 410  IF (IPUP .EQ. 0) GO TO 430
      IPUP = MITER
      GO TO 220
 430  TN = TOLD
      NCF = NCF + 1
      RMAX = 2.0E0
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (ABS(H) .LE. HMIN*1.00001E0) GO TO 670
      IF (NCF .EQ. 10) GO TO 670
      RH = 0.25E0
      IPUP = MITER
      IREDO = 1
      GO TO 170
C-----------------------------------------------------------------------
C THE CORRECTOR HAS CONVERGED.  IPUP IS SET TO -1 IF MITER .NE. 0,
C TO SIGNAL THAT THE JACOBIAN INVOLVED MAY NEED UPDATING LATER.
C THE LOCAL ERROR TEST IS MADE AND CONTROL PASSES TO STATEMENT 500
C IF IT FAILS.
C-----------------------------------------------------------------------
 450  IF (MITER .NE. 0) IPUP = -1
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = VNWRMS (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0E0) GO TO 500
C-----------------------------------------------------------------------
C AFTER A SUCCESSFUL STEP, UPDATE THE YH ARRAY.
C CONSIDER CHANGING H IF IALTH = 1.  OTHERWISE DECREASE IALTH BY 1.
C IF IALTH IS THEN 1 AND NQ .LT. MAXORD, THEN ACOR IS SAVED FOR
C USE IN A POSSIBLE ORDER INCREASE ON THE NEXT STEP.
C IF A CHANGE IN H IS CONSIDERED, AN INCREASE OR DECREASE IN ORDER
C BY ONE IS CONSIDERED ALSO.  A CHANGE IN H IS MADE ONLY IF IT IS BY A
C FACTOR OF AT LEAST 1.1.  IF NOT, IALTH IS SET TO 3 TO PREVENT
C TESTING FOR THAT MANY STEPS.
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 J = 1,L
        DO 470 I = 1,N
 470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 690
      IF (L .EQ. LMAX) GO TO 690
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 690
C-----------------------------------------------------------------------
C THE ERROR TEST FAILED.  KFLAG KEEPS TRACK OF MULTIPLE FAILURES.
C RESTORE TN AND THE YH ARRAY TO THEIR PREVIOUS VALUES, AND PREPARE
C TO TRY THE STEP AGAIN.  COMPUTE THE OPTIMUM STEP SIZE FOR THIS OR
C ONE LOWER ORDER.  AFTER 2 OR MORE FAILURES, H IS FORCED TO DECREASE
C BY A FACTOR OF 0.2 OR LESS.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0E0
      IF (ABS(H) .LE. HMIN*1.00001E0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0E0
      GO TO 540
C-----------------------------------------------------------------------
C REGARDLESS OF THE SUCCESS OR FAILURE OF THE STEP, FACTORS
C RHDN, RHSM, AND RHUP ARE COMPUTED, BY WHICH H COULD BE MULTIPLIED
C AT ORDER NQ - 1, ORDER NQ, OR ORDER NQ + 1, RESPECTIVELY.
C IN THE CASE OF FAILURE, RHUP = 0.0 TO AVOID AN ORDER INCREASE.
C THE LARGEST OF THESE IS DETERMINED AND THE NEW ORDER CHOSEN
C ACCORDINGLY.  IF THE ORDER IS TO BE INCREASED, WE COMPUTE ONE
C ADDITIONAL SCALED DERIVATIVE.
C-----------------------------------------------------------------------
 520  RHUP = 0.0E0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = VNWRMS (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0E0/(L+1)
      RHUP = 1.0E0/(1.4E0*DUP**EXUP + 0.0000014E0)
 540  EXSM = 1.0E0/L
      RHSM = 1.0E0/(1.2E0*DSM**EXSM + 0.0000012E0)
      RHDN = 0.0E0
      IF (NQ .EQ. 1) GO TO 560
      DDN = VNWRMS (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0E0/NQ
      RHDN = 1.0E0/(1.3E0*DDN**EXDN + 0.0000013E0)
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0E0) RH = 1.0E0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1E0) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 690
 620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1E0)) GO TO 610
      IF (KFLAG .LE. -2) RH = MIN(RH,0.2E0)
C-----------------------------------------------------------------------
C IF THERE IS A CHANGE OF ORDER, RESET NQ, L, AND THE COEFFICIENTS.
C IN ANY CASE H IS RESET ACCORDING TO RH AND THE YH ARRAY IS RESCALED.
C THEN EXIT FROM 680 IF THE STEP WAS OK, OR REDO THE STEP OTHERWISE.
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
C CONTROL REACHES THIS SECTION IF 3 OR MORE FAILURES HAVE OCCURRED.
C IF 10 FAILURES HAVE OCCURRED, EXIT WITH KFLAG = -1.
C IT IS ASSUMED THAT THE DERIVATIVES THAT HAVE ACCUMULATED IN THE
C YH ARRAY HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST
C DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO 1.  THEN
C H IS REDUCED BY A FACTOR OF 10, AND THE STEP IS RETRIED,
C UNTIL IT SUCCEEDS OR H REACHES HMIN.
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1E0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
 645    Y(I) = YH(I,1)
      CALL F (TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
C-----------------------------------------------------------------------
C ALL RETURNS ARE MADE THROUGH THIS SECTION.  H IS SAVED IN HOLD
C TO ALLOW THE CALLER TO CHANGE H ON THE NEXT STEP.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 700
 670  KFLAG = -2
      GO TO 700
 680  RMAX = 10.0E0
 690  R = 1.0E0/TESCO(2,NQU)
      DO 695 I = 1,N
 695    ACOR(I) = ACOR(I)*R
 700  HOLD = H
      JSTART = 1
      RETURN
C----------------------- END OF SUBROUTINE STOD  -----------------------
      END
