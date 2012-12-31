*DECK DERKFS
      SUBROUTINE DERKFS (F, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID, H,
     +   TOLFAC, YP, F1, F2, F3, F4, F5, YS, TOLD, DTSIGN, U26, RER,
     +   INIT, KSTEPS, KOP, IQUIT, STIFF, NONSTF, NTSTEP, NSTIFS, RPAR,
     +   IPAR)
C***BEGIN PROLOGUE  DERKFS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DERKF
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (DERKFS-S, DRKFS-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     Fehlberg Fourth-Fifth order Runge-Kutta Method
C **********************************************************************
C
C     DERKFS integrates a system of first order ordinary differential
C     equations as described in the comments for DERKF .
C
C     The arrays YP,F1,F2,F3,F4,F5,and YS  (of length at least NEQ)
C     appear in the call list for variable dimensioning purposes.
C
C     The variables H,TOLFAC,TOLD,DTSIGN,U26,RER,INIT,KSTEPS,KOP,IQUIT,
C     STIFF,NONSTF,NTSTEP, and NSTIFS are used internally by the code
C     and appear in the call list to eliminate local retention of
C     variables between calls. Accordingly, these variables and the
C     array YP should not be altered.
C     Items of possible interest are
C         H  - An appropriate step size to be used for the next step
C         TOLFAC - Factor of change in the tolerances
C         YP - Derivative of solution vector at T
C         KSTEPS - Counter on the number of steps attempted
C
C **********************************************************************
C
C***SEE ALSO  DERKF
C***ROUTINES CALLED  DEFEHL, HSTART, HVNRM, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800501  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891024  Changed references from VNORM to HVNRM.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls, replace GOTOs with
C           IF-THEN-ELSEs.  (RWC)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DERKFS
C
      LOGICAL HFAILD,OUTPUT,STIFF,NONSTF
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3, XERN4
C
      DIMENSION Y(*),YP(*),F1(*),F2(*),F3(*),F4(*),F5(*),
     1          YS(*),INFO(15),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
C
      EXTERNAL F
C
C.......................................................................
C
C  A FIFTH ORDER METHOD WILL GENERALLY NOT BE CAPABLE OF DELIVERING
C  ACCURACIES NEAR LIMITING PRECISION ON COMPUTERS WITH LONG
C  WORDLENGTHS. TO PROTECT AGAINST LIMITING PRECISION DIFFICULTIES
C  ARISING FROM UNREASONABLE ACCURACY REQUESTS, AN APPROPRIATE
C  TOLERANCE THRESHOLD REMIN IS ASSIGNED FOR THIS METHOD. THIS VALUE
C  SHOULD NOT BE CHANGED ACROSS DIFFERENT MACHINES.
C
      SAVE REMIN, MXSTEP, MXKOP
      DATA REMIN/1.E-12/
C
C.......................................................................
C
C  THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
C  NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MXSTEP, THE COUNTER
C  IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE
C  WORK.
C
      DATA MXSTEP/500/
C
C.......................................................................
C
C  INEFFICIENCY CAUSED BY TOO FREQUENT OUTPUT IS MONITORED BY COUNTING
C  THE NUMBER OF STEP SIZES WHICH ARE SEVERELY SHORTENED DUE SOLELY TO
C  THE CHOICE OF OUTPUT POINTS. WHEN THE NUMBER OF ABUSES EXCEED MXKOP,
C  THE COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE
C  MISUSE OF THE CODE.
C
      DATA MXKOP/100/
C
C.......................................................................
C
C***FIRST EXECUTABLE STATEMENT  DERKFS
      IF (INFO(1) .EQ. 0) THEN
C
C ON THE FIRST CALL , PERFORM INITIALIZATION --
C        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
C        FUNCTION ROUTINE  R1MACH. THE USER MUST MAKE SURE THAT THE
C        VALUES SET IN R1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
C
         U = R1MACH(4)
C                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS
         U26 = 26.*U
         RER = 2.*U+REMIN
C                       -- SET TERMINATION FLAG
         IQUIT = 0
C                       -- SET INITIALIZATION INDICATOR
         INIT = 0
C                       -- SET COUNTER FOR IMPACT OF OUTPUT POINTS
         KOP = 0
C                       -- SET COUNTER FOR ATTEMPTED STEPS
         KSTEPS = 0
C                       -- SET INDICATORS FOR STIFFNESS DETECTION
         STIFF = .FALSE.
         NONSTF = .FALSE.
C                       -- SET STEP COUNTERS FOR STIFFNESS DETECTION
         NTSTEP = 0
         NSTIFS = 0
C                       -- RESET INFO(1) FOR SUBSEQUENT CALLS
         INFO(1) = 1
      ENDIF
C
C.......................................................................
C
C        CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
C
      IF (INFO(1) .NE. 0 .AND. INFO(1) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(1)
         CALL XERMSG ('SLATEC', 'DERKFS',
     *      'IN DERKF, INFO(1) MUST BE SET TO 0 ' //
     *      'FOR THE START OF A NEW PROBLEM, AND MUST BE SET TO 1 ' //
     *      'FOLLOWING AN INTERRUPTED TASK.  YOU ARE ATTEMPTING TO ' //
     *      'CONTINUE THE INTEGRATION ILLEGALLY BY CALLING THE CODE ' //
     *      'WITH  INFO(1) = ' // XERN1, 3, 1)
         IDID = -33
      ENDIF
C
      IF (INFO(2) .NE. 0 .AND. INFO(2) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(2)
         CALL XERMSG ('SLATEC', 'DERKFS',
     *      'IN DERKF, INFO(2) MUST BE 0 OR 1 INDICATING SCALAR ' //
     *      'AND VECTOR ERROR TOLERANCES, RESPECTIVELY.  YOU HAVE ' //
     *      'CALLED THE CODE WITH INFO(2) = ' // XERN1, 4, 1)
         IDID = -33
      ENDIF
C
      IF (INFO(3) .NE. 0 .AND. INFO(3) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(3)
         CALL XERMSG ('SLATEC', 'DERKFS',
     *      'IN DERKF, INFO(3) MUST BE 0 OR 1 INDICATING THE ' //
     *      'OR INTERMEDIATE-OUTPUT MODE OF INTEGRATION, ' //
     *      'RESPECTIVELY.  YOU HAVE CALLED THE CODE ' //
     *      'WITH  INFO(3) = ' // XERN1, 5, 1)
         IDID = -33
      ENDIF
C
      IF (NEQ .LT. 1) THEN
         WRITE (XERN1, '(I8)') NEQ
         CALL XERMSG ('SLATEC', 'DERKFS',
     *      'IN DERKF, THE NUMBER OF EQUATIONS NEQ MUST BE A ' //
     *      'POSITIVE INTEGER.  YOU HAVE CALLED THE ' //
     *      'CODE WITH NEQ = ' // XERN1, 6, 1)
         IDID = -33
      ENDIF
C
      NRTOLP = 0
      NATOLP = 0
      DO 10 K=1,NEQ
         IF (NRTOLP .EQ. 0 .AND. RTOL(K) .LT. 0.D0) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') RTOL(K)
            CALL XERMSG ('SLATEC', 'DERKFS',
     *         'IN DERKF, THE RELATIVE ERROR ' //
     *         'TOLERANCES RTOL MUST BE NON-NEGATIVE.  YOU HAVE ' //
     *         'CALLED THE CODE WITH  RTOL(' // XERN1 // ') = ' //
     *         XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' //
     *         'NO FURTHER CHECKING OF RTOL COMPONENTS IS DONE.', 7, 1)
            IDID = -33
            NRTOLP = 1
         ENDIF
C
         IF (NATOLP .EQ. 0 .AND. ATOL(K) .LT. 0.D0) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') ATOL(K)
            CALL XERMSG ('SLATEC', 'DERKFS',
     *         'IN DERKF, THE ABSOLUTE ERROR ' //
     *         'TOLERANCES ATOL MUST BE NON-NEGATIVE.  YOU HAVE ' //
     *         'CALLED THE CODE WITH  ATOL(' // XERN1 // ') = ' //
     *         XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' //
     *         'NO FURTHER CHECKING OF ATOL COMPONENTS IS DONE.', 8, 1)
            IDID = -33
            NATOLP = 1
         ENDIF
C
         IF (INFO(2) .EQ. 0) GO TO 20
         IF (NATOLP.GT.0 .AND. NRTOLP.GT.0) GO TO 20
   10 CONTINUE
C
C
C     CHECK SOME CONTINUATION POSSIBILITIES
C
   20 IF (INIT .NE. 0) THEN
         IF (T .EQ. TOUT) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DERKFS',
     *         'IN DERKF, YOU HAVE CALLED THE ' //
     *         'CODE WITH  T = TOUT = ' // XERN3 // '$$THIS IS NOT ' //
     *         'ALLOWED ON CONTINUATION CALLS.', 9, 1)
            IDID=-33
         ENDIF
C
         IF (T .NE. TOLD) THEN
            WRITE (XERN3, '(1PE15.6)') TOLD
            WRITE (XERN4, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DERKFS',
     *         'IN DERKF, YOU HAVE CHANGED THE ' //
     *         'VALUE OF T FROM ' // XERN3 // ' TO ' // XERN4 //
     *         '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
            IDID=-33
         ENDIF
C
         IF (INIT .NE. 1) THEN
            IF (DTSIGN*(TOUT-T) .LT. 0.D0) THEN
               WRITE (XERN3, '(1PE15.6)') TOUT
               CALL XERMSG ('SLATEC', 'DERKFS',
     *            'IN DERKF, BY CALLING THE CODE ' //
     *            'WITH TOUT = ' // XERN3 // ' YOU ARE ATTEMPTING ' //
     *            'TO CHANGE THE DIRECTION OF INTEGRATION.$$THIS IS ' //
     *            'NOT ALLOWED WITHOUT RESTARTING.', 11, 1)
               IDID=-33
            ENDIF
         ENDIF
      ENDIF
C
C     INVALID INPUT DETECTED
C
      IF (IDID .EQ. (-33)) THEN
         IF (IQUIT .NE. (-33)) THEN
            IQUIT = -33
            GOTO 909
         ELSE
            CALL XERMSG ('SLATEC', 'DERKFS',
     *         'IN DERKF, INVALID INPUT WAS ' //
     *         'DETECTED ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE ' //
     *         'TO PROCEED BECAUSE YOU HAVE NOT CORRECTED THE ' //
     *         'PROBLEM, SO EXECUTION IS BEING TERMINATED.', 12, 2)
            RETURN
         ENDIF
      ENDIF
C
C.......................................................................
C
C     RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS
C     ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE,
C     THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE
C     RER WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE.
C
      DO 50 K=1,NEQ
        IF (RTOL(K)+ATOL(K) .GT. 0.) GO TO 45
        RTOL(K)=RER
        IDID=-2
   45   IF (INFO(2) .EQ. 0) GO TO 55
   50   CONTINUE
C
   55 IF (IDID .NE. (-2)) GO TO 60
C
C                       RTOL=ATOL=0 ON INPUT, SO RTOL WAS CHANGED TO A
C                                                SMALL POSITIVE VALUE
      TOLFAC=1.
      GO TO 909
C
C     BRANCH ON STATUS OF INITIALIZATION INDICATOR
C            INIT=0 MEANS INITIAL DERIVATIVES AND STARTING STEP SIZE
C                   NOT YET COMPUTED
C            INIT=1 MEANS STARTING STEP SIZE NOT YET COMPUTED
C            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED
C
   60 IF (INIT .EQ. 0) GO TO 65
      IF (INIT .EQ. 1) GO TO 70
      GO TO 80
C
C.......................................................................
C
C     MORE INITIALIZATION --
C                         -- EVALUATE INITIAL DERIVATIVES
C
   65 INIT=1
      A=T
      CALL F(A,Y,YP,RPAR,IPAR)
      IF (T .EQ. TOUT) GO TO 666
C
C                         -- SET SIGN OF INTEGRATION DIRECTION  AND
C                         -- ESTIMATE STARTING STEP SIZE
C
   70 INIT=2
      DTSIGN=SIGN(1.,TOUT-T)
      U=R1MACH(4)
      BIG=SQRT(R1MACH(2))
      UTE=U**0.375
      DY=UTE*HVNRM(Y,NEQ)
      IF (DY .EQ. 0.) DY=UTE
      KTOL=1
      DO 75 K=1,NEQ
        IF (INFO(2) .EQ. 1)  KTOL=K
        TOL=RTOL(KTOL)*ABS(Y(K))+ATOL(KTOL)
        IF (TOL .EQ. 0.) TOL=DY*RTOL(KTOL)
   75   F1(K)=TOL
C
      CALL HSTART (F,NEQ,T,TOUT,Y,YP,F1,4,U,BIG,F2,F3,F4,F5,RPAR,IPAR,H)
C
C.......................................................................
C
C     SET STEP SIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
C     AND SET OUTPUT POINT INDICATOR
C
   80 DT=TOUT-T
      H=SIGN(H,DT)
      OUTPUT= .FALSE.
C
C     TEST TO SEE IF DERKF IS BEING SEVERELY IMPACTED BY TOO MANY
C     OUTPUT POINTS
C
      IF (ABS(H) .GE. 2.*ABS(DT)) KOP=KOP+1
      IF (KOP .LE. MXKOP) GO TO 85
C
C                       UNNECESSARY FREQUENCY OF OUTPUT IS RESTRICTING
C                                                 THE STEP SIZE CHOICE
      IDID=-5
      KOP=0
      GO TO 909
C
   85 IF (ABS(DT) .GT. U26*ABS(T)) GO TO 100
C
C     IF TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND RETURN
C
      DO 90 K=1,NEQ
   90   Y(K)=Y(K)+DT*YP(K)
      A=TOUT
      CALL F(A,Y,YP,RPAR,IPAR)
      KSTEPS=KSTEPS+1
      GO TO 666
C
C **********************************************************************
C **********************************************************************
C     STEP BY STEP INTEGRATION
C
  100 HFAILD= .FALSE.
C
C     TO PROTECT AGAINST IMPOSSIBLE ACCURACY REQUESTS, COMPUTE A
C     TOLERANCE FACTOR BASED ON THE REQUESTED ERROR TOLERANCE AND A
C     LEVEL OF ACCURACY ACHIEVABLE AT LIMITING PRECISION
C
      TOLFAC=0.
      KTOL=1
      DO 125 K=1,NEQ
        IF (INFO(2) .EQ. 1) KTOL=K
        ET=RTOL(KTOL)*ABS(Y(K))+ATOL(KTOL)
        IF (ET .GT. 0.) GO TO 120
        TOLFAC=MAX(TOLFAC,RER/RTOL(KTOL))
        GO TO 125
  120   TOLFAC=MAX(TOLFAC,ABS(Y(K))*(RER/ET))
  125   CONTINUE
      IF (TOLFAC .LE. 1.) GO TO 150
C
C                       REQUESTED ERROR UNATTAINABLE DUE TO LIMITED
C                                               PRECISION AVAILABLE
      TOLFAC=2.*TOLFAC
      IDID=-2
      GO TO 909
C
C     SET SMALLEST ALLOWABLE STEP SIZE
C
  150 HMIN=U26*ABS(T)
C
C     ADJUST STEP SIZE IF NECESSARY TO HIT THE OUTPUT POINT --
C     LOOK AHEAD TWO STEPS TO AVOID DRASTIC CHANGES IN THE STEP SIZE AND
C     THUS LESSEN THE IMPACT OF OUTPUT POINTS ON THE CODE.
C     STRETCH THE STEP SIZE BY, AT MOST, AN AMOUNT EQUAL TO THE
C     SAFETY FACTOR OF 9/10.
C
      DT=TOUT-T
      IF (ABS(DT) .GE. 2.*ABS(H)) GO TO 200
      IF (ABS(DT) .GT. ABS(H)/0.9) GO TO 175
C
C     THE NEXT STEP, IF SUCCESSFUL, WILL COMPLETE THE INTEGRATION TO
C     THE OUTPUT POINT
C
      OUTPUT= .TRUE.
      H=DT
      GO TO 200
C
  175 H=0.5*DT
C
C
C **********************************************************************
C     CORE INTEGRATOR FOR TAKING A SINGLE STEP
C **********************************************************************
C     TO AVOID PROBLEMS WITH ZERO CROSSINGS, RELATIVE ERROR IS MEASURED
C     USING THE AVERAGE OF THE MAGNITUDES OF THE SOLUTION AT THE
C     BEGINNING AND END OF A STEP.
C     THE ERROR ESTIMATE FORMULA HAS BEEN GROUPED TO CONTROL LOSS OF
C     SIGNIFICANCE.
C     LOCAL ERROR ESTIMATES FOR A FIRST ORDER METHOD USING THE SAME
C     STEP SIZE AS THE FEHLBERG METHOD ARE CALCULATED AS PART OF THE
C     TEST FOR STIFFNESS.
C     TO DISTINGUISH THE VARIOUS ARGUMENTS, H IS NOT PERMITTED
C     TO BECOME SMALLER THAN 26 UNITS OF ROUNDOFF IN T.
C     PRACTICAL LIMITS ON THE CHANGE IN THE STEP SIZE ARE ENFORCED TO
C     SMOOTH THE STEP SIZE SELECTION PROCESS AND TO AVOID EXCESSIVE
C     CHATTERING ON PROBLEMS HAVING DISCONTINUITIES.
C     TO PREVENT UNNECESSARY FAILURES, THE CODE USES 9/10 THE STEP SIZE
C     IT ESTIMATES WILL SUCCEED.
C     AFTER A STEP FAILURE, THE STEP SIZE IS NOT ALLOWED TO INCREASE FOR
C     THE NEXT ATTEMPTED STEP. THIS MAKES THE CODE MORE EFFICIENT ON
C     PROBLEMS HAVING DISCONTINUITIES AND MORE EFFECTIVE IN GENERAL
C     SINCE LOCAL EXTRAPOLATION IS BEING USED AND EXTRA CAUTION SEEMS
C     WARRANTED.
C.......................................................................
C
C     MONITOR NUMBER OF STEPS ATTEMPTED
C
  200 IF (KSTEPS .LE. MXSTEP) GO TO 222
C
C                       A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED
      IDID=-1
      KSTEPS=0
      IF (.NOT. STIFF) GO TO 909
C
C                       PROBLEM APPEARS TO BE STIFF
      IDID=-4
      STIFF= .FALSE.
      NONSTF= .FALSE.
      NTSTEP=0
      NSTIFS=0
      GO TO 909
C
C     ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H
C
  222 CALL DEFEHL(F,NEQ,T,Y,H,YP,F1,F2,F3,F4,F5,YS,RPAR,IPAR)
      KSTEPS=KSTEPS+1
C
C.......................................................................
C
C     COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR
C     ESTIMATES.  NOTE THAT RELATIVE ERROR IS MEASURED WITH RESPECT TO
C     THE AVERAGE OF THE MAGNITUDES OF THE SOLUTION AT THE BEGINNING
C     AND END OF THE STEP.
C     LOCAL ERROR ESTIMATES FOR A SPECIAL FIRST ORDER METHOD ARE
C     CALCULATED ONLY WHEN THE STIFFNESS DETECTION IS TURNED ON.
C
      EEOET=0.
      ESTIFF=0.
      KTOL=1
      DO 350 K=1,NEQ
        YAVG=0.5*(ABS(Y(K))+ABS(YS(K)))
        IF (INFO(2) .EQ. 1) KTOL=K
        ET=RTOL(KTOL)*YAVG+ATOL(KTOL)
        IF (ET .GT. 0.) GO TO 325
C
C                       PURE RELATIVE ERROR INAPPROPRIATE WHEN SOLUTION
C                                                              VANISHES
        IDID=-3
        GO TO 909
C
  325   EE=ABS((-2090.*YP(K)+(21970.*F3(K)-15048.*F4(K)))+
     1                        (22528.*F2(K)-27360.*F5(K)))
        IF (STIFF .OR. NONSTF) GO TO 350
        ES=ABS(H*(0.055455*YP(K)-0.035493*F1(K)-0.036571*F2(K)+
     1            0.023107*F3(K)-0.009515*F4(K)+0.003017*F5(K)))
        ESTIFF=MAX(ESTIFF,ES/ET)
  350   EEOET=MAX(EEOET,EE/ET)
C
      ESTTOL=ABS(H)*EEOET/752400.
C
      IF (ESTTOL .LE. 1.) GO TO 500
C
C.......................................................................
C
C     UNSUCCESSFUL STEP
C
      IF (ABS(H) .GT. HMIN) GO TO 400
C
C                       REQUESTED ERROR UNATTAINABLE AT SMALLEST
C                                            ALLOWABLE STEP SIZE
      TOLFAC=1.69*ESTTOL
      IDID=-2
      GO TO 909
C
C                       REDUCE THE STEP SIZE , TRY AGAIN
C                       THE DECREASE IS LIMITED TO A FACTOR OF 1/10
C
  400 HFAILD= .TRUE.
      OUTPUT= .FALSE.
      S=0.1
      IF (ESTTOL .LT. 59049.) S=0.9/ESTTOL**0.2
      H=SIGN(MAX(S*ABS(H),HMIN),H)
      GO TO 200
C
C.......................................................................
C
C     SUCCESSFUL STEP
C                       STORE SOLUTION AT T+H
C                       AND EVALUATE DERIVATIVES THERE
C
  500 T=T+H
      DO 525 K=1,NEQ
  525   Y(K)=YS(K)
      A=T
      CALL F(A,Y,YP,RPAR,IPAR)
C
C                       CHOOSE NEXT STEP SIZE
C                       THE INCREASE IS LIMITED TO A FACTOR OF 5
C                       IF STEP FAILURE HAS JUST OCCURRED, NEXT
C                          STEP SIZE IS NOT ALLOWED TO INCREASE
C
      S=5.
      IF (ESTTOL .GT. 1.889568E-4) S=0.9/ESTTOL**0.2
      IF (HFAILD) S=MIN(S,1.)
      H=SIGN(MAX(S*ABS(H),HMIN),H)
C
C.......................................................................
C
C     CHECK FOR STIFFNESS (IF NOT ALREADY DETECTED)
C
C     IN A SEQUENCE OF 50 SUCCESSFUL STEPS BY THE FEHLBERG METHOD, 25
C     SUCCESSFUL STEPS BY THE FIRST ORDER METHOD INDICATES STIFFNESS
C     AND TURNS THE TEST OFF. IF 26 FAILURES BY THE FIRST ORDER METHOD
C     OCCUR, THE TEST IS TURNED OFF UNTIL THIS SEQUENCE OF 50 STEPS
C     BY THE FEHLBERG METHOD IS COMPLETED.
C
      IF (STIFF) GO TO 600
      NTSTEP=MOD(NTSTEP+1,50)
      IF (NTSTEP .EQ. 1) NONSTF= .FALSE.
      IF (NONSTF) GO TO 600
      IF (ESTIFF .GT. 1.) GO TO 550
C
C                       SUCCESSFUL STEP WITH FIRST ORDER METHOD
      NSTIFS=NSTIFS+1
C                       TURN TEST OFF AFTER 25 INDICATIONS OF STIFFNESS
      IF (NSTIFS .EQ. 25) STIFF= .TRUE.
      GO TO 600
C
C                       UNSUCCESSFUL STEP WITH FIRST ORDER METHOD
  550 IF (NTSTEP-NSTIFS .LE. 25) GO TO 600
C                       TURN STIFFNESS DETECTION OFF FOR THIS BLOCK OF
C                                                          FIFTY STEPS
      NONSTF= .TRUE.
C                       RESET STIFF STEP COUNTER
      NSTIFS=0
C
C **********************************************************************
C     END OF CORE INTEGRATOR
C **********************************************************************
C
C
C     SHOULD WE TAKE ANOTHER STEP
C
  600 IF (OUTPUT) GO TO 666
      IF (INFO(3) .EQ. 0) GO TO 100
C
C **********************************************************************
C **********************************************************************
C
C     INTEGRATION SUCCESSFULLY COMPLETED
C
C                 ONE-STEP MODE
      IDID=1
      TOLD=T
      RETURN
C
C                 INTERVAL MODE
  666 IDID=2
      T=TOUT
      TOLD=T
      RETURN
C
C     INTEGRATION TASK INTERRUPTED
C
  909 INFO(1)=-1
      TOLD=T
      IF (IDID .NE. (-2)) RETURN
C
C                       THE ERROR TOLERANCES ARE INCREASED TO VALUES
C                               WHICH ARE APPROPRIATE FOR CONTINUING
      RTOL(1)=TOLFAC*RTOL(1)
      ATOL(1)=TOLFAC*ATOL(1)
      IF (INFO(2) .EQ. 0) RETURN
      DO 939 K=2,NEQ
        RTOL(K)=TOLFAC*RTOL(K)
  939   ATOL(K)=TOLFAC*ATOL(K)
      RETURN
      END
