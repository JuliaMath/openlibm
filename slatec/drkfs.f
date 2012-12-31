*DECK DRKFS
      SUBROUTINE DRKFS (DF, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID, H,
     +   TOLFAC, YP, F1, F2, F3, F4, F5, YS, TOLD, DTSIGN, U26, RER,
     +   INIT, KSTEPS, KOP, IQUIT, STIFF, NONSTF, NTSTEP, NSTIFS, RPAR,
     +   IPAR)
C***BEGIN PROLOGUE  DRKFS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDERKF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (DERKFS-S, DRKFS-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     Fehlberg Fourth-Fifth Order Runge-Kutta Method
C **********************************************************************
C
C     DRKFS integrates a system of first order ordinary differential
C     equations as described in the comments for DDERKF .
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
C***SEE ALSO  DDERKF
C***ROUTINES CALLED  D1MACH, DFEHL, DHSTRT, DHVNRM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891024  Changed references from DVNORM to DHVNRM.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls, change GOTOs to
C           IF-THEN-ELSEs.  (RWC)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DRKFS
C
      INTEGER IDID, INFO, INIT, IPAR, IQUIT, K, KOP, KSTEPS, KTOL,
     1      MXKOP, MXSTEP, NATOLP, NEQ, NRTOLP, NSTIFS, NTSTEP
      DOUBLE PRECISION A, ATOL, BIG, D1MACH,
     1      DT, DTSIGN, DHVNRM, DY, EE, EEOET, ES, ESTIFF,
     2      ESTTOL, ET, F1, F2, F3, F4, F5, H, HMIN, REMIN, RER, RPAR,
     3      RTOL, S, T, TOL, TOLD, TOLFAC, TOUT, U, U26, UTE, Y, YAVG,
     4      YP, YS
      LOGICAL HFAILD,OUTPUT,STIFF,NONSTF
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3, XERN4
C
      DIMENSION Y(*),YP(*),F1(*),F2(*),F3(*),F4(*),F5(*),
     1          YS(*),INFO(15),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
C
      EXTERNAL DF
C
C     ..................................................................
C
C       A FIFTH ORDER METHOD WILL GENERALLY NOT BE CAPABLE OF DELIVERING
C       ACCURACIES NEAR LIMITING PRECISION ON COMPUTERS WITH LONG
C       WORDLENGTHS. TO PROTECT AGAINST LIMITING PRECISION DIFFICULTIES
C       ARISING FROM UNREASONABLE ACCURACY REQUESTS, AN APPROPRIATE
C       TOLERANCE THRESHOLD REMIN IS ASSIGNED FOR THIS METHOD. THIS
C       VALUE SHOULD NOT BE CHANGED ACROSS DIFFERENT MACHINES.
C
      SAVE REMIN, MXSTEP, MXKOP
      DATA REMIN /1.0D-12/
C
C     ..................................................................
C
C       THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
C       NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MXSTEP, THE
C       COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE
C       EXCESSIVE WORK.
C
      DATA MXSTEP /500/
C
C     ..................................................................
C
C       INEFFICIENCY CAUSED BY TOO FREQUENT OUTPUT IS MONITORED BY
C       COUNTING THE NUMBER OF STEP SIZES WHICH ARE SEVERELY SHORTENED
C       DUE SOLELY TO THE CHOICE OF OUTPUT POINTS. WHEN THE NUMBER OF
C       ABUSES EXCEED MXKOP, THE COUNTER IS RESET TO ZERO AND THE USER
C       IS INFORMED ABOUT POSSIBLE MISUSE OF THE CODE.
C
      DATA MXKOP /100/
C
C     ..................................................................
C
C***FIRST EXECUTABLE STATEMENT  DRKFS
      IF (INFO(1) .EQ. 0) THEN
C
C ON THE FIRST CALL , PERFORM INITIALIZATION --
C        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
C        FUNCTION ROUTINE  D1MACH. THE USER MUST MAKE SURE THAT THE
C        VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
C
         U = D1MACH(4)
C                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS
         U26 = 26.0D0*U
         RER = 2.0D0*U + REMIN
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
         CALL XERMSG ('SLATEC', 'DRKFS',
     *      'IN DDERKF, INFO(1) MUST BE SET TO 0 ' //
     *      'FOR THE START OF A NEW PROBLEM, AND MUST BE SET TO 1 ' //
     *      'FOLLOWING AN INTERRUPTED TASK.  YOU ARE ATTEMPTING TO ' //
     *      'CONTINUE THE INTEGRATION ILLEGALLY BY CALLING THE CODE ' //
     *      'WITH  INFO(1) = ' // XERN1, 3, 1)
         IDID = -33
      ENDIF
C
      IF (INFO(2) .NE. 0 .AND. INFO(2) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(2)
         CALL XERMSG ('SLATEC', 'DRKFS',
     *      'IN DDERKF, INFO(2) MUST BE 0 OR 1 ' //
     *      'INDICATING SCALAR AND VECTOR ERROR TOLERANCES, ' //
     *      'RESPECTIVELY.  YOU HAVE CALLED THE CODE WITH INFO(2) = ' //
     *      XERN1, 4, 1)
         IDID = -33
      ENDIF
C
      IF (INFO(3) .NE. 0 .AND. INFO(3) .NE. 1) THEN
         WRITE (XERN1, '(I8)') INFO(3)
         CALL XERMSG ('SLATEC', 'DRKFS',
     *      'IN DDERKF, INFO(3) MUST BE 0 OR 1 ' //
     *      'INDICATING THE INTERVAL OR INTERMEDIATE-OUTPUT MODE OF ' //
     *      'INTEGRATION, RESPECTIVELY.  YOU HAVE CALLED THE CODE ' //
     *      'WITH  INFO(3) = ' // XERN1, 5, 1)
         IDID = -33
      ENDIF
C
      IF (NEQ .LT. 1) THEN
         WRITE (XERN1, '(I8)') NEQ
         CALL XERMSG ('SLATEC', 'DRKFS',
     *      'IN DDERKF, THE NUMBER OF EQUATIONS ' //
     *      'NEQ MUST BE A POSITIVE INTEGER.  YOU HAVE CALLED THE ' //
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
            CALL XERMSG ('SLATEC', 'DRKFS',
     *         'IN DDERKF, THE RELATIVE ERROR ' //
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
            CALL XERMSG ('SLATEC', 'DRKFS',
     *         'IN DDERKF, THE ABSOLUTE ERROR ' //
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
            CALL XERMSG ('SLATEC', 'DRKFS',
     *         'IN DDERKF, YOU HAVE CALLED THE ' //
     *         'CODE WITH  T = TOUT = ' // XERN3 // '$$THIS IS NOT ' //
     *         'ALLOWED ON CONTINUATION CALLS.', 9, 1)
            IDID=-33
         ENDIF
C
         IF (T .NE. TOLD) THEN
            WRITE (XERN3, '(1PE15.6)') TOLD
            WRITE (XERN4, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DRKFS',
     *         'IN DDERKF, YOU HAVE CHANGED THE ' //
     *         'VALUE OF T FROM ' // XERN3 // ' TO ' // XERN4 //
     *         '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
            IDID=-33
         ENDIF
C
         IF (INIT .NE. 1) THEN
            IF (DTSIGN*(TOUT-T) .LT. 0.D0) THEN
               WRITE (XERN3, '(1PE15.6)') TOUT
               CALL XERMSG ('SLATEC', 'DRKFS',
     *            'IN DDERKF, BY CALLING THE CODE WITH TOUT = ' //
     *            XERN3 // ' YOU ARE ATTEMPTING TO CHANGE THE ' //
     *            'DIRECTION OF INTEGRATION.$$THIS IS NOT ALLOWED ' //
     *            'WITHOUT RESTARTING.', 11, 1)
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
            GOTO 540
         ELSE
            CALL XERMSG ('SLATEC', 'DRKFS',
     *         'IN DDERKF, INVALID INPUT WAS ' //
     *         'DETECTED ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE ' //
     *         'TO PROCEED BECAUSE YOU HAVE NOT CORRECTED THE ' //
     *         'PROBLEM, SO EXECUTION IS BEING TERMINATED.', 12, 2)
            RETURN
         ENDIF
      ENDIF
C
C           ............................................................
C
C                RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND
C                INTERPRETED AS ASKING FOR THE MOST ACCURATE SOLUTION
C                POSSIBLE. IN THIS CASE, THE RELATIVE ERROR TOLERANCE
C                RTOL IS RESET TO THE SMALLEST VALUE RER WHICH IS LIKELY
C                TO BE REASONABLE FOR THIS METHOD AND MACHINE.
C
            DO 190 K = 1, NEQ
               IF (RTOL(K) + ATOL(K) .GT. 0.0D0) GO TO 180
                  RTOL(K) = RER
                  IDID = -2
  180          CONTINUE
C           ...EXIT
               IF (INFO(2) .EQ. 0) GO TO 200
  190       CONTINUE
  200       CONTINUE
C
            IF (IDID .NE. (-2)) GO TO 210
C
C              RTOL=ATOL=0 ON INPUT, SO RTOL WAS CHANGED TO A
C                                       SMALL POSITIVE VALUE
               TOLFAC = 1.0D0
            GO TO 530
  210       CONTINUE
C
C                       BRANCH ON STATUS OF INITIALIZATION INDICATOR
C                              INIT=0 MEANS INITIAL DERIVATIVES AND
C                              STARTING STEP SIZE
C                                     NOT YET COMPUTED
C                              INIT=1 MEANS STARTING STEP SIZE NOT YET
C                              COMPUTED INIT=2 MEANS NO FURTHER
C                              INITIALIZATION REQUIRED
C
                        IF (INIT .EQ. 0) GO TO 220
C                    ......EXIT
                           IF (INIT .EQ. 1) GO TO 240
C                 .........EXIT
                           GO TO 260
  220                   CONTINUE
C
C                       ................................................
C
C                            MORE INITIALIZATION --
C                                                -- EVALUATE INITIAL
C                                                DERIVATIVES
C
                        INIT = 1
                        A = T
                        CALL DF(A,Y,YP,RPAR,IPAR)
                        IF (T .NE. TOUT) GO TO 230
C
C                          INTERVAL MODE
                           IDID = 2
                           T = TOUT
                           TOLD = T
C     .....................EXIT
                           GO TO 560
  230                   CONTINUE
  240                CONTINUE
C
C                    -- SET SIGN OF INTEGRATION DIRECTION  AND
C                    -- ESTIMATE STARTING STEP SIZE
C
                     INIT = 2
                     DTSIGN = SIGN(1.0D0,TOUT-T)
                     U = D1MACH(4)
                     BIG = SQRT(D1MACH(2))
                     UTE = U**0.375D0
                     DY = UTE*DHVNRM(Y,NEQ)
                     IF (DY .EQ. 0.0D0) DY = UTE
                     KTOL = 1
                     DO 250 K = 1, NEQ
                        IF (INFO(2) .EQ. 1) KTOL = K
                        TOL = RTOL(KTOL)*ABS(Y(K)) + ATOL(KTOL)
                        IF (TOL .EQ. 0.0D0) TOL = DY*RTOL(KTOL)
                        F1(K) = TOL
  250                CONTINUE
C
                     CALL DHSTRT(DF,NEQ,T,TOUT,Y,YP,F1,4,U,BIG,F2,F3,F4,
     1                           F5,RPAR,IPAR,H)
  260             CONTINUE
C
C                 ......................................................
C
C                      SET STEP SIZE FOR INTEGRATION IN THE DIRECTION
C                      FROM T TO TOUT AND SET OUTPUT POINT INDICATOR
C
                  DT = TOUT - T
                  H = SIGN(H,DT)
                  OUTPUT = .FALSE.
C
C                 TEST TO SEE IF DDERKF IS BEING SEVERELY IMPACTED BY
C                 TOO MANY OUTPUT POINTS
C
                  IF (ABS(H) .GE. 2.0D0*ABS(DT)) KOP = KOP + 1
                  IF (KOP .LE. MXKOP) GO TO 270
C
C                    UNNECESSARY FREQUENCY OF OUTPUT IS RESTRICTING
C                                              THE STEP SIZE CHOICE
                     IDID = -5
                     KOP = 0
                  GO TO 510
  270             CONTINUE
C
                     IF (ABS(DT) .GT. U26*ABS(T)) GO TO 290
C
C                       IF TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND
C                       RETURN
C
                        DO 280 K = 1, NEQ
                           Y(K) = Y(K) + DT*YP(K)
  280                   CONTINUE
                        A = TOUT
                        CALL DF(A,Y,YP,RPAR,IPAR)
                        KSTEPS = KSTEPS + 1
                     GO TO 500
  290                CONTINUE
C                       BEGIN BLOCK PERMITTING ...EXITS TO 490
C
C                          *********************************************
C                          *********************************************
C                               STEP BY STEP INTEGRATION
C
  300                      CONTINUE
C                             BEGIN BLOCK PERMITTING ...EXITS TO 480
                                 HFAILD = .FALSE.
C
C                                TO PROTECT AGAINST IMPOSSIBLE ACCURACY
C                                REQUESTS, COMPUTE A TOLERANCE FACTOR
C                                BASED ON THE REQUESTED ERROR TOLERANCE
C                                AND A LEVEL OF ACCURACY ACHIEVABLE AT
C                                LIMITING PRECISION
C
                                 TOLFAC = 0.0D0
                                 KTOL = 1
                                 DO 330 K = 1, NEQ
                                    IF (INFO(2) .EQ. 1) KTOL = K
                                    ET = RTOL(KTOL)*ABS(Y(K))
     1                                   + ATOL(KTOL)
                                    IF (ET .GT. 0.0D0) GO TO 310
                                       TOLFAC = MAX(TOLFAC,
     1                                                RER/RTOL(KTOL))
                                    GO TO 320
  310                               CONTINUE
                                       TOLFAC = MAX(TOLFAC,
     1                                                ABS(Y(K))
     2                                                *(RER/ET))
  320                               CONTINUE
  330                            CONTINUE
                                 IF (TOLFAC .LE. 1.0D0) GO TO 340
C
C                          REQUESTED ERROR UNATTAINABLE DUE TO LIMITED
C                                                  PRECISION AVAILABLE
                                    TOLFAC = 2.0D0*TOLFAC
                                    IDID = -2
C              .....................EXIT
                                    GO TO 520
  340                            CONTINUE
C
C                                SET SMALLEST ALLOWABLE STEP SIZE
C
                                 HMIN = U26*ABS(T)
C
C                                ADJUST STEP SIZE IF NECESSARY TO HIT
C                                THE OUTPUT POINT -- LOOK AHEAD TWO
C                                STEPS TO AVOID DRASTIC CHANGES IN THE
C                                STEP SIZE AND THUS LESSEN THE IMPACT OF
C                                OUTPUT POINTS ON THE CODE.  STRETCH THE
C                                STEP SIZE BY, AT MOST, AN AMOUNT EQUAL
C                                TO THE SAFETY FACTOR OF 9/10.
C
                                 DT = TOUT - T
                                 IF (ABS(DT) .GE. 2.0D0*ABS(H))
     1                              GO TO 370
                                    IF (ABS(DT) .GT. ABS(H)/0.9D0)
     1                                 GO TO 350
C
C                                      THE NEXT STEP, IF SUCCESSFUL,
C                                      WILL COMPLETE THE INTEGRATION TO
C                                      THE OUTPUT POINT
C
                                       OUTPUT = .TRUE.
                                       H = DT
                                    GO TO 360
  350                               CONTINUE
C
                                       H = 0.5D0*DT
  360                               CONTINUE
  370                            CONTINUE
C
C
C                                ***************************************
C                                     CORE INTEGRATOR FOR TAKING A
C                                     SINGLE STEP
C                                ***************************************
C                                     TO AVOID PROBLEMS WITH ZERO
C                                     CROSSINGS, RELATIVE ERROR IS
C                                     MEASURED USING THE AVERAGE OF THE
C                                     MAGNITUDES OF THE SOLUTION AT THE
C                                     BEGINNING AND END OF A STEP.
C                                     THE ERROR ESTIMATE FORMULA HAS
C                                     BEEN GROUPED TO CONTROL LOSS OF
C                                     SIGNIFICANCE.
C                                     LOCAL ERROR ESTIMATES FOR A FIRST
C                                     ORDER METHOD USING THE SAME
C                                     STEP SIZE AS THE FEHLBERG METHOD
C                                     ARE CALCULATED AS PART OF THE
C                                     TEST FOR STIFFNESS.
C                                     TO DISTINGUISH THE VARIOUS
C                                     ARGUMENTS, H IS NOT PERMITTED
C                                     TO BECOME SMALLER THAN 26 UNITS OF
C                                     ROUNDOFF IN T.  PRACTICAL LIMITS
C                                     ON THE CHANGE IN THE STEP SIZE ARE
C                                     ENFORCED TO SMOOTH THE STEP SIZE
C                                     SELECTION PROCESS AND TO AVOID
C                                     EXCESSIVE CHATTERING ON PROBLEMS
C                                     HAVING DISCONTINUITIES.  TO
C                                     PREVENT UNNECESSARY FAILURES, THE
C                                     CODE USES 9/10 THE STEP SIZE
C                                     IT ESTIMATES WILL SUCCEED.
C                                     AFTER A STEP FAILURE, THE STEP
C                                     SIZE IS NOT ALLOWED TO INCREASE
C                                     FOR THE NEXT ATTEMPTED STEP. THIS
C                                     MAKES THE CODE MORE EFFICIENT ON
C                                     PROBLEMS HAVING DISCONTINUITIES
C                                     AND MORE EFFECTIVE IN GENERAL
C                                     SINCE LOCAL EXTRAPOLATION IS BEING
C                                     USED AND EXTRA CAUTION SEEMS
C                                     WARRANTED.
C                                .......................................
C
C                                     MONITOR NUMBER OF STEPS ATTEMPTED
C
  380                            CONTINUE
                                    IF (KSTEPS .LE. MXSTEP) GO TO 390
C
C                                      A SIGNIFICANT AMOUNT OF WORK HAS
C                                      BEEN EXPENDED
                                       IDID = -1
                                       KSTEPS = 0
C              ........................EXIT
                                       IF (.NOT.STIFF) GO TO 520
C
C                                      PROBLEM APPEARS TO BE STIFF
                                       IDID = -4
                                       STIFF = .FALSE.
                                       NONSTF = .FALSE.
                                       NTSTEP = 0
                                       NSTIFS = 0
C              ........................EXIT
                                       GO TO 520
  390                               CONTINUE
C
C                                   ADVANCE AN APPROXIMATE SOLUTION OVER
C                                   ONE STEP OF LENGTH H
C
                                    CALL DFEHL(DF,NEQ,T,Y,H,YP,F1,F2,F3,
     1                                         F4,F5,YS,RPAR,IPAR)
                                    KSTEPS = KSTEPS + 1
C
C                                   ....................................
C
C                                        COMPUTE AND TEST ALLOWABLE
C                                        TOLERANCES VERSUS LOCAL ERROR
C                                        ESTIMATES.  NOTE THAT RELATIVE
C                                        ERROR IS MEASURED WITH RESPECT
C                                        TO THE AVERAGE OF THE
C                                        MAGNITUDES OF THE SOLUTION AT
C                                        THE BEGINNING AND END OF THE
C                                        STEP.  LOCAL ERROR ESTIMATES
C                                        FOR A SPECIAL FIRST ORDER
C                                        METHOD ARE CALCULATED ONLY WHEN
C                                        THE STIFFNESS DETECTION IS
C                                        TURNED ON.
C
                                    EEOET = 0.0D0
                                    ESTIFF = 0.0D0
                                    KTOL = 1
                                    DO 420 K = 1, NEQ
                                       YAVG = 0.5D0
     1                                        *(ABS(Y(K))
     2                                          + ABS(YS(K)))
                                       IF (INFO(2) .EQ. 1) KTOL = K
                                       ET = RTOL(KTOL)*YAVG + ATOL(KTOL)
                                       IF (ET .GT. 0.0D0) GO TO 400
C
C           PURE RELATIVE ERROR INAPPROPRIATE WHEN SOLUTION
C                                                  VANISHES
                                          IDID = -3
C              ...........................EXIT
                                          GO TO 520
  400                                  CONTINUE
C
                                       EE = ABS((-2090.0D0*YP(K)
     1                                            +(21970.0D0*F3(K)
     2                                              -15048.0D0*F4(K)))
     3                                           +(22528.0D0*F2(K)
     4                                             -27360.0D0*F5(K)))
                                       IF (STIFF .OR. NONSTF) GO TO 410
                                          ES = ABS(H
     1                                              *(0.055455D0*YP(K)
     2                                                -0.035493D0*F1(K)
     3                                                -0.036571D0*F2(K)
     4                                                +0.023107D0*F3(K)
     5                                                -0.009515D0*F4(K)
     6                                                +0.003017D0*F5(K))
     7                                                )
                                          ESTIFF = MAX(ESTIFF,ES/ET)
  410                                  CONTINUE
                                       EEOET = MAX(EEOET,EE/ET)
  420                               CONTINUE
C
                                    ESTTOL = ABS(H)*EEOET/752400.0D0
C
C                                ...EXIT
                                    IF (ESTTOL .LE. 1.0D0) GO TO 440
C
C                                   ....................................
C
C                                        UNSUCCESSFUL STEP
C
                                    IF (ABS(H) .GT. HMIN) GO TO 430
C
C                             REQUESTED ERROR UNATTAINABLE AT SMALLEST
C                                                  ALLOWABLE STEP SIZE
                                       TOLFAC = 1.69D0*ESTTOL
                                       IDID = -2
C              ........................EXIT
                                       GO TO 520
  430                               CONTINUE
C
C                                   REDUCE THE STEP SIZE , TRY AGAIN
C                                   THE DECREASE IS LIMITED TO A FACTOR
C                                   OF 1/10
C
                                    HFAILD = .TRUE.
                                    OUTPUT = .FALSE.
                                    S = 0.1D0
                                    IF (ESTTOL .LT. 59049.0D0)
     1                                 S = 0.9D0/ESTTOL**0.2D0
                                    H = SIGN(MAX(S*ABS(H),HMIN),H)
                                 GO TO 380
  440                            CONTINUE
C
C                                .......................................
C
C                                SUCCESSFUL STEP
C                                                  STORE SOLUTION AT T+H
C                                                  AND EVALUATE
C                                                  DERIVATIVES THERE
C
                                 T = T + H
                                 DO 450 K = 1, NEQ
                                    Y(K) = YS(K)
  450                            CONTINUE
                                 A = T
                                 CALL DF(A,Y,YP,RPAR,IPAR)
C
C                                CHOOSE NEXT STEP SIZE
C                                THE INCREASE IS LIMITED TO A FACTOR OF
C                                5 IF STEP FAILURE HAS JUST OCCURRED,
C                                NEXT
C                                   STEP SIZE IS NOT ALLOWED TO INCREASE
C
                                 S = 5.0D0
                                 IF (ESTTOL .GT. 1.889568D-4)
     1                              S = 0.9D0/ESTTOL**0.2D0
                                 IF (HFAILD) S = MIN(S,1.0D0)
                                 H = SIGN(MAX(S*ABS(H),HMIN),H)
C
C                                .......................................
C
C                                     CHECK FOR STIFFNESS (IF NOT
C                                     ALREADY DETECTED)
C
C                                     IN A SEQUENCE OF 50 SUCCESSFUL
C                                     STEPS BY THE FEHLBERG METHOD, 25
C                                     SUCCESSFUL STEPS BY THE FIRST
C                                     ORDER METHOD INDICATES STIFFNESS
C                                     AND TURNS THE TEST OFF. IF 26
C                                     FAILURES BY THE FIRST ORDER METHOD
C                                     OCCUR, THE TEST IS TURNED OFF
C                                     UNTIL THIS SEQUENCE OF 50 STEPS BY
C                                     THE FEHLBERG METHOD IS COMPLETED.
C
C                             ...EXIT
                                 IF (STIFF) GO TO 480
                                 NTSTEP = MOD(NTSTEP+1,50)
                                 IF (NTSTEP .EQ. 1) NONSTF = .FALSE.
C                             ...EXIT
                                 IF (NONSTF) GO TO 480
                                 IF (ESTIFF .GT. 1.0D0) GO TO 460
C
C                                   SUCCESSFUL STEP WITH FIRST ORDER
C                                   METHOD
                                    NSTIFS = NSTIFS + 1
C                                   TURN TEST OFF AFTER 25 INDICATIONS
C                                   OF STIFFNESS
                                    IF (NSTIFS .EQ. 25) STIFF = .TRUE.
                                 GO TO 470
  460                            CONTINUE
C
C                                UNSUCCESSFUL STEP WITH FIRST ORDER
C                                METHOD
                                 IF (NTSTEP - NSTIFS .LE. 25) GO TO 470
C               TURN STIFFNESS DETECTION OFF FOR THIS BLOCK OF
C                                                  FIFTY STEPS
                                    NONSTF = .TRUE.
C                                   RESET STIFF STEP COUNTER
                                    NSTIFS = 0
  470                            CONTINUE
  480                         CONTINUE
C
C                             ******************************************
C                                  END OF CORE INTEGRATOR
C                             ******************************************
C
C
C                                  SHOULD WE TAKE ANOTHER STEP
C
C                       ......EXIT
                              IF (OUTPUT) GO TO 490
                           IF (INFO(3) .EQ. 0) GO TO 300
C
C                          *********************************************
C                          *********************************************
C
C                               INTEGRATION SUCCESSFULLY COMPLETED
C
C                                           ONE-STEP MODE
                           IDID = 1
                           TOLD = T
C     .....................EXIT
                           GO TO 560
  490                   CONTINUE
  500                CONTINUE
C
C                    INTERVAL MODE
                     IDID = 2
                     T = TOUT
                     TOLD = T
C     ...............EXIT
                     GO TO 560
  510             CONTINUE
  520          CONTINUE
  530       CONTINUE
  540    CONTINUE
C
C        INTEGRATION TASK INTERRUPTED
C
         INFO(1) = -1
         TOLD = T
C     ...EXIT
         IF (IDID .NE. (-2)) GO TO 560
C
C        THE ERROR TOLERANCES ARE INCREASED TO VALUES
C                WHICH ARE APPROPRIATE FOR CONTINUING
         RTOL(1) = TOLFAC*RTOL(1)
         ATOL(1) = TOLFAC*ATOL(1)
C     ...EXIT
         IF (INFO(2) .EQ. 0) GO TO 560
         DO 550 K = 2, NEQ
            RTOL(K) = TOLFAC*RTOL(K)
            ATOL(K) = TOLFAC*ATOL(K)
  550    CONTINUE
  560 CONTINUE
      RETURN
      END
