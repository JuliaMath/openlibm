*DECK LSOD
      SUBROUTINE LSOD (F, NEQ, T, Y, TOUT, RTOL, ATOL, IDID, YPOUT, YH,
     +   YH1, EWT, SAVF, ACOR, WM, IWM, JAC, INTOUT, TSTOP, TOLFAC,
     +   DELSGN, RPAR, IPAR)
C***BEGIN PROLOGUE  LSOD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DEBDF
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LSOD-S, DLSOD-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   DEBDF  merely allocates storage for  LSOD  to relieve the user of
C   the inconvenience of a long call list.  Consequently  LSOD  is used
C   as described in the comments for  DEBDF .
C
C***SEE ALSO  DEBDF
C***ROUTINES CALLED  HSTART, INTYD, R1MACH, STOD, VNWRMS, XERMSG
C***COMMON BLOCKS    DEBDF1
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C***END PROLOGUE  LSOD
C
C
      LOGICAL INTOUT
C
      DIMENSION Y(*),YPOUT(*),YH(NEQ,6),YH1(*),EWT(*),SAVF(*),
     1          ACOR(*),WM(*),IWM(*),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3, XERN4
C
      COMMON /DEBDF1/ TOLD, ROWNS(210),
     1   EL0, H, HMIN, HMXI, HU, X, U,
     2   IQUIT, INIT, LYH, LEWT, LACOR, LSAVF, LWM, KSTEPS,
     3   IBEGIN, ITOL, IINTEG, ITSTOP, IJAC, IBAND, IOWNS(6),
     4   IER, JSTART, KFLAG, LDUM, METH, MITER, MAXORD, N, NQ, NST,
     5   NFE, NJE, NQU
C
      EXTERNAL F, JAC
C
C.......................................................................
C
C  THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
C  NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE COUNTER
C  IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE
C  WORK.
C
      SAVE MAXNUM
      DATA MAXNUM/500/
C
C.......................................................................
C
C***FIRST EXECUTABLE STATEMENT  LSOD
      IF (IBEGIN .EQ. 0) THEN
C
C        ON THE FIRST CALL , PERFORM INITIALIZATION --
C        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
C        FUNCTION ROUTINE R1MACH. THE USER MUST MAKE SURE THAT THE
C        VALUES SET IN R1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
C
         U = R1MACH(4)
C                          -- SET ASSOCIATED MACHINE DEPENDENT PARAMETER
         WM(1) = SQRT(U)
C                          -- SET TERMINATION FLAG
         IQUIT = 0
C                          -- SET INITIALIZATION INDICATOR
         INIT = 0
C                          -- SET COUNTER FOR ATTEMPTED STEPS
         KSTEPS = 0
C                          -- SET INDICATOR FOR INTERMEDIATE-OUTPUT
         INTOUT = .FALSE.
C                          -- SET START INDICATOR FOR STOD CODE
         JSTART = 0
C                          -- SET BDF METHOD INDICATOR
         METH = 2
C                          -- SET MAXIMUM ORDER FOR BDF METHOD
         MAXORD = 5
C                          -- SET ITERATION MATRIX INDICATOR
C
         IF (IJAC .EQ. 0 .AND. IBAND .EQ. 0) MITER = 2
         IF (IJAC .EQ. 1 .AND. IBAND .EQ. 0) MITER = 1
         IF (IJAC .EQ. 0 .AND. IBAND .EQ. 1) MITER = 5
         IF (IJAC .EQ. 1 .AND. IBAND .EQ. 1) MITER = 4
C
C                          -- SET OTHER NECESSARY ITEMS IN COMMON BLOCK
         N = NEQ
         NST = 0
         NJE = 0
         HMXI = 0.
         NQ = 1
         H = 1.
C                          -- RESET IBEGIN FOR SUBSEQUENT CALLS
         IBEGIN=1
      ENDIF
C
C.......................................................................
C
C      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
C
      IF (NEQ .LT. 1) THEN
         WRITE (XERN1, '(I8)') NEQ
         CALL XERMSG ('SLATEC', 'LSOD',
     *      'IN DEBDF, THE NUMBER OF EQUATIONS MUST BE A POSITIVE ' //
     *      'INTEGER.$$YOU HAVE CALLED THE CODE WITH NEQ = ' // XERN1,
     *      6, 1)
         IDID=-33
      ENDIF
C
      NRTOLP = 0
      NATOLP = 0
      DO 60 K = 1,NEQ
         IF (NRTOLP .LE. 0) THEN
            IF (RTOL(K) .LT. 0.) THEN
               WRITE (XERN1, '(I8)') K
               WRITE (XERN3, '(1PE15.6)') RTOL(K)
               CALL XERMSG ('SLATEC', 'LSOD',
     *            'IN DEBDF, THE RELATIVE ERROR TOLERANCES MUST ' //
     *            'BE NON-NEGATIVE.$$YOU HAVE CALLED THE CODE WITH ' //
     *            'RTOL(' // XERN1 // ') = ' // XERN3 // '$$IN THE ' //
     *            'CASE OF VECTOR ERROR TOLERANCES, NO FURTHER ' //
     *            'CHECKING OF RTOL COMPONENTS IS DONE.', 7, 1)
               IDID = -33
               IF (NATOLP .GT. 0) GO TO 70
               NRTOLP = 1
            ELSEIF (NATOLP .GT. 0) THEN
               GO TO 50
            ENDIF
         ENDIF
C
         IF (ATOL(K) .LT. 0.) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') ATOL(K)
            CALL XERMSG ('SLATEC', 'LSOD',
     *         'IN DEBDF, THE ABSOLUTE ERROR ' //
     *         'TOLERANCES MUST BE NON-NEGATIVE.$$YOU HAVE CALLED ' //
     *         'THE CODE WITH ATOL(' // XERN1 // ') = ' // XERN3 //
     *         '$$IN THE CASE OF VECTOR ERROR TOLERANCES, NO FURTHER '
     *         // 'CHECKING OF ATOL COMPONENTS IS DONE.', 8, 1)
            IDID=-33
            IF (NRTOLP .GT. 0) GO TO 70
            NATOLP=1
         ENDIF
   50    IF (ITOL .EQ. 0) GO TO 70
   60 CONTINUE
C
   70 IF (ITSTOP .EQ. 1) THEN
         IF (SIGN(1.,TOUT-T) .NE. SIGN(1.,TSTOP-T) .OR.
     1      ABS(TOUT-T) .GT. ABS(TSTOP-T)) THEN
            WRITE (XERN3, '(1PE15.6)') TOUT
            WRITE (XERN4, '(1PE15.6)') TSTOP
            CALL XERMSG ('SLATEC', 'LSOD',
     *         'IN DEBDF, YOU HAVE CALLED THE ' //
     *         'CODE WITH TOUT = ' // XERN3 // '$$BUT YOU HAVE ' //
     *         'ALSO TOLD THE CODE NOT TO INTEGRATE PAST THE POINT ' //
     *         'TSTOP = ' // XERN4 // ' BY SETTING INFO(4) = 1.  ' //
     *         'THESE INSTRUCTIONS CONFLICT.', 14, 1)
            IDID=-33
         ENDIF
      ENDIF
C
C        CHECK SOME CONTINUATION POSSIBILITIES
C
      IF (INIT .NE. 0) THEN
         IF (T .EQ. TOUT) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'LSOD',
     *         'IN DEBDF, YOU HAVE CALLED THE CODE WITH T = TOUT = ' //
     *         XERN3 // '  THIS IS NOT ALLOWED ON CONTINUATION CALLS.',
     *         9, 1)
            IDID=-33
         ENDIF
C
         IF (T .NE. TOLD) THEN
            WRITE (XERN3, '(1PE15.6)') TOLD
            WRITE (XERN4, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'LSOD',
     *         'IN DEBDF, YOU HAVE CHANGED THE VALUE OF T FROM ' //
     *         XERN3 // ' TO ' // XERN4 //
     *         '  THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
            IDID=-33
         ENDIF
C
         IF (INIT .NE. 1) THEN
            IF (DELSGN*(TOUT-T) .LT. 0.) THEN
               WRITE (XERN3, '(1PE15.6)') TOUT
               CALL XERMSG ('SLATEC', 'LSOD',
     *            'IN DEBDF, BY CALLING THE CODE WITH TOUT = ' //
     *            XERN3 // ' YOU ARE ATTEMPTING TO CHANGE THE ' //
     *            'DIRECTION OF INTEGRATION.$$' //
     *            'THIS IS NOT ALLOWED WITHOUT RESTARTING.', 11, 1)
               IDID=-33
            ENDIF
         ENDIF
      ENDIF
C
      IF (IDID .EQ. (-33)) THEN
         IF (IQUIT .NE. (-33)) THEN
C                       INVALID INPUT DETECTED
            IQUIT=-33
            IBEGIN=-1
         ELSE
            CALL XERMSG ('SLATEC', 'LSOD',
     *         'IN DEBDF, INVALID INPUT WAS ' //
     *         'DETECTED ON SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE ' //
     *         'TO PROCEED BECAUSE YOU HAVE NOT CORRECTED THE ' //
     *         'PROBLEM, SO EXECUTION IS BEING TERMINATED.', 12, 2)
         ENDIF
         RETURN
      ENDIF
C
C.......................................................................
C
C     RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS
C     ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE,
C     THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE
C     100*U WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE
C
      DO 170 K=1,NEQ
        IF (RTOL(K)+ATOL(K) .GT. 0.) GO TO 160
        RTOL(K)=100.*U
        IDID=-2
  160   IF (ITOL .EQ. 0) GO TO 180
  170   CONTINUE
C
  180 IF (IDID .NE. (-2)) GO TO 190
C                       RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
C                                                SMALL POSITIVE VALUE
      IBEGIN=-1
      RETURN
C
C     BRANCH ON STATUS OF INITIALIZATION INDICATOR
C            INIT=0 MEANS INITIAL DERIVATIVES AND NOMINAL STEP SIZE
C                   AND DIRECTION NOT YET SET
C            INIT=1 MEANS NOMINAL STEP SIZE AND DIRECTION NOT YET SET
C            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED
C
  190 IF (INIT .EQ. 0) GO TO 200
      IF (INIT .EQ. 1) GO TO 220
      GO TO 240
C
C.......................................................................
C
C     MORE INITIALIZATION --
C                         -- EVALUATE INITIAL DERIVATIVES
C
  200 INIT=1
      CALL F(T,Y,YH(1,2),RPAR,IPAR)
      NFE=1
      IF (T .NE. TOUT) GO TO 220
      IDID=2
      DO 210 L = 1,NEQ
  210    YPOUT(L) = YH(L,2)
      TOLD=T
      RETURN
C
C                         -- COMPUTE INITIAL STEP SIZE
C                         -- SAVE SIGN OF INTEGRATION DIRECTION
C                         -- SET INDEPENDENT AND DEPENDENT VARIABLES
C                                              X AND YH(*) FOR STOD
C
  220 LTOL = 1
      DO 225 L=1,NEQ
        IF (ITOL .EQ. 1) LTOL = L
        TOL = RTOL(LTOL)*ABS(Y(L)) + ATOL(LTOL)
        IF (TOL .EQ. 0.) GO TO 380
  225   EWT(L) = TOL
C
      BIG = SQRT(R1MACH(2))
      CALL HSTART (F,NEQ,T,TOUT,Y,YH(1,2),EWT,1,U,BIG,
     1             YH(1,3),YH(1,4),YH(1,5),YH(1,6),RPAR,IPAR,H)
C
      DELSGN = SIGN(1.0,TOUT-T)
      X = T
      DO 230 L = 1,NEQ
        YH(L,1) = Y(L)
  230   YH(L,2) = H*YH(L,2)
      INIT = 2
C
C.......................................................................
C
C   ON EACH CALL SET INFORMATION WHICH DETERMINES THE ALLOWED INTERVAL
C   OF INTEGRATION BEFORE RETURNING WITH AN ANSWER AT TOUT
C
  240 DEL = TOUT - T
      ABSDEL = ABS(DEL)
C
C.......................................................................
C
C   IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN
C
  250 IF (ABS(X-T) .LT. ABSDEL) GO TO 270
      CALL INTYD(TOUT,0,YH,NEQ,Y,INTFLG)
      CALL INTYD(TOUT,1,YH,NEQ,YPOUT,INTFLG)
      IDID = 3
      IF (X .NE. TOUT) GO TO 260
      IDID = 2
      INTOUT = .FALSE.
  260 T = TOUT
      TOLD = T
      RETURN
C
C   IF CANNOT GO PAST TSTOP AND SUFFICIENTLY CLOSE,
C   EXTRAPOLATE AND RETURN
C
  270 IF (ITSTOP .NE. 1) GO TO 290
      IF (ABS(TSTOP-X) .GE. 100.*U*ABS(X)) GO TO 290
      DT = TOUT - X
      DO 280 L = 1,NEQ
  280   Y(L) = YH(L,1) + (DT/H)*YH(L,2)
      CALL F(TOUT,Y,YPOUT,RPAR,IPAR)
      NFE = NFE + 1
      IDID = 3
      T = TOUT
      TOLD = T
      RETURN
C
  290 IF (IINTEG .EQ. 0  .OR.  .NOT.INTOUT) GO TO 300
C
C   INTERMEDIATE-OUTPUT MODE
C
      IDID = 1
      GO TO 500
C
C.......................................................................
C
C     MONITOR NUMBER OF STEPS ATTEMPTED
C
  300 IF (KSTEPS .LE. MAXNUM) GO TO 330
C
C                       A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED
      IDID=-1
      KSTEPS=0
      IBEGIN = -1
      GO TO 500
C
C.......................................................................
C
C   LIMIT STEP SIZE AND SET WEIGHT VECTOR
C
  330 HMIN = 100.*U*ABS(X)
      HA = MAX(ABS(H),HMIN)
      IF (ITSTOP .NE. 1) GO TO 340
      HA = MIN(HA,ABS(TSTOP-X))
  340 H = SIGN(HA,H)
      LTOL = 1
      DO 350 L = 1,NEQ
        IF (ITOL .EQ. 1) LTOL = L
        EWT(L) = RTOL(LTOL)*ABS(YH(L,1)) + ATOL(LTOL)
        IF (EWT(L) .LE. 0.0) GO TO 380
  350   CONTINUE
      TOLFAC = U*VNWRMS(NEQ,YH,EWT)
      IF (TOLFAC .LE. 1.) GO TO 400
C
C                       TOLERANCES TOO SMALL
      IDID = -2
      TOLFAC = 2.*TOLFAC
      RTOL(1) = TOLFAC*RTOL(1)
      ATOL(1) = TOLFAC*ATOL(1)
      IF (ITOL .EQ. 0) GO TO 370
      DO 360 L = 2,NEQ
        RTOL(L) = TOLFAC*RTOL(L)
  360   ATOL(L) = TOLFAC*ATOL(L)
  370 IBEGIN = -1
      GO TO 500
C
C                       RELATIVE ERROR CRITERION INAPPROPRIATE
  380 IDID = -3
      IBEGIN = -1
      GO TO 500
C
C.......................................................................
C
C     TAKE A STEP
C
  400 CALL STOD(NEQ,Y,YH,NEQ,YH1,EWT,SAVF,ACOR,WM,IWM,F,JAC,RPAR,IPAR)
C
      JSTART = -2
      INTOUT = .TRUE.
      IF (KFLAG .EQ. 0) GO TO 250
C
C.......................................................................
C
      IF (KFLAG .EQ. -1) GO TO 450
C
C                       REPEATED CORRECTOR CONVERGENCE FAILURES
      IDID = -6
      IBEGIN = -1
      GO TO 500
C
C                       REPEATED ERROR TEST FAILURES
  450 IDID = -7
      IBEGIN = -1
C
C.......................................................................
C
C                       STORE VALUES BEFORE RETURNING TO DEBDF
  500 DO 555 L = 1,NEQ
        Y(L) = YH(L,1)
  555   YPOUT(L) = YH(L,2)/H
      T = X
      TOLD = T
      INTOUT = .FALSE.
      RETURN
      END
