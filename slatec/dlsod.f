*DECK DLSOD
      SUBROUTINE DLSOD (DF, NEQ, T, Y, TOUT, RTOL, ATOL, IDID, YPOUT,
     +   YH, YH1, EWT, SAVF, ACOR, WM, IWM, DJAC, INTOUT, TSTOP, TOLFAC,
     +   DELSGN, RPAR, IPAR)
C***BEGIN PROLOGUE  DLSOD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (LSOD-S, DLSOD-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   DDEBDF  merely allocates storage for  DLSOD  to relieve the user of
C   the inconvenience of a long call list.  Consequently  DLSOD  is used
C   as described in the comments for  DDEBDF .
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  D1MACH, DHSTRT, DINTYD, DSTOD, DVNRMS, XERMSG
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C***END PROLOGUE  DLSOD
C
      INTEGER IBAND, IBEGIN, IDID, IER, IINTEG, IJAC, INIT, INTFLG,
     1      IOWNS, IPAR, IQUIT, ITOL, ITSTOP, IWM, JSTART, K, KFLAG,
     2      KSTEPS, L, LACOR, LDUM, LEWT, LSAVF, LTOL, LWM, LYH, MAXNUM,
     3      MAXORD, METH, MITER, N, NATOLP, NEQ, NFE, NJE, NQ, NQU,
     4      NRTOLP, NST
      DOUBLE PRECISION ABSDEL, ACOR, ATOL, BIG, D1MACH, DEL,
     1      DELSGN, DT, DVNRMS, EL0, EWT,
     2      H, HA, HMIN, HMXI, HU, ROWNS, RPAR, RTOL, SAVF, T, TOL,
     3      TOLD, TOLFAC, TOUT, TSTOP, U, WM, X, Y, YH, YH1, YPOUT
      LOGICAL INTOUT
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3, XERN4
C
      DIMENSION Y(*),YPOUT(*),YH(NEQ,6),YH1(*),EWT(*),SAVF(*),
     1          ACOR(*),WM(*),IWM(*),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
C
C
      COMMON /DDEBD1/ TOLD,ROWNS(210),EL0,H,HMIN,HMXI,HU,X,U,IQUIT,INIT,
     1                LYH,LEWT,LACOR,LSAVF,LWM,KSTEPS,IBEGIN,ITOL,
     2                IINTEG,ITSTOP,IJAC,IBAND,IOWNS(6),IER,JSTART,
     3                KFLAG,LDUM,METH,MITER,MAXORD,N,NQ,NST,NFE,NJE,NQU
C
      EXTERNAL DF, DJAC
C
C     ..................................................................
C
C       THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
C       NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE
C       COUNTER IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE
C       EXCESSIVE WORK.
      SAVE MAXNUM
C
      DATA MAXNUM /500/
C
C     ..................................................................
C
C***FIRST EXECUTABLE STATEMENT  DLSOD
      IF (IBEGIN .EQ. 0) THEN
C
C        ON THE FIRST CALL , PERFORM INITIALIZATION --
C        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
C        FUNCTION ROUTINE D1MACH. THE USER MUST MAKE SURE THAT THE
C        VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
C
         U = D1MACH(4)
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
C                          -- SET START INDICATOR FOR DSTOD CODE
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
         HMXI = 0.0D0
         NQ = 1
         H = 1.0D0
C                          -- RESET IBEGIN FOR SUBSEQUENT CALLS
         IBEGIN = 1
      ENDIF
C
C     ..................................................................
C
C      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
C
      IF (NEQ .LT. 1) THEN
         WRITE (XERN1, '(I8)') NEQ
         CALL XERMSG ('SLATEC', 'DLSOD',
     *      'IN DDEBDF, THE NUMBER OF EQUATIONS MUST BE A ' //
     *      'POSITIVE INTEGER.$$YOU HAVE CALLED THE CODE WITH NEQ = ' //
     *      XERN1, 6, 1)
         IDID=-33
      ENDIF
C
      NRTOLP = 0
      NATOLP = 0
      DO 60 K = 1, NEQ
         IF (NRTOLP .LE. 0) THEN
            IF (RTOL(K) .LT. 0.) THEN
               WRITE (XERN1, '(I8)') K
               WRITE (XERN3, '(1PE15.6)') RTOL(K)
               CALL XERMSG ('SLATEC', 'DLSOD',
     *            'IN DDEBDF, THE RELATIVE ERROR TOLERANCES MUST ' //
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
            CALL XERMSG ('SLATEC', 'DLSOD',
     *         'IN DDEBDF, THE ABSOLUTE ERROR ' //
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
         IF (SIGN(1.0D0,TOUT-T) .NE. SIGN(1.0D0,TSTOP-T) .OR.
     1      ABS(TOUT-T) .GT. ABS(TSTOP-T)) THEN
            WRITE (XERN3, '(1PE15.6)') TOUT
            WRITE (XERN4, '(1PE15.6)') TSTOP
            CALL XERMSG ('SLATEC', 'DLSOD',
     *         'IN DDEBDF, YOU HAVE CALLED THE ' //
     *         'CODE WITH TOUT = ' // XERN3 // '$$BUT YOU HAVE ' //
     *         'ALSO TOLD THE CODE NOT TO INTEGRATE PAST THE POINT ' //
     *         'TSTOP = ' // XERN4 // ' BY SETTING INFO(4) = 1.$$' //
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
            CALL XERMSG ('SLATEC', 'DLSOD',
     *         'IN DDEBDF, YOU HAVE CALLED THE CODE WITH T = TOUT = ' //
     *         XERN3 // '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.',
     *         9, 1)
            IDID=-33
         ENDIF
C
         IF (T .NE. TOLD) THEN
            WRITE (XERN3, '(1PE15.6)') TOLD
            WRITE (XERN4, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DLSOD',
     *         'IN DDEBDF, YOU HAVE CHANGED THE VALUE OF T FROM ' //
     *         XERN3 // ' TO ' // XERN4 //
     *         '  THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 10, 1)
            IDID=-33
         ENDIF
C
         IF (INIT .NE. 1) THEN
            IF (DELSGN*(TOUT-T) .LT. 0.0D0) THEN
               WRITE (XERN3, '(1PE15.6)') TOUT
               CALL XERMSG ('SLATEC', 'DLSOD',
     *            'IN DDEBDF, BY CALLING THE CODE WITH TOUT = ' //
     *            XERN3 // ' YOU ARE ATTEMPTING TO CHANGE THE ' //
     *            'DIRECTION OF INTEGRATION.$$THIS IS NOT ALLOWED ' //
     *            'WITHOUT RESTARTING.', 11, 1)
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
            CALL XERMSG ('SLATEC', 'DLSOD',
     *         'IN DDEBDF, INVALID INPUT WAS DETECTED ON ' //
     *         'SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED ' //
     *         'BECAUSE YOU HAVE NOT CORRECTED THE PROBLEM, ' //
     *         'SO EXECUTION IS BEING TERMINATED.', 12, 2)
         ENDIF
         RETURN
      ENDIF
C
C        ...............................................................
C
C             RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED
C             AS ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS
C             CASE, THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE
C             SMALLEST VALUE 100*U WHICH IS LIKELY TO BE REASONABLE FOR
C             THIS METHOD AND MACHINE
C
      DO 180 K = 1, NEQ
         IF (RTOL(K) + ATOL(K) .GT. 0.0D0) GO TO 170
            RTOL(K) = 100.0D0*U
            IDID = -2
  170    CONTINUE
C     ...EXIT
         IF (ITOL .EQ. 0) GO TO 190
  180 CONTINUE
  190 CONTINUE
C
      IF (IDID .NE. (-2)) GO TO 200
C        RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
C                                 SMALL POSITIVE VALUE
         IBEGIN = -1
      GO TO 460
  200 CONTINUE
C        BEGIN BLOCK PERMITTING ...EXITS TO 450
C           BEGIN BLOCK PERMITTING ...EXITS TO 430
C              BEGIN BLOCK PERMITTING ...EXITS TO 260
C                 BEGIN BLOCK PERMITTING ...EXITS TO 230
C
C                    BRANCH ON STATUS OF INITIALIZATION INDICATOR
C                           INIT=0 MEANS INITIAL DERIVATIVES AND
C                           NOMINAL STEP SIZE
C                                  AND DIRECTION NOT YET SET
C                           INIT=1 MEANS NOMINAL STEP SIZE AND
C                           DIRECTION NOT YET SET INIT=2 MEANS NO
C                           FURTHER INITIALIZATION REQUIRED
C
                     IF (INIT .EQ. 0) GO TO 210
C                 ......EXIT
                        IF (INIT .EQ. 1) GO TO 230
C              .........EXIT
                        GO TO 260
  210                CONTINUE
C
C                    ................................................
C
C                         MORE INITIALIZATION --
C                                             -- EVALUATE INITIAL
C                                             DERIVATIVES
C
                     INIT = 1
                     CALL DF(T,Y,YH(1,2),RPAR,IPAR)
                     NFE = 1
C                 ...EXIT
                     IF (T .NE. TOUT) GO TO 230
                     IDID = 2
                     DO 220 L = 1, NEQ
                        YPOUT(L) = YH(L,2)
  220                CONTINUE
                     TOLD = T
C        ............EXIT
                     GO TO 450
  230             CONTINUE
C
C                 -- COMPUTE INITIAL STEP SIZE
C                 -- SAVE SIGN OF INTEGRATION DIRECTION
C                 -- SET INDEPENDENT AND DEPENDENT VARIABLES
C                                      X AND YH(*) FOR DSTOD
C
                  LTOL = 1
                  DO 240 L = 1, NEQ
                     IF (ITOL .EQ. 1) LTOL = L
                     TOL = RTOL(LTOL)*ABS(Y(L)) + ATOL(LTOL)
                     IF (TOL .EQ. 0.0D0) GO TO 390
                     EWT(L) = TOL
  240             CONTINUE
C
                  BIG = SQRT(D1MACH(2))
                  CALL DHSTRT(DF,NEQ,T,TOUT,Y,YH(1,2),EWT,1,U,BIG,
     1                        YH(1,3),YH(1,4),YH(1,5),YH(1,6),RPAR,
     2                        IPAR,H)
C
                  DELSGN = SIGN(1.0D0,TOUT-T)
                  X = T
                  DO 250 L = 1, NEQ
                     YH(L,1) = Y(L)
                     YH(L,2) = H*YH(L,2)
  250             CONTINUE
                  INIT = 2
  260          CONTINUE
C
C              ......................................................
C
C                 ON EACH CALL SET INFORMATION WHICH DETERMINES THE
C                 ALLOWED INTERVAL OF INTEGRATION BEFORE RETURNING
C                 WITH AN ANSWER AT TOUT
C
               DEL = TOUT - T
               ABSDEL = ABS(DEL)
C
C              ......................................................
C
C                 IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND
C                 RETURN
C
  270          CONTINUE
C                 BEGIN BLOCK PERMITTING ...EXITS TO 400
C                    BEGIN BLOCK PERMITTING ...EXITS TO 380
                        IF (ABS(X-T) .LT. ABSDEL) GO TO 290
                           CALL DINTYD(TOUT,0,YH,NEQ,Y,INTFLG)
                           CALL DINTYD(TOUT,1,YH,NEQ,YPOUT,INTFLG)
                           IDID = 3
                           IF (X .NE. TOUT) GO TO 280
                              IDID = 2
                              INTOUT = .FALSE.
  280                      CONTINUE
                           T = TOUT
                           TOLD = T
C        ..................EXIT
                           GO TO 450
  290                   CONTINUE
C
C                       IF CANNOT GO PAST TSTOP AND SUFFICIENTLY
C                       CLOSE, EXTRAPOLATE AND RETURN
C
                        IF (ITSTOP .NE. 1) GO TO 310
                        IF (ABS(TSTOP-X) .GE. 100.0D0*U*ABS(X))
     1                     GO TO 310
                           DT = TOUT - X
                           DO 300 L = 1, NEQ
                              Y(L) = YH(L,1) + (DT/H)*YH(L,2)
  300                      CONTINUE
                           CALL DF(TOUT,Y,YPOUT,RPAR,IPAR)
                           NFE = NFE + 1
                           IDID = 3
                           T = TOUT
                           TOLD = T
C        ..................EXIT
                           GO TO 450
  310                   CONTINUE
C
                        IF (IINTEG .EQ. 0 .OR. .NOT.INTOUT) GO TO 320
C
C                          INTERMEDIATE-OUTPUT MODE
C
                           IDID = 1
                        GO TO 370
  320                   CONTINUE
C
C                       .............................................
C
C                            MONITOR NUMBER OF STEPS ATTEMPTED
C
                        IF (KSTEPS .LE. MAXNUM) GO TO 330
C
C                          A SIGNIFICANT AMOUNT OF WORK HAS BEEN
C                          EXPENDED
                           IDID = -1
                           KSTEPS = 0
                           IBEGIN = -1
                        GO TO 370
  330                   CONTINUE
C
C                          ..........................................
C
C                             LIMIT STEP SIZE AND SET WEIGHT VECTOR
C
                           HMIN = 100.0D0*U*ABS(X)
                           HA = MAX(ABS(H),HMIN)
                           IF (ITSTOP .EQ. 1)
     1                        HA = MIN(HA,ABS(TSTOP-X))
                           H = SIGN(HA,H)
                           LTOL = 1
                           DO 340 L = 1, NEQ
                              IF (ITOL .EQ. 1) LTOL = L
                              EWT(L) = RTOL(LTOL)*ABS(YH(L,1))
     1                                 + ATOL(LTOL)
C                    .........EXIT
                              IF (EWT(L) .LE. 0.0D0) GO TO 380
  340                      CONTINUE
                           TOLFAC = U*DVNRMS(NEQ,YH,EWT)
C                 .........EXIT
                           IF (TOLFAC .LE. 1.0D0) GO TO 400
C
C                          TOLERANCES TOO SMALL
                           IDID = -2
                           TOLFAC = 2.0D0*TOLFAC
                           RTOL(1) = TOLFAC*RTOL(1)
                           ATOL(1) = TOLFAC*ATOL(1)
                           IF (ITOL .EQ. 0) GO TO 360
                              DO 350 L = 2, NEQ
                                 RTOL(L) = TOLFAC*RTOL(L)
                                 ATOL(L) = TOLFAC*ATOL(L)
  350                         CONTINUE
  360                      CONTINUE
                           IBEGIN = -1
  370                   CONTINUE
C           ............EXIT
                        GO TO 430
  380                CONTINUE
C
C                    RELATIVE ERROR CRITERION INAPPROPRIATE
  390                CONTINUE
                     IDID = -3
                     IBEGIN = -1
C           .........EXIT
                     GO TO 430
  400             CONTINUE
C
C                 ...................................................
C
C                      TAKE A STEP
C
                  CALL DSTOD(NEQ,Y,YH,NEQ,YH1,EWT,SAVF,ACOR,WM,IWM,
     1                       DF,DJAC,RPAR,IPAR)
C
                  JSTART = -2
                  INTOUT = .TRUE.
               IF (KFLAG .EQ. 0) GO TO 270
C
C              ......................................................
C
               IF (KFLAG .EQ. -1) GO TO 410
C
C                 REPEATED CORRECTOR CONVERGENCE FAILURES
                  IDID = -6
                  IBEGIN = -1
               GO TO 420
  410          CONTINUE
C
C                 REPEATED ERROR TEST FAILURES
                  IDID = -7
                  IBEGIN = -1
  420          CONTINUE
  430       CONTINUE
C
C           .........................................................
C
C                                  STORE VALUES BEFORE RETURNING TO
C                                  DDEBDF
            DO 440 L = 1, NEQ
               Y(L) = YH(L,1)
               YPOUT(L) = YH(L,2)/H
  440       CONTINUE
            T = X
            TOLD = T
            INTOUT = .FALSE.
  450    CONTINUE
  460 CONTINUE
      RETURN
      END
