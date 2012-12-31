*DECK LA05AD
      SUBROUTINE LA05AD (A, IND, NZ, IA, N, IP, IW, W, G, U)
C***BEGIN PROLOGUE  LA05AD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DSPLP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (LA05AS-S, LA05AD-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
C     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
C     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
C     THE FINAL LETTER =D= IN THE NAMES USED HERE.
C     REVISIONS MADE BY R J HANSON, SNLA, AUGUST, 1979.
C     REVISED SEP. 13, 1979.
C
C     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
C     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
C     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
C     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
C     DSPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
C
C IP(I,1),IP(I,2) POINT TO THE START OF ROW/COL I.
C IW(I,1),IW(I,2) HOLD THE NUMBER OF NON-ZEROS IN ROW/COL I.
C DURING THE MAIN BODY OF THIS SUBROUTINE THE VECTORS IW(.,3),IW(.,5),
C     IW(.,7) ARE USED TO HOLD DOUBLY LINKED LISTS OF ROWS THAT HAVE
C     NOT BEEN PIVOTAL AND HAVE EQUAL NUMBERS OF NON-ZEROS.
C IW(.,4),IW(.,6),IW(.,8) HOLD SIMILAR LISTS FOR THE COLUMNS.
C IW(I,3),IW(I,4) HOLD FIRST ROW/COLUMN TO HAVE I NON-ZEROS
C     OR ZERO IF THERE ARE NONE.
C IW(I,5), IW(I,6) HOLD ROW/COL NUMBER OF ROW/COL PRIOR TO ROW/COL I
C     IN ITS LIST, OR ZERO IF NONE.
C IW(I,7), IW(I,8) HOLD ROW/COL NUMBER OF ROW/COL AFTER ROW/COL I
C     IN ITS LIST, OR ZERO IF NONE.
C FOR ROWS/COLS THAT HAVE BEEN PIVOTAL IW(I,5),IW(I,6) HOLD NEGATION OF
C     POSITION OF ROW/COL I IN THE PIVOTAL ORDERING.
C
C***SEE ALSO  DSPLP
C***ROUTINES CALLED  D1MACH, LA05ED, MC20AD, XERMSG, XSETUN
C***COMMON BLOCKS    LA05DD
C***REVISION HISTORY  (YYMMDD)
C   811215  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890605  Added D1MACH to list of DOUBLE PRECISION variables.
C   890605  Corrected references to XERRWV.  (WRB)
C           (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900402  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C***END PROLOGUE  LA05AD
      INTEGER IP(N,2)
      INTEGER IND(IA,2), IW(N,8)
      DOUBLE PRECISION A(*), AMAX, AU, AM, D1MACH, EPS, G, U, SMALL,
     *                 W(*)
      LOGICAL FIRST
      CHARACTER*8 XERN0, XERN1, XERN2
C
      COMMON /LA05DD/ SMALL, LP, LENL, LENU, NCP, LROW, LCOL
C EPS IS THE RELATIVE ACCURACY OF FLOATING-POINT COMPUTATION
      SAVE EPS, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  LA05AD
      IF (FIRST) THEN
         EPS = 2.0D0 * D1MACH(4)
      ENDIF
      FIRST = .FALSE.
C
C     SET THE OUTPUT UNIT NUMBER FOR THE ERROR PROCESSOR.
C     THE USAGE OF THIS ERROR PROCESSOR IS DOCUMENTED IN THE
C     SANDIA LABS. TECH. REPT. SAND78-1189, BY R E JONES.
      CALL XSETUN(LP)
      IF (U.GT.1.0D0) U = 1.0D0
      IF (U.LT.EPS) U = EPS
      IF (N.LT.1) GO TO 670
      G = 0.
      DO 50 I=1,N
         W(I) = 0.
         DO 40 J=1,5
            IW(I,J) = 0
   40    CONTINUE
   50 CONTINUE
C
C FLUSH OUT SMALL ENTRIES, COUNT ELEMENTS IN ROWS AND COLUMNS
      L = 1
      LENU = NZ
      DO 80 IDUMMY=1,NZ
         IF (L.GT.LENU) GO TO 90
         DO 60 K=L,LENU
            IF (ABS(A(K)).LE.SMALL) GO TO 70
            I = IND(K,1)
            J = IND(K,2)
            G = MAX(ABS(A(K)),G)
            IF (I.LT.1 .OR. I.GT.N) GO TO 680
            IF (J.LT.1 .OR. J.GT.N) GO TO 680
            IW(I,1) = IW(I,1) + 1
            IW(J,2) = IW(J,2) + 1
   60    CONTINUE
         GO TO 90
   70    L = K
         A(L) = A(LENU)
         IND(L,1) = IND(LENU,1)
         IND(L,2) = IND(LENU,2)
         LENU = LENU - 1
   80 CONTINUE
C
   90 LENL = 0
      LROW = LENU
      LCOL = LROW
C MCP IS THE MAXIMUM NUMBER OF COMPRESSES PERMITTED BEFORE AN
C     ERROR RETURN RESULTS.
      MCP = MAX(N/10,20)
      NCP = 0
C CHECK FOR NULL ROW OR COLUMN AND INITIALIZE IP(I,2) TO POINT
C     JUST BEYOND WHERE THE LAST COMPONENT OF COLUMN I OF A WILL
C     BE STORED.
      K = 1
      DO 110 IR=1,N
         K = K + IW(IR,2)
         IP(IR,2) = K
         DO 100 L=1,2
            IF (IW(IR,L).LE.0) GO TO 700
  100    CONTINUE
  110 CONTINUE
C REORDER BY ROWS
C CHECK FOR DOUBLE ENTRIES WHILE USING THE NEWLY CONSTRUCTED
C     ROW FILE TO CONSTRUCT THE COLUMN FILE. NOTE THAT BY PUTTING
C    THE ENTRIES IN BACKWARDS AND DECREASING IP(J,2) EACH TIME IT
C     IS USED WE AUTOMATICALLY LEAVE IT POINTING TO THE FIRST ELEMENT.
      CALL MC20AD(N, LENU, A, IND(1,2), IP, IND(1,1), 0)
      KL = LENU
      DO 130 II=1,N
         IR = N + 1 - II
         KP = IP(IR,1)
         DO 120 K=KP,KL
            J = IND(K,2)
            IF (IW(J,5).EQ.IR) GO TO 660
            IW(J,5) = IR
            KR = IP(J,2) - 1
            IP(J,2) = KR
            IND(KR,1) = IR
  120    CONTINUE
         KL = KP - 1
  130 CONTINUE
C
C SET UP LINKED LISTS OF ROWS AND COLS WITH EQUAL NUMBERS OF NON-ZEROS.
      DO 150 L=1,2
         DO 140 I=1,N
            NZ = IW(I,L)
            IN = IW(NZ,L+2)
            IW(NZ,L+2) = I
            IW(I,L+6) = IN
            IW(I,L+4) = 0
            IF (IN.NE.0) IW(IN,L+4) = I
  140    CONTINUE
  150 CONTINUE
C
C
C START OF MAIN ELIMINATION LOOP.
      DO 590 IPV=1,N
C FIND PIVOT. JCOST IS MARKOWITZ COST OF CHEAPEST PIVOT FOUND SO FAR,
C     WHICH IS IN ROW IPP AND COLUMN JP.
         JCOST = N*N
C LOOP ON LENGTH OF COLUMN TO BE SEARCHED
         DO 240 NZ=1,N
            IF (JCOST.LE.(NZ-1)**2) GO TO 250
            J = IW(NZ,4)
C SEARCH COLUMNS WITH NZ NON-ZEROS.
            DO 190 IDUMMY=1,N
               IF (J.LE.0) GO TO 200
               KP = IP(J,2)
               KL = KP + IW(J,2) - 1
               DO 180 K=KP,KL
                  I = IND(K,1)
                  KCOST = (NZ-1)*(IW(I,1)-1)
                  IF (KCOST.GE.JCOST) GO TO 180
                  IF (NZ.EQ.1) GO TO 170
C FIND LARGEST ELEMENT IN ROW OF POTENTIAL PIVOT.
                  AMAX = 0.
                  K1 = IP(I,1)
                  K2 = IW(I,1) + K1 - 1
                  DO 160 KK=K1,K2
                     AMAX = MAX(AMAX,ABS(A(KK)))
                     IF (IND(KK,2).EQ.J) KJ = KK
  160             CONTINUE
C PERFORM STABILITY TEST.
                  IF (ABS(A(KJ)).LT.AMAX*U) GO TO 180
  170             JCOST = KCOST
                  IPP = I
                  JP = J
                  IF (JCOST.LE.(NZ-1)**2) GO TO 250
  180          CONTINUE
               J = IW(J,8)
  190       CONTINUE
C SEARCH ROWS WITH NZ NON-ZEROS.
  200       I = IW(NZ,3)
            DO 230 IDUMMY=1,N
               IF (I.LE.0) GO TO 240
               AMAX = 0.
               KP = IP(I,1)
               KL = KP + IW(I,1) - 1
C FIND LARGEST ELEMENT IN THE ROW
               DO 210 K=KP,KL
                  AMAX = MAX(ABS(A(K)),AMAX)
  210          CONTINUE
               AU = AMAX*U
               DO 220 K=KP,KL
C PERFORM STABILITY TEST.
                  IF (ABS(A(K)).LT.AU) GO TO 220
                  J = IND(K,2)
                  KCOST = (NZ-1)*(IW(J,2)-1)
                  IF (KCOST.GE.JCOST) GO TO 220
                  JCOST = KCOST
                  IPP = I
                  JP = J
                  IF (JCOST.LE.(NZ-1)**2) GO TO 250
  220          CONTINUE
               I = IW(I,7)
  230       CONTINUE
  240    CONTINUE
C
C PIVOT FOUND.
C REMOVE ROWS AND COLUMNS INVOLVED IN ELIMINATION FROM ORDERING VECTORS.
  250    KP = IP(JP,2)
         KL = IW(JP,2) + KP - 1
         DO 290 L=1,2
            DO 280 K=KP,KL
               I = IND(K,L)
               IL = IW(I,L+4)
               IN = IW(I,L+6)
               IF (IL.EQ.0) GO TO 260
               IW(IL,L+6) = IN
               GO TO 270
  260          NZ = IW(I,L)
               IW(NZ,L+2) = IN
  270          IF (IN.GT.0) IW(IN,L+4) = IL
  280       CONTINUE
            KP = IP(IPP,1)
            KL = KP + IW(IPP,1) - 1
  290    CONTINUE
C STORE PIVOT
         IW(IPP,5) = -IPV
         IW(JP,6) = -IPV
C ELIMINATE PIVOTAL ROW FROM COLUMN FILE AND FIND PIVOT IN ROW FILE.
         DO 320 K=KP,KL
            J = IND(K,2)
            KPC = IP(J,2)
            IW(J,2) = IW(J,2) - 1
            KLC = KPC + IW(J,2)
            DO 300 KC=KPC,KLC
               IF (IPP.EQ.IND(KC,1)) GO TO 310
  300       CONTINUE
  310       IND(KC,1) = IND(KLC,1)
            IND(KLC,1) = 0
            IF (J.EQ.JP) KR = K
  320    CONTINUE
C BRING PIVOT TO FRONT OF PIVOTAL ROW.
         AU = A(KR)
         A(KR) = A(KP)
         A(KP) = AU
         IND(KR,2) = IND(KP,2)
         IND(KP,2) = JP
C
C PERFORM ELIMINATION ITSELF, LOOPING ON NON-ZEROS IN PIVOT COLUMN.
         NZC = IW(JP,2)
         IF (NZC.EQ.0) GO TO 550
         DO 540 NC=1,NZC
            KC = IP(JP,2) + NC - 1
            IR = IND(KC,1)
C SEARCH NON-PIVOT ROW FOR ELEMENT TO BE ELIMINATED.
            KR = IP(IR,1)
            KRL = KR + IW(IR,1) - 1
            DO 330 KNP=KR,KRL
               IF (JP.EQ.IND(KNP,2)) GO TO 340
  330       CONTINUE
C BRING ELEMENT TO BE ELIMINATED TO FRONT OF ITS ROW.
  340       AM = A(KNP)
            A(KNP) = A(KR)
            A(KR) = AM
            IND(KNP,2) = IND(KR,2)
            IND(KR,2) = JP
            AM = -A(KR)/A(KP)
C COMPRESS ROW FILE UNLESS IT IS CERTAIN THAT THERE IS ROOM FOR NEW ROW.
            IF (LROW+IW(IR,1)+IW(IPP,1)+LENL.LE.IA) GO TO 350
            IF (NCP.GE.MCP .OR. LENU+IW(IR,1)+IW(IPP,1)+LENL.GT.IA) GO
     *       TO 710
            CALL LA05ED(A, IND(1,2), IP, N, IW, IA, .TRUE.)
            KP = IP(IPP,1)
            KR = IP(IR,1)
  350       KRL = KR + IW(IR,1) - 1
            KQ = KP + 1
            KPL = KP + IW(IPP,1) - 1
C PLACE PIVOT ROW (EXCLUDING PIVOT ITSELF) IN W.
            IF (KQ.GT.KPL) GO TO 370
            DO 360 K=KQ,KPL
               J = IND(K,2)
               W(J) = A(K)
  360       CONTINUE
  370       IP(IR,1) = LROW + 1
C
C TRANSFER MODIFIED ELEMENTS.
            IND(KR,2) = 0
            KR = KR + 1
            IF (KR.GT.KRL) GO TO 430
            DO 420 KS=KR,KRL
               J = IND(KS,2)
               AU = A(KS) + AM*W(J)
               IND(KS,2) = 0
C IF ELEMENT IS VERY SMALL REMOVE IT FROM U.
               IF (ABS(AU).LE.SMALL) GO TO 380
               G = MAX(G,ABS(AU))
               LROW = LROW + 1
               A(LROW) = AU
               IND(LROW,2) = J
               GO TO 410
  380          LENU = LENU - 1
C REMOVE ELEMENT FROM COL FILE.
               K = IP(J,2)
               KL = K + IW(J,2) - 1
               IW(J,2) = KL - K
               DO 390 KK=K,KL
                  IF (IND(KK,1).EQ.IR) GO TO 400
  390          CONTINUE
  400          IND(KK,1) = IND(KL,1)
               IND(KL,1) = 0
  410          W(J) = 0.
  420       CONTINUE
C
C SCAN PIVOT ROW FOR FILLS.
  430       IF (KQ.GT.KPL) GO TO 520
            DO 510 KS=KQ,KPL
               J = IND(KS,2)
               AU = AM*W(J)
               IF (ABS(AU).LE.SMALL) GO TO 500
               LROW = LROW + 1
               A(LROW) = AU
               IND(LROW,2) = J
               LENU = LENU + 1
C
C CREATE FILL IN COLUMN FILE.
               NZ = IW(J,2)
               K = IP(J,2)
               KL = K + NZ - 1
               IF (NZ .EQ. 0) GO TO 460
C IF POSSIBLE PLACE NEW ELEMENT AT END OF PRESENT ENTRY.
               IF (KL.NE.LCOL) GO TO 440
               IF (LCOL+LENL.GE.IA) GO TO 460
               LCOL = LCOL + 1
               GO TO 450
  440          IF (IND(KL+1,1).NE.0) GO TO 460
  450          IND(KL+1,1) = IR
               GO TO 490
C NEW ENTRY HAS TO BE CREATED.
  460          IF (LCOL+LENL+NZ+1.LT.IA) GO TO 470
C COMPRESS COLUMN FILE IF THERE IS NOT ROOM FOR NEW ENTRY.
               IF (NCP.GE.MCP .OR. LENU+LENL+NZ+1.GE.IA) GO TO 710
               CALL LA05ED(A, IND, IP(1,2), N, IW(1,2), IA, .FALSE.)
               K = IP(J,2)
               KL = K + NZ - 1
C TRANSFER OLD ENTRY INTO NEW.
  470          IP(J,2) = LCOL + 1
               IF (KL .LT. K) GO TO 485
               DO 480 KK=K,KL
                  LCOL = LCOL + 1
                  IND(LCOL,1) = IND(KK,1)
                  IND(KK,1) = 0
  480          CONTINUE
  485          CONTINUE
C ADD NEW ELEMENT.
               LCOL = LCOL + 1
               IND(LCOL,1) = IR
  490          G = MAX(G,ABS(AU))
               IW(J,2) = NZ + 1
  500          W(J) = 0.
  510       CONTINUE
  520       IW(IR,1) = LROW + 1 - IP(IR,1)
C
C STORE MULTIPLIER
            IF (LENL+LCOL+1.LE.IA) GO TO 530
C COMPRESS COL FILE IF NECESSARY.
            IF (NCP.GE.MCP) GO TO 710
            CALL LA05ED(A, IND, IP(1,2), N, IW(1,2), IA, .FALSE.)
  530       K = IA - LENL
            LENL = LENL + 1
            A(K) = AM
            IND(K,1) = IPP
            IND(K,2) = IR
            LENU = LENU - 1
  540    CONTINUE
C
C INSERT ROWS AND COLUMNS INVOLVED IN ELIMINATION IN LINKED LISTS
C     OF EQUAL NUMBERS OF NON-ZEROS.
  550    K1 = IP(JP,2)
         K2 = IW(JP,2) + K1 - 1
         IW(JP,2) = 0
         DO 580 L=1,2
            IF (K2.LT.K1) GO TO 570
            DO 560 K=K1,K2
               IR = IND(K,L)
               IF (L.EQ.1) IND(K,L) = 0
               NZ = IW(IR,L)
               IF (NZ.LE.0) GO TO 720
               IN = IW(NZ,L+2)
               IW(IR,L+6) = IN
               IW(IR,L+4) = 0
               IW(NZ,L+2) = IR
               IF (IN.NE.0) IW(IN,L+4) = IR
  560       CONTINUE
  570       K1 = IP(IPP,1) + 1
            K2 = IW(IPP,1) + K1 - 2
  580    CONTINUE
  590 CONTINUE
C
C RESET COLUMN FILE TO REFER TO U AND STORE ROW/COL NUMBERS IN
C     PIVOTAL ORDER IN IW(.,3),IW(.,4)
      DO 600 I=1,N
         J = -IW(I,5)
         IW(J,3) = I
         J = -IW(I,6)
         IW(J,4) = I
         IW(I,2) = 0
  600 CONTINUE
      DO 620 I=1,N
         KP = IP(I,1)
         KL = IW(I,1) + KP - 1
         DO 610 K=KP,KL
            J = IND(K,2)
            IW(J,2) = IW(J,2) + 1
  610    CONTINUE
  620 CONTINUE
      K = 1
      DO 630 I=1,N
         K = K + IW(I,2)
         IP(I,2) = K
  630 CONTINUE
      LCOL = K - 1
      DO 650 II=1,N
         I = IW(II,3)
         KP = IP(I,1)
         KL = IW(I,1) + KP - 1
         DO 640 K=KP,KL
            J = IND(K,2)
            KN = IP(J,2) - 1
            IP(J,2) = KN
            IND(KN,1) = I
  640    CONTINUE
  650 CONTINUE
      RETURN
C
C     THE FOLLOWING INSTRUCTIONS IMPLEMENT THE FAILURE EXITS.
C
  660 IF (LP.GT.0) THEN
         WRITE (XERN1, '(I8)') IR
         WRITE (XERN2, '(I8)') J
         CALL XERMSG ('SLATEC', 'LA05AD', 'MORE THAN ONE MATRIX ' //
     *      'ENTRY.  HERE ROW = ' // XERN1 // ' AND COL = ' // XERN2,
     *      -4, 1)
      ENDIF
      G = -4.
      RETURN
C
  670 IF (LP.GT.0) CALL XERMSG ('SLATEC', 'LA05AD',
     *   'THE ORDER OF THE SYSTEM, N, IS NOT POSITIVE.', -1, 1)
      G = -1.0D0
      RETURN
C
  680 IF (LP.GT.0) THEN
         WRITE (XERN0, '(I8)') K
         WRITE (XERN1, '(I8)') I
         WRITE (XERN2, '(I8)') J
         CALL XERMSG ('SLATEC', 'LA05AD', 'ELEMENT K = ' // XERN0 //
     *      ' IS OUT OF BOUNDS.$$HERE ROW = ' // XERN1 //
     *      ' AND COL = ' // XERN2, -3, 1)
      ENDIF
      G = -3.
      RETURN
C
  700 IF (LP.GT.0) THEN
         WRITE (XERN1, '(I8)') L
         CALL XERMSG ('SLATEC', 'LA05AD', 'ROW OR COLUMN HAS NO ' //
     *      'ELEMENTS.  HERE INDEX = ' // XERN1, -2, 1)
      ENDIF
      G = -2.
      RETURN
C
  710 IF (LP.GT.0) CALL XERMSG ('SLATEC', 'LA05AD',
     *   'LENGTHS OF ARRAYS A(*) AND IND(*,2) ARE TOO SMALL.', -7, 1)
      G = -7.
      RETURN
C
  720 IPV = IPV + 1
      IW(IPV,1) = IR
      DO 730 I=1,N
         II = -IW(I,L+4)
         IF (II.GT.0) IW(II,1) = I
  730 CONTINUE
C
      IF (LP.GT.0) THEN
         XERN1 = 'ROWS'
         IF (L.EQ.2) XERN1 = 'COLUMNS'
         CALL XERMSG ('SLATEC', 'LA05AD', 'DEPENDANT ' // XERN1, -5, 1)
C
  740    WRITE (XERN1, '(I8)') IW(I,1)
         XERN2 = ' '
         IF (I+1.LE.IPV) WRITE (XERN2, '(I8)') IW(I+1,1)
         CALL XERMSG ('SLATEC', 'LA05AD',
     *      'DEPENDENT VECTOR INDICES ARE ' // XERN1 // ' AND ' //
     *      XERN2, -5, 1)
         I = I + 2
         IF (I.LE.IPV) GO TO 740
      ENDIF
      G = -5.
      RETURN
      END
