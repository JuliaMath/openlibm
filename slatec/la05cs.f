*DECK LA05CS
      SUBROUTINE LA05CS (A, IND, IA, N, IP, IW, W, G, U, MM)
C***BEGIN PROLOGUE  LA05CS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SPLP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LA05CS-S, LA05CD-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
C     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
C     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
C     THE FINAL LETTER =S= IN THE NAMES USED HERE.
C     REVISED SEP. 13, 1979.
C
C     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
C     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
C     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
C     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
C     SPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
C
C***SEE ALSO  SPLP
C***ROUTINES CALLED  LA05ES, XERMSG, XSETUN
C***COMMON BLOCKS    LA05DS
C***REVISION HISTORY  (YYMMDD)
C   811215  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890605  Corrected references to XERRWV.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900402  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920410  Corrected second dimension on IW declaration.  (WRB)
C   920422  Changed upper limit on DO from LAST to LAST-1.  (WRB)
C***END PROLOGUE  LA05CS
      REAL A(*), G, U, AM, W(*), SMALL, AU
      INTEGER IND(IA,2), IW(N,8)
      INTEGER IP(N,2)
      CHARACTER*8 XERN1
C
      COMMON /LA05DS/ SMALL, LP, LENL, LENU, NCP, LROW, LCOL
C***FIRST EXECUTABLE STATEMENT  LA05CS
      CALL XSETUN(LP)
      IF (G.LT.0.0E0) GO TO 620
      JM = MM
C MCP LIMITS THE VALUE OF NCP PERMITTED BEFORE AN ERROR RETURN RESULTS.
      MCP = NCP + 20
C REMOVE OLD COLUMN
      LENU = LENU - IW(JM,2)
      KP = IP(JM,2)
      IM = IND(KP,1)
      KL = KP + IW(JM,2) - 1
      IW(JM,2) = 0
      DO 30 K=KP,KL
         I = IND(K,1)
         IND(K,1) = 0
         KR = IP(I,1)
         NZ = IW(I,1) - 1
         IW(I,1) = NZ
         KRL = KR + NZ
         DO 10 KM=KR,KRL
            IF (IND(KM,2).EQ.JM) GO TO 20
   10    CONTINUE
   20    A(KM) = A(KRL)
         IND(KM,2) = IND(KRL,2)
         IND(KRL,2) = 0
   30 CONTINUE
C
C INSERT NEW COLUMN
      DO 110 II=1,N
         I = IW(II,3)
         IF (I.EQ.IM) M = II
         IF (ABS(W(I)).LE.SMALL) GO TO 100
         LENU = LENU + 1
         LAST = II
         IF (LCOL+LENL.LT.IA) GO TO 40
C COMPRESS COLUMN FILE IF NECESSARY.
         IF (NCP.GE.MCP .OR. LENL+LENU.GE.IA) GO TO 610
         CALL LA05ES(A, IND, IP(1,2), N, IW(1,2), IA, .FALSE.)
   40    LCOL = LCOL + 1
         NZ = IW(JM,2)
         IF (NZ.EQ.0) IP(JM,2) = LCOL
         IW(JM,2) = NZ + 1
         IND(LCOL,1) = I
         NZ = IW(I,1)
         KPL = IP(I,1) + NZ
         IF (KPL.GT.LROW) GO TO 50
         IF (IND(KPL,2).EQ.0) GO TO 90
C NEW ENTRY HAS TO BE CREATED.
   50    IF (LENL+LROW+NZ.LT.IA) GO TO 60
         IF (NCP.GE.MCP .OR. LENL+LENU+NZ.GE.IA) GO TO 610
C COMPRESS ROW FILE IF NECESSARY.
         CALL LA05ES(A, IND(1,2), IP, N, IW, IA, .TRUE.)
   60    KP = IP(I,1)
         IP(I,1) = LROW + 1
         IF (NZ.EQ.0) GO TO 80
         KPL = KP + NZ - 1
         DO 70 K=KP,KPL
            LROW = LROW + 1
            A(LROW) = A(K)
            IND(LROW,2) = IND(K,2)
            IND(K,2) = 0
   70    CONTINUE
   80    LROW = LROW + 1
         KPL = LROW
C PLACE NEW ELEMENT AT END OF ROW.
   90    IW(I,1) = NZ + 1
         A(KPL) = W(I)
         IND(KPL,2) = JM
  100    W(I) = 0.0E0
  110 CONTINUE
      IF (IW(IM,1).EQ.0 .OR. IW(JM,2).EQ.0 .OR. M.GT.LAST) GO TO 590
C
C FIND COLUMN SINGLETONS, OTHER THAN THE SPIKE. NON-SINGLETONS ARE
C     MARKED WITH W(J)=1. ONLY IW(.,3) IS REVISED AND IW(.,4) IS USED
C     FOR WORKSPACE.
      INS = M
      M1 = M
      W(JM) = 1.0E0
      DO 140 II=M,LAST
         I = IW(II,3)
         J = IW(II,4)
         IF (W(J).EQ.0.0E0) GO TO 130
         KP = IP(I,1)
         KL = KP + IW(I,1) - 1
         DO 120 K=KP,KL
            J = IND(K,2)
            W(J) = 1.0E0
  120    CONTINUE
         IW(INS,4) = I
         INS = INS + 1
         GO TO 140
C PLACE SINGLETONS IN NEW POSITION.
  130    IW(M1,3) = I
         M1 = M1 + 1
  140 CONTINUE
C PLACE NON-SINGLETONS IN NEW POSITION.
      IJ = M + 1
      DO 150 II=M1,LAST-1
         IW(II,3) = IW(IJ,4)
         IJ = IJ + 1
  150 CONTINUE
C PLACE SPIKE AT END.
      IW(LAST,3) = IM
C
C FIND ROW SINGLETONS, APART FROM SPIKE ROW. NON-SINGLETONS ARE MARKED
C     WITH W(I)=2. AGAIN ONLY IW(.,3) IS REVISED AND IW(.,4) IS USED
C     FOR WORKSPACE.
      LAST1 = LAST
      JNS = LAST
      W(IM) = 2.0E0
      J = JM
      DO 180 IJ=M1,LAST
         II = LAST + M1 - IJ
         I = IW(II,3)
         IF (W(I).NE.2.0E0) GO TO 170
         K = IP(I,1)
         IF (II.NE.LAST) J = IND(K,2)
         KP = IP(J,2)
         KL = KP + IW(J,2) - 1
         IW(JNS,4) = I
         JNS = JNS - 1
         DO 160 K=KP,KL
            I = IND(K,1)
            W(I) = 2.0E0
  160    CONTINUE
         GO TO 180
  170    IW(LAST1,3) = I
         LAST1 = LAST1 - 1
  180 CONTINUE
      DO 190 II=M1,LAST1
         JNS = JNS + 1
         I = IW(JNS,4)
         W(I) = 3.0E0
         IW(II,3) = I
  190 CONTINUE
C
C DEAL WITH SINGLETON SPIKE COLUMN. NOTE THAT BUMP ROWS ARE MARKED BY
C    W(I)=3.0E0
      DO 230 II=M1,LAST1
         KP = IP(JM,2)
         KL = KP + IW(JM,2) - 1
         IS = 0
         DO 200 K=KP,KL
            L = IND(K,1)
            IF (W(L).NE.3.0E0) GO TO 200
            IF (IS.NE.0) GO TO 240
            I = L
            KNP = K
            IS = 1
  200    CONTINUE
         IF (IS.EQ.0) GO TO 590
C MAKE A(I,JM) A PIVOT.
         IND(KNP,1) = IND(KP,1)
         IND(KP,1) = I
         KP = IP(I,1)
         DO 210 K=KP,IA
            IF (IND(K,2).EQ.JM) GO TO 220
  210    CONTINUE
  220    AM = A(KP)
         A(KP) = A(K)
         A(K) = AM
         IND(K,2) = IND(KP,2)
         IND(KP,2) = JM
         JM = IND(K,2)
         IW(II,4) = I
         W(I) = 2.0E0
  230 CONTINUE
      II = LAST1
      GO TO 260
  240 IN = M1
      DO 250 IJ=II,LAST1
         IW(IJ,4) = IW(IN,3)
         IN = IN + 1
  250 CONTINUE
  260 LAST2 = LAST1 - 1
      IF (M1.EQ.LAST1) GO TO 570
      DO 270 I=M1,LAST2
         IW(I,3) = IW(I,4)
  270 CONTINUE
      M1 = II
      IF (M1.EQ.LAST1) GO TO 570
C
C CLEAR W
      DO 280 I=1,N
         W(I) = 0.0E0
  280 CONTINUE
C
C PERFORM ELIMINATION
      IR = IW(LAST1,3)
      DO 560 II=M1,LAST1
         IPP = IW(II,3)
         KP = IP(IPP,1)
         KR = IP(IR,1)
         JP = IND(KP,2)
         IF (II.EQ.LAST1) JP = JM
C SEARCH NON-PIVOT ROW FOR ELEMENT TO BE ELIMINATED.
C  AND BRING IT TO FRONT OF ITS ROW
         KRL = KR + IW(IR,1) - 1
         DO 290 KNP=KR,KRL
            IF (JP.EQ.IND(KNP,2)) GO TO 300
  290    CONTINUE
         IF (II-LAST1) 560, 590, 560
C BRING ELEMENT TO BE ELIMINATED TO FRONT OF ITS ROW.
  300    AM = A(KNP)
         A(KNP) = A(KR)
         A(KR) = AM
         IND(KNP,2) = IND(KR,2)
         IND(KR,2) = JP
         IF (II.EQ.LAST1) GO TO 310
         IF (ABS(A(KP)).LT.U*ABS(AM)) GO TO 310
         IF (ABS(AM).LT.U*ABS(A(KP))) GO TO 340
         IF (IW(IPP,1).LE.IW(IR,1)) GO TO 340
C PERFORM INTERCHANGE
  310    IW(LAST1,3) = IPP
         IW(II,3) = IR
         IR = IPP
         IPP = IW(II,3)
         K = KR
         KR = KP
         KP = K
         KJ = IP(JP,2)
         DO 320 K=KJ,IA
            IF (IND(K,1).EQ.IPP) GO TO 330
  320    CONTINUE
  330    IND(K,1) = IND(KJ,1)
         IND(KJ,1) = IPP
  340    IF (A(KP).EQ.0.0E0) GO TO 590
         IF (II.EQ.LAST1) GO TO 560
         AM = -A(KR)/A(KP)
C COMPRESS ROW FILE UNLESS IT IS CERTAIN THAT THERE IS ROOM FOR NEW ROW.
         IF (LROW+IW(IR,1)+IW(IPP,1)+LENL.LE.IA) GO TO 350
         IF (NCP.GE.MCP .OR. LENU+IW(IR,1)+IW(IPP,1)+LENL.GT.IA) GO TO
     *    610
         CALL LA05ES(A, IND(1,2), IP, N, IW, IA, .TRUE.)
         KP = IP(IPP,1)
         KR = IP(IR,1)
  350    KRL = KR + IW(IR,1) - 1
         KQ = KP + 1
         KPL = KP + IW(IPP,1) - 1
C PLACE PIVOT ROW (EXCLUDING PIVOT ITSELF) IN W.
         IF (KQ.GT.KPL) GO TO 370
         DO 360 K=KQ,KPL
            J = IND(K,2)
            W(J) = A(K)
  360    CONTINUE
  370    IP(IR,1) = LROW + 1
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
  380       LENU = LENU - 1
C REMOVE ELEMENT FROM COL FILE.
            K = IP(J,2)
            KL = K + IW(J,2) - 1
            IW(J,2) = KL - K
            DO 390 KK=K,KL
               IF (IND(KK,1).EQ.IR) GO TO 400
  390       CONTINUE
  400       IND(KK,1) = IND(KL,1)
            IND(KL,1) = 0
  410       W(J) = 0.0E0
  420    CONTINUE
C
C SCAN PIVOT ROW FOR FILLS.
  430    IF (KQ.GT.KPL) GO TO 520
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
C IF POSSIBLE PLACE NEW ELEMENT AT END OF PRESENT ENTRY.
            IF (KL.NE.LCOL) GO TO 440
            IF (LCOL+LENL.GE.IA) GO TO 460
            LCOL = LCOL + 1
            GO TO 450
  440       IF (IND(KL+1,1).NE.0) GO TO 460
  450       IND(KL+1,1) = IR
            GO TO 490
C NEW ENTRY HAS TO BE CREATED.
  460       IF (LCOL+LENL+NZ+1.LT.IA) GO TO 470
C COMPRESS COLUMN FILE IF THERE IS NOT ROOM FOR NEW ENTRY.
            IF (NCP.GE.MCP .OR. LENU+LENL+NZ+1.GE.IA) GO TO 610
            CALL LA05ES(A, IND, IP(1,2), N, IW(1,2), IA, .FALSE.)
            K = IP(J,2)
            KL = K + NZ - 1
C TRANSFER OLD ENTRY INTO NEW.
  470       IP(J,2) = LCOL + 1
            DO 480 KK=K,KL
               LCOL = LCOL + 1
               IND(LCOL,1) = IND(KK,1)
               IND(KK,1) = 0
  480       CONTINUE
C ADD NEW ELEMENT.
            LCOL = LCOL + 1
            IND(LCOL,1) = IR
  490       G = MAX(G,ABS(AU))
            IW(J,2) = NZ + 1
  500       W(J) = 0.0E0
  510    CONTINUE
  520    IW(IR,1) = LROW + 1 - IP(IR,1)
C
C STORE MULTIPLIER
         IF (LENL+LCOL+1.LE.IA) GO TO 530
C COMPRESS COL FILE IF NECESSARY.
         IF (NCP.GE.MCP) GO TO 610
         CALL LA05ES(A, IND, IP(1,2), N, IW(1,2), IA, .FALSE.)
  530    K = IA - LENL
         LENL = LENL + 1
         A(K) = AM
         IND(K,1) = IPP
         IND(K,2) = IR
C CREATE BLANK IN PIVOTAL COLUMN.
         KP = IP(JP,2)
         NZ = IW(JP,2) - 1
         KL = KP + NZ
         DO 540 K=KP,KL
            IF (IND(K,1).EQ.IR) GO TO 550
  540    CONTINUE
  550    IND(K,1) = IND(KL,1)
         IW(JP,2) = NZ
         IND(KL,1) = 0
         LENU = LENU - 1
  560 CONTINUE
C
C CONSTRUCT COLUMN PERMUTATION AND STORE IT IN IW(.,4)
  570 DO 580 II=M,LAST
         I = IW(II,3)
         K = IP(I,1)
         J = IND(K,2)
         IW(II,4) = J
  580 CONTINUE
      RETURN
C
C     THE FOLLOWING INSTRUCTIONS IMPLEMENT THE FAILURE EXITS.
C
  590 IF (LP.GT.0) THEN
         WRITE (XERN1, '(I8)') MM
         CALL XERMSG ('SLATEC', 'LA05CS', 'SINGULAR MATRIX AFTER ' //
     *      'REPLACEMENT OF COLUMN.  INDEX = ' // XERN1, -6, 1)
      ENDIF
      G = -6.0E0
      RETURN
C
  610 IF (LP.GT.0) CALL XERMSG ('SLATEC', 'LA05CS',
     *   'LENGTHS OF ARRAYS A(*) AND IND(*,2) ARE TOO SMALL.', -7, 1)
      G = -7.0E0
      RETURN
C
  620 IF (LP.GT.0) CALL XERMSG ('SLATEC', 'LA05CS',
     *   'EARLIER ENTRY GAVE ERROR RETURN.', -8, 2)
      G = -8.0E0
      RETURN
      END
