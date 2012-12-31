*DECK DU11LS
      SUBROUTINE DU11LS (A, MDA, M, N, UB, DB, MODE, NP, KRANK, KSURE,
     +   H, W, EB, IC, IR)
C***BEGIN PROLOGUE  DU11LS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DLLSIA
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (U11LS-S, DU11LS-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   **** Double Precision version of U11LS ****
C
C       This routine performs a QR factorization of A
C       using Householder transformations. Row and
C       column pivots are chosen to reduce the growth
C       of round-off and to help detect possible rank
C       deficiency.
C
C***SEE ALSO  DLLSIA
C***ROUTINES CALLED  DAXPY, DDOT, DNRM2, DSCAL, DSWAP, IDAMAX, ISWAP,
C                    XERMSG
C***REVISION HISTORY  (YYMMDD)
C   810801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced variable.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DU11LS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DDOT,DNRM2
      DIMENSION A(MDA,*),UB(*),DB(*),H(*),W(*),EB(*)
      INTEGER IC(*),IR(*)
C
C        INITIALIZATION
C
C***FIRST EXECUTABLE STATEMENT  DU11LS
      J=0
      KRANK=N
      DO 10 I=1,N
      IC(I)=I
   10 CONTINUE
      DO 12 I=1,M
      IR(I)=I
   12 CONTINUE
C
C        DETERMINE REL AND ABS ERROR VECTORS
C
C
C
C        CALCULATE COL LENGTH
C
      DO 30 I=1,N
      H(I)=DNRM2(M,A(1,I),1)
      W(I)=H(I)
   30 CONTINUE
C
C         INITIALIZE ERROR BOUNDS
C
      DO  40 I=1,N
      EB(I)=MAX(DB(I),UB(I)*H(I))
      UB(I)=EB(I)
      DB(I)=0.0D0
   40 CONTINUE
C
C          DISCARD SELF DEPENDENT COLUMNS
C
      I=1
   50 IF(EB(I).GE.H(I)) GO TO 60
      IF(I.EQ.KRANK) GO TO 70
      I=I+1
      GO TO 50
C
C          MATRIX REDUCTION
C
   60 CONTINUE
      KK=KRANK
      KRANK=KRANK-1
      IF(MODE.EQ.0) RETURN
      IF(I.GT.NP) GO TO  64
      CALL XERMSG ('SLATEC', 'DU11LS',
     +   'FIRST NP COLUMNS ARE LINEARLY DEPENDENT', 8, 0)
      KRANK=I-1
      RETURN
   64 CONTINUE
      IF(I.GT.KRANK) GO TO 70
      CALL DSWAP(1,EB(I),1,EB(KK),1)
      CALL DSWAP(1,UB(I),1,UB(KK),1)
      CALL DSWAP(1,W(I),1,W(KK),1)
      CALL DSWAP(1,H(I),1,H(KK),1)
      CALL ISWAP(1,IC(I),1,IC(KK),1)
      CALL DSWAP(M,A(1,I),1,A(1,KK),1)
      GO TO 50
C
C           TEST FOR ZERO RANK
C
   70 IF(KRANK.GT.0) GO TO 80
      KRANK=0
      KSURE=0
      RETURN
   80 CONTINUE
C
C        M A I N    L O O P
C
  110 CONTINUE
      J=J+1
      JP1=J+1
      JM1=J-1
      KZ=KRANK
      IF(J.LE.NP) KZ=J
C
C        EACH COL HAS MM=M-J+1 COMPONENTS
C
      MM=M-J+1
C
C         UB DETERMINES COLUMN PIVOT
C
  115 IMIN=J
      IF(H(J).EQ.0.D0) GO TO 170
      RMIN=UB(J)/H(J)
      DO 120 I=J,KZ
      IF(UB(I).GE.H(I)*RMIN) GO TO 120
      RMIN=UB(I)/H(I)
      IMIN=I
  120 CONTINUE
C
C     TEST FOR RANK DEFICIENCY
C
      IF(RMIN.LT.1.0D0) GO TO 200
      TT=(EB(IMIN)+ABS(DB(IMIN)))/H(IMIN)
      IF(TT.GE.1.0D0) GO TO 170
C     COMPUTE EXACT UB
      DO 125 I=1,JM1
      W(I)=A(I,IMIN)
  125 CONTINUE
      L=JM1
  130 W(L)=W(L)/A(L,L)
      IF(L.EQ.1) GO TO 150
      LM1=L-1
      DO 140 I=L,JM1
      W(LM1)=W(LM1)-A(LM1,I)*W(I)
  140 CONTINUE
      L=LM1
      GO TO 130
  150 TT=EB(IMIN)
      DO 160 I=1,JM1
      TT=TT+ABS(W(I))*EB(I)
  160 CONTINUE
      UB(IMIN)=TT
      IF(UB(IMIN)/H(IMIN).GE.1.0D0) GO TO 170
      GO TO 200
C
C        MATRIX REDUCTION
C
  170 CONTINUE
      KK=KRANK
      KRANK=KRANK-1
      KZ=KRANK
      IF(MODE.EQ.0) RETURN
      IF(J.GT.NP) GO TO 172
      CALL XERMSG ('SLATEC', 'DU11LS',
     +   'FIRST NP COLUMNS ARE LINEARLY DEPENDENT', 8, 0)
      KRANK=J-1
      RETURN
  172 CONTINUE
      IF(IMIN.GT.KRANK) GO TO 180
      CALL ISWAP(1,IC(IMIN),1,IC(KK),1)
      CALL DSWAP(M,A(1,IMIN),1,A(1,KK),1)
      CALL DSWAP(1,EB(IMIN),1,EB(KK),1)
      CALL DSWAP(1,UB(IMIN),1,UB(KK),1)
      CALL DSWAP(1,DB(IMIN),1,DB(KK),1)
      CALL DSWAP(1,W(IMIN),1,W(KK),1)
      CALL DSWAP(1,H(IMIN),1,H(KK),1)
  180 IF(J.GT.KRANK) GO TO 300
      GO TO 115
C
C        COLUMN PIVOT
C
  200 IF(IMIN.EQ.J) GO TO 230
      CALL DSWAP(1,H(J),1,H(IMIN),1)
      CALL DSWAP(M,A(1,J),1,A(1,IMIN),1)
      CALL DSWAP(1,EB(J),1,EB(IMIN),1)
      CALL DSWAP(1,UB(J),1,UB(IMIN),1)
      CALL DSWAP(1,DB(J),1,DB(IMIN),1)
      CALL DSWAP(1,W(J),1,W(IMIN),1)
      CALL ISWAP(1,IC(J),1,IC(IMIN),1)
C
C        ROW PIVOT
C
  230 CONTINUE
      JMAX=IDAMAX(MM,A(J,J),1)
      JMAX=JMAX+J-1
      IF(JMAX.EQ.J) GO TO 240
      CALL DSWAP(N,A(J,1),MDA,A(JMAX,1),MDA)
      CALL ISWAP(1,IR(J),1,IR(JMAX),1)
  240 CONTINUE
C
C     APPLY HOUSEHOLDER TRANSFORMATION
C
      TN=DNRM2(MM,A(J,J),1)
      IF(TN.EQ.0.0D0) GO TO 170
      IF(A(J,J).NE.0.0D0) TN=SIGN(TN,A(J,J))
      CALL DSCAL(MM,1.0D0/TN,A(J,J),1)
      A(J,J)=A(J,J)+1.0D0
      IF(J.EQ.N) GO TO 250
      DO 248 I=JP1,N
      BB=-DDOT(MM,A(J,J),1,A(J,I),1)/A(J,J)
      CALL DAXPY(MM,BB,A(J,J),1,A(J,I),1)
      IF(I.LE.NP) GO TO 248
      IF(H(I).EQ.0.0D0) GO TO 248
      TT=1.0D0-(ABS(A(J,I))/H(I))**2
      TT=MAX(TT,0.0D0)
      T=TT
      TT=1.0D0+.05D0*TT*(H(I)/W(I))**2
      IF(TT.EQ.1.0D0) GO TO 244
      H(I)=H(I)*SQRT(T)
      GO TO 246
  244 CONTINUE
      H(I)=DNRM2(M-J,A(J+1,I),1)
      W(I)=H(I)
  246 CONTINUE
  248 CONTINUE
  250 CONTINUE
      H(J)=A(J,J)
      A(J,J)=-TN
C
C
C          UPDATE UB, DB
C
      UB(J)=UB(J)/ABS(A(J,J))
      DB(J)=(SIGN(EB(J),DB(J))+DB(J))/A(J,J)
      IF(J.EQ.KRANK) GO TO 300
      DO 260 I=JP1,KRANK
      UB(I)=UB(I)+ABS(A(J,I))*UB(J)
      DB(I)=DB(I)-A(J,I)*DB(J)
  260 CONTINUE
      GO TO 110
C
C        E N D    M A I N    L O O P
C
  300 CONTINUE
C
C        COMPUTE KSURE
C
      KM1=KRANK-1
      DO 318 I=1,KM1
      IS=0
      KMI=KRANK-I
      DO 315 II=1,KMI
      IF(UB(II).LE.UB(II+1)) GO TO 315
      IS=1
      TEMP=UB(II)
      UB(II)=UB(II+1)
      UB(II+1)=TEMP
  315 CONTINUE
      IF(IS.EQ.0) GO TO 320
  318 CONTINUE
  320 CONTINUE
      KSURE=0
      SUM=0.0D0
      DO 328 I=1,KRANK
      R2=UB(I)*UB(I)
      IF(R2+SUM.GE.1.0D0) GO TO 330
      SUM=SUM+R2
      KSURE=KSURE+1
  328 CONTINUE
  330 CONTINUE
C
C     IF SYSTEM IS OF REDUCED RANK AND MODE = 2
C     COMPLETE THE DECOMPOSITION FOR SHORTEST LEAST SQUARES SOLUTION
C
      IF(KRANK.EQ.N .OR. MODE.LT.2) GO TO 360
      NMK=N-KRANK
      KP1=KRANK+1
      I=KRANK
  340 TN=DNRM2(NMK,A(I,KP1),MDA)/A(I,I)
      TN=A(I,I)*SQRT(1.0D0+TN*TN)
      CALL DSCAL(NMK,1.0D0/TN,A(I,KP1),MDA)
      W(I)=A(I,I)/TN+1.0D0
      A(I,I)=-TN
      IF(I.EQ.1) GO TO 350
      IM1=I-1
      DO 345 II=1,IM1
      TT=-DDOT(NMK,A(II,KP1),MDA,A(I,KP1),MDA)/W(I)
      TT=TT-A(II,I)
      CALL DAXPY(NMK,TT,A(I,KP1),MDA,A(II,KP1),MDA)
      A(II,I)=A(II,I)+TT*W(I)
  345 CONTINUE
      I=I-1
      GO TO 340
  350 CONTINUE
  360 CONTINUE
      RETURN
      END
