*DECK LA05BS
      SUBROUTINE LA05BS (A, IND, IA, N, IP, IW, W, G, B, TRANS)
C***BEGIN PROLOGUE  LA05BS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SPLP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LA05BS-S, LA05BD-D)
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
C IP(I,1),IP(I,2) POINT TO START OF ROW/COLUMN I OF U.
C IW(I,1),IW(I,2) ARE LENGTHS OF ROW/COL I OF U.
C IW(.,3),IW(.,4) HOLD ROW/COL NUMBERS IN PIVOTAL ORDER.
C
C***SEE ALSO  SPLP
C***ROUTINES CALLED  XERMSG, XSETUN
C***COMMON BLOCKS    LA05DS
C***REVISION HISTORY  (YYMMDD)
C   811215  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900402  Added TYPE section.  (WRB)
C   920410  Corrected second dimension on IW declaration.  (WRB)
C***END PROLOGUE  LA05BS
      REAL A(IA), B(*), AM, W(*), G, SMALL
      LOGICAL TRANS
      INTEGER IND(IA,2), IW(N,8)
      INTEGER IP(N,2)
      COMMON /LA05DS/ SMALL, LP, LENL, LENU, NCP, LROW, LCOL
C***FIRST EXECUTABLE STATEMENT  LA05BS
      IF (G.LT.0.) GO TO 130
      KLL = IA - LENL + 1
      IF (TRANS) GO TO 80
C
C     MULTIPLY VECTOR BY INVERSE OF L
      IF (LENL.LE.0) GO TO 20
      L1 = IA + 1
      DO 10 KK=1,LENL
         K = L1 - KK
         I = IND(K,1)
         IF (B(I).EQ.0.) GO TO 10
         J = IND(K,2)
         B(J) = B(J) + A(K)*B(I)
   10 CONTINUE
   20 DO 30 I=1,N
         W(I) = B(I)
         B(I) = 0.
   30 CONTINUE
C
C     MULTIPLY VECTOR BY INVERSE OF U
      N1 = N + 1
      DO 70 II=1,N
         I = N1 - II
         I = IW(I,3)
         AM = W(I)
         KP = IP(I,1)
         IF (KP.GT.0) GO TO 50
         KP = -KP
         IP(I,1) = KP
         NZ = IW(I,1)
         KL = KP - 1 + NZ
         K2 = KP + 1
         DO 40 K=K2,KL
            J = IND(K,2)
            AM = AM - A(K)*B(J)
   40    CONTINUE
   50    IF (AM.EQ.0.) GO TO 70
         J = IND(KP,2)
         B(J) = AM/A(KP)
         KPC = IP(J,2)
         KL = IW(J,2) + KPC - 1
         IF (KL.EQ.KPC) GO TO 70
         K2 = KPC + 1
         DO 60 K=K2,KL
            I = IND(K,1)
            IP(I,1) = -ABS(IP(I,1))
   60    CONTINUE
   70 CONTINUE
      GO TO 140
C
C     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF U
   80 DO 90 I=1,N
         W(I) = B(I)
         B(I) = 0.
   90 CONTINUE
      DO 110 II=1,N
         I = IW(II,4)
         AM = W(I)
         IF (AM.EQ.0.) GO TO 110
         J = IW(II,3)
         KP = IP(J,1)
         AM = AM/A(KP)
         B(J) = AM
         KL = IW(J,1) + KP - 1
         IF (KP.EQ.KL) GO TO 110
         K2 = KP + 1
         DO 100 K=K2,KL
            I = IND(K,2)
            W(I) = W(I) - AM*A(K)
  100    CONTINUE
  110 CONTINUE
C
C     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF L
      IF (KLL.GT.IA) RETURN
      DO 120 K=KLL,IA
         J = IND(K,2)
         IF (B(J).EQ.0.) GO TO 120
         I = IND(K,1)
         B(I) = B(I) + A(K)*B(J)
  120 CONTINUE
      GO TO 140
C
  130 CALL XSETUN(LP)
      IF (LP .GT. 0) CALL XERMSG ('SLATEC', 'LA05BS',
     +   'EARLIER ENTRY GAVE ERROR RETURN.', -8, 2)
  140 RETURN
      END
