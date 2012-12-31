*DECK PROCP
      SUBROUTINE PROCP (ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, Y, M, A,
     +   B, C, D, U, W)
C***BEGIN PROLOGUE  PROCP
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CBLKTR
C***LIBRARY   SLATEC
C***TYPE      COMPLEX (PRODP-C, PROCP-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C PROCP applies a sequence of matrix operations to the vector X and
C stores the result in Y (periodic boundary conditions).
C
C BD,BM1,BM2 are arrays containing roots of certain B polynomials.
C ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively.
C AA         Array containing scalar multipliers of the vector X.
C NA         is the length of the array AA.
C X,Y        The matrix operations are applied to X and the result is Y.
C A,B,C      are arrays which contain the tridiagonal matrix.
C M          is the order of the matrix.
C D,U,W      are working arrays.
C IS         determines whether or not a change in sign is made.
C
C***SEE ALSO  CBLKTR
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  PROCP
C
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       ,
     1                Y(*)       ,D(*)       ,U(*)       ,BD(*)      ,
     2                BM1(*)     ,BM2(*)     ,AA(*)      ,W(*)
      COMPLEX         X          ,Y          ,A          ,B          ,
     1                C          ,D          ,U          ,W          ,
     2                DEN        ,YM         ,V          ,BH         ,AM
C***FIRST EXECUTABLE STATEMENT  PROCP
      DO 101 J=1,M
         Y(J) = X(J)
         W(J) = Y(J)
  101 CONTINUE
      MM = M-1
      MM2 = M-2
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 IF (IA) 105,105,103
  103 RT = AA(IA)
      IF (ND .EQ. 0) RT = -RT
      IA = IA-1
      DO 104 J=1,M
         Y(J) = RT*W(J)
  104 CONTINUE
  105 IF (ID) 128,128,106
  106 RT = BD(ID)
      ID = ID-1
      IF (ID .EQ. 0) IBR = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      BH = B(M)-RT
      YM = Y(M)
      DEN = B(1)-RT
      D(1) = C(1)/DEN
      U(1) = A(1)/DEN
      W(1) = Y(1)/DEN
      V = C(M)
      IF (MM2-2) 109,107,107
  107 DO 108 J=2,MM2
         DEN = B(J)-RT-A(J)*D(J-1)
         D(J) = C(J)/DEN
         U(J) = -A(J)*U(J-1)/DEN
         W(J) = (Y(J)-A(J)*W(J-1))/DEN
         BH = BH-V*U(J-1)
         YM = YM-V*W(J-1)
         V = -V*D(J-1)
  108 CONTINUE
  109 DEN = B(M-1)-RT-A(M-1)*D(M-2)
      D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
      W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/DEN
      AM = A(M)-V*D(M-2)
      BH = BH-V*U(M-2)
      YM = YM-V*W(M-2)
      DEN = BH-AM*D(M-1)
      IF (ABS(DEN)) 110,111,110
  110 W(M) = (YM-AM*W(M-1))/DEN
      GO TO 112
  111 W(M) = (1.,0.)
  112 W(M-1) = W(M-1)-D(M-1)*W(M)
      DO 113 J=2,MM
         K = M-J
         W(K) = W(K)-D(K)*W(K+1)-U(K)*W(M)
  113 CONTINUE
      IF (NA) 116,116,102
  114 DO 115 J=1,M
         Y(J) = W(J)
  115 CONTINUE
      IBR = 1
      GO TO 102
  116 IF (M1) 117,117,118
  117 IF (M2) 114,114,123
  118 IF (M2) 120,120,119
  119 IF (ABS(BM1(M1))-ABS(BM2(M2))) 123,123,120
  120 IF (IBR) 121,121,122
  121 IF (ABS(BM1(M1)-BD(ID))-ABS(BM1(M1)-RT)) 114,122,122
  122 RT = RT-BM1(M1)
      M1 = M1-1
      GO TO 126
  123 IF (IBR) 124,124,125
  124 IF (ABS(BM2(M2)-BD(ID))-ABS(BM2(M2)-RT)) 114,125,125
  125 RT = RT-BM2(M2)
      M2 = M2-1
  126 DO 127 J=1,M
         Y(J) = Y(J)+RT*W(J)
  127 CONTINUE
      GO TO 102
  128 RETURN
      END
