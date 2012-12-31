*DECK CPRODP
      SUBROUTINE CPRODP (ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, YY, M,
     +   A, B, C, D, U, Y)
C***BEGIN PROLOGUE  CPRODP
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BLKTRI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (CPRODP-S, CPROCP-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C PRODP applies a sequence of matrix operations to the vector X and
C stores the result in YY. (Periodic boundary conditions and COMPLEX
C case)
C
C BD,BM1,BM2     are arrays containing roots of certain B polynomials.
C ND,NM1,NM2     are the lengths of the arrays BD,BM1,BM2 respectively.
C AA             Array containing scalar multipliers of the vector X.
C NA             is the length of the array AA.
C X,YY      The matrix operations are applied to X and the result is YY.
C A,B,C          are arrays which contain the tridiagonal matrix.
C M              is the order of the matrix.
C D,U,Y          are working arrays.
C ISGN           determines whether or not a change in sign is made.
C
C***SEE ALSO  BLKTRI
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CPRODP
C
      COMPLEX         Y          ,D          ,U          ,V          ,
     1                DEN        ,BH         ,YM         ,AM         ,
     2                Y1         ,Y2         ,YH         ,BD         ,
     3                CRT
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       ,
     1                Y(*)       ,D(*)       ,U(*)       ,BD(*)      ,
     2                BM1(*)     ,BM2(*)     ,AA(*)      ,YY(*)
C***FIRST EXECUTABLE STATEMENT  CPRODP
      DO 101 J=1,M
         Y(J) = CMPLX(X(J),0.)
  101 CONTINUE
      MM = M-1
      MM2 = M-2
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 IFLG = 0
      IF (ID) 111,111,103
  103 CRT = BD(ID)
      ID = ID-1
      IFLG = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      BH = B(M)-CRT
      YM = Y(M)
      DEN = B(1)-CRT
      D(1) = C(1)/DEN
      U(1) = A(1)/DEN
      Y(1) = Y(1)/DEN
      V = CMPLX(C(M),0.)
      IF (MM2-2) 106,104,104
  104 DO 105 J=2,MM2
         DEN = B(J)-CRT-A(J)*D(J-1)
         D(J) = C(J)/DEN
         U(J) = -A(J)*U(J-1)/DEN
         Y(J) = (Y(J)-A(J)*Y(J-1))/DEN
         BH = BH-V*U(J-1)
         YM = YM-V*Y(J-1)
         V = -V*D(J-1)
  105 CONTINUE
  106 DEN = B(M-1)-CRT-A(M-1)*D(M-2)
      D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
      Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/DEN
      AM = A(M)-V*D(M-2)
      BH = BH-V*U(M-2)
      YM = YM-V*Y(M-2)
      DEN = BH-AM*D(M-1)
      IF (ABS(DEN)) 107,108,107
  107 Y(M) = (YM-AM*Y(M-1))/DEN
      GO TO 109
  108 Y(M) = (1.,0.)
  109 Y(M-1) = Y(M-1)-D(M-1)*Y(M)
      DO 110 J=2,MM
         K = M-J
         Y(K) = Y(K)-D(K)*Y(K+1)-U(K)*Y(M)
  110 CONTINUE
  111 IF (M1) 112,112,114
  112 IF (M2) 123,123,113
  113 RT = BM2(M2)
      M2 = M2-1
      GO TO 119
  114 IF (M2) 115,115,116
  115 RT = BM1(M1)
      M1 = M1-1
      GO TO 119
  116 IF (ABS(BM1(M1))-ABS(BM2(M2))) 118,118,117
  117 RT = BM1(M1)
      M1 = M1-1
      GO TO 119
  118 RT = BM2(M2)
      M2 = M2-1
C
C MATRIX MULTIPLICATION
C
  119 YH = Y(1)
      Y1 = (B(1)-RT)*Y(1)+C(1)*Y(2)+A(1)*Y(M)
      IF (MM-2) 122,120,120
  120 DO 121 J=2,MM
         Y2 = A(J)*Y(J-1)+(B(J)-RT)*Y(J)+C(J)*Y(J+1)
         Y(J-1) = Y1
         Y1 = Y2
  121 CONTINUE
  122 Y(M) = A(M)*Y(M-1)+(B(M)-RT)*Y(M)+C(M)*YH
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  123 IF (IA) 126,126,124
  124 RT = AA(IA)
      IA = IA-1
      IFLG = 1
C
C SCALAR MULTIPLICATION
C
      DO 125 J=1,M
         Y(J) = RT*Y(J)
  125 CONTINUE
  126 IF (IFLG) 127,127,102
  127 DO 128 J=1,M
         YY(J) = REAL(Y(J))
  128 CONTINUE
      RETURN
      END
