*DECK PROD
      SUBROUTINE PROD (ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, Y, M, A,
     +   B, C, D, W, U)
C***BEGIN PROLOGUE  PROD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BLKTRI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (PROD-S, PROC-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C PROD applies a sequence of matrix operations to the vector X and
C stores the result in Y.
C
C BD,BM1,BM2 are arrays containing roots of certain B polynomials.
C ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively.
C AA         Array containing scalar multipliers of the vector X.
C NA         is the length of the array AA.
C X,Y        The matrix operations are applied to X and the result is Y.
C A,B,C      are arrays which contain the tridiagonal matrix.
C M          is the order of the matrix.
C D,W,U      are working arrays.
C IS         determines whether or not a change in sign is made.
C
C***SEE ALSO  BLKTRI
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  PROD
C
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       ,
     1                Y(*)       ,D(*)       ,W(*)       ,BD(*)      ,
     2                BM1(*)     ,BM2(*)     ,AA(*)      ,U(*)
C***FIRST EXECUTABLE STATEMENT  PROD
      DO 101 J=1,M
         W(J) = X(J)
         Y(J) = W(J)
  101 CONTINUE
      MM = M-1
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 IF (IA) 105,105,103
  103 RT = AA(IA)
      IF (ND .EQ. 0) RT = -RT
      IA = IA-1
C
C SCALAR MULTIPLICATION
C
      DO 104 J=1,M
         Y(J) = RT*W(J)
  104 CONTINUE
  105 IF (ID) 125,125,106
  106 RT = BD(ID)
      ID = ID-1
      IF (ID .EQ. 0) IBR = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      D(M) = A(M)/(B(M)-RT)
      W(M) = Y(M)/(B(M)-RT)
      DO 107 J=2,MM
         K = M-J
         DEN = B(K+1)-RT-C(K+1)*D(K+2)
         D(K+1) = A(K+1)/DEN
         W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
  107 CONTINUE
      DEN = B(1)-RT-C(1)*D(2)
      W(1) = 1.
      IF (DEN) 108,109,108
  108 W(1) = (Y(1)-C(1)*W(2))/DEN
  109 DO 110 J=2,M
         W(J) = W(J)-D(J)*W(J-1)
  110 CONTINUE
      IF (NA) 113,113,102
  111 DO 112 J=1,M
         Y(J) = W(J)
  112 CONTINUE
      IBR = 1
      GO TO 102
  113 IF (M1) 114,114,115
  114 IF (M2) 111,111,120
  115 IF (M2) 117,117,116
  116 IF (ABS(BM1(M1))-ABS(BM2(M2))) 120,120,117
  117 IF (IBR) 118,118,119
  118 IF (ABS(BM1(M1)-BD(ID))-ABS(BM1(M1)-RT)) 111,119,119
  119 RT = RT-BM1(M1)
      M1 = M1-1
      GO TO 123
  120 IF (IBR) 121,121,122
  121 IF (ABS(BM2(M2)-BD(ID))-ABS(BM2(M2)-RT)) 111,122,122
  122 RT = RT-BM2(M2)
      M2 = M2-1
  123 DO 124 J=1,M
         Y(J) = Y(J)+RT*W(J)
  124 CONTINUE
      GO TO 102
  125 RETURN
      END
