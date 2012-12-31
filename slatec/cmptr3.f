*DECK CMPTR3
      SUBROUTINE CMPTR3 (M, A, B, C, K, Y1, Y2, Y3, TCOS, D, W1, W2, W3)
C***BEGIN PROLOGUE  CMPTR3
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CMGNBN
C***LIBRARY   SLATEC
C***TYPE      COMPLEX (TRI3-S, CMPTR3-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve tridiagonal systems.
C
C***SEE ALSO  CMGNBN
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890206  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CMPTR3
      COMPLEX         A          ,B          ,C          ,Y1         ,
     1                Y2         ,Y3         ,TCOS       ,D          ,
     2                W1         ,W2         ,W3         ,X          ,
     3                XX         ,Z
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,K(4)       ,
     1                TCOS(*)    ,Y1(*)      ,Y2(*)      ,Y3(*)      ,
     2                D(*)       ,W1(*)      ,W2(*)      ,W3(*)
      INTEGER K1P1, K2P1, K3P1, K4P1
C
C***FIRST EXECUTABLE STATEMENT  CMPTR3
      MM1 = M-1
      K1 = K(1)
      K2 = K(2)
      K3 = K(3)
      K4 = K(4)
      K1P1 = K1+1
      K2P1 = K2+1
      K3P1 = K3+1
      K4P1 = K4+1
      K2K3K4 = K2+K3+K4
      IF (K2K3K4 .EQ. 0) GO TO 101
      L1 = K1P1/K2P1
      L2 = K1P1/K3P1
      L3 = K1P1/K4P1
      LINT1 = 1
      LINT2 = 1
      LINT3 = 1
      KINT1 = K1
      KINT2 = KINT1+K2
      KINT3 = KINT2+K3
  101 CONTINUE
      DO 115 N=1,K1
         X = TCOS(N)
         IF (K2K3K4 .EQ. 0) GO TO 107
         IF (N .NE. L1) GO TO 103
         DO 102 I=1,M
            W1(I) = Y1(I)
  102    CONTINUE
  103    IF (N .NE. L2) GO TO 105
         DO 104 I=1,M
            W2(I) = Y2(I)
  104    CONTINUE
  105    IF (N .NE. L3) GO TO 107
         DO 106 I=1,M
            W3(I) = Y3(I)
  106    CONTINUE
  107    CONTINUE
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y1(1) = Y1(1)*Z
         Y2(1) = Y2(1)*Z
         Y3(1) = Y3(1)*Z
         DO 108 I=2,M
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y1(I) = (Y1(I)-A(I)*Y1(I-1))*Z
            Y2(I) = (Y2(I)-A(I)*Y2(I-1))*Z
            Y3(I) = (Y3(I)-A(I)*Y3(I-1))*Z
  108    CONTINUE
         DO 109 IP=1,MM1
            I = M-IP
            Y1(I) = Y1(I)-D(I)*Y1(I+1)
            Y2(I) = Y2(I)-D(I)*Y2(I+1)
            Y3(I) = Y3(I)-D(I)*Y3(I+1)
  109    CONTINUE
         IF (K2K3K4 .EQ. 0) GO TO 115
         IF (N .NE. L1) GO TO 111
         I = LINT1+KINT1
         XX = X-TCOS(I)
         DO 110 I=1,M
            Y1(I) = XX*Y1(I)+W1(I)
  110    CONTINUE
         LINT1 = LINT1+1
         L1 = (LINT1*K1P1)/K2P1
  111    IF (N .NE. L2) GO TO 113
         I = LINT2+KINT2
         XX = X-TCOS(I)
         DO 112 I=1,M
            Y2(I) = XX*Y2(I)+W2(I)
  112    CONTINUE
         LINT2 = LINT2+1
         L2 = (LINT2*K1P1)/K3P1
  113    IF (N .NE. L3) GO TO 115
         I = LINT3+KINT3
         XX = X-TCOS(I)
         DO 114 I=1,M
            Y3(I) = XX*Y3(I)+W3(I)
  114    CONTINUE
         LINT3 = LINT3+1
         L3 = (LINT3*K1P1)/K4P1
  115 CONTINUE
      RETURN
      END
