*DECK POSTG2
      SUBROUTINE POSTG2 (NPEROD, N, M, A, BB, C, IDIMQ, Q, B, B2, B3, W,
     +   W2, W3, D, TCOS, P)
C***BEGIN PROLOGUE  POSTG2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to POISTG
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (POSTG2-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve Poisson's equation on a staggered grid.
C
C***SEE ALSO  POISTG
C***ROUTINES CALLED  COSGEN, S1MERG, TRI3, TRIX
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   920130  Modified to use merge routine S1MERG rather than deleted
C           routine MERGE.  (WRB)
C***END PROLOGUE  POSTG2
C
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                K(4)       ,P(*)
      EQUIVALENCE     (K(1),K1)  ,(K(2),K2)  ,(K(3),K3)  ,(K(4),K4)
C***FIRST EXECUTABLE STATEMENT  POSTG2
      NP = NPEROD
      FNUM = 0.5*(NP/3)
      FNUM2 = 0.5*(NP/2)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      IF (NR .LE. 3) GO TO 142
  101 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 .EQ. NR) NROD = 0
      JSTART = 1
      JSTOP = NLAST-JR
      IF (NROD .EQ. 0) JSTOP = JSTOP-I2R
      I2RBY2 = I2R/2
      IF (JSTOP .GE. JSTART) GO TO 102
      J = JR
      GO TO 115
  102 CONTINUE
C
C     REGULAR REDUCTION.
C
      IJUMP = 1
      DO 114 J=JSTART,JSTOP,JR
         JP1 = J+I2RBY2
         JP2 = J+I2R
         JP3 = JP2+I2RBY2
         JM1 = J-I2RBY2
         JM2 = J-I2R
         JM3 = JM2-I2RBY2
         IF (J .NE. 1) GO TO 106
         CALL COSGEN (I2R,1,FNUM,0.5,TCOS)
         IF (I2R .NE. 1) GO TO 104
         DO 103 I=1,MR
            B(I) = Q(I,1)
            Q(I,1) = Q(I,2)
  103    CONTINUE
         GO TO 112
  104    DO 105 I=1,MR
            B(I) = Q(I,1)+0.5*(Q(I,JP2)-Q(I,JP1)-Q(I,JP3))
            Q(I,1) = Q(I,JP2)+Q(I,1)-Q(I,JP1)
  105    CONTINUE
         GO TO 112
  106    CONTINUE
         GO TO (107,108),IJUMP
  107    CONTINUE
         IJUMP = 2
         CALL COSGEN (I2R,1,0.5,0.0,TCOS)
  108    CONTINUE
         IF (I2R .NE. 1) GO TO 110
         DO 109 I=1,MR
            B(I) = 2.*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  109    CONTINUE
         GO TO 112
  110    DO 111 I=1,MR
            FI = Q(I,J)
            Q(I,J) = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = FI+Q(I,J)-Q(I,JM3)-Q(I,JP3)
  111    CONTINUE
  112    CONTINUE
         CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
         DO 113 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  113    CONTINUE
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
  114 CONTINUE
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
      J = JSTOP+JR
  115 NLAST = J
      JM1 = J-I2RBY2
      JM2 = J-I2R
      JM3 = JM2-I2RBY2
      IF (NROD .EQ. 0) GO TO 125
C
C     ODD NUMBER OF UNKNOWNS
C
      IF (I2R .NE. 1) GO TO 117
      DO 116 I=1,MR
         B(I) = Q(I,J)
         Q(I,J) = Q(I,JM2)
  116 CONTINUE
      GO TO 123
  117 DO 118 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  118 CONTINUE
      IF (NRODPR .NE. 0) GO TO 120
      DO 119 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)
  119 CONTINUE
      IP = IP-MR
      GO TO 122
  120 CONTINUE
      DO 121 I=1,MR
         Q(I,J) = Q(I,J)-Q(I,JM1)+Q(I,JM2)
  121 CONTINUE
  122 IF (LR .EQ. 0) GO TO 123
      CALL COSGEN (LR,1,FNUM2,0.5,TCOS(KR+1))
  123 CONTINUE
      CALL COSGEN (KR,1,FNUM2,0.5,TCOS)
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 124 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
  124 CONTINUE
      KR = KR+I2R
      GO TO 141
  125 CONTINUE
C
C     EVEN NUMBER OF UNKNOWNS
C
      JP1 = J+I2RBY2
      JP2 = J+I2R
      IF (I2R .NE. 1) GO TO 129
      DO 126 I=1,MR
         B(I) = Q(I,J)
  126 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      IP = 0
      IPSTOR = MR
      DO 127 I=1,MR
         P(I) = B(I)
         B(I) = B(I)+Q(I,N)
  127 CONTINUE
      TCOS(1) = -1.+2*(NP/2)
      TCOS(2) = 0.
      CALL TRIX (1,1,MR,A,BB,C,B,TCOS,D,W)
      DO 128 I=1,MR
         Q(I,J) = Q(I,JM2)+P(I)+B(I)
  128 CONTINUE
      GO TO 140
  129 CONTINUE
      DO 130 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  130 CONTINUE
      IF (NRODPR .NE. 0) GO TO 132
      DO 131 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  131 CONTINUE
      GO TO 134
  132 CONTINUE
      DO 133 I=1,MR
         B(I) = B(I)+Q(I,JP2)-Q(I,JP1)
  133 CONTINUE
  134 CONTINUE
      CALL COSGEN (I2R,1,0.5,0.0,TCOS)
      CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
      IP = IP+MR
      IPSTOR = MAX(IPSTOR,IP+MR)
      DO 135 I=1,MR
         II = IP+I
         P(II) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         B(I) = P(II)+Q(I,JP2)
  135 CONTINUE
      IF (LR .EQ. 0) GO TO 136
      CALL COSGEN (LR,1,FNUM2,0.5,TCOS(I2R+1))
      CALL S1MERG (TCOS,0,I2R,I2R,LR,KR)
      GO TO 138
  136 DO 137 I=1,I2R
         II = KR+I
         TCOS(II) = TCOS(I)
  137 CONTINUE
  138 CALL COSGEN (KR,1,FNUM2,0.5,TCOS)
      CALL TRIX (KR,KR,MR,A,BB,C,B,TCOS,D,W)
      DO 139 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)+B(I)
  139 CONTINUE
  140 CONTINUE
      LR = KR
      KR = KR+JR
  141 CONTINUE
      NR = (NLAST-1)/JR+1
      IF (NR .LE. 3) GO TO 142
      I2R = JR
      NRODPR = NROD
      GO TO 101
  142 CONTINUE
C
C      BEGIN SOLUTION
C
      J = 1+JR
      JM1 = J-I2R
      JP1 = J+I2R
      JM2 = NLAST-I2R
      IF (NR .EQ. 2) GO TO 180
      IF (LR .NE. 0) GO TO 167
      IF (N .NE. 3) GO TO 156
C
C     CASE N = 3.
C
      GO TO (143,148,143),NP
  143 DO 144 I=1,MR
         B(I) = Q(I,2)
         B2(I) = Q(I,1)+Q(I,3)
         B3(I) = 0.
  144 CONTINUE
      GO TO (146,146,145),NP
  145 TCOS(1) = -1.
      TCOS(2) = 1.
      K1 = 1
      GO TO 147
  146 TCOS(1) = -2.
      TCOS(2) = 1.
      TCOS(3) = -1.
      K1 = 2
  147 K2 = 1
      K3 = 0
      K4 = 0
      GO TO 150
  148 DO 149 I=1,MR
         B(I) = Q(I,2)
         B2(I) = Q(I,3)
         B3(I) = Q(I,1)
  149 CONTINUE
      CALL COSGEN (3,1,0.5,0.0,TCOS)
      TCOS(4) = -1.
      TCOS(5) = 1.
      TCOS(6) = -1.
      TCOS(7) = 1.
      K1 = 3
      K2 = 2
      K3 = 1
      K4 = 1
  150 CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 151 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  151 CONTINUE
      GO TO (153,153,152),NP
  152 TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  153 DO 154 I=1,MR
         Q(I,2) = B(I)
         B(I) = Q(I,1)+B(I)
  154 CONTINUE
      TCOS(1) = -1.+4.*FNUM
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 155 I=1,MR
         Q(I,1) = B(I)
  155 CONTINUE
      JR = 1
      I2R = 0
      GO TO 188
C
C     CASE N = 2**P+1
C
  156 CONTINUE
      DO 157 I=1,MR
         B(I) = Q(I,J)+Q(I,1)-Q(I,JM1)+Q(I,NLAST)-Q(I,JM2)
  157 CONTINUE
      GO TO (158,160,158),NP
  158 DO 159 I=1,MR
         B2(I) = Q(I,1)+Q(I,NLAST)+Q(I,J)-Q(I,JM1)-Q(I,JP1)
         B3(I) = 0.
  159 CONTINUE
      K1 = NLAST-1
      K2 = NLAST+JR-1
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(NLAST))
      TCOS(K2) = 2*NP-4
      CALL COSGEN (JR,1,0.5-FNUM,0.5,TCOS(K2+1))
      K3 = (3-NP)/2
      CALL S1MERG (TCOS,K1,JR-K3,K2-K3,JR+K3,0)
      K1 = K1-1+K3
      CALL COSGEN (JR,1,FNUM,0.5,TCOS(K1+1))
      K2 = JR
      K3 = 0
      K4 = 0
      GO TO 162
  160 DO 161 I=1,MR
         FI = (Q(I,J)-Q(I,JM1)-Q(I,JP1))/2.
         B2(I) = Q(I,1)+FI
         B3(I) = Q(I,NLAST)+FI
  161 CONTINUE
      K1 = NLAST+JR-1
      K2 = K1+JR-1
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K1+1))
      CALL COSGEN (NLAST,1,0.5,0.0,TCOS(K2+1))
      CALL S1MERG (TCOS,K1,JR-1,K2,NLAST,0)
      K3 = K1+NLAST-1
      K4 = K3+JR
      CALL COSGEN (JR,1,0.5,0.5,TCOS(K3+1))
      CALL COSGEN (JR,1,0.0,0.5,TCOS(K4+1))
      CALL S1MERG (TCOS,K3,JR,K4,JR,K1)
      K2 = NLAST-1
      K3 = JR
      K4 = JR
  162 CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 163 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  163 CONTINUE
      IF (NP .NE. 3) GO TO 164
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  164 DO 165 I=1,MR
         Q(I,J) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         B(I) = Q(I,J)+Q(I,1)
  165 CONTINUE
      CALL COSGEN (JR,1,FNUM,0.5,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 166 I=1,MR
         Q(I,1) = Q(I,1)-Q(I,JM1)+B(I)
  166 CONTINUE
      GO TO 188
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  167 CONTINUE
      DO 168 I=1,MR
         B(I) = Q(I,1)-Q(I,JM1)+Q(I,J)
  168 CONTINUE
      IF (NROD .NE. 0) GO TO 170
      DO 169 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  169 CONTINUE
      GO TO 172
  170 DO 171 I=1,MR
         B(I) = B(I)+Q(I,NLAST)-Q(I,JM2)
  171 CONTINUE
  172 CONTINUE
      DO 173 I=1,MR
         T = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         Q(I,J) = T
         B2(I) = Q(I,NLAST)+T
         B3(I) = Q(I,1)+T
  173 CONTINUE
      K1 = KR+2*JR
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K1+1))
      K2 = K1+JR
      TCOS(K2) = 2*NP-4
      K4 = (NP-1)*(3-NP)
      K3 = K2+1-K4
      CALL COSGEN (KR+JR+K4,1,K4/2.,1.-K4,TCOS(K3))
      K4 = 1-NP/3
      CALL S1MERG (TCOS,K1,JR-K4,K2-K4,KR+JR+K4,0)
      IF (NP .EQ. 3) K1 = K1-1
      K2 = KR+JR
      K4 = K1+K2
      CALL COSGEN (KR,1,FNUM2,0.5,TCOS(K4+1))
      K3 = K4+KR
      CALL COSGEN (JR,1,FNUM,0.5,TCOS(K3+1))
      CALL S1MERG (TCOS,K4,KR,K3,JR,K1)
      K4 = K3+JR
      CALL COSGEN (LR,1,FNUM2,0.5,TCOS(K4+1))
      CALL S1MERG (TCOS,K3,JR,K4,LR,K1+K2)
      CALL COSGEN (KR,1,FNUM2,0.5,TCOS(K3+1))
      K3 = KR
      K4 = KR
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 174 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  174 CONTINUE
      IF (NP .NE. 3) GO TO 175
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  175 DO 176 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+Q(I,J)
  176 CONTINUE
      CALL COSGEN (JR,1,FNUM,0.5,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      IF (JR .NE. 1) GO TO 178
      DO 177 I=1,MR
         Q(I,1) = B(I)
  177 CONTINUE
      GO TO 188
  178 CONTINUE
      DO 179 I=1,MR
         Q(I,1) = Q(I,1)-Q(I,JM1)+B(I)
  179 CONTINUE
      GO TO 188
  180 CONTINUE
C
C     CASE OF GENERAL N AND NR = 2 .
C
      DO 181 I=1,MR
         II = IP+I
         B3(I) = 0.
         B(I) = Q(I,1)+P(II)
         Q(I,1) = Q(I,1)-Q(I,JM1)
         B2(I) = Q(I,1)+Q(I,NLAST)
  181 CONTINUE
      K1 = KR+JR
      K2 = K1+JR
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K1+1))
      GO TO (182,183,182),NP
  182 TCOS(K2) = 2*NP-4
      CALL COSGEN (KR,1,0.0,1.0,TCOS(K2+1))
      GO TO 184
  183 CALL COSGEN (KR+1,1,0.5,0.0,TCOS(K2))
  184 K4 = 1-NP/3
      CALL S1MERG (TCOS,K1,JR-K4,K2-K4,KR+K4,0)
      IF (NP .EQ. 3) K1 = K1-1
      K2 = KR
      CALL COSGEN (KR,1,FNUM2,0.5,TCOS(K1+1))
      K4 = K1+KR
      CALL COSGEN (LR,1,FNUM2,0.5,TCOS(K4+1))
      K3 = LR
      K4 = 0
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 185 I=1,MR
         B(I) = B(I)+B2(I)
  185 CONTINUE
      IF (NP .NE. 3) GO TO 186
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  186 DO 187 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  187 CONTINUE
  188 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      J = NLAST-JR
      DO 189 I=1,MR
         B(I) = Q(I,NLAST)+Q(I,J)
  189 CONTINUE
      JM2 = NLAST-I2R
      IF (JR .NE. 1) GO TO 191
      DO 190 I=1,MR
         Q(I,NLAST) = 0.
  190 CONTINUE
      GO TO 195
  191 CONTINUE
      IF (NROD .NE. 0) GO TO 193
      DO 192 I=1,MR
         II = IP+I
         Q(I,NLAST) = P(II)
  192 CONTINUE
      IP = IP-MR
      GO TO 195
  193 DO 194 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)-Q(I,JM2)
  194 CONTINUE
  195 CONTINUE
      CALL COSGEN (KR,1,FNUM2,0.5,TCOS)
      CALL COSGEN (LR,1,FNUM2,0.5,TCOS(KR+1))
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 196 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)+B(I)
  196 CONTINUE
      NLASTP = NLAST
  197 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR .EQ. 0) GO TO 210
      JSTART = 1+JR
      KR = KR-JR
      IF (NLAST+JR .GT. N) GO TO 198
      KR = KR-JR
      NLAST = NLAST+JR
      JSTOP = NLAST-JSTEP
      GO TO 199
  198 CONTINUE
      JSTOP = NLAST-JR
  199 CONTINUE
      LR = KR-JR
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      DO 209 J=JSTART,JSTOP,JSTEP
         JM2 = J-JR
         JP2 = J+JR
         IF (J .NE. JR) GO TO 201
         DO 200 I=1,MR
            B(I) = Q(I,J)+Q(I,JP2)
  200    CONTINUE
         GO TO 203
  201    CONTINUE
         DO 202 I=1,MR
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  202    CONTINUE
  203    CONTINUE
         IF (JR .NE. 1) GO TO 205
         DO 204 I=1,MR
            Q(I,J) = 0.
  204    CONTINUE
         GO TO 207
  205    CONTINUE
         JM1 = J-I2R
         JP1 = J+I2R
         DO 206 I=1,MR
            Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  206    CONTINUE
  207    CONTINUE
         CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
         DO 208 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  208    CONTINUE
  209 CONTINUE
      NROD = 1
      IF (NLAST+I2R .LE. N) NROD = 0
      IF (NLASTP .NE. NLAST) GO TO 188
      GO TO 197
  210 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
