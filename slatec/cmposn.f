*DECK CMPOSN
      SUBROUTINE CMPOSN (M, N, ISTAG, MIXBND, A, BB, C, Q, IDIMQ, B, B2,
     +   B3, W, W2, W3, D, TCOS, P)
C***BEGIN PROLOGUE  CMPOSN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CMGNBN
C***LIBRARY   SLATEC
C***TYPE      COMPLEX (POISN2-S, CMPOSN-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve Poisson's equation with Neumann boundary
C     conditions.
C
C     ISTAG = 1 if the last diagonal block is A.
C     ISTAG = 2 if the last diagonal block is A-I.
C     MIXBND = 1 if have Neumann boundary conditions at both boundaries.
C     MIXBND = 2 if have Neumann boundary conditions at bottom and
C     Dirichlet condition at top.  (For this case, must have ISTAG = 1)
C
C***SEE ALSO  CMGNBN
C***ROUTINES CALLED  C1MERG, CMPCSG, CMPTR3, CMPTRX
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   920130  Modified to use merge routine C1MERG rather than deleted
C           routine CMPMRG.  (WRB)
C***END PROLOGUE  CMPOSN
C
      COMPLEX         A          ,BB         ,C          ,Q          ,
     1                B          ,B2         ,B3         ,W          ,
     2                W2         ,W3         ,D          ,TCOS       ,
     3                P          ,FI         ,T
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                K(4)       ,P(*)
      EQUIVALENCE     (K(1),K1)  ,(K(2),K2)  ,(K(3),K3)  ,(K(4),K4)
C***FIRST EXECUTABLE STATEMENT  CMPOSN
      FISTAG = 3-ISTAG
      FNUM = 1./ISTAG
      FDEN = 0.5*(ISTAG-1)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      GO TO (101,103),ISTAG
  101 CONTINUE
      DO 102 I=1,MR
         Q(I,N) = .5*Q(I,N)
  102 CONTINUE
      GO TO (103,104),MIXBND
  103 IF (N .LE. 3) GO TO 155
  104 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 .EQ. NR) NROD = 0
      GO TO (105,106),MIXBND
  105 JSTART = 1
      GO TO 107
  106 JSTART = JR
      NROD = 1-NROD
  107 CONTINUE
      JSTOP = NLAST-JR
      IF (NROD .EQ. 0) JSTOP = JSTOP-I2R
      CALL CMPCSG (I2R,1,0.5,0.0,TCOS)
      I2RBY2 = I2R/2
      IF (JSTOP .GE. JSTART) GO TO 108
      J = JR
      GO TO 116
  108 CONTINUE
C
C     REGULAR REDUCTION.
C
      DO 115 J=JSTART,JSTOP,JR
         JP1 = J+I2RBY2
         JP2 = J+I2R
         JP3 = JP2+I2RBY2
         JM1 = J-I2RBY2
         JM2 = J-I2R
         JM3 = JM2-I2RBY2
         IF (J .NE. 1) GO TO 109
         JM1 = JP1
         JM2 = JP2
         JM3 = JP3
  109    CONTINUE
         IF (I2R .NE. 1) GO TO 111
         IF (J .EQ. 1) JM2 = JP2
         DO 110 I=1,MR
            B(I) = 2.*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  110    CONTINUE
         GO TO 113
  111    CONTINUE
         DO 112 I=1,MR
            FI = Q(I,J)
            Q(I,J) = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = FI+Q(I,J)-Q(I,JM3)-Q(I,JP3)
  112    CONTINUE
  113    CONTINUE
         CALL CMPTRX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
         DO 114 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  114    CONTINUE
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
  115 CONTINUE
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
      J = JSTOP+JR
  116 NLAST = J
      JM1 = J-I2RBY2
      JM2 = J-I2R
      JM3 = JM2-I2RBY2
      IF (NROD .EQ. 0) GO TO 128
C
C     ODD NUMBER OF UNKNOWNS
C
      IF (I2R .NE. 1) GO TO 118
      DO 117 I=1,MR
         B(I) = FISTAG*Q(I,J)
         Q(I,J) = Q(I,JM2)
  117 CONTINUE
      GO TO 126
  118 DO 119 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  119 CONTINUE
      IF (NRODPR .NE. 0) GO TO 121
      DO 120 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)
  120 CONTINUE
      IP = IP-MR
      GO TO 123
  121 CONTINUE
      DO 122 I=1,MR
         Q(I,J) = Q(I,J)-Q(I,JM1)+Q(I,JM2)
  122 CONTINUE
  123 IF (LR .EQ. 0) GO TO 124
      CALL CMPCSG (LR,1,0.5,FDEN,TCOS(KR+1))
      GO TO 126
  124 CONTINUE
      DO 125 I=1,MR
         B(I) = FISTAG*B(I)
  125 CONTINUE
  126 CONTINUE
      CALL CMPCSG (KR,1,0.5,FDEN,TCOS)
      CALL CMPTRX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 127 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
  127 CONTINUE
      KR = KR+I2R
      GO TO 151
  128 CONTINUE
C
C     EVEN NUMBER OF UNKNOWNS
C
      JP1 = J+I2RBY2
      JP2 = J+I2R
      IF (I2R .NE. 1) GO TO 135
      DO 129 I=1,MR
         B(I) = Q(I,J)
  129 CONTINUE
      CALL CMPTRX (1,0,MR,A,BB,C,B,TCOS,D,W)
      IP = 0
      IPSTOR = MR
      GO TO (133,130),ISTAG
  130 DO 131 I=1,MR
         P(I) = B(I)
         B(I) = B(I)+Q(I,N)
  131 CONTINUE
      TCOS(1) = CMPLX(1.,0.)
      TCOS(2) = CMPLX(0.,0.)
      CALL CMPTRX (1,1,MR,A,BB,C,B,TCOS,D,W)
      DO 132 I=1,MR
         Q(I,J) = Q(I,JM2)+P(I)+B(I)
  132 CONTINUE
      GO TO 150
  133 CONTINUE
      DO 134 I=1,MR
         P(I) = B(I)
         Q(I,J) = Q(I,JM2)+2.*Q(I,JP2)+3.*B(I)
  134 CONTINUE
      GO TO 150
  135 CONTINUE
      DO 136 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  136 CONTINUE
      IF (NRODPR .NE. 0) GO TO 138
      DO 137 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  137 CONTINUE
      GO TO 140
  138 CONTINUE
      DO 139 I=1,MR
         B(I) = B(I)+Q(I,JP2)-Q(I,JP1)
  139 CONTINUE
  140 CONTINUE
      CALL CMPTRX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
      IP = IP+MR
      IPSTOR = MAX(IPSTOR,IP+MR)
      DO 141 I=1,MR
         II = IP+I
         P(II) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         B(I) = P(II)+Q(I,JP2)
  141 CONTINUE
      IF (LR .EQ. 0) GO TO 142
      CALL CMPCSG (LR,1,0.5,FDEN,TCOS(I2R+1))
      CALL C1MERG (TCOS,0,I2R,I2R,LR,KR)
      GO TO 144
  142 DO 143 I=1,I2R
         II = KR+I
         TCOS(II) = TCOS(I)
  143 CONTINUE
  144 CALL CMPCSG (KR,1,0.5,FDEN,TCOS)
      IF (LR .NE. 0) GO TO 145
      GO TO (146,145),ISTAG
  145 CONTINUE
      CALL CMPTRX (KR,KR,MR,A,BB,C,B,TCOS,D,W)
      GO TO 148
  146 CONTINUE
      DO 147 I=1,MR
         B(I) = FISTAG*B(I)
  147 CONTINUE
  148 CONTINUE
      DO 149 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)+B(I)
  149 CONTINUE
  150 CONTINUE
      LR = KR
      KR = KR+JR
  151 CONTINUE
      GO TO (152,153),MIXBND
  152 NR = (NLAST-1)/JR+1
      IF (NR .LE. 3) GO TO 155
      GO TO 154
  153 NR = NLAST/JR
      IF (NR .LE. 1) GO TO 192
  154 I2R = JR
      NRODPR = NROD
      GO TO 104
  155 CONTINUE
C
C      BEGIN SOLUTION
C
      J = 1+JR
      JM1 = J-I2R
      JP1 = J+I2R
      JM2 = NLAST-I2R
      IF (NR .EQ. 2) GO TO 184
      IF (LR .NE. 0) GO TO 170
      IF (N .NE. 3) GO TO 161
C
C     CASE N = 3.
C
      GO TO (156,168),ISTAG
  156 CONTINUE
      DO 157 I=1,MR
         B(I) = Q(I,2)
  157 CONTINUE
      TCOS(1) = CMPLX(0.,0.)
      CALL CMPTRX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 158 I=1,MR
         Q(I,2) = B(I)
         B(I) = 4.*B(I)+Q(I,1)+2.*Q(I,3)
  158 CONTINUE
      TCOS(1) = CMPLX(-2.,0.)
      TCOS(2) = CMPLX(2.,0.)
      I1 = 2
      I2 = 0
      CALL CMPTRX (I1,I2,MR,A,BB,C,B,TCOS,D,W)
      DO 159 I=1,MR
         Q(I,2) = Q(I,2)+B(I)
         B(I) = Q(I,1)+2.*Q(I,2)
  159 CONTINUE
      TCOS(1) = (0.,0.)
      CALL CMPTRX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 160 I=1,MR
         Q(I,1) = B(I)
  160 CONTINUE
      JR = 1
      I2R = 0
      GO TO 194
C
C     CASE N = 2**P+1
C
  161 CONTINUE
      GO TO (162,170),ISTAG
  162 CONTINUE
      DO 163 I=1,MR
         B(I) = Q(I,J)+.5*Q(I,1)-Q(I,JM1)+Q(I,NLAST)-Q(I,JM2)
  163 CONTINUE
      CALL CMPCSG (JR,1,0.5,0.0,TCOS)
      CALL CMPTRX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 164 I=1,MR
         Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
         B(I) = Q(I,1)+2.*Q(I,NLAST)+4.*Q(I,J)
  164 CONTINUE
      JR2 = 2*JR
      CALL CMPCSG (JR,1,0.0,0.0,TCOS)
      DO 165 I=1,JR
         I1 = JR+I
         I2 = JR+1-I
         TCOS(I1) = -TCOS(I2)
  165 CONTINUE
      CALL CMPTRX (JR2,0,MR,A,BB,C,B,TCOS,D,W)
      DO 166 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+2.*Q(I,J)
  166 CONTINUE
      CALL CMPCSG (JR,1,0.5,0.0,TCOS)
      CALL CMPTRX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 167 I=1,MR
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)+B(I)
  167 CONTINUE
      GO TO 194
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  168 DO 169 I=1,MR
         B(I) = Q(I,2)
         Q(I,2) = (0.,0.)
         B2(I) = Q(I,3)
         B3(I) = Q(I,1)
  169 CONTINUE
      JR = 1
      I2R = 0
      J = 2
      GO TO 177
  170 CONTINUE
      DO 171 I=1,MR
         B(I) = .5*Q(I,1)-Q(I,JM1)+Q(I,J)
  171 CONTINUE
      IF (NROD .NE. 0) GO TO 173
      DO 172 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  172 CONTINUE
      GO TO 175
  173 DO 174 I=1,MR
         B(I) = B(I)+Q(I,NLAST)-Q(I,JM2)
  174 CONTINUE
  175 CONTINUE
      DO 176 I=1,MR
         T = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         Q(I,J) = T
         B2(I) = Q(I,NLAST)+T
         B3(I) = Q(I,1)+2.*T
  176 CONTINUE
  177 CONTINUE
      K1 = KR+2*JR-1
      K2 = KR+JR
      TCOS(K1+1) = (-2.,0.)
      K4 = K1+3-ISTAG
      CALL CMPCSG (K2+ISTAG-2,1,0.0,FNUM,TCOS(K4))
      K4 = K1+K2+1
      CALL CMPCSG (JR-1,1,0.0,1.0,TCOS(K4))
      CALL C1MERG (TCOS,K1,K2,K1+K2,JR-1,0)
      K3 = K1+K2+LR
      CALL CMPCSG (JR,1,0.5,0.0,TCOS(K3+1))
      K4 = K3+JR+1
      CALL CMPCSG (KR,1,0.5,FDEN,TCOS(K4))
      CALL C1MERG (TCOS,K3,JR,K3+JR,KR,K1)
      IF (LR .EQ. 0) GO TO 178
      CALL CMPCSG (LR,1,0.5,FDEN,TCOS(K4))
      CALL C1MERG (TCOS,K3,JR,K3+JR,LR,K3-LR)
      CALL CMPCSG (KR,1,0.5,FDEN,TCOS(K4))
  178 K3 = KR
      K4 = KR
      CALL CMPTR3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 179 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  179 CONTINUE
      TCOS(1) = (2.,0.)
      CALL CMPTRX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 180 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+2.*Q(I,J)
  180 CONTINUE
      CALL CMPCSG (JR,1,0.5,0.0,TCOS)
      CALL CMPTRX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      IF (JR .NE. 1) GO TO 182
      DO 181 I=1,MR
         Q(I,1) = B(I)
  181 CONTINUE
      GO TO 194
  182 CONTINUE
      DO 183 I=1,MR
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)+B(I)
  183 CONTINUE
      GO TO 194
  184 CONTINUE
      IF (N .NE. 2) GO TO 188
C
C     CASE  N = 2
C
      DO 185 I=1,MR
         B(I) = Q(I,1)
  185 CONTINUE
      TCOS(1) = (0.,0.)
      CALL CMPTRX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 186 I=1,MR
         Q(I,1) = B(I)
         B(I) = 2.*(Q(I,2)+B(I))*FISTAG
  186 CONTINUE
      TCOS(1) = CMPLX(-FISTAG,0.)
      TCOS(2) = CMPLX(2.,0.)
      CALL CMPTRX (2,0,MR,A,BB,C,B,TCOS,D,W)
      DO 187 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  187 CONTINUE
      JR = 1
      I2R = 0
      GO TO 194
  188 CONTINUE
C
C     CASE OF GENERAL N AND NR = 2 .
C
      DO 189 I=1,MR
         II = IP+I
         B3(I) = (0.,0.)
         B(I) = Q(I,1)+2.*P(II)
         Q(I,1) = .5*Q(I,1)-Q(I,JM1)
         B2(I) = 2.*(Q(I,1)+Q(I,NLAST))
  189 CONTINUE
      K1 = KR+JR-1
      TCOS(K1+1) = (-2.,0.)
      K4 = K1+3-ISTAG
      CALL CMPCSG (KR+ISTAG-2,1,0.0,FNUM,TCOS(K4))
      K4 = K1+KR+1
      CALL CMPCSG (JR-1,1,0.0,1.0,TCOS(K4))
      CALL C1MERG (TCOS,K1,KR,K1+KR,JR-1,0)
      CALL CMPCSG (KR,1,0.5,FDEN,TCOS(K1+1))
      K2 = KR
      K4 = K1+K2+1
      CALL CMPCSG (LR,1,0.5,FDEN,TCOS(K4))
      K3 = LR
      K4 = 0
      CALL CMPTR3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 190 I=1,MR
         B(I) = B(I)+B2(I)
  190 CONTINUE
      TCOS(1) = (2.,0.)
      CALL CMPTRX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 191 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  191 CONTINUE
      GO TO 194
  192 DO 193 I=1,MR
         B(I) = Q(I,NLAST)
  193 CONTINUE
      GO TO 196
  194 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      J = NLAST-JR
      DO 195 I=1,MR
         B(I) = Q(I,NLAST)+Q(I,J)
  195 CONTINUE
  196 JM2 = NLAST-I2R
      IF (JR .NE. 1) GO TO 198
      DO 197 I=1,MR
         Q(I,NLAST) = (0.,0.)
  197 CONTINUE
      GO TO 202
  198 CONTINUE
      IF (NROD .NE. 0) GO TO 200
      DO 199 I=1,MR
         II = IP+I
         Q(I,NLAST) = P(II)
  199 CONTINUE
      IP = IP-MR
      GO TO 202
  200 DO 201 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)-Q(I,JM2)
  201 CONTINUE
  202 CONTINUE
      CALL CMPCSG (KR,1,0.5,FDEN,TCOS)
      CALL CMPCSG (LR,1,0.5,FDEN,TCOS(KR+1))
      IF (LR .NE. 0) GO TO 204
      DO 203 I=1,MR
         B(I) = FISTAG*B(I)
  203 CONTINUE
  204 CONTINUE
      CALL CMPTRX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 205 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)+B(I)
  205 CONTINUE
      NLASTP = NLAST
  206 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR .EQ. 0) GO TO 222
      GO TO (207,208),MIXBND
  207 JSTART = 1+JR
      GO TO 209
  208 JSTART = JR
  209 CONTINUE
      KR = KR-JR
      IF (NLAST+JR .GT. N) GO TO 210
      KR = KR-JR
      NLAST = NLAST+JR
      JSTOP = NLAST-JSTEP
      GO TO 211
  210 CONTINUE
      JSTOP = NLAST-JR
  211 CONTINUE
      LR = KR-JR
      CALL CMPCSG (JR,1,0.5,0.0,TCOS)
      DO 221 J=JSTART,JSTOP,JSTEP
         JM2 = J-JR
         JP2 = J+JR
         IF (J .NE. JR) GO TO 213
         DO 212 I=1,MR
            B(I) = Q(I,J)+Q(I,JP2)
  212    CONTINUE
         GO TO 215
  213    CONTINUE
         DO 214 I=1,MR
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  214    CONTINUE
  215    CONTINUE
         IF (JR .NE. 1) GO TO 217
         DO 216 I=1,MR
            Q(I,J) = (0.,0.)
  216    CONTINUE
         GO TO 219
  217    CONTINUE
         JM1 = J-I2R
         JP1 = J+I2R
         DO 218 I=1,MR
            Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  218    CONTINUE
  219    CONTINUE
         CALL CMPTRX (JR,0,MR,A,BB,C,B,TCOS,D,W)
         DO 220 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  220    CONTINUE
  221 CONTINUE
      NROD = 1
      IF (NLAST+I2R .LE. N) NROD = 0
      IF (NLASTP .NE. NLAST) GO TO 194
      GO TO 206
  222 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = CMPLX(REAL(IPSTOR),0.)
      RETURN
      END
