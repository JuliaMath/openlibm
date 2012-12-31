*DECK POS3D1
      SUBROUTINE POS3D1 (LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, XRT,
     +   YRT, T, D, WX, WY, C1, C2, BB)
C***BEGIN PROLOGUE  POS3D1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to POIS3D
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (POS3D1-S)
C***AUTHOR  (UNKNOWN)
C***SEE ALSO  POIS3D
C***ROUTINES CALLED  COSQB, COSQF, COSQI, COST, COSTI, PIMACH, RFFTB,
C                    RFFTF, RFFTI, SINQB, SINQF, SINQI, SINT, SINTI,
C                    TRIDQ
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891009  Removed unreferenced variable.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900308  Changed call to TRID to call to TRIDQ.  (WRB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  POS3D1
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,
     1                F(LDIMF,MDIMF,*)       ,XRT(*)     ,YRT(*)     ,
     2                T(*)       ,D(*)       ,WX(*)      ,WY(*)      ,
     3                BB(*)
C***FIRST EXECUTABLE STATEMENT  POS3D1
      PI = PIMACH(DUM)
      LR = L
      MR = M
      NR = N
C
C     GENERATE TRANSFORM ROOTS
C
      LRDEL = ((LP-1)*(LP-3)*(LP-5))/3
      SCALX = LR+LRDEL
      DX = PI/(2.*SCALX)
      GO TO (108,103,101,102,101),LP
  101 DI = 0.5
      SCALX = 2.*SCALX
      GO TO 104
  102 DI = 1.0
      GO TO 104
  103 DI = 0.0
  104 DO 105 I=1,LR
         XRT(I) = -4.*C1*(SIN((I-DI)*DX))**2
  105 CONTINUE
      SCALX = 2.*SCALX
      GO TO (112,106,110,107,111),LP
  106 CALL SINTI (LR,WX)
      GO TO 112
  107 CALL COSTI (LR,WX)
      GO TO 112
  108 XRT(1) = 0.
      XRT(LR) = -4.*C1
      DO 109 I=3,LR,2
         XRT(I-1) = -4.*C1*(SIN((I-1)*DX))**2
         XRT(I) = XRT(I-1)
  109 CONTINUE
      CALL RFFTI (LR,WX)
      GO TO 112
  110 CALL SINQI (LR,WX)
      GO TO 112
  111 CALL COSQI (LR,WX)
  112 CONTINUE
      MRDEL = ((MP-1)*(MP-3)*(MP-5))/3
      SCALY = MR+MRDEL
      DY = PI/(2.*SCALY)
      GO TO (120,115,113,114,113),MP
  113 DJ = 0.5
      SCALY = 2.*SCALY
      GO TO 116
  114 DJ = 1.0
      GO TO 116
  115 DJ = 0.0
  116 DO 117 J=1,MR
         YRT(J) = -4.*C2*(SIN((J-DJ)*DY))**2
  117 CONTINUE
      SCALY = 2.*SCALY
      GO TO (124,118,122,119,123),MP
  118 CALL SINTI (MR,WY)
      GO TO 124
  119 CALL COSTI (MR,WY)
      GO TO 124
  120 YRT(1) = 0.
      YRT(MR) = -4.*C2
      DO 121 J=3,MR,2
         YRT(J-1) = -4.*C2*(SIN((J-1)*DY))**2
         YRT(J) = YRT(J-1)
  121 CONTINUE
      CALL RFFTI (MR,WY)
      GO TO 124
  122 CALL SINQI (MR,WY)
      GO TO 124
  123 CALL COSQI (MR,WY)
  124 CONTINUE
      IFWRD = 1
  125 CONTINUE
C
C     TRANSFORM X
C
      DO 141 J=1,MR
         DO 140 K=1,NR
            DO 126 I=1,LR
               T(I) = F(I,J,K)
  126       CONTINUE
            GO TO (127,130,131,134,135),LP
  127       GO TO (128,129),IFWRD
  128       CALL RFFTF (LR,T,WX)
            GO TO 138
  129       CALL RFFTB (LR,T,WX)
            GO TO 138
  130       CALL SINT (LR,T,WX)
            GO TO 138
  131       GO TO (132,133),IFWRD
  132       CALL SINQF (LR,T,WX)
            GO TO 138
  133       CALL SINQB (LR,T,WX)
            GO TO 138
  134       CALL COST (LR,T,WX)
            GO TO 138
  135       GO TO (136,137),IFWRD
  136       CALL COSQF (LR,T,WX)
            GO TO 138
  137       CALL COSQB (LR,T,WX)
  138       CONTINUE
            DO 139 I=1,LR
               F(I,J,K) = T(I)
  139       CONTINUE
  140    CONTINUE
  141 CONTINUE
      GO TO (142,164),IFWRD
C
C     TRANSFORM Y
C
  142 CONTINUE
      DO 158 I=1,LR
         DO 157 K=1,NR
            DO 143 J=1,MR
               T(J) = F(I,J,K)
  143       CONTINUE
            GO TO (144,147,148,151,152),MP
  144       GO TO (145,146),IFWRD
  145       CALL RFFTF (MR,T,WY)
            GO TO 155
  146       CALL RFFTB (MR,T,WY)
            GO TO 155
  147       CALL SINT (MR,T,WY)
            GO TO 155
  148       GO TO (149,150),IFWRD
  149       CALL SINQF (MR,T,WY)
            GO TO 155
  150       CALL SINQB (MR,T,WY)
            GO TO 155
  151       CALL COST (MR,T,WY)
            GO TO 155
  152       GO TO (153,154),IFWRD
  153       CALL COSQF (MR,T,WY)
            GO TO 155
  154       CALL COSQB (MR,T,WY)
  155       CONTINUE
            DO 156 J=1,MR
               F(I,J,K) = T(J)
  156       CONTINUE
  157    CONTINUE
  158 CONTINUE
      GO TO (159,125),IFWRD
  159 CONTINUE
C
C     SOLVE TRIDIAGONAL SYSTEMS IN Z
C
      DO 163 I=1,LR
         DO 162 J=1,MR
            DO 160 K=1,NR
               BB(K) = B(K)+XRT(I)+YRT(J)
               T(K) = F(I,J,K)
  160       CONTINUE
            CALL TRIDQ (NR,A,BB,C,T,D)
            DO 161 K=1,NR
               F(I,J,K) = T(K)
  161       CONTINUE
  162    CONTINUE
  163 CONTINUE
      IFWRD = 2
      GO TO 142
  164 CONTINUE
      DO 167 I=1,LR
         DO 166 J=1,MR
            DO 165 K=1,NR
               F(I,J,K) = F(I,J,K)/(SCALX*SCALY)
  165       CONTINUE
  166    CONTINUE
  167 CONTINUE
      RETURN
      END
