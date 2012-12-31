*DECK CBLKT1
      SUBROUTINE CBLKT1 (N, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, B, W1,
     +   W2, W3, WD, WW, WU, PRDCT, CPRDCT)
C***BEGIN PROLOGUE  CBLKT1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CBLKTR
C***LIBRARY   SLATEC
C***TYPE      COMPLEX (BLKTR1-S, CBLKT1-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C CBLKT1 solves the linear system of routine CBLKTR.
C
C B  contains the roots of all the B polynomials.
C W1,W2,W3,WD,WW,WU  are all working arrays.
C PRDCT is either PROCP or PROC depending on whether the boundary
C conditions in the M direction are periodic or not.
C CPRDCT is either CPROCP or CPROC which are called if some of the zeros
C of the B polynomials are complex.
C
C***SEE ALSO  CBLKTR
C***ROUTINES CALLED  INXCA, INXCB, INXCC
C***COMMON BLOCKS    CCBLK
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CBLKT1
C
      DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,AM(*)      ,
     1                BM(*)      ,CM(*)      ,B(*)       ,W1(*)      ,
     2                W2(*)      ,W3(*)      ,WD(*)      ,WW(*)      ,
     3                WU(*)      ,Y(IDIMY,*)
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
      COMPLEX         AM         ,BM         ,CM         ,Y          ,
     1                W1         ,W2         ,W3         ,WD         ,
     2                WW         ,WU
C***FIRST EXECUTABLE STATEMENT  CBLKT1
      KDO = K-1
      DO 109 L=1,KDO
         IR = L-1
         I2 = 2**IR
         I1 = I2/2
         I3 = I2+I1
         I4 = I2+I2
         IRM1 = IR-1
         CALL INXCB (I2,IR,IM2,NM2)
         CALL INXCB (I1,IRM1,IM3,NM3)
         CALL INXCB (I3,IRM1,IM1,NM1)
         CALL PRDCT (NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,Y(1,I2),W3,
     1               M,AM,BM,CM,WD,WW,WU)
         IF = 2**K
         DO 108 I=I4,IF,I4
            IF (I-NM) 101,101,108
  101       IPI1 = I+I1
            IPI2 = I+I2
            IPI3 = I+I3
            CALL INXCC (I,IR,IDXC,NC)
            IF (I-IF) 102,108,108
  102       CALL INXCA (I,IR,IDXA,NA)
            CALL INXCB (I-I1,IRM1,IM1,NM1)
            CALL INXCB (IPI2,IR,IP2,NP2)
            CALL INXCB (IPI1,IRM1,IP1,NP1)
            CALL INXCB (IPI3,IRM1,IP3,NP3)
            CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W3,W1,M,AM,
     1                  BM,CM,WD,WW,WU)
            IF (IPI2-NM) 105,105,103
  103       DO 104 J=1,M
               W3(J) = (0.,0.)
               W2(J) = (0.,0.)
  104       CONTINUE
            GO TO 106
  105       CALL PRDCT (NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,
     1                  Y(1,IPI2),W3,M,AM,BM,CM,WD,WW,WU)
            CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W3,W2,M,AM,
     1                  BM,CM,WD,WW,WU)
  106       DO 107 J=1,M
               Y(J,I) = W1(J)+W2(J)+Y(J,I)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      IF (NPP) 132,110,132
C
C     THE PERIODIC CASE IS TREATED USING THE CAPACITANCE MATRIX METHOD
C
  110 IF = 2**K
      I = IF/2
      I1 = I/2
      CALL INXCB (I-I1,K-2,IM1,NM1)
      CALL INXCB (I+I1,K-2,IP1,NP1)
      CALL INXCB (I,K-1,IZ,NZ)
      CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,Y(1,I),W1,M,AM,
     1            BM,CM,WD,WW,WU)
      IZR = I
      DO 111 J=1,M
         W2(J) = W1(J)
  111 CONTINUE
      DO 113 LL=2,K
         L = K-LL+1
         IR = L-1
         I2 = 2**IR
         I1 = I2/2
         I = I2
         CALL INXCC (I,IR,IDXC,NC)
         CALL INXCB (I,IR,IZ,NZ)
         CALL INXCB (I-I1,IR-1,IM1,NM1)
         CALL INXCB (I+I1,IR-1,IP1,NP1)
         CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W1,W1,M,AM,BM,
     1               CM,WD,WW,WU)
         DO 112 J=1,M
            W1(J) = Y(J,I)+W1(J)
  112    CONTINUE
         CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,W1,M,AM,
     1               BM,CM,WD,WW,WU)
  113 CONTINUE
      DO 118 LL=2,K
         L = K-LL+1
         IR = L-1
         I2 = 2**IR
         I1 = I2/2
         I4 = I2+I2
         IFD = IF-I2
         DO 117 I=I2,IFD,I4
            IF (I-I2-IZR) 117,114,117
  114       IF (I-NM) 115,115,118
  115       CALL INXCA (I,IR,IDXA,NA)
            CALL INXCB (I,IR,IZ,NZ)
            CALL INXCB (I-I1,IR-1,IM1,NM1)
            CALL INXCB (I+I1,IR-1,IP1,NP1)
            CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W2,W2,M,AM,
     1                  BM,CM,WD,WW,WU)
            DO 116 J=1,M
               W2(J) = Y(J,I)+W2(J)
  116       CONTINUE
            CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W2,W2,M,
     1                  AM,BM,CM,WD,WW,WU)
            IZR = I
            IF (I-NM) 117,119,117
  117    CONTINUE
  118 CONTINUE
  119 DO 120 J=1,M
         Y(J,NM+1) = Y(J,NM+1)-CN(NM+1)*W1(J)-AN(NM+1)*W2(J)
  120 CONTINUE
      CALL INXCB (IF/2,K-1,IM1,NM1)
      CALL INXCB (IF,K-1,IP,NP)
      IF (NCMPLX) 121,122,121
  121 CALL CPRDCT (NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1),
     1             Y(1,NM+1),M,AM,BM,CM,W1,W3,WW)
      GO TO 123
  122 CALL PRDCT (NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1),
     1            Y(1,NM+1),M,AM,BM,CM,WD,WW,WU)
  123 DO 124 J=1,M
         W1(J) = AN(1)*Y(J,NM+1)
         W2(J) = CN(NM)*Y(J,NM+1)
         Y(J,1) = Y(J,1)-W1(J)
         Y(J,NM) = Y(J,NM)-W2(J)
  124 CONTINUE
      DO 126 L=1,KDO
         IR = L-1
         I2 = 2**IR
         I4 = I2+I2
         I1 = I2/2
         I = I4
         CALL INXCA (I,IR,IDXA,NA)
         CALL INXCB (I-I2,IR,IM2,NM2)
         CALL INXCB (I-I2-I1,IR-1,IM3,NM3)
         CALL INXCB (I-I1,IR-1,IM1,NM1)
         CALL PRDCT (NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,W1,W1,M,AM,
     1               BM,CM,WD,WW,WU)
         CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W1,W1,M,AM,BM,
     1               CM,WD,WW,WU)
         DO 125 J=1,M
            Y(J,I) = Y(J,I)-W1(J)
  125    CONTINUE
  126 CONTINUE
C
      IZR = NM
      DO 131 L=1,KDO
         IR = L-1
         I2 = 2**IR
         I1 = I2/2
         I3 = I2+I1
         I4 = I2+I2
         IRM1 = IR-1
         DO 130 I=I4,IF,I4
            IPI1 = I+I1
            IPI2 = I+I2
            IPI3 = I+I3
            IF (IPI2-IZR) 127,128,127
  127       IF (I-IZR) 130,131,130
  128       CALL INXCC (I,IR,IDXC,NC)
            CALL INXCB (IPI2,IR,IP2,NP2)
            CALL INXCB (IPI1,IRM1,IP1,NP1)
            CALL INXCB (IPI3,IRM1,IP3,NP3)
            CALL PRDCT (NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,W2,W2,M,
     1                  AM,BM,CM,WD,WW,WU)
            CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W2,W2,M,AM,
     1                  BM,CM,WD,WW,WU)
            DO 129 J=1,M
               Y(J,I) = Y(J,I)-W2(J)
  129       CONTINUE
            IZR = I
            GO TO 131
  130    CONTINUE
  131 CONTINUE
C
C BEGIN BACK SUBSTITUTION PHASE
C
  132 DO 144 LL=1,K
         L = K-LL+1
         IR = L-1
         IRM1 = IR-1
         I2 = 2**IR
         I1 = I2/2
         I4 = I2+I2
         IFD = IF-I2
         DO 143 I=I2,IFD,I4
            IF (I-NM) 133,133,143
  133       IMI1 = I-I1
            IMI2 = I-I2
            IPI1 = I+I1
            IPI2 = I+I2
            CALL INXCA (I,IR,IDXA,NA)
            CALL INXCC (I,IR,IDXC,NC)
            CALL INXCB (I,IR,IZ,NZ)
            CALL INXCB (IMI1,IRM1,IM1,NM1)
            CALL INXCB (IPI1,IRM1,IP1,NP1)
            IF (I-I2) 134,134,136
  134       DO 135 J=1,M
               W1(J) = (0.,0.)
  135       CONTINUE
            GO TO 137
  136       CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),Y(1,IMI2),
     1                  W1,M,AM,BM,CM,WD,WW,WU)
  137       IF (IPI2-NM) 140,140,138
  138       DO 139 J=1,M
               W2(J) = (0.,0.)
  139       CONTINUE
            GO TO 141
  140       CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),Y(1,IPI2),
     1                  W2,M,AM,BM,CM,WD,WW,WU)
  141       DO 142 J=1,M
               W1(J) = Y(J,I)+W1(J)+W2(J)
  142       CONTINUE
            CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,Y(1,I),
     1                  M,AM,BM,CM,WD,WW,WU)
  143    CONTINUE
  144 CONTINUE
      RETURN
      END
