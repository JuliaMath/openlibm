*DECK POISP2
      SUBROUTINE POISP2 (M, N, A, BB, C, Q, IDIMQ, B, B2, B3, W, W2, W3,
     +   D, TCOS, P)
C***BEGIN PROLOGUE  POISP2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (POISP2-S, CMPOSP-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Subroutine to solve Poisson equation with periodic boundary
C     conditions.
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  POISD2, POISN2
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  POISP2
C
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                P(*)
C***FIRST EXECUTABLE STATEMENT  POISP2
      MR = M
      NR = (N+1)/2
      NRM1 = NR-1
      IF (2*NR .NE. N) GO TO 107
C
C     EVEN NUMBER OF UNKNOWNS
C
      DO 102 J=1,NRM1
         NRMJ = NR-J
         NRPJ = NR+J
         DO 101 I=1,MR
            S = Q(I,NRMJ)-Q(I,NRPJ)
            T = Q(I,NRMJ)+Q(I,NRPJ)
            Q(I,NRMJ) = S
            Q(I,NRPJ) = T
  101    CONTINUE
  102 CONTINUE
      DO 103 I=1,MR
         Q(I,NR) = 2.*Q(I,NR)
         Q(I,N) = 2.*Q(I,N)
  103 CONTINUE
      CALL POISD2 (MR,NRM1,1,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
      IPSTOR = W(1)
      CALL POISN2 (MR,NR+1,1,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,
     1             TCOS,P)
      IPSTOR = MAX(IPSTOR,INT(W(1)))
      DO 105 J=1,NRM1
         NRMJ = NR-J
         NRPJ = NR+J
         DO 104 I=1,MR
            S = .5*(Q(I,NRPJ)+Q(I,NRMJ))
            T = .5*(Q(I,NRPJ)-Q(I,NRMJ))
            Q(I,NRMJ) = S
            Q(I,NRPJ) = T
  104    CONTINUE
  105 CONTINUE
      DO 106 I=1,MR
         Q(I,NR) = .5*Q(I,NR)
         Q(I,N) = .5*Q(I,N)
  106 CONTINUE
      GO TO 118
  107 CONTINUE
C
C     ODD  NUMBER OF UNKNOWNS
C
      DO 109 J=1,NRM1
         NRPJ = N+1-J
         DO 108 I=1,MR
            S = Q(I,J)-Q(I,NRPJ)
            T = Q(I,J)+Q(I,NRPJ)
            Q(I,J) = S
            Q(I,NRPJ) = T
  108    CONTINUE
  109 CONTINUE
      DO 110 I=1,MR
         Q(I,NR) = 2.*Q(I,NR)
  110 CONTINUE
      LH = NRM1/2
      DO 112 J=1,LH
         NRMJ = NR-J
         DO 111 I=1,MR
            S = Q(I,J)
            Q(I,J) = Q(I,NRMJ)
            Q(I,NRMJ) = S
  111    CONTINUE
  112 CONTINUE
      CALL POISD2 (MR,NRM1,2,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
      IPSTOR = W(1)
      CALL POISN2 (MR,NR,2,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,
     1             TCOS,P)
      IPSTOR = MAX(IPSTOR,INT(W(1)))
      DO 114 J=1,NRM1
         NRPJ = NR+J
         DO 113 I=1,MR
            S = .5*(Q(I,NRPJ)+Q(I,J))
            T = .5*(Q(I,NRPJ)-Q(I,J))
            Q(I,NRPJ) = T
            Q(I,J) = S
  113    CONTINUE
  114 CONTINUE
      DO 115 I=1,MR
         Q(I,NR) = .5*Q(I,NR)
  115 CONTINUE
      DO 117 J=1,LH
         NRMJ = NR-J
         DO 116 I=1,MR
            S = Q(I,J)
            Q(I,J) = Q(I,NRMJ)
            Q(I,NRMJ) = S
  116    CONTINUE
  117 CONTINUE
  118 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
