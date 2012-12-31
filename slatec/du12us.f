*DECK DU12US
      SUBROUTINE DU12US (A, MDA, M, N, B, MDB, NB, MODE, KRANK, RNORM,
     +   H, W, IR, IC)
C***BEGIN PROLOGUE  DU12US
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DULSIA
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (U12US-S, DU12US-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C        Given the Householder LQ factorization of A, this
C        subroutine solves the system AX=B. If the system
C        is of reduced rank, this routine returns a solution
C        according to the selected mode.
C
C       Note - If MODE.NE.2, W is never accessed.
C
C***SEE ALSO  DULSIA
C***ROUTINES CALLED  DAXPY, DDOT, DNRM2, DSWAP
C***REVISION HISTORY  (YYMMDD)
C   810801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DU12US
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DDOT,DNRM2
      DIMENSION A(MDA,*),B(MDB,*),RNORM(*),H(*),W(*)
      INTEGER IC(*),IR(*)
C***FIRST EXECUTABLE STATEMENT  DU12US
      K=KRANK
      KP1=K+1
C
C        RANK=0
C
      IF(K.GT.0) GO TO 410
      DO 404 JB=1,NB
      RNORM(JB)=DNRM2(M,B(1,JB),1)
  404 CONTINUE
      DO 406 JB=1,NB
      DO 406 I=1,N
      B(I,JB)=0.0D0
  406 CONTINUE
      RETURN
C
C     REORDER B TO REFLECT ROW INTERCHANGES
C
  410 CONTINUE
      I=0
  412 I=I+1
      IF(I.EQ.M) GO TO 418
      J=IR(I)
      IF(J.EQ.I) GO TO 412
      IF(J.LT.0) GO TO 412
      IR(I)=-IR(I)
      DO 413 JB=1,NB
      RNORM(JB)=B(I,JB)
  413 CONTINUE
      IJ=I
  414 DO 415 JB=1,NB
      B(IJ,JB)=B(J,JB)
  415 CONTINUE
      IJ=J
      J=IR(IJ)
      IR(IJ)=-IR(IJ)
      IF(J.NE.I) GO TO 414
      DO 416 JB=1,NB
      B(IJ,JB)=RNORM(JB)
  416 CONTINUE
      GO TO 412
  418 CONTINUE
      DO 420 I=1,M
      IR(I)=ABS(IR(I))
  420 CONTINUE
C
C     IF A IS OF REDUCED RANK AND MODE=2,
C     APPLY HOUSEHOLDER TRANSFORMATIONS TO B
C
      IF(MODE.LT.2 .OR. K.EQ.M) GO TO 440
      MMK=M-K
      DO 430 JB=1,NB
      DO 425 J=1,K
      I=KP1-J
      TT=-DDOT(MMK,A(KP1,I),1,B(KP1,JB),1)/W(I)
      TT=TT-B(I,JB)
      CALL DAXPY(MMK,TT,A(KP1,I),1,B(KP1,JB),1)
      B(I,JB)=B(I,JB)+TT*W(I)
  425 CONTINUE
  430 CONTINUE
C
C     FIND NORMS OF RESIDUAL VECTOR(S)..(BEFORE OVERWRITE B)
C
  440 DO 442 JB=1,NB
      RNORM(JB)=DNRM2((M-K),B(KP1,JB),1)
  442 CONTINUE
C
C     BACK SOLVE LOWER TRIANGULAR L
C
      DO 450 JB=1,NB
      DO 448 I=1,K
      B(I,JB)=B(I,JB)/A(I,I)
      IF(I.EQ.K) GO TO 450
      IP1=I+1
      CALL DAXPY(K-I,-B(I,JB),A(IP1,I),1,B(IP1,JB),1)
  448 CONTINUE
  450 CONTINUE
C
C
C      TRUNCATED SOLUTION
C
      IF(K.EQ.N) GO TO 462
      DO 460 JB=1,NB
      DO 460 I=KP1,N
      B(I,JB)=0.0D0
  460 CONTINUE
C
C     APPLY HOUSEHOLDER TRANSFORMATIONS TO B
C
  462 DO 470 I=1,K
      J=KP1-I
      TT=A(J,J)
      A(J,J)=H(J)
      DO 465 JB=1,NB
      BB=-DDOT(N-J+1,A(J,J),MDA,B(J,JB),1)/H(J)
      CALL DAXPY(N-J+1,BB,A(J,J),MDA,B(J,JB),1)
  465 CONTINUE
      A(J,J)=TT
  470 CONTINUE
C
C
C     REORDER B TO REFLECT COLUMN INTERCHANGES
C
      I=0
  482 I=I+1
      IF(I.EQ.N) GO TO 488
      J=IC(I)
      IF(J.EQ.I) GO TO 482
      IF(J.LT.0) GO TO 482
      IC(I)=-IC(I)
  484 CALL DSWAP(NB,B(J,1),MDB,B(I,1),MDB)
      IJ=IC(J)
      IC(J)=-IC(J)
      J=IJ
      IF(J.EQ.I) GO TO 482
      GO TO 484
  488 CONTINUE
      DO 490 I=1,N
      IC(I)=ABS(IC(I))
  490 CONTINUE
C
C        SOLUTION VECTORS ARE IN FIRST N ROWS OF B(,)
C
      RETURN
      END
