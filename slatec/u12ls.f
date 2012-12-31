*DECK U12LS
      SUBROUTINE U12LS (A, MDA, M, N, B, MDB, NB, MODE, KRANK, RNORM, H,
     +   W, IC, IR)
C***BEGIN PROLOGUE  U12LS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to LLSIA
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (U12LS-S, DU12LS-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C        Given the Householder QR factorization of A, this
C        subroutine solves the system AX=B. If the system
C        is of reduced rank, this routine returns a solution
C        according to the selected mode.
C
C       Note - If MODE.NE.2, W is never accessed.
C
C***SEE ALSO  LLSIA
C***ROUTINES CALLED  SAXPY, SDOT, SNRM2, SSWAP
C***REVISION HISTORY  (YYMMDD)
C   810801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  U12LS
      DIMENSION A(MDA,*),B(MDB,*),RNORM(*),H(*),W(*)
      INTEGER IC(*),IR(*)
C***FIRST EXECUTABLE STATEMENT  U12LS
      K=KRANK
      KP1=K+1
C
C        RANK=0
C
      IF(K.GT.0) GO TO 410
      DO 404 JB=1,NB
      RNORM(JB)=SNRM2(M,B(1,JB),1)
  404 CONTINUE
      DO 406 JB=1,NB
      DO 406 I=1,N
      B(I,JB)=0.0
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
C     APPLY HOUSEHOLDER TRANSFORMATIONS TO B
C
      DO 430 J=1,K
      TT=A(J,J)
      A(J,J)=H(J)
      DO 425 I=1,NB
      BB=-SDOT(M-J+1,A(J,J),1,B(J,I),1)/H(J)
      CALL SAXPY(M-J+1,BB,A(J,J),1,B(J,I),1)
  425 CONTINUE
      A(J,J)=TT
  430 CONTINUE
C
C        FIND NORMS OF RESIDUAL VECTOR(S)..(BEFORE OVERWRITE B)
C
      DO 440 JB=1,NB
      RNORM(JB)=SNRM2((M-K),B(KP1,JB),1)
  440 CONTINUE
C
C     BACK SOLVE UPPER TRIANGULAR R
C
      I=K
  442 DO 444 JB=1,NB
      B(I,JB)=B(I,JB)/A(I,I)
  444 CONTINUE
      IF(I.EQ.1) GO TO 450
      IM1=I-1
      DO 448 JB=1,NB
      CALL SAXPY(IM1,-B(I,JB),A(1,I),1,B(1,JB),1)
  448 CONTINUE
      I=IM1
      GO TO 442
  450 CONTINUE
C
C     RANK LT N
C
C      TRUNCATED SOLUTION
C
      IF(K.EQ.N) GO TO 480
      DO 460 JB=1,NB
      DO 460 I=KP1,N
      B(I,JB)=0.0
  460 CONTINUE
      IF(MODE.EQ.1) GO TO 480
C
C      MINIMAL LENGTH SOLUTION
C
      NMK=N-K
      DO 470 JB=1,NB
      DO 465 I=1,K
      TT=-SDOT(NMK,A(I,KP1),MDA,B(KP1,JB),1)/W(I)
      TT=TT-B(I,JB)
      CALL SAXPY(NMK,TT,A(I,KP1),MDA,B(KP1,JB),1)
      B(I,JB)=B(I,JB)+TT*W(I)
  465 CONTINUE
  470 CONTINUE
C
C
C     REORDER B TO REFLECT COLUMN INTERCHANGES
C
  480 CONTINUE
      I=0
  482 I=I+1
      IF(I.EQ.N) GO TO 488
      J=IC(I)
      IF(J.EQ.I) GO TO 482
      IF(J.LT.0) GO TO 482
      IC(I)=-IC(I)
  484 CALL SSWAP(NB,B(J,1),MDB,B(I,1),MDB)
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
