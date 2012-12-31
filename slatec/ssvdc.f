*DECK SSVDC
      SUBROUTINE SSVDC (X, LDX, N, P, S, E, U, LDU, V, LDV, WORK, JOB,
     +   INFO)
C***BEGIN PROLOGUE  SSVDC
C***PURPOSE  Perform the singular value decomposition of a rectangular
C            matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D6
C***TYPE      SINGLE PRECISION (SSVDC-S, DSVDC-D, CSVDC-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX,
C             SINGULAR VALUE DECOMPOSITION
C***AUTHOR  Stewart, G. W., (U. of Maryland)
C***DESCRIPTION
C
C     SSVDC is a subroutine to reduce a real NxP matrix X by orthogonal
C     transformations U and V to diagonal form.  The elements S(I) are
C     the singular values of X.  The columns of U are the corresponding
C     left singular vectors, and the columns of V the right singular
C     vectors.
C
C     On Entry
C
C         X         REAL(LDX,P), where LDX .GE. N.
C                   X contains the matrix whose singular value
C                   decomposition is to be computed.  X is
C                   destroyed by SSVDC.
C
C         LDX       INTEGER
C                   LDX is the leading dimension of the array X.
C
C         N         INTEGER
C                   N is the number of rows of the matrix X.
C
C         P         INTEGER
C                   P is the number of columns of the matrix X.
C
C         LDU       INTEGER
C                   LDU is the leading dimension of the array U.
C                   (See below).
C
C         LDV       INTEGER
C                   LDV is the leading dimension of the array V.
C                   (See below).
C
C         WORK      REAL(N)
C                   work is a scratch array.
C
C         JOB       INTEGER
C                   JOB controls the computation of the singular
C                   vectors.  It has the decimal expansion AB
C                   with the following meaning
C
C                        A .EQ. 0  Do not compute the left singular
C                                  vectors.
C                        A .EQ. 1  Return the N left singular vectors
C                                  in U.
C                        A .GE. 2  Return the first MIN(N,P) singular
C                                  vectors in U.
C                        B .EQ. 0  Do not compute the right singular
C                                  vectors.
C                        B .EQ. 1  Return the right singular vectors
C                                  in V.
C
C     On Return
C
C         S         REAL(MM), where MM=MIN(N+1,P).
C                   The first MIN(N,P) entries of S contain the
C                   singular values of X arranged in descending
C                   order of magnitude.
C
C         E         REAL(P).
C                   E ordinarily contains zeros.  However, see the
C                   discussion of INFO for exceptions.
C
C         U         REAL(LDU,K), where LDU .GE. N.  If JOBA .EQ. 1, then
C                                   K .EQ. N.  If JOBA .GE. 2 , then
C                                   K .EQ. MIN(N,P).
C                   U contains the matrix of right singular vectors.
C                   U is not referenced if JOBA .EQ. 0.  If N .LE. P
C                   or if JOBA .EQ. 2, then U may be identified with X
C                   in the subroutine call.
C
C         V         REAL(LDV,P), where LDV .GE. P.
C                   V contains the matrix of right singular vectors.
C                   V is not referenced if JOB .EQ. 0.  If P .LE. N,
C                   then V may be identified with X in the
C                   subroutine call.
C
C         INFO      INTEGER.
C                   the singular values (and their corresponding
C                   singular vectors) S(INFO+1),S(INFO+2),...,S(M)
C                   are correct (here M=MIN(N,P)).  Thus if
C                   INFO .EQ. 0, all the singular values and their
C                   vectors are correct.  In any event, the matrix
C                   B = TRANS(U)*X*V is the bidiagonal matrix
C                   with the elements of S on its diagonal and the
C                   elements of E on its super-diagonal (TRANS(U)
C                   is the transpose of U).  Thus the singular
C                   values of X and B are the same.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  SAXPY, SDOT, SNRM2, SROT, SROTG, SSCAL, SSWAP
C***REVISION HISTORY  (YYMMDD)
C   790319  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SSVDC
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      REAL X(LDX,*),S(*),E(*),U(LDU,*),V(LDV,*),WORK(*)
C
C
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,
     1        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      REAL SDOT,T
      REAL B,C,CS,EL,EMM1,F,G,SNRM2,SCALE,SHIFT,SL,SM,SN,SMM1,T1,TEST,
     1     ZTEST
      LOGICAL WANTU,WANTV
C***FIRST EXECUTABLE STATEMENT  SSVDC
C
C     SET THE MAXIMUM NUMBER OF ITERATIONS.
C
      MAXIT = 30
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
C
C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
C
      INFO = 0
      NCT = MIN(N-1,P)
      NRT = MAX(0,MIN(P-2,N))
      LU = MAX(NCT,NRT)
      IF (LU .LT. 1) GO TO 170
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .GT. NCT) GO TO 20
C
C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
C           PLACE THE L-TH DIAGONAL IN S(L).
C
            S(L) = SNRM2(N-L+1,X(L,L),1)
            IF (S(L) .EQ. 0.0E0) GO TO 10
               IF (X(L,L) .NE. 0.0E0) S(L) = SIGN(S(L),X(L,L))
               CALL SSCAL(N-L+1,1.0E0/S(L),X(L,L),1)
               X(L,L) = 1.0E0 + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (S(L) .EQ. 0.0E0) GO TO 30
C
C              APPLY THE TRANSFORMATION.
C
               T = -SDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
               CALL SAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
            E(J) = X(L,J)
   40    CONTINUE
   50    CONTINUE
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70
C
C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
C           MULTIPLICATION.
C
            DO 60 I = L, N
               U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
         IF (L .GT. NRT) GO TO 150
C
C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
C           L-TH SUPER-DIAGONAL IN E(L).
C
            E(L) = SNRM2(P-L,E(LP1),1)
            IF (E(L) .EQ. 0.0E0) GO TO 80
               IF (E(LP1) .NE. 0.0E0) E(L) = SIGN(E(L),E(LP1))
               CALL SSCAL(P-L,1.0E0/E(L),E(LP1),1)
               E(LP1) = 1.0E0 + E(LP1)
   80       CONTINUE
            E(L) = -E(L)
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0E0) GO TO 120
C
C              APPLY THE TRANSFORMATION.
C
               DO 90 I = LP1, N
                  WORK(I) = 0.0E0
   90          CONTINUE
               DO 100 J = LP1, P
                  CALL SAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
               DO 110 J = LP1, P
                  CALL SAXPY(N-L,-E(J)/E(LP1),WORK(LP1),1,X(LP1,J),1)
  110          CONTINUE
  120       CONTINUE
            IF (.NOT.WANTV) GO TO 140
C
C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
C              BACK MULTIPLICATION.
C
               DO 130 I = LP1, P
                  V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
C
C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
C
      M = MIN(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = 0.0E0
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0E0
C
C     IF REQUIRED, GENERATE U.
C
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = 0.0E0
  180       CONTINUE
            U(J,J) = 1.0E0
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (S(L) .EQ. 0.0E0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
               DO 210 J = LP1, NCU
                  T = -SDOT(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
                  CALL SAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
               CALL SSCAL(N-L+1,-1.0E0,U(L,L),1)
               U(L,L) = 1.0E0 + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = 0.0E0
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = 0.0E0
  260          CONTINUE
               U(L,L) = 1.0E0
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
C
C     IF IT IS REQUIRED, GENERATE V.
C
      IF (.NOT.WANTV) GO TO 350
         DO 340 LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (E(L) .EQ. 0.0E0) GO TO 320
               DO 310 J = LP1, P
                  T = -SDOT(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
                  CALL SAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = 0.0E0
  330       CONTINUE
            V(L,L) = 1.0E0
  340    CONTINUE
  350 CONTINUE
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
      MM = M
      ITER = 0
  360 CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
         IF (M .EQ. 0) GO TO 620
C
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
C        FLAG AND RETURN.
C
         IF (ITER .LT. MAXIT) GO TO 370
            INFO = M
            GO TO 620
  370    CONTINUE
C
C        THIS SECTION OF THE PROGRAM INSPECTS FOR
C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
C        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
C
C           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
C           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
C           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
C                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
C           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
C
         DO 390 LL = 1, M
            L = M - LL
            IF (L .EQ. 0) GO TO 400
            TEST = ABS(S(L)) + ABS(S(L+1))
            ZTEST = TEST + ABS(E(L))
            IF (ZTEST .NE. TEST) GO TO 380
               E(L) = 0.0E0
               GO TO 400
  380       CONTINUE
  390    CONTINUE
  400    CONTINUE
         IF (L .NE. M - 1) GO TO 410
            KASE = 4
         GO TO 480
  410    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 430 LLS = LP1, MP1
               LS = M - LLS + LP1
               IF (LS .EQ. L) GO TO 440
               TEST = 0.0E0
               IF (LS .NE. M) TEST = TEST + ABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + ABS(E(LS-1))
               ZTEST = TEST + ABS(S(LS))
               IF (ZTEST .NE. TEST) GO TO 420
                  S(LS) = 0.0E0
                  GO TO 440
  420          CONTINUE
  430       CONTINUE
  440       CONTINUE
            IF (LS .NE. L) GO TO 450
               KASE = 3
            GO TO 470
  450       CONTINUE
            IF (LS .NE. M) GO TO 460
               KASE = 1
            GO TO 470
  460       CONTINUE
               KASE = 2
               L = LS
  470       CONTINUE
  480    CONTINUE
         L = L + 1
C
C        PERFORM THE TASK INDICATED BY KASE.
C
         GO TO (490,520,540,570), KASE
C
C        DEFLATE NEGLIGIBLE S(M).
C
  490    CONTINUE
            MM1 = M - 1
            F = E(M-1)
            E(M-1) = 0.0E0
            DO 510 KK = L, MM1
               K = MM1 - KK + L
               T1 = S(K)
               CALL SROTG(T1,F,CS,SN)
               S(K) = T1
               IF (K .EQ. L) GO TO 500
                  F = -SN*E(K-1)
                  E(K-1) = CS*E(K-1)
  500          CONTINUE
               IF (WANTV) CALL SROT(P,V(1,K),1,V(1,M),1,CS,SN)
  510       CONTINUE
         GO TO 610
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  520    CONTINUE
            F = E(L-1)
            E(L-1) = 0.0E0
            DO 530 K = L, M
               T1 = S(K)
               CALL SROTG(T1,F,CS,SN)
               S(K) = T1
               F = -SN*E(K)
               E(K) = CS*E(K)
               IF (WANTU) CALL SROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  530       CONTINUE
         GO TO 610
C
C        PERFORM ONE QR STEP.
C
  540    CONTINUE
C
C           CALCULATE THE SHIFT.
C
            SCALE = MAX(ABS(S(M)),ABS(S(M-1)),ABS(E(M-1)),ABS(S(L)),
     1                    ABS(E(L)))
            SM = S(M)/SCALE
            SMM1 = S(M-1)/SCALE
            EMM1 = E(M-1)/SCALE
            SL = S(L)/SCALE
            EL = E(L)/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0E0
            C = (SM*EMM1)**2
            SHIFT = 0.0E0
            IF (B .EQ. 0.0E0 .AND. C .EQ. 0.0E0) GO TO 550
               SHIFT = SQRT(B**2+C)
               IF (B .LT. 0.0E0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  550       CONTINUE
            F = (SL + SM)*(SL - SM) - SHIFT
            G = SL*EL
C
C           CHASE ZEROS.
C
            MM1 = M - 1
            DO 560 K = L, MM1
               CALL SROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = F
               F = CS*S(K) + SN*E(K)
               E(K) = CS*E(K) - SN*S(K)
               G = SN*S(K+1)
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL SROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL SROTG(F,G,CS,SN)
               S(K) = F
               F = CS*E(K) + SN*S(K+1)
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*E(K+1)
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N)
     1            CALL SROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  560       CONTINUE
            E(M-1) = F
            ITER = ITER + 1
         GO TO 610
C
C        CONVERGENCE.
C
  570    CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE.
C
            IF (S(L) .GE. 0.0E0) GO TO 580
               S(L) = -S(L)
               IF (WANTV) CALL SSCAL(P,-1.0E0,V(1,L),1)
  580       CONTINUE
C
C           ORDER THE SINGULAR VALUE.
C
  590       IF (L .EQ. MM) GO TO 600
               IF (S(L) .GE. S(L+1)) GO TO 600
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P)
     1            CALL SSWAP(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N)
     1            CALL SSWAP(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 590
  600       CONTINUE
            ITER = 0
            M = M - 1
  610    CONTINUE
      GO TO 360
  620 CONTINUE
      RETURN
      END
