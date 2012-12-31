*DECK DLPDP
      SUBROUTINE DLPDP (A, MDA, M, N1, N2, PRGOPT, X, WNORM, MODE, WS,
     +   IS)
C***BEGIN PROLOGUE  DLPDP
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DLSEI
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (LPDP-S, DLPDP-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C  **** Double Precision version of LPDP ****
C     DIMENSION A(MDA,N+1),PRGOPT(*),X(N),WS((M+2)*(N+7)),IS(M+N+1),
C     where N=N1+N2.  This is a slight overestimate for WS(*).
C
C     Determine an N1-vector W, and
C               an N2-vector Z
C     which minimizes the Euclidean length of W
C     subject to G*W+H*Z .GE. Y.
C     This is the least projected distance problem, LPDP.
C     The matrices G and H are of respective
C     dimensions M by N1 and M by N2.
C
C     Called by subprogram DLSI( ).
C
C     The matrix
C                (G H Y)
C
C     occupies rows 1,...,M and cols 1,...,N1+N2+1 of A(*,*).
C
C     The solution (W) is returned in X(*).
C                  (Z)
C
C     The value of MODE indicates the status of
C     the computation after returning to the user.
C
C          MODE=1  The solution was successfully obtained.
C
C          MODE=2  The inequalities are inconsistent.
C
C***SEE ALSO  DLSEI
C***ROUTINES CALLED  DCOPY, DDOT, DNRM2, DSCAL, DWNNLS
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR section.  (WRB)
C***END PROLOGUE  DLPDP
C
      INTEGER I, IS(*), IW, IX, J, L, M, MDA, MODE, MODEW, N, N1, N2,
     *     NP1
      DOUBLE PRECISION A(MDA,*), DDOT, DNRM2, FAC, ONE,
     *     PRGOPT(*), RNORM, SC, WNORM, WS(*), X(*), YNORM, ZERO
      SAVE ZERO, ONE, FAC
      DATA ZERO,ONE /0.0D0,1.0D0/, FAC /0.1D0/
C***FIRST EXECUTABLE STATEMENT  DLPDP
      N = N1 + N2
      MODE = 1
      IF (M .GT. 0) GO TO 20
         IF (N .LE. 0) GO TO 10
            X(1) = ZERO
            CALL DCOPY(N,X,0,X,1)
   10    CONTINUE
         WNORM = ZERO
      GO TO 200
   20 CONTINUE
C        BEGIN BLOCK PERMITTING ...EXITS TO 190
            NP1 = N + 1
C
C           SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE.
            DO 40 I = 1, M
               SC = DNRM2(N,A(I,1),MDA)
               IF (SC .EQ. ZERO) GO TO 30
                  SC = ONE/SC
                  CALL DSCAL(NP1,SC,A(I,1),MDA)
   30          CONTINUE
   40       CONTINUE
C
C           SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO).
            YNORM = DNRM2(M,A(1,NP1),1)
            IF (YNORM .EQ. ZERO) GO TO 50
               SC = ONE/YNORM
               CALL DSCAL(M,SC,A(1,NP1),1)
   50       CONTINUE
C
C           SCALE COLS OF MATRIX H.
            J = N1 + 1
   60       IF (J .GT. N) GO TO 70
               SC = DNRM2(M,A(1,J),1)
               IF (SC .NE. ZERO) SC = ONE/SC
               CALL DSCAL(M,SC,A(1,J),1)
               X(J) = SC
               J = J + 1
            GO TO 60
   70       CONTINUE
            IF (N1 .LE. 0) GO TO 130
C
C              COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*).
               IW = 0
               DO 80 I = 1, M
C
C                 MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY.
                  CALL DCOPY(N2,A(I,N1+1),MDA,WS(IW+1),1)
                  IW = IW + N2
C
C                 MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY.
                  CALL DCOPY(N1,A(I,1),MDA,WS(IW+1),1)
                  IW = IW + N1
C
C                 MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY.
                  WS(IW+1) = A(I,NP1)
                  IW = IW + 1
   80          CONTINUE
               WS(IW+1) = ZERO
               CALL DCOPY(N,WS(IW+1),0,WS(IW+1),1)
               IW = IW + N
               WS(IW+1) = ONE
               IW = IW + 1
C
C              SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U.GE.0.  THE
C              MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR
C              F = TRANSPOSE OF (0,...,0,1).
               IX = IW + 1
               IW = IW + M
C
C              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
C              DWNNLS( ).
               IS(1) = 0
               IS(2) = 0
               CALL DWNNLS(WS,NP1,N2,NP1-N2,M,0,PRGOPT,WS(IX),RNORM,
     *                     MODEW,IS,WS(IW+1))
C
C              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W.
               SC = ONE - DDOT(M,A(1,NP1),1,WS(IX),1)
               IF (ONE + FAC*ABS(SC) .EQ. ONE .OR. RNORM .LE. ZERO)
     *            GO TO 110
                  SC = ONE/SC
                  DO 90 J = 1, N1
                     X(J) = SC*DDOT(M,A(1,J),1,WS(IX),1)
   90             CONTINUE
C
C                 COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS
C                 VECTOR.
                  DO 100 I = 1, M
                     A(I,NP1) = A(I,NP1) - DDOT(N1,A(I,1),MDA,X,1)
  100             CONTINUE
               GO TO 120
  110          CONTINUE
                  MODE = 2
C        .........EXIT
                  GO TO 190
  120          CONTINUE
  130       CONTINUE
            IF (N2 .LE. 0) GO TO 180
C
C              COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*).
               IW = 0
               DO 140 I = 1, M
                  CALL DCOPY(N2,A(I,N1+1),MDA,WS(IW+1),1)
                  IW = IW + N2
                  WS(IW+1) = A(I,NP1)
                  IW = IW + 1
  140          CONTINUE
               WS(IW+1) = ZERO
               CALL DCOPY(N2,WS(IW+1),0,WS(IW+1),1)
               IW = IW + N2
               WS(IW+1) = ONE
               IW = IW + 1
               IX = IW + 1
               IW = IW + M
C
C              SOLVE RV=S SUBJECT TO V.GE.0.  THE MATRIX R =(TRANSPOSE
C              OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE
C              OF (0,...,0,1)).
C
C              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
C              DWNNLS( ).
               IS(1) = 0
               IS(2) = 0
               CALL DWNNLS(WS,N2+1,0,N2+1,M,0,PRGOPT,WS(IX),RNORM,MODEW,
     *                     IS,WS(IW+1))
C
C              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z.
               SC = ONE - DDOT(M,A(1,NP1),1,WS(IX),1)
               IF (ONE + FAC*ABS(SC) .EQ. ONE .OR. RNORM .LE. ZERO)
     *            GO TO 160
                  SC = ONE/SC
                  DO 150 J = 1, N2
                     L = N1 + J
                     X(L) = SC*DDOT(M,A(1,L),1,WS(IX),1)*X(L)
  150             CONTINUE
               GO TO 170
  160          CONTINUE
                  MODE = 2
C        .........EXIT
                  GO TO 190
  170          CONTINUE
  180       CONTINUE
C
C           ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION.
            CALL DSCAL(N,YNORM,X,1)
            WNORM = DNRM2(N1,X,1)
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
