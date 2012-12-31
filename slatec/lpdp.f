*DECK LPDP
      SUBROUTINE LPDP (A, MDA, M, N1, N2, PRGOPT, X, WNORM, MODE, WS,
     +   IS)
C***BEGIN PROLOGUE  LPDP
C***SUBSIDIARY
C***PURPOSE  Subsidiary to LSEI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LPDP-S, DLPDP-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
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
C     Called by subprogram LSI( ).
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
C***SEE ALSO  LSEI
C***ROUTINES CALLED  SCOPY, SDOT, SNRM2, SSCAL, WNNLS
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR section.  (WRB)
C***END PROLOGUE  LPDP
C
C     SUBROUTINES CALLED
C
C     WNNLS         SOLVES A NONNEGATIVELY CONSTRAINED LINEAR LEAST
C                   SQUARES PROBLEM WITH LINEAR EQUALITY CONSTRAINTS.
C                   PART OF THIS PACKAGE.
C
C++
C     SDOT,         SUBROUTINES FROM THE BLAS PACKAGE.
C     SSCAL,SNRM2,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308.
C     SCOPY
C
      REAL             A(MDA,*), PRGOPT(*), WS(*), WNORM, X(*)
      INTEGER IS(*)
      REAL             FAC, ONE, RNORM, SC, YNORM, ZERO
      REAL             SDOT, SNRM2
      SAVE ZERO, ONE, FAC
      DATA ZERO, ONE /0.E0,1.E0/, FAC /0.1E0/
C***FIRST EXECUTABLE STATEMENT  LPDP
      N = N1 + N2
      MODE = 1
      IF (.NOT.(M.LE.0)) GO TO 20
      IF (.NOT.(N.GT.0)) GO TO 10
      X(1) = ZERO
      CALL SCOPY(N, X, 0, X, 1)
   10 WNORM = ZERO
      RETURN
   20 NP1 = N + 1
C
C     SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE.
      DO 40 I=1,M
        SC = SNRM2(N,A(I,1),MDA)
        IF (.NOT.(SC.NE.ZERO)) GO TO 30
        SC = ONE/SC
        CALL SSCAL(NP1, SC, A(I,1), MDA)
   30   CONTINUE
   40 CONTINUE
C
C     SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO).
      YNORM = SNRM2(M,A(1,NP1),1)
      IF (.NOT.(YNORM.NE.ZERO)) GO TO 50
      SC = ONE/YNORM
      CALL SSCAL(M, SC, A(1,NP1), 1)
C
C     SCALE COLS OF MATRIX H.
   50 J = N1 + 1
   60 IF (.NOT.(J.LE.N)) GO TO 70
      SC = SNRM2(M,A(1,J),1)
      IF (SC.NE.ZERO) SC = ONE/SC
      CALL SSCAL(M, SC, A(1,J), 1)
      X(J) = SC
      J = J + 1
      GO TO 60
   70 IF (.NOT.(N1.GT.0)) GO TO 130
C
C     COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*).
      IW = 0
      DO 80 I=1,M
C
C     MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY.
        CALL SCOPY(N2, A(I,N1+1), MDA, WS(IW+1), 1)
        IW = IW + N2
C
C     MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY.
        CALL SCOPY(N1, A(I,1), MDA, WS(IW+1), 1)
        IW = IW + N1
C
C     MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY.
        WS(IW+1) = A(I,NP1)
        IW = IW + 1
   80 CONTINUE
      WS(IW+1) = ZERO
      CALL SCOPY(N, WS(IW+1), 0, WS(IW+1), 1)
      IW = IW + N
      WS(IW+1) = ONE
      IW = IW + 1
C
C     SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U.GE.0.  THE
C     MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR
C     F = TRANSPOSE OF (0,...,0,1).
      IX = IW + 1
      IW = IW + M
C
C     DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF WNNLS( ).
      IS(1) = 0
      IS(2) = 0
      CALL WNNLS(WS, NP1, N2, NP1-N2, M, 0, PRGOPT, WS(IX), RNORM,
     1 MODEW, IS, WS(IW+1))
C
C     COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W.
      SC = ONE - SDOT(M,A(1,NP1),1,WS(IX),1)
      IF (.NOT.(ONE+FAC*ABS(SC).NE.ONE .AND. RNORM.GT.ZERO)) GO TO 110
      SC = ONE/SC
      DO 90 J=1,N1
        X(J) = SC*SDOT(M,A(1,J),1,WS(IX),1)
   90 CONTINUE
C
C     COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS VECTOR.
      DO 100 I=1,M
        A(I,NP1) = A(I,NP1) - SDOT(N1,A(I,1),MDA,X,1)
  100 CONTINUE
      GO TO 120
  110 MODE = 2
      RETURN
  120 CONTINUE
  130 IF (.NOT.(N2.GT.0)) GO TO 180
C
C     COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*).
      IW = 0
      DO 140 I=1,M
        CALL SCOPY(N2, A(I,N1+1), MDA, WS(IW+1), 1)
        IW = IW + N2
        WS(IW+1) = A(I,NP1)
        IW = IW + 1
  140 CONTINUE
      WS(IW+1) = ZERO
      CALL SCOPY(N2, WS(IW+1), 0, WS(IW+1), 1)
      IW = IW + N2
      WS(IW+1) = ONE
      IW = IW + 1
      IX = IW + 1
      IW = IW + M
C
C     SOLVE RV=S SUBJECT TO V.GE.0.  THE MATRIX R =(TRANSPOSE
C     OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE
C     OF (0,...,0,1)).
C
C     DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF WNNLS( ).
      IS(1) = 0
      IS(2) = 0
      CALL WNNLS(WS, N2+1, 0, N2+1, M, 0, PRGOPT, WS(IX), RNORM, MODEW,
     1 IS, WS(IW+1))
C
C     COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z.
      SC = ONE - SDOT(M,A(1,NP1),1,WS(IX),1)
      IF (.NOT.(ONE+FAC*ABS(SC).NE.ONE .AND. RNORM.GT.ZERO)) GO TO 160
      SC = ONE/SC
      DO 150 J=1,N2
        L = N1 + J
        X(L) = SC*SDOT(M,A(1,L),1,WS(IX),1)*X(L)
  150 CONTINUE
      GO TO 170
  160 MODE = 2
      RETURN
  170 CONTINUE
C
C     ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION.
  180 CALL SSCAL(N, YNORM, X, 1)
      WNORM = SNRM2(N1,X,1)
      RETURN
      END
