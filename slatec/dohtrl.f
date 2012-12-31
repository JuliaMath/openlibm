*DECK DOHTRL
      SUBROUTINE DOHTRL (Q, N, NRDA, DIAG, IRANK, DIV, TD)
C***BEGIN PROLOGUE  DOHTRL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP and DSUDS
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (OHTROL-S, DOHTRL-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     For a rank deficient problem, additional orthogonal
C     HOUSEHOLDER transformations are applied to the left side
C     of Q to further reduce the triangular form.
C     Thus, after application of the routines DORTHR and DOHTRL
C     to the original matrix, the result is a nonsingular
C     triangular matrix while the remainder of the matrix
C     has been zeroed out.
C
C***SEE ALSO  DBVSUP, DSUDS
C***ROUTINES CALLED  DDOT
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DOHTRL
      DOUBLE PRECISION DDOT
      INTEGER IRANK, IRP, J, K, KIR, KIRM, L, N, NMIR, NRDA
      DOUBLE PRECISION DD, DIAG(*), DIAGK, DIV(*), Q(NRDA,*), QS, SIG,
     1     SQD, TD(*), TDV
C***FIRST EXECUTABLE STATEMENT  DOHTRL
      NMIR = N - IRANK
      IRP = IRANK + 1
      DO 40 K = 1, IRANK
         KIR = IRP - K
         DIAGK = DIAG(KIR)
         SIG = (DIAGK*DIAGK) + DDOT(NMIR,Q(IRP,KIR),1,Q(IRP,KIR),1)
         DD = SIGN(SQRT(SIG),-DIAGK)
         DIV(KIR) = DD
         TDV = DIAGK - DD
         TD(KIR) = TDV
         IF (K .EQ. IRANK) GO TO 30
            KIRM = KIR - 1
            SQD = DD*DIAGK - SIG
            DO 20 J = 1, KIRM
               QS = ((TDV*Q(KIR,J))
     1               + DDOT(NMIR,Q(IRP,J),1,Q(IRP,KIR),1))/SQD
               Q(KIR,J) = Q(KIR,J) + QS*TDV
               DO 10 L = IRP, N
                  Q(L,J) = Q(L,J) + QS*Q(L,KIR)
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
