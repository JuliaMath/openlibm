*DECK DWNLT1
      SUBROUTINE DWNLT1 (I, LEND, MEND, IR, MDW, RECALC, IMAX, HBAR, H,
     +   SCALE, W)
C***BEGIN PROLOGUE  DWNLT1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIT
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (WNLT1-S, DWNLT1-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     To update the column Sum Of Squares and find the pivot column.
C     The column Sum of Squares Vector will be updated at each step.
C     When numerically necessary, these values will be recomputed.
C
C***SEE ALSO  DWNLIT
C***ROUTINES CALLED  IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
C   900604  DP version created from SP version.  (RWC)
C***END PROLOGUE  DWNLT1
      INTEGER I, IMAX, IR, LEND, MDW, MEND
      DOUBLE PRECISION H(*), HBAR, SCALE(*), W(MDW,*)
      LOGICAL RECALC
C
      EXTERNAL IDAMAX
      INTEGER IDAMAX
C
      INTEGER J, K
C
C***FIRST EXECUTABLE STATEMENT  DWNLT1
      IF (IR.NE.1 .AND. (.NOT.RECALC)) THEN
C
C        Update column SS=sum of squares.
C
         DO 10 J=I,LEND
            H(J) = H(J) - SCALE(IR-1)*W(IR-1,J)**2
   10    CONTINUE
C
C        Test for numerical accuracy.
C
         IMAX = IDAMAX(LEND-I+1, H(I), 1) + I - 1
         RECALC = (HBAR+1.E-3*H(IMAX)) .EQ. HBAR
      ENDIF
C
C     If required, recalculate column SS, using rows IR through MEND.
C
      IF (RECALC) THEN
         DO 30 J=I,LEND
            H(J) = 0.D0
            DO 20 K=IR,MEND
               H(J) = H(J) + SCALE(K)*W(K,J)**2
   20       CONTINUE
   30    CONTINUE
C
C        Find column with largest SS.
C
         IMAX = IDAMAX(LEND-I+1, H(I), 1) + I - 1
         HBAR = H(IMAX)
      ENDIF
      RETURN
      END
