*DECK BSPLVN
      SUBROUTINE BSPLVN (T, JHIGH, INDEX, X, ILEFT, VNIKX)
C***BEGIN PROLOGUE  BSPLVN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to FC
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BSPLVN-S, DFSPVN-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C Calculates the value of all possibly nonzero B-splines at *X* of
C  order MAX(JHIGH,(J+1)(INDEX-1)) on *T*.
C
C***SEE ALSO  FC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  BSPLVN
      DIMENSION T(*),VNIKX(*)
      DIMENSION DELTAM(20),DELTAP(20)
      SAVE J, DELTAM, DELTAP
      DATA J/1/,(DELTAM(I),I=1,20),(DELTAP(I),I=1,20)/40*0./
C***FIRST EXECUTABLE STATEMENT  BSPLVN
                                       GO TO (10,20),INDEX
   10 J = 1
      VNIKX(1) = 1.
      IF (J .GE. JHIGH)                GO TO 99
C
   20    IPJ = ILEFT+J
         DELTAP(J) = T(IPJ) - X
         IMJP1 = ILEFT-J+1
         DELTAM(J) = X - T(IMJP1)
         VMPREV = 0.
         JP1 = J+1
         DO 26 L=1,J
            JP1ML = JP1-L
            VM = VNIKX(L)/(DELTAP(L) + DELTAM(JP1ML))
            VNIKX(L) = VM*DELTAP(L) + VMPREV
   26       VMPREV = VM*DELTAM(JP1ML)
         VNIKX(JP1) = VMPREV
         J = JP1
         IF (J .LT. JHIGH)             GO TO 20
C
   99                                  RETURN
      END
