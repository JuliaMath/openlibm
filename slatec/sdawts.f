*DECK SDAWTS
      SUBROUTINE SDAWTS (NEQ, IWT, RTOL, ATOL, Y, WT, RPAR, IPAR)
C***BEGIN PROLOGUE  SDAWTS
C***SUBSIDIARY
C***PURPOSE  Set error weight vector for SDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      SINGLE PRECISION (SDAWTS-S, DDAWTS-D)
C***AUTHOR  Petzold, Linda R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------------
C     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR
C     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
C     I=1,-,N.
C     RTOL AND ATOL ARE SCALARS IF IWT = 0,
C     AND VECTORS IF IWT = 1.
C-----------------------------------------------------------------------
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  SDAWTS
C
      INTEGER  NEQ, IWT, IPAR(*)
      REAL  RTOL(*), ATOL(*), Y(*), WT(*), RPAR(*)
C
      INTEGER  I
      REAL  ATOLI, RTOLI
C
C***FIRST EXECUTABLE STATEMENT  SDAWTS
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 20 I=1,NEQ
         IF (IWT .EQ.0) GO TO 10
           RTOLI=RTOL(I)
           ATOLI=ATOL(I)
10         WT(I)=RTOLI*ABS(Y(I))+ATOLI
20         CONTINUE
      RETURN
C-----------END OF SUBROUTINE SDAWTS------------------------------------
      END
