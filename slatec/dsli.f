*DECK DSLI
      SUBROUTINE DSLI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C***BEGIN PROLOGUE  DSLI
C***PURPOSE  SLAP MSOLVE for Lower Triangle Matrix.
C            This routine acts as an interface between the SLAP generic
C            MSOLVE calling convention and the routine that actually
C                      -1
C            computes L  B = X.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A3
C***TYPE      DOUBLE PRECISION (SSLI-S, DSLI-D)
C***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C       It is assumed that RWORK and IWORK have initialized with
C       the information required for DSLI2:
C          IWORK(1) = NEL
C          IWORK(2) = Starting location of IEL in IWORK.
C          IWORK(3) = Starting location of JEL in IWORK.
C          IWORK(4) = Starting location of EL in RWORK.
C       See the DESCRIPTION of DSLI2 for details.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DSLI2
C***REVISION HISTORY  (YYMMDD)
C   871119  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   920511  Added complete declaration section.  (WRB)
C   921113  Corrected C***CATEGORY line.  (FNF)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  DSLI
C     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RWORK(*), X(N)
      INTEGER IA(NELT), IWORK(10), JA(NELT)
C     .. Local Scalars ..
      INTEGER LOCEL, LOCIEL, LOCJEL, NEL
C     .. External Subroutines ..
      EXTERNAL DSLI2
C***FIRST EXECUTABLE STATEMENT  DSLI
C
      NEL = IWORK(1)
      LOCIEL = IWORK(2)
      LOCJEL = IWORK(3)
      LOCEL = IWORK(4)
      CALL DSLI2(N, B, X, NEL, IWORK(LOCIEL), IWORK(LOCJEL),
     $     RWORK(LOCEL))
C
      RETURN
C------------- LAST LINE OF DSLI FOLLOWS ----------------------------
      END
