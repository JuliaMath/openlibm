*DECK DSDI
      SUBROUTINE DSDI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C***BEGIN PROLOGUE  DSDI
C***PURPOSE  Diagonal Matrix Vector Multiply.
C            Routine to calculate the product  X = DIAG*B, where DIAG
C            is a diagonal matrix.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D1B4
C***TYPE      DOUBLE PRECISION (SSDI-S, DSDI-D)
C***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(10)
C     DOUBLE PRECISION B(N), X(N), A(NELT), RWORK(USER DEFINED)
C
C     CALL DSDI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Vector to multiply the diagonal by.
C X      :OUT      Double Precision X(N).
C         Result of DIAG*B.
C NELT   :DUMMY    Integer.
C IA     :DUMMY    Integer IA(NELT).
C JA     :DUMMY    Integer JA(NELT).
C A      :DUMMY    Double Precision A(NELT).
C ISYM   :DUMMY    Integer.
C         These are for compatibility with SLAP MSOLVE calling sequence.
C RWORK  :IN       Double Precision RWORK(USER DEFINED).
C         Work array holding the diagonal of some matrix to scale
C         B by.  This array must be set by the user or by a call
C         to the SLAP routine DSDS or DSD2S.  The length of RWORK
C         must be >= IWORK(4)+N.
C IWORK  :IN       Integer IWORK(10).
C         IWORK(4) holds the offset into RWORK for the diagonal matrix
C         to scale B by.  This is usually set up by the SLAP pre-
C         conditioner setup routines DSDS or DSD2S.
C
C *Description:
C         This routine is supplied with the SLAP package to perform
C         the  MSOLVE  operation for iterative drivers that require
C         diagonal  Scaling  (e.g., DSDCG, DSDBCG).   It  conforms
C         to the SLAP MSOLVE CALLING CONVENTION  and hence does not
C         require an interface routine as do some of the other pre-
C         conditioners supplied with SLAP.
C
C***SEE ALSO  DSDS, DSD2S
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   871119  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   920511  Added complete declaration section.  (WRB)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  DSDI
C     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RWORK(*), X(N)
      INTEGER IA(NELT), IWORK(10), JA(NELT)
C     .. Local Scalars ..
      INTEGER I, LOCD
C***FIRST EXECUTABLE STATEMENT  DSDI
C
C         Determine where the inverse of the diagonal
C         is in the work array and then scale by it.
C
      LOCD = IWORK(4) - 1
      DO 10 I = 1, N
         X(I) = RWORK(LOCD+I)*B(I)
 10   CONTINUE
      RETURN
C------------- LAST LINE OF DSDI FOLLOWS ----------------------------
      END
