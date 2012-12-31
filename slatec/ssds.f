*DECK SSDS
      SUBROUTINE SSDS (N, NELT, IA, JA, A, ISYM, DINV)
C***BEGIN PROLOGUE  SSDS
C***PURPOSE  Diagonal Scaling Preconditioner SLAP Set Up.
C            Routine to compute the inverse of the diagonal of a matrix
C            stored in the SLAP Column format.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2E
C***TYPE      SINGLE PRECISION (SSDS-S, DSDS-D)
C***KEYWORDS  DIAGONAL, SLAP SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     REAL    A(NELT), DINV(N)
C
C     CALL SSDS( N, NELT, IA, JA, A, ISYM, DINV )
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Real A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C DINV   :OUT      Real DINV(N).
C         Upon return this array holds 1./DIAG(A).
C
C *Description
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       With the SLAP  format  all  of  the   "inner  loops" of this
C       routine should vectorize  on  machines with hardware support
C       for vector   gather/scatter  operations.  Your compiler  may
C       require a compiler directive to  convince it that  there are
C       no  implicit  vector  dependencies.  Compiler directives for
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
C       supplied with the standard SLAP distribution.
C
C
C *Cautions:
C       This routine assumes that the diagonal of A is all  non-zero
C       and that the operation DINV = 1.0/DIAG(A) will not underflow
C       or overflow.    This  is done so that the  loop  vectorizes.
C       Matrices  with zero or near zero or very  large entries will
C       have numerical difficulties  and  must  be fixed before this
C       routine is called.
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
C***END PROLOGUE  SSDS
C     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
C     .. Array Arguments ..
      REAL A(NELT), DINV(N)
      INTEGER IA(NELT), JA(NELT)
C     .. Local Scalars ..
      INTEGER ICOL
C***FIRST EXECUTABLE STATEMENT  SSDS
C
C         Assume the Diagonal elements are the first in each column.
C         This loop should *VECTORIZE*.  If it does not you may have
C         to add a compiler directive.  We do not check for a zero
C         (or near zero) diagonal element since this would interfere
C         with vectorization.  If this makes you nervous put a check
C         in!  It will run much slower.
C
      DO 10 ICOL = 1, N
         DINV(ICOL) = 1.0E0/A(JA(ICOL))
 10   CONTINUE
C
      RETURN
C------------- LAST LINE OF SSDS FOLLOWS ----------------------------
      END
