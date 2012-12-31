*DECK SSMTV
      SUBROUTINE SSMTV (N, X, Y, NELT, IA, JA, A, ISYM)
C***BEGIN PROLOGUE  SSMTV
C***PURPOSE  SLAP Column Format Sparse Matrix Transpose Vector Product.
C            Routine to calculate the sparse matrix vector product:
C            Y = A'*X, where ' denotes transpose.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D1B4
C***TYPE      SINGLE PRECISION (SSMTV-S, DSMTV-D)
C***KEYWORDS  MATRIX TRANSPOSE VECTOR MULTIPLY, SLAP, SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM
C     REAL     X(N), Y(N), A(NELT)
C
C     CALL SSMTV(N, X, Y, NELT, IA, JA, A, ISYM )
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C X      :IN       Real X(N).
C         The vector that should be multiplied by the transpose of
C         the matrix.
C Y      :OUT      Real Y(N).
C         The product of the transpose of the matrix and the vector.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Real A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
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
C       With  the SLAP  format  the "inner  loops" of  this  routine
C       should vectorize   on machines with   hardware  support  for
C       vector gather/scatter operations.  Your compiler may require
C       a  compiler directive  to  convince   it that there  are  no
C       implicit vector  dependencies.  Compiler directives  for the
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
C       with the standard SLAP distribution.
C
C *Cautions:
C     This   routine   assumes  that  the matrix A is stored in SLAP
C     Column format.  It does not check  for  this (for  speed)  and
C     evil, ugly, ornery and nasty things  will happen if the matrix
C     data  structure  is,  in fact, not SLAP Column.  Beware of the
C     wrong data structure!!!
C
C***SEE ALSO  SSMV
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
C***END PROLOGUE  SSMTV
C     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
C     .. Array Arguments ..
      REAL A(NELT), X(N), Y(N)
      INTEGER IA(NELT), JA(NELT)
C     .. Local Scalars ..
      INTEGER I, IBGN, ICOL, IEND, IROW, J, JBGN, JEND
C***FIRST EXECUTABLE STATEMENT  SSMTV
C
C         Zero out the result vector.
C
      DO 10 I = 1, N
         Y(I) = 0
 10   CONTINUE
C
C         Multiply by A-Transpose.
C         A-Transpose is stored by rows...
CVD$R NOCONCUR
      DO 30 IROW = 1, N
         IBGN = JA(IROW)
         IEND = JA(IROW+1)-1
CVD$ ASSOC
         DO 20 I = IBGN, IEND
            Y(IROW) = Y(IROW) + A(I)*X(IA(I))
 20      CONTINUE
 30   CONTINUE
C
      IF( ISYM.EQ.1 ) THEN
C
C         The matrix is non-symmetric.  Need to get the other half in...
C         This loops assumes that the diagonal is the first entry in
C         each column.
C
         DO 50 ICOL = 1, N
            JBGN = JA(ICOL)+1
            JEND = JA(ICOL+1)-1
            IF( JBGN.GT.JEND ) GOTO 50
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NODEPCHK
            DO 40 J = JBGN, JEND
               Y(IA(J)) = Y(IA(J)) + A(J)*X(ICOL)
 40         CONTINUE
 50      CONTINUE
      ENDIF
      RETURN
C------------- LAST LINE OF SSMTV FOLLOWS ----------------------------
      END
