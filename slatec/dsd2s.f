*DECK DSD2S
      SUBROUTINE DSD2S (N, NELT, IA, JA, A, ISYM, DINV)
C***BEGIN PROLOGUE  DSD2S
C***PURPOSE  Diagonal Scaling Preconditioner SLAP Normal Eqns Set Up.
C            Routine to compute the inverse of the diagonal of the
C            matrix A*A', where A is stored in SLAP-Column format.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2E
C***TYPE      DOUBLE PRECISION (SSD2S-S, DSD2S-D)
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
C     DOUBLE PRECISION A(NELT), DINV(N)
C
C     CALL DSD2S( N, NELT, IA, JA, A, ISYM, DINV )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C DINV   :OUT      Double Precision DINV(N).
C         Upon return this array holds 1./DIAG(A*A').
C
C *Description
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
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
C       and that the operation DINV = 1.0/DIAG(A*A') will not under-
C       flow or overflow. This is done so that the loop  vectorizes.
C       Matrices  with zero or near zero or very  large entries will
C       have numerical difficulties  and  must  be fixed before this
C       routine is called.
C
C***SEE ALSO  DSDCGN
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   920511  Added complete declaration section.  (WRB)
C   921113  Corrected C***CATEGORY line.  (FNF)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  DSD2S
C     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), DINV(N)
      INTEGER IA(NELT), JA(NELT)
C     .. Local Scalars ..
      INTEGER I, K, KBGN, KEND
C***FIRST EXECUTABLE STATEMENT  DSD2S
      DO 10 I = 1, N
         DINV(I) = 0
 10   CONTINUE
C
C         Loop over each column.
CVD$R NOCONCUR
      DO 40 I = 1, N
         KBGN = JA(I)
         KEND = JA(I+1) - 1
C
C         Add in the contributions for each row that has a non-zero
C         in this column.
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NODEPCHK
         DO 20 K = KBGN, KEND
            DINV(IA(K)) = DINV(IA(K)) + A(K)**2
 20      CONTINUE
         IF( ISYM.EQ.1 ) THEN
C
C         Lower triangle stored by columns => upper triangle stored by
C         rows with Diagonal being the first entry.  Loop across the
C         rest of the row.
            KBGN = KBGN + 1
            IF( KBGN.LE.KEND ) THEN
               DO 30 K = KBGN, KEND
                  DINV(I) = DINV(I) + A(K)**2
 30            CONTINUE
            ENDIF
         ENDIF
 40   CONTINUE
      DO 50 I=1,N
         DINV(I) = 1.0D0/DINV(I)
 50   CONTINUE
C
      RETURN
C------------- LAST LINE OF DSD2S FOLLOWS ----------------------------
      END
