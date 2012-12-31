*DECK SSDSCL
      SUBROUTINE SSDSCL (N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB,
     +   ITOL)
C***BEGIN PROLOGUE  SSDSCL
C***PURPOSE  Diagonal Scaling of system Ax = b.
C            This routine scales (and unscales) the system  Ax = b
C            by symmetric diagonal scaling.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2E
C***TYPE      SINGLE PRECISION (SSDSCL-S, DSDSCL-D)
C***KEYWORDS  DIAGONAL, SLAP SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C    This routine scales (and unscales) the system Ax = b by symmetric
C    diagonal scaling.  The new system is:
C         -1/2  -1/2  1/2      -1/2
C        D    AD    (D   x) = D    b
C    when scaling is selected with the JOB parameter.  When unscaling
C    is selected this process is reversed.  The true solution is also
C    scaled or unscaled if ITOL is set appropriately, see below.
C
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB, ITOL
C     REAL    A(NELT), X(N), B(N), DINV(N)
C
C     CALL SSDSCL( N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB, ITOL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
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
C X      :INOUT    Real X(N).
C         Initial guess that will be later used in the iterative
C         solution.
C         of the scaled system.
C B      :INOUT    Real B(N).
C         Right hand side vector.
C DINV   :INOUT    Real DINV(N).
C         Upon return this array holds 1./DIAG(A).
C         This is an input if JOB = 0.
C JOB    :IN       Integer.
C         Flag indicating whether to scale or not.
C         JOB non-zero means do scaling.
C         JOB = 0 means do unscaling.
C ITOL   :IN       Integer.
C         Flag indicating what type of error estimation to do in the
C         iterative method.  When ITOL = 11 the exact solution from
C         common block SSLBLK will be used.  When the system is scaled
C         then the true solution must also be scaled.  If ITOL is not
C         11 then this vector is not referenced.
C
C *Common Blocks:
C SOLN    :INOUT   Real SOLN(N).  COMMON BLOCK /SSLBLK/
C         The true solution, SOLN, is scaled (or unscaled) if ITOL is
C         set to 11, see above.
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
C       and that the operation DINV = 1.0/DIAG(A)  will  not  under-
C       flow or overflow. This is done so that the loop  vectorizes.
C       Matrices  with zero or near zero or very  large entries will
C       have numerical difficulties  and  must  be fixed before this
C       routine is called.
C
C***SEE ALSO  SSDCG
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    SSLBLK
C***REVISION HISTORY  (YYMMDD)
C   871119  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
C   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   921113  Corrected C***CATEGORY line.  (FNF)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  SSDSCL
C     .. Scalar Arguments ..
      INTEGER ISYM, ITOL, JOB, N, NELT
C     .. Array Arguments ..
      REAL A(NELT), B(N), DINV(N), X(N)
      INTEGER IA(NELT), JA(NELT)
C     .. Arrays in Common ..
      REAL SOLN(1)
C     .. Local Scalars ..
      REAL DI
      INTEGER ICOL, J, JBGN, JEND
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     .. Common blocks ..
      COMMON /SSLBLK/ SOLN
C***FIRST EXECUTABLE STATEMENT  SSDSCL
C
C         SCALING...
C
      IF( JOB.NE.0 ) THEN
         DO 10 ICOL = 1, N
            DINV(ICOL) = 1.0E0/SQRT( A(JA(ICOL)) )
 10      CONTINUE
      ELSE
C
C         UNSCALING...
C
         DO 15 ICOL = 1, N
            DINV(ICOL) = 1.0E0/DINV(ICOL)
 15      CONTINUE
      ENDIF
C
      DO 30 ICOL = 1, N
         JBGN = JA(ICOL)
         JEND = JA(ICOL+1)-1
         DI = DINV(ICOL)
         DO 20 J = JBGN, JEND
            A(J) = DINV(IA(J))*A(J)*DI
 20      CONTINUE
 30   CONTINUE
C
      DO 40 ICOL = 1, N
         B(ICOL) = B(ICOL)*DINV(ICOL)
         X(ICOL) = X(ICOL)/DINV(ICOL)
 40   CONTINUE
C
C         Check to see if we need to scale the "true solution" as well.
C
      IF( ITOL.EQ.11 ) THEN
         DO 50 ICOL = 1, N
            SOLN(ICOL) = SOLN(ICOL)/DINV(ICOL)
 50      CONTINUE
      ENDIF
C
      RETURN
C------------- LAST LINE OF SSDSCL FOLLOWS ----------------------------
      END
