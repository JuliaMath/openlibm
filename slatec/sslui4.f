*DECK SSLUI4
      SUBROUTINE SSLUI4 (N, B, X, IL, JL, L, DINV, IU, JU, U)
C***BEGIN PROLOGUE  SSLUI4
C***PURPOSE  SLAP Backsolve for LDU Factorization.
C            Routine to solve a system of the form  (L*D*U)' X = B,
C            where L is a unit lower triangular matrix, D is a diagonal
C            matrix, and U is a unit upper triangular matrix and '
C            denotes transpose.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2E
C***TYPE      SINGLE PRECISION (SSLUI4-S, DSLUI4-D)
C***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE,
C             SLAP, SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER N, IL(NL), JL(NL), IU(NU), JU(NU)
C     REAL    B(N), X(N), L(NL), DINV(N), U(NU)
C
C     CALL SSLUI4( N, B, X, IL, JL, L, DINV, IU, JU, U )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right hand side.
C X      :OUT      Real X(N).
C         Solution of (L*D*U)trans x = b.
C IL     :IN       Integer IL(NL).
C JL     :IN       Integer JL(NL).
C L      :IN       Real     L(NL).
C         IL, JL, L contain the unit lower triangular  factor of the
C         incomplete decomposition of some matrix stored in SLAP Row
C         format.  The diagonal of ones *IS* stored.  This structure
C         can    be  set  up  by   the  SSILUS  routine.   See   the
C         "Description",  below for  more  details about  the   SLAP
C         format.  (NL is the number of non-zeros in the L array.)
C DINV   :IN       Real DINV(N).
C         Inverse of the diagonal matrix D.
C IU     :IN       Integer IU(NU).
C JU     :IN       Integer JU(NU).
C U      :IN       Real     U(NU).
C         IU, JU, U contain the  unit upper triangular factor of the
C         incomplete  decomposition of some  matrix stored  in  SLAP
C         Column  format.   The diagonal of  ones *IS* stored.  This
C         structure can be set up by the  SSILUS routine.  See   the
C         "Description",  below for  more  details  about  the  SLAP
C         format.  (NU is the number of non-zeros in the U array.)
C
C *Description:
C       This routine is supplied with the SLAP package as  a routine
C       to  perform  the  MTSOLV  operation  in  the SBCG  iteration
C       routine for the  driver SSLUBC.   It must  be called via the
C       SLAP  MTSOLV calling  sequence convention interface  routine
C       SSLUTI.
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
C               **** SLAP MSOLVE CALLING CONVENTION ****
C
C       IL, JL, L should contain the unit lower triangular factor of
C       the incomplete decomposition of the A matrix  stored in SLAP
C       Row format.  IU, JU, U should contain  the unit upper factor
C       of the  incomplete decomposition of  the A matrix  stored in
C       SLAP Column format This ILU factorization can be computed by
C       the SSILUS routine. The diagonals (which are all one's) are
C       stored.
C
C       =================== S L A P Column format ==================
C
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
C       ==================== S L A P Row format ====================
C
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the real
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
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
C***SEE ALSO  SSILUS
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
C   921113  Corrected C***CATEGORY line.  (FNF)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  SSLUI4
C     .. Scalar Arguments ..
      INTEGER N
C     .. Array Arguments ..
      REAL B(N), DINV(N), L(*), U(*), X(N)
      INTEGER IL(*), IU(*), JL(*), JU(*)
C     .. Local Scalars ..
      INTEGER I, ICOL, IROW, J, JBGN, JEND
C***FIRST EXECUTABLE STATEMENT  SSLUI4
      DO 10 I=1,N
         X(I) = B(I)
 10   CONTINUE
C
C         Solve  U'*Y = X,  storing result in X, U stored by columns.
      DO 80 IROW = 2, N
         JBGN = JU(IROW)
         JEND = JU(IROW+1) - 1
         IF( JBGN.LE.JEND ) THEN
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ ASSOC
CVD$ NODEPCHK
            DO 70 J = JBGN, JEND
               X(IROW) = X(IROW) - U(J)*X(IU(J))
 70         CONTINUE
         ENDIF
 80   CONTINUE
C
C         Solve  D*Z = Y,  storing result in X.
      DO 90 I = 1, N
         X(I) = X(I)*DINV(I)
 90   CONTINUE
C
C         Solve  L'*X = Z, L stored by rows.
      DO 110 ICOL = N, 2, -1
         JBGN = IL(ICOL)
         JEND = IL(ICOL+1) - 1
         IF( JBGN.LE.JEND ) THEN
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NODEPCHK
            DO 100 J = JBGN, JEND
               X(JL(J)) = X(JL(J)) - L(J)*X(ICOL)
 100        CONTINUE
         ENDIF
 110  CONTINUE
      RETURN
C------------- LAST LINE OF SSLUI4 FOLLOWS ----------------------------
      END
