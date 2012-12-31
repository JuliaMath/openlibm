*DECK DS2LT
      SUBROUTINE DS2LT (N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL)
C***BEGIN PROLOGUE  DS2LT
C***PURPOSE  Lower Triangle Preconditioner SLAP Set Up.
C            Routine to store the lower triangle of a matrix stored
C            in the SLAP Column format.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2E
C***TYPE      DOUBLE PRECISION (SS2LT-S, DS2LT-D)
C***KEYWORDS  LINEAR SYSTEM, LOWER TRIANGLE, SLAP SPARSE
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
C     INTEGER NEL, IEL(NEL), JEL(NEL)
C     DOUBLE PRECISION A(NELT), EL(NEL)
C
C     CALL DS2LT( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of non-zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C NEL    :OUT      Integer.
C         Number of non-zeros in the lower triangle of A.   Also
C         corresponds to the length of the IEL, JEL, EL arrays.
C IEL    :OUT      Integer IEL(NEL).
C JEL    :OUT      Integer JEL(NEL).
C EL     :OUT      Double Precision     EL(NEL).
C         IEL, JEL, EL contain the lower triangle of the A matrix
C         stored in SLAP Column format.  See "Description", below,
C         for more details bout the SLAP Column format.
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
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  DS2LT
C     .. Scalar Arguments ..
      INTEGER ISYM, N, NEL, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), EL(NELT)
      INTEGER IA(NELT), IEL(NEL), JA(NELT), JEL(NEL)
C     .. Local Scalars ..
      INTEGER I, ICOL, J, JBGN, JEND
C***FIRST EXECUTABLE STATEMENT  DS2LT
      IF( ISYM.EQ.0 ) THEN
C
C         The matrix is stored non-symmetricly.  Pick out the lower
C         triangle.
C
         NEL = 0
         DO 20 ICOL = 1, N
            JEL(ICOL) = NEL+1
            JBGN = JA(ICOL)
            JEND = JA(ICOL+1)-1
CVD$ NOVECTOR
            DO 10 J = JBGN, JEND
               IF( IA(J).GE.ICOL ) THEN
                  NEL = NEL + 1
                  IEL(NEL) = IA(J)
                  EL(NEL)  = A(J)
               ENDIF
 10         CONTINUE
 20      CONTINUE
         JEL(N+1) = NEL+1
      ELSE
C
C         The matrix is symmetric and only the lower triangle is
C         stored.  Copy it to IEL, JEL, EL.
C
         NEL = NELT
         DO 30 I = 1, NELT
            IEL(I) = IA(I)
            EL(I) = A(I)
 30      CONTINUE
         DO 40 I = 1, N+1
            JEL(I) = JA(I)
 40      CONTINUE
      ENDIF
      RETURN
C------------- LAST LINE OF DS2LT FOLLOWS ----------------------------
      END
